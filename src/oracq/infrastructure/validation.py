"RIR 的封闭操作集、类型、别名、资源和无环调用检查。"

from __future__ import annotations

import json
import math
import re
from typing import cast

from oracq.infrastructure.ir import (
    QRAM,
    VERSION,
    Adjoint,
    Call,
    Control,
    Instruction,
    Load,
    Module,
    Primitive,
    Program,
    Ref,
    RegType,
    Repeat,
    Store,
    ValidationError,
)

KINDS = {"bits", "uint", "sint", "rational"}
"""寄存器类型允许的种类集合。"""
UNARY = {"h", "x", "y", "z", "s", "t", "rx", "ry", "rz", "phase"}
"""允许的一元量子基元操作名集合。"""
ROTATIONS = {"rx", "ry", "rz", "phase"}
"""带角度参数的旋转类基元名集合，是 ``UNARY`` 的子集。"""
BINARY = {"xor", "swap"}
"""允许的二元量子基元操作名集合。"""
IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")
"""合法标识符的正则模式：以字母或下划线开头，后随字母、数字或下划线。"""


def require(condition: object, message: str) -> None:
    """断言结构检查条件成立，否则抛出 ``ValidationError``。

    Args:
        condition: 按布尔语义解释的检查结果。
        message: 失败时写入异常的错误说明。

    Raises:
        ValidationError: ``condition`` 为假。
    """
    if not condition:
        raise ValidationError(message)


def integer(value: object, lo: int, hi: int, message: str) -> None:
    """断言 ``value`` 是位于闭区间 ``lo..hi`` 内的整数。

    Args:
        value: 待检查的数值。
        lo: 允许的最小值（含端点）。
        hi: 允许的最大值（含端点）。
        message: 失败时写入异常的错误说明。

    Raises:
        ValidationError: ``value`` 不是 ``int`` 或越出区间。
    """
    require(type(value) is int and lo <= value <= hi, message)


def name(value: object) -> None:
    """断言 ``value`` 是匹配 ``IDENTIFIER`` 模式的合法标识符字符串。

    Args:
        value: 待检查的标识符，如寄存器名、资源名或模块名。

    Raises:
        ValidationError: ``value`` 不是字符串或包含非法字符。
    """
    require(
        isinstance(value, str) and IDENTIFIER.fullmatch(value) is not None, f"非法标识符：{value!r}"
    )


def reg_type(dtype: RegType) -> None:
    """断言 ``dtype`` 是种类已知、宽度为 0..64 的寄存器类型。

    Args:
        dtype: 待检查的存储类型对象。

    Raises:
        ValidationError: ``dtype`` 不是 ``RegType``、种类未知或宽度越界。
    """
    require(isinstance(dtype, RegType) and dtype.kind in KINDS, "未知寄存器类型")
    integer(dtype.width, 0, 64, "每个寄存器或视图的宽度必须为 0..64")


def locations(ref: Ref) -> tuple[tuple[str, int], ...]:
    """展开视图引用覆盖的全部量子位坐标。

    Args:
        ref: 待展开的寄存器视图。

    Returns:
        tuple: ``(寄存器名, 位序号)`` 二元组按视图分段顺序组成的元组。
    """
    return tuple((s.register, bit) for s in ref.parts for bit in range(s.start, s.start + s.width))


def capture_map(module: Module) -> dict[str, str]:
    """解析链接器的资源来源属性；保证映射引用存在且一一对应。"""
    raw = dict(module.attributes).get("binding_captures", "{}")
    require(isinstance(raw, str), "binding_captures 必须是 JSON 字符串")
    try:
        result = json.loads(cast(str, raw))
    except ValueError as exc:
        raise ValidationError("binding_captures 不是有效 JSON") from exc
    require(isinstance(result, dict), "binding_captures 必须是对象")
    existing = {r.name for r in module.resources}
    for logical, local in result.items():
        name(logical)
        name(local)
        require(local in existing, "binding_captures 引用了不存在的资源")
    require(len(set(result.values())) == len(result), "binding_captures 局部资源名重复")
    return cast(dict[str, str], result)


def _validate(program: Program) -> Program:
    """对程序执行全部跨节点结构检查，通过后原样返回。"""
    require(isinstance(program, Program) and type(program.modules) is tuple, "需要不可变 Program")
    require(program.version in {"0.1", "0.2", VERSION}, f"不支持 RIR 版本：{program.version}")
    name(program.entry)
    modules = program.module_map
    require(len(modules) == len(program.modules), "模块名重复")
    require(program.entry in modules, "入口模块不存在")
    graph: dict[str, set[str]] = {key: set() for key in modules}

    for module in program.modules:
        name(module.name)
        require(
            type(module.registers) is tuple and type(module.resources) is tuple,
            "模块签名必须不可变",
        )
        require(type(module.attributes) is tuple, "模块属性必须不可变")
        require(
            all(type(pair) is tuple and len(pair) == 2 for pair in module.attributes),
            "属性条目必须是二元组",
        )
        attributes = dict(module.attributes)
        require(len(attributes) == len(module.attributes), "模块属性键重复")
        for key, value in attributes.items():
            name(key)
            require(type(value) in (str, int, float, bool), "属性必须是标量")
            if type(value) is float:
                require(math.isfinite(value), "属性浮点数必须有限")
        for capability in ("supports_adjoint", "supports_controlled"):
            if capability in attributes:
                require(type(attributes[capability]) is bool, "能力字段必须为布尔值")
        require(type(module.locals) is tuple, "局部寄存器列表必须不可变")
        require(program.version == VERSION or not module.locals, "旧 RIR 不支持局部寄存器")
        if module.body is None:
            require(not module.locals, "开放声明不能定义私有工作区")
        all_registers = module.registers + module.locals
        regs = {reg.name: reg.type for reg in all_registers}
        resources = {res.name: res.type for res in module.resources}
        require(len(regs) == len(all_registers), "寄存器名重复")
        require(len(resources) == len(module.resources), "资源名重复")
        require(not regs.keys() & resources.keys(), "寄存器与资源名字冲突")
        for reg in all_registers:
            name(reg.name)
            reg_type(reg.type)
        require(sum(reg.type.width for reg in module.registers) > 0, "模块必须有非空量子接口")
        for resource in module.resources:
            name(resource.name)
            integer(resource.type.address_width, 1, 64, "QRAM 地址宽度必须为 1..64")
            integer(resource.type.data_width, 1, 64, "QRAM 数据宽度必须为 1..64")
        if "binding_captures" in attributes:
            capture_map(module)

        def check_ref(ref: Ref, regs: dict[str, RegType] = regs) -> set[tuple[str, int]]:
            """校验单个视图的合法性与无重叠，返回其覆盖的量子位集合。"""
            require(isinstance(ref, Ref), "需要寄存器视图")
            reg_type(ref.type)
            require(type(ref.parts) is tuple, "视图必须不可变")
            require(sum(s.width for s in ref.parts) == ref.width, "视图宽度与分段不符")
            for span in ref.parts:
                require(span.register in regs, f"未知寄存器：{span.register}")
                width = regs[span.register].width
                integer(span.start, 0, width, "视图起点越界")
                integer(span.width, 0, width - span.start, "视图终点越界")
            locs = locations(ref)
            require(len(set(locs)) == len(locs), "视图内存在重叠量子位")
            return set(locs)

        def distinct(refs: tuple[Ref, ...], protected: frozenset[tuple[str, int]]) -> None:
            """断言各操作数互不重叠，且不修改受保护的控制量子位。"""
            used: set[tuple[str, int]] = set()
            for ref in refs:
                current = check_ref(ref)
                require(not current & used, "操作数存在别名或重叠")
                require(not current & protected, "操作数修改了受保护的控制寄存器")
                used |= current

        def body(
            nodes: tuple[Instruction, ...],
            protected: frozenset[tuple[str, int]] = frozenset(),
            depth: int = 0,
            resources: dict[str, QRAM] = resources,
            module: Module = module,
            unitary: bool = True,
        ) -> None:
            """递归校验指令体：基元元数、别名、控制保护与调用匹配。"""
            require(depth < 128, "嵌套深度超过 127")
            require(type(nodes) is tuple, "指令体必须不可变")
            for node in nodes:
                if isinstance(node, Primitive):
                    require(type(node.operands) is tuple, "基元操作数必须不可变")
                    require(node.op in UNARY | BINARY | {"add_const", "gphase"}, "未知基元")
                    expected = 2 if node.op in BINARY else 0 if node.op == "gphase" else 1
                    require(len(node.operands) == expected, "基元参数数量不匹配")
                    distinct(node.operands, protected)
                    if node.op in BINARY:
                        require(
                            node.operands[0].width == node.operands[1].width, "二元操作宽度不符"
                        )
                    if node.op in ROTATIONS | {"gphase"}:
                        require(
                            # and 短路保证进入 isfinite 时 angle 已是有限实数类型。
                            type(node.angle) in (float, int) and math.isfinite(cast(float, node.angle)),
                            "旋转角必须是有限实数",
                        )
                    else:
                        require(node.angle is None, "此基元不接受角度")
                    if node.op == "add_const":
                        require(node.operands[0].type.kind == "uint", "add_const 需要 uint")
                        integer(
                            node.value,
                            0,
                            (1 << node.operands[0].width) - 1,
                            "加法常量超出寄存器范围",
                        )
                    else:
                        require(node.value is None, "此基元不接受整数参数")
                elif isinstance(node, Load):
                    require(node.resource in resources, "QRAM 资源没有声明")
                    spec = resources[node.resource]
                    distinct((node.address, node.data), protected)
                    require(node.address.width == spec.address_width, "QRAM 地址宽度不符")
                    require(node.data.width == spec.data_width, "QRAM 数据宽度不符")
                elif isinstance(node, Store):
                    require(unitary, "Store 是非酉副作用，不能出现在 Control 或 Adjoint 体内")
                    require(program.version == VERSION, "旧 RIR 不支持 Store")
                    require(node.resource in resources, "QRAM 资源没有声明")
                    spec = resources[node.resource]
                    distinct((node.address, node.data), protected)
                    require(node.address.width == spec.address_width, "QRAM 地址宽度不符")
                    require(node.data.width == spec.data_width, "QRAM 数据宽度不符")
                elif isinstance(node, Call):
                    require(
                        type(node.arguments) is tuple and type(node.resources) is tuple,
                        "调用参数必须不可变",
                    )
                    require(node.module in modules, f"未知模块：{node.module}")
                    target = modules[node.module]
                    require(len(node.arguments) == len(target.registers), "模块量子参数数量不符")
                    require(len(node.resources) == len(target.resources), "模块资源参数数量不符")
                    distinct(node.arguments, protected)
                    for actual, formal in zip(node.arguments, target.registers, strict=True):
                        require(actual.type == formal.type, f"模块参数类型不符：{formal.name}")
                    # actual/formal 在上一循环绑定为 Ref/Register，本循环承载 str/Resource。
                    for actual, formal in zip(  # type: ignore[assignment]
                        node.resources, target.resources, strict=True
                    ):
                        require(
                            actual in resources and resources[cast(str, actual)] == formal.type,
                            "模块 QRAM 实参类型不符",
                        )
                    graph[module.name].add(node.module)
                elif isinstance(node, Repeat):
                    integer(node.count, 0, 2**63 - 1, "重复次数必须为 0..2^63-1")
                    body(node.body, protected, depth + 1, unitary=unitary)
                elif isinstance(node, Control):
                    locs = check_ref(node.register)
                    require(bool(locs), "控制寄存器不能为空")
                    integer(node.value, 0, (1 << node.register.width) - 1, "控制值越界")
                    require(not locs & protected, "嵌套控制寄存器重叠")
                    body(node.body, protected | locs, depth + 1, unitary=False)
                elif isinstance(node, Adjoint):
                    body(node.body, protected, depth + 1, unitary=False)
                else:
                    raise ValidationError(f"未知指令类型：{type(node).__name__}")

        if module.body is None:
            require(program.version != "0.1", "RIR 0.1 不支持开放声明")
            require(
                isinstance(attributes.get("oracle_paradigm"), str),
                "开放声明必须指定 oracle_paradigm",
            )
        else:
            body(module.body)

    visited: set[str]
    active: set[str]
    visited, active = set(), set()

    def visit(key: str) -> None:
        """沿调用图深度优先检测递归与过深的模块调用。"""
        require(key not in active, "模块调用图存在递归")
        require(len(active) < 128, "模块调用深度超过 127")
        if key in visited:
            return
        active.add(key)
        for child in sorted(graph[key]):
            visit(child)
        active.remove(key)
        visited.add(key)

    for key in modules:
        visit(key)
    from oracq.infrastructure.linking import capability_table

    inferred = capability_table(program)

    def check_demands(
        nodes: tuple[Instruction, ...] | None, controlled: bool = False, inverse: bool = False
    ) -> None:
        """检查控制或伴随语境下调用的模块具备相应能力。"""
        for node in nodes or ():
            if isinstance(node, Call):
                if controlled:
                    require(
                        inferred[node.module]["supports_controlled"],
                        f"{node.module} 缺少 supports_controlled",
                    )
                if inverse:
                    require(
                        inferred[node.module]["supports_adjoint"],
                        f"{node.module} 缺少 supports_adjoint",
                    )
            elif isinstance(node, Control):
                check_demands(node.body, True, inverse)
            elif isinstance(node, Adjoint):
                check_demands(node.body, controlled, not inverse)
            elif isinstance(node, Repeat):
                check_demands(node.body, controlled, inverse)

    for module in program.modules:
        check_demands(module.body)
    return program


def validate(program: Program, *, require_closed: bool = False) -> Program:
    """对 RIR 程序执行完整结构验证并原样返回。

    覆盖模块签名与属性、寄存器和视图规则、基元元数与别名约束、
    控制寄存器保护、QRAM 资源匹配、无环调用图等跨节点规则，并检查
    控制或伴随语境下调用的模块具备相应能力；规则与
    docs/reference/rir.md 保持一致。

    Args:
        program: 待验证的 ``Program``。
        require_closed: 为真时还要求程序不含未绑定的开放 oracle 声明。

    Returns:
        Program: 通过验证的同一 ``program`` 对象。

    Raises:
        ValidationError: 程序违反结构或语义约束；或捕获到表明输入不是
            合法 RIR 对象的标准异常。
    """
    try:
        result = _validate(program)
        if require_closed:
            from oracq.infrastructure.linking import unresolved

            missing = unresolved(result)
            if missing:
                raise ValidationError(
                    "未绑定 oracle："
                    + "; ".join(f"{r.name} ({r.paradigm}; {' -> '.join(r.path)})" for r in missing)
                )
        return result
    except ValidationError:
        raise
    except (AttributeError, TypeError, KeyError, ValueError, OverflowError, RecursionError) as exc:
        raise ValidationError(f"非法 RIR 对象：{exc}") from exc
