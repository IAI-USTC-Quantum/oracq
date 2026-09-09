"RIR 的封闭操作集、类型、别名、资源和无环调用检查。"

from __future__ import annotations

import math
import re

from pyqecclang.infrastructure.ir import (
    VERSION,
    Adjoint,
    Call,
    Control,
    Load,
    Primitive,
    Program,
    Ref,
    RegType,
    Repeat,
    ValidationError,
)

KINDS = {"bits", "uint", "sint", "rational"}
UNARY = {"h", "x", "y", "z", "s", "t", "rx", "ry", "rz", "phase"}
ROTATIONS = {"rx", "ry", "rz", "phase"}
BINARY = {"xor", "swap"}
IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")


def require(condition, message):
    if not condition:
        raise ValidationError(message)


def integer(value, lo, hi, message):
    require(type(value) is int and lo <= value <= hi, message)


def name(value):
    require(
        isinstance(value, str) and IDENTIFIER.fullmatch(value) is not None, f"非法标识符：{value!r}"
    )


def reg_type(dtype: RegType):
    require(isinstance(dtype, RegType) and dtype.kind in KINDS, "未知寄存器类型")
    integer(dtype.width, 0, 64, "每个寄存器或视图的宽度必须为 0..64")


def locations(ref: Ref):
    return tuple((s.register, bit) for s in ref.parts for bit in range(s.start, s.start + s.width))


def _validate(program: Program) -> Program:
    require(isinstance(program, Program) and type(program.modules) is tuple, "需要不可变 Program")
    require(program.version in {"0.1", "0.2", VERSION}, f"不支持 RIR 版本：{program.version}")
    name(program.entry)
    modules = program.module_map
    require(len(modules) == len(program.modules), "模块名重复")
    require(program.entry in modules, "入口模块不存在")
    graph = {key: set() for key in modules}

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

        def check_ref(ref, regs=regs):
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

        def distinct(refs, protected):
            used = set()
            for ref in refs:
                current = check_ref(ref)
                require(not current & used, "操作数存在别名或重叠")
                require(not current & protected, "操作数修改了受保护的控制寄存器")
                used |= current

        def body(nodes, protected=frozenset(), depth=0, resources=resources, module=module):
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
                            type(node.angle) in (float, int) and math.isfinite(node.angle),
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
                    for actual, formal in zip(node.resources, target.resources, strict=True):
                        require(
                            actual in resources and resources[actual] == formal.type,
                            "模块 QRAM 实参类型不符",
                        )
                    graph[module.name].add(node.module)
                elif isinstance(node, Repeat):
                    integer(node.count, 0, 2**63 - 1, "重复次数必须为 0..2^63-1")
                    body(node.body, protected, depth + 1)
                elif isinstance(node, Control):
                    locs = check_ref(node.register)
                    require(bool(locs), "控制寄存器不能为空")
                    integer(node.value, 0, (1 << node.register.width) - 1, "控制值越界")
                    require(not locs & protected, "嵌套控制寄存器重叠")
                    body(node.body, protected | locs, depth + 1)
                elif isinstance(node, Adjoint):
                    body(node.body, protected, depth + 1)
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

    visited, active = set(), set()

    def visit(key):
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
    from pyqecclang.infrastructure.linking import capability_table

    inferred = capability_table(program)

    def check_demands(nodes, controlled=False, inverse=False):
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
    try:
        result = _validate(program)
        if require_closed:
            from pyqecclang.infrastructure.linking import unresolved

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
