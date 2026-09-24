"开放 oracle 的缺口分析、部分绑定和 QRAM 捕获资源提升。"

from __future__ import annotations

import hashlib
import json
from collections import deque
from collections.abc import Iterator
from dataclasses import asdict, dataclass, replace
from typing import NoReturn, cast

from oracq.infrastructure.builder import Operation
from oracq.infrastructure.ir import (
    VERSION,
    Adjoint,
    Call,
    Control,
    Instruction,
    Module,
    Program,
    Ref,
    Repeat,
    Resource,
    Span,
    Store,
    ValidationError,
)
from oracq.infrastructure.validation import capture_map, validate


@dataclass(frozen=True)
class OracleRequirement:
    """描述一个未实现的开放 oracle 声明及其调用位置。

    Attributes:
        name: 开放声明的模块名，即绑定时使用的槽名。
        paradigm: 声明指定的 ``oracle_paradigm`` 范式。
        path: 从入口模块到该声明的最短调用路径上的模块名序列。
        registers: 声明的形式寄存器签名元组。
        attributes: 声明的属性条目元组，每项为键值二元组。
    """
    name: str
    paradigm: str
    path: tuple[str, ...]
    registers: tuple
    attributes: tuple


@dataclass(frozen=True)
class Binding:
    """把一个开放声明槽绑定到实现操作，并可携带 QRAM 资源映射。

    Attributes:
        operation: 提供实现的 ``Operation``；其模块名必须不同于声明槽名。
        resources: 实现的资源形式参数到入口逻辑资源名的映射；未显式
            覆盖的资源按 ``槽名__资源名`` 自动捕获并提升为入口资源。
    """
    operation: Operation
    resources: dict[str, str] | None = None


@dataclass(frozen=True)
class BindingIssue:
    """绑定失败的结构化原因；不声称证明 oracle 的数学语义。"""

    code: str
    path: tuple[str, ...]
    expected: object
    actual: object
    message: str


class BindingError(ValidationError):
    """仍兼容 ValidationError 的绑定诊断。"""

    def __init__(self, issues: tuple[BindingIssue, ...]) -> None:
        self.issues = issues
        super().__init__("; ".join(f"{i.code} {' -> '.join(i.path)}: {i.message}" for i in issues))


@dataclass(frozen=True)
class BindingReport:
    """绑定清单、程序指纹与剩余依赖；只保存可交换的数据。"""

    source_digest: str
    result_digest: str | None
    bindings: tuple[tuple[str, str], ...]
    resource_mappings: tuple[tuple[str, tuple[tuple[str, str], ...]], ...]
    remaining: tuple[OracleRequirement, ...]
    resources: tuple[Resource, ...]
    issues: tuple[BindingIssue, ...] = ()

    @property
    def ok(self) -> bool:
        """本次绑定是否通过结构检查；不要求所有槽位均已闭合。"""
        return not self.issues

    def to_dict(self) -> dict[str, object]:
        """输出 JSON 友好的报告，独立于可执行 RIR 保存。"""
        return {"ok": self.ok, **asdict(self)}


@dataclass(frozen=True)
class BindingResult:
    """一次绑定的程序与报告，避免先检查后再次链接。"""

    program: Program | None
    report: BindingReport

    def require(self) -> Program:
        """绑定成功时取得程序，否则抛出结构化问题。"""
        if self.report.issues:
            raise BindingError(self.report.issues)
        if self.program is None:
            raise ValidationError("绑定结果缺少程序")
        return self.program


def bind_with_report(
    program: Program | Operation, bindings: dict[str, Binding | Operation]
) -> BindingResult:
    """执行一次纯函数绑定，返回可直接使用的程序及可追溯报告。

    无效输入程序仍由 validate 拒绝；候选实现的绑定错误记录在报告中。
    数据内容不属于 RIR 指纹，运行记录须另行标识 QRAM 数据快照。
    """
    from oracq.infrastructure.serialization import dumps

    source = program.program() if isinstance(program, Operation) else validate(program)
    digest = hashlib.sha256(dumps(source).encode()).hexdigest()
    entries = tuple(sorted(
        (slot, (item.operation if isinstance(item, Binding) else item).module.name)
        for slot, item in bindings.items()
    ))
    mappings = tuple(sorted(
        (slot, tuple(sorted((item.resources or {}).items())) if isinstance(item, Binding) else ())
        for slot, item in bindings.items()
    ))
    try:
        linked = bind(source, bindings)
    except ValidationError as exc:
        issues = exc.issues if isinstance(exc, BindingError) else (
            BindingIssue("BIND_STRUCTURE", (source.entry,), "valid linked program", None, str(exc)),
        )
        return BindingResult(None, BindingReport(
            digest, None, entries, mappings, unresolved(source), source.main.resources, issues,
        ))
    connections = []
    for slot, _ in entries:
        wrapper = linked.module_map[slot]
        if wrapper.body is None or not wrapper.body or not isinstance(wrapper.body[0], Call):
            continue
        call = wrapper.body[0]
        logical_names = {local: logical for logical, local in capture_map(wrapper).items()}
        connections.append((slot, tuple(
            (formal.name, logical_names.get(actual, actual))
            for formal, actual in zip(linked.module_map[call.module].resources, call.resources, strict=True)
        )))
    return BindingResult(linked, BindingReport(
        digest, hashlib.sha256(dumps(linked).encode()).hexdigest(), entries, tuple(connections),
        unresolved(linked), linked.main.resources,
    ))


def calls(nodes: tuple[Instruction, ...] | None) -> Iterator[Call]:
    """按出现顺序产出指令体（含嵌套结构块）中的全部模块调用节点。

    Args:
        nodes: 待扫描的指令体；``None`` 或空元组视为无内容。

    Returns:
        Iterator[Call]: 按出现顺序惰性产出的 ``Call`` 节点，含嵌套块内。
    """
    for node in nodes or ():
        if isinstance(node, Call):
            yield node
        elif isinstance(node, (Repeat, Control, Adjoint)):
            yield from calls(node.body)


def stores(nodes: tuple[Instruction, ...] | None) -> Iterator[Store]:
    """列出指令体（含嵌套结构块）中的全部 Store 副作用。

    Args:
        nodes: 待扫描的指令体；``None`` 或空元组视为无内容。

    Returns:
        Iterator[Store]: 按出现顺序惰性产出的 ``Store`` 节点，含嵌套块内。
    """
    for node in nodes or ():
        if isinstance(node, Store):
            yield node
        elif isinstance(node, (Repeat, Control, Adjoint)):
            yield from stores(node.body)


def uses_store(program: Program) -> bool:
    """入口可达的模块中是否存在 QRAM 随机写。

    Args:
        program: 待检查的 RIR 程序，从入口沿调用图遍历。

    Returns:
        bool: 任一入口可达模块体内存在 ``Store`` 节点时为 True。
    """
    modules = program.module_map
    pending, seen = [program.entry], set()
    while pending:
        key = pending.pop()
        if key in seen:
            continue
        seen.add(key)
        body = modules[key].body
        if body is not None and any(True for _ in stores(body)):
            return True
        pending.extend(call.module for call in calls(body))
    return False


def unresolved(program: Program) -> tuple[OracleRequirement, ...]:
    """列出入口结构可达的未实现声明，并给出首条最短调用路径。

    Args:
        program: 待分析的 RIR 程序；内部先做结构校验。

    Returns:
        tuple[OracleRequirement, ...]: 按声明名排序的开放声明需求列表，每项
        记录从入口出发的首条最短调用路径。
    """
    validate(program)
    modules = program.module_map
    pending: deque[tuple[str, tuple[str, ...]]] = deque([(program.entry, (program.entry,))])
    visited, result = set(), []
    while pending:
        key, path = pending.popleft()
        if key in visited:
            continue
        visited.add(key)
        module = modules[key]
        if module.body is None:
            result.append(
                OracleRequirement(
                    key,
                    # 校验保证开放声明的 oracle_paradigm 属性为 str。
                    cast(str, dict(module.attributes)["oracle_paradigm"]),
                    path,
                    module.registers,
                    module.attributes,
                )
            )
        else:
            for call in calls(module.body):
                pending.append((call.module, path + (call.module,)))
    return tuple(sorted(result, key=lambda item: item.name))


def capability_table(program: Program) -> dict[str, dict[str, bool]]:
    """一次遍历推导全部模块的变换能力。

    Args:
        program: 待分析的 RIR 程序，遍历其全部模块。

    Returns:
        dict[str, dict[str, bool]]: 模块名到能力字典的映射；含 ``Store``
        副作用的模块各项能力均为 False。
    """
    cache: dict[str, dict[str, bool]]
    modules, cache = program.module_map, {}

    def infer(name: str) -> dict[str, bool]:
        """按声明与被调模块能力保守合取，推导单个模块的能力字典。"""
        if name in cache:
            return cache[name]
        module = modules[name]
        attrs = dict(module.attributes)
        # 校验保证 supports_adjoint/supports_controlled 属性存在时必为布尔值。
        result: dict[str, bool] = {
            cap: cast(bool, attrs.get(cap, True)) for cap in ("supports_adjoint", "supports_controlled")
        }
        if module.body is not None and any(True for _ in stores(module.body)):
            result = {cap: False for cap in result}
        for call in calls(module.body):
            child = infer(call.module)
            result = {cap: result[cap] and child[cap] for cap in result}
        cache[name] = result
        return result

    for key in modules:
        infer(key)
    return cache


def capabilities(program: Program, key: str | None = None) -> dict[str, bool]:
    """返回指定模块（默认入口）的变换能力字典。

    能力取值由模块声明与被调用模块的能力保守合取推导，见
    ``capability_table``。

    Args:
        program: 待分析的 ``Program``。
        key: 目标模块名；省略时使用入口模块。

    Returns:
        dict: 以 ``supports_adjoint`` 和 ``supports_controlled`` 为键的布尔字典。
    """
    return capability_table(program)[key or program.entry]


def bind(program: Program | Operation, bindings: dict[str, Binding | Operation]) -> Program:
    """绑定已声明的槽；新增的 QRAM 资源沿模块图显式提升，其他槽可以继续开放。

    Args:
        program: 含开放声明槽的 ``Program`` 或 ``Operation``。
        bindings: 声明槽名到 ``Binding``（或裸 ``Operation``）的映射。

    Returns:
        Program: 绑定并提升捕获资源后重新校验的程序；未涉及的槽保持开放。
    """
    if isinstance(program, Operation):
        program = program.program()
    validate(program)
    modules = program.module_map.copy()
    global_types = {r.name: r.type for r in program.main.resources}
    captures = {}

    def reject(code: str, slot: str, expected: object, actual: object, message: str) -> NoReturn:
        """诊断记录从入口到槽位的调用路径。"""
        paths = {item.name: item.path for item in unresolved(program)}
        raise BindingError((BindingIssue(code, paths.get(slot, (program.entry, slot)), expected, actual, message),))

    def add(module: Module) -> None:
        """把实现模块并入模块表，同名且不同定义时抛错。"""
        existing = modules.get(module.name)
        if existing is not None and existing != module:
            raise ValidationError(f"绑定实现的模块名冲突：{module.name}")
        modules[module.name] = module

    for slot, item in bindings.items():
        if slot not in modules or modules[slot].body is not None:
            reject("BIND_TARGET", slot, "open declaration", slot, f"绑定目标不是开放声明：{slot}")
        declaration = modules[slot]
        binding = item if isinstance(item, Binding) else Binding(item)
        implementation = binding.operation
        implementation.program()
        target = implementation.module
        if target.name == slot:
            reject("BIND_NAME", slot, "distinct implementation name", target.name, "实现必须使用不同于声明槽的模块名")
        if tuple(r.type for r in target.registers) != tuple(r.type for r in declaration.registers):
            reject("BIND_SIGNATURE", slot, [asdict(r) for r in declaration.registers], [asdict(r) for r in target.registers], f"{slot} 的寄存器类型/宽度不匹配；形状变化请重新生成算法")
        expected, offered = dict(declaration.attributes), dict(target.attributes)
        for key in ("be_alpha", "fixed_width", "fixed_fraction", "fixed_signed", "rounding"):
            if key in expected and offered.get(key) != expected[key]:
                reject("BIND_ATTRIBUTE", slot, {key: expected[key]}, {key: offered.get(key)}, f"{slot} 的 {key} 不匹配；请按新的常量重新生成算法")
        for key in ("zero_input", "clean_work"):
            # 缺省的历史注解仍可绑定；明确矛盾的声明不能由包装模块覆盖。
            if expected.get(key) is True and offered.get(key) is False:
                reject("BIND_PROMISE", slot, {key: True}, {key: False}, f"{slot} 的 {key} 声明冲突")
        role = offered.get("oracle_paradigm")
        if (
            role
            and role != expected["oracle_paradigm"]
            and expected["oracle_paradigm"] != "unitary"
        ):
            reject("BIND_ROLE", slot, expected["oracle_paradigm"], role, f"{slot} 的 oracle paradigm 不匹配")
        provided_caps = capabilities(implementation.program())
        for capability in ("supports_adjoint", "supports_controlled"):
            if expected.get(capability, True) and not provided_caps[capability]:
                reject("BIND_CAPABILITY", slot, capability, False, f"{slot} 缺少要求的能力 {capability}")
        for module in (*implementation.dependencies, target):
            add(module)
        explicit = binding.resources or {}
        if set(explicit) - {r.name for r in target.resources}:
            raise ValidationError("绑定包含未知资源形式参数")
        formal_resources = {r.name: r.type for r in declaration.resources}
        resource_arguments = []
        # resource 先后承载 Resource 形参对象与捕获资源名字符串，按联合类型注解。
        resource: Resource | str
        for resource in target.resources:
            if resource.name not in explicit and resource.name in formal_resources:
                if formal_resources[resource.name] != resource.type:
                    raise ValidationError("声明和实现的资源类型不符")
                resource_arguments.append(resource.name)
            else:
                actual = explicit.get(resource.name, f"{slot}__{resource.name}")
                from oracq.infrastructure.validation import name

                name(actual)
                previous = global_types.get(actual)
                if previous is not None and previous != resource.type:
                    raise ValidationError(f"捕获资源类型冲突：{actual}")
                global_types[actual] = resource.type
                captures[actual] = resource.type
                resource_arguments.append("@" + actual)
        attributes = dict(declaration.attributes)
        attributes["oracle_bound_to"] = target.name
        attributes["implementation_status"] = "bound"
        arguments = tuple(
            Ref((Span(r.name, 0, r.type.width),), r.type) for r in declaration.registers
        )
        modules[slot] = replace(
            declaration,
            body=(Call(target.name, arguments, tuple(resource_arguments)),),
            attributes=tuple(sorted(attributes.items())),
        )

    requirements: dict[str, tuple[str, ...]]
    active: set[str]
    requirements, active = {}, set()

    def need(key: str) -> tuple[str, ...]:
        """收集模块及其被调链所需的全部捕获资源名，并检测循环调用。"""
        if key in active:
            raise ValidationError("绑定引入了循环调用")
        if key in requirements:
            return requirements[key]
        active.add(key)
        result: set[str] = set()
        for call in calls(modules[key].body):
            result.update(r[1:] for r in call.resources if r.startswith("@"))
            if call.module not in modules:
                raise ValidationError(f"绑定后存在未知模块：{call.module}")
            result.update(need(call.module))
        active.remove(key)
        requirements[key] = tuple(sorted(result))
        return requirements[key]

    for key in modules:
        need(key)
    local_names, additions, capture_maps = {}, {}, {}
    for key, module in modules.items():
        occupied = {r.name for r in module.registers + module.locals} | {
            r.name for r in module.resources
        }
        existing = {r.name: r.type for r in module.resources}
        previous_captures = capture_map(module)
        names, extra = {}, []
        for resource in requirements[key]:
            if resource in previous_captures:
                local = previous_captures[resource]
                if existing[local] != captures[resource]:
                    raise ValidationError("已捕获资源类型冲突")
                names[resource] = local
                continue
            if key == program.entry:
                local = resource
                if local in existing:
                    if existing[local] != captures[resource]:
                        raise ValidationError("入口资源类型冲突")
                    names[resource] = local
                    continue
                if local in occupied:
                    raise ValidationError("入口捕获资源与寄存器名字冲突")
            else:
                prefix = "capture_" + hashlib.sha256(resource.encode()).hexdigest()[:16]
                local = prefix
                while local in occupied:
                    local += "_"
            occupied.add(local)
            names[resource] = local
            extra.append((resource, Resource(local, captures[resource])))
        local_names[key] = names
        additions[key] = tuple(extra)
        capture_maps[key] = previous_captures | {logical: r.name for logical, r in extra}

    ordered_resources: dict[str, tuple[Resource, ...]] = {}
    permutations: dict[str, tuple[int, ...]] = {}
    for key, module in modules.items():
        unsorted_resources = module.resources + tuple(r for _, r in additions[key])
        captured_names = set(capture_maps[key].values())
        by_name = {r.name: r for r in unsorted_resources}
        ordered_resources[key] = (
            tuple(r for r in unsorted_resources if r.name not in captured_names)
            + tuple(by_name[local] for _, local in sorted(capture_maps[key].items()))
        )
        positions = {r.name: i for i, r in enumerate(unsorted_resources)}
        permutations[key] = tuple(positions[r.name] for r in ordered_resources[key])

    def rewrite(
        nodes: tuple[Instruction, ...] | None, owner: str
    ) -> tuple[Instruction, ...] | None:
        """改写体内调用的捕获资源绑定，并把提升资源追加到被调实参。"""
        if nodes is None:
            return None
        result: list[Instruction] = []
        for node in nodes:
            if isinstance(node, Call):
                values = tuple(
                    local_names[owner][r[1:]] if r.startswith("@") else r for r in node.resources
                )
                values += tuple(local_names[owner][key] for key, _ in additions[node.module])
                values = tuple(values[i] for i in permutations[node.module])
                result.append(replace(node, resources=values))
            elif isinstance(node, (Repeat, Control, Adjoint)):
                # 传入的 node.body 为具体元组时 rewrite 必返回元组（None 仅来自 None 入参）。
                result.append(
                    replace(node, body=cast("tuple[Instruction, ...]", rewrite(node.body, owner)))
                )
            else:
                result.append(node)
        return tuple(result)

    linked = tuple(
        replace(
            module,
            resources=ordered_resources[key],
            body=rewrite(module.body, key),
            attributes=tuple(sorted({
                **dict(module.attributes),
                **({"binding_captures": json.dumps(capture_maps[key], sort_keys=True)} if capture_maps[key] else {}),
            }.items())),
        )
        for key, module in sorted(modules.items())
    )
    return validate(Program(program.entry, linked, VERSION))
