"开放 oracle 的缺口分析、部分绑定和 QRAM 捕获资源提升。"

from __future__ import annotations

import hashlib
from collections import deque
from dataclasses import dataclass, replace

from pyqecclang.infrastructure.builder import Operation
from pyqecclang.infrastructure.ir import (
    VERSION,
    Adjoint,
    Call,
    Control,
    Program,
    Ref,
    Repeat,
    Resource,
    Span,
    Store,
    ValidationError,
)
from pyqecclang.infrastructure.validation import validate


@dataclass(frozen=True)
class OracleRequirement:
    name: str
    paradigm: str
    path: tuple[str, ...]
    registers: tuple
    attributes: tuple


@dataclass(frozen=True)
class Binding:
    operation: Operation
    resources: dict[str, str] | None = None


def calls(nodes):
    for node in nodes or ():
        if isinstance(node, Call):
            yield node
        elif isinstance(node, (Repeat, Control, Adjoint)):
            yield from calls(node.body)


def stores(nodes):
    """列出指令体（含嵌套结构块）中的全部 Store 副作用。"""
    for node in nodes or ():
        if isinstance(node, Store):
            yield node
        elif isinstance(node, (Repeat, Control, Adjoint)):
            yield from stores(node.body)


def uses_store(program: Program) -> bool:
    """入口可达的模块中是否存在 QRAM 随机写。"""
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


def unresolved(program: Program) -> tuple[OracleRequirement, ...]:
    """列出入口结构可达的未实现声明，并给出首条最短调用路径。"""
    validate(program)
    modules = program.module_map
    pending = deque([(program.entry, (program.entry,))])
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
                    dict(module.attributes)["oracle_paradigm"],
                    path,
                    module.registers,
                    module.attributes,
                )
            )
        else:
            for call in calls(module.body):
                pending.append((call.module, path + (call.module,)))
    return tuple(sorted(result, key=lambda item: item.name))


def capability_table(program: Program):
    """一次遍历推导全部模块的变换能力。"""
    modules, cache = program.module_map, {}

    def infer(name):
        if name in cache:
            return cache[name]
        module = modules[name]
        attrs = dict(module.attributes)
        result = {cap: attrs.get(cap, True) for cap in ("supports_adjoint", "supports_controlled")}
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


def capabilities(program: Program, key: str | None = None):
    return capability_table(program)[key or program.entry]


def bind(program: Program | Operation, bindings: dict[str, Binding | Operation]) -> Program:
    """绑定已声明的槽；新增的 QRAM 资源沿模块图显式提升，其他槽可以继续开放。"""
    if isinstance(program, Operation):
        program = program.program()
    validate(program)
    modules = program.module_map.copy()
    global_types = {r.name: r.type for r in program.main.resources}
    captures = {}

    def add(module):
        existing = modules.get(module.name)
        if existing is not None and existing != module:
            raise ValidationError(f"绑定实现的模块名冲突：{module.name}")
        modules[module.name] = module

    for slot, item in bindings.items():
        if slot not in modules or modules[slot].body is not None:
            raise ValidationError(f"绑定目标不是开放声明：{slot}")
        declaration = modules[slot]
        binding = item if isinstance(item, Binding) else Binding(item)
        implementation = binding.operation
        implementation.program()
        target = implementation.module
        if target.name == slot:
            raise ValidationError("实现必须使用不同于声明槽的模块名")
        if tuple(r.type for r in target.registers) != tuple(r.type for r in declaration.registers):
            raise ValidationError(f"{slot} 的寄存器类型/宽度不匹配；形状变化请重新生成算法")
        expected, offered = dict(declaration.attributes), dict(target.attributes)
        if "be_alpha" in expected and offered.get("be_alpha") != expected["be_alpha"]:
            raise ValidationError(f"{slot} 的 be_alpha 不匹配；请按新的常量重新生成算法")
        role = offered.get("oracle_paradigm")
        if (
            role
            and role != expected["oracle_paradigm"]
            and expected["oracle_paradigm"] != "unitary"
        ):
            raise ValidationError(f"{slot} 的 oracle paradigm 不匹配")
        provided_caps = capabilities(implementation.program())
        for capability in ("supports_adjoint", "supports_controlled"):
            if expected.get(capability, True) and not provided_caps[capability]:
                raise ValidationError(f"{slot} 缺少要求的能力 {capability}")
        for module in (*implementation.dependencies, target):
            add(module)
        explicit = binding.resources or {}
        if set(explicit) - {r.name for r in target.resources}:
            raise ValidationError("绑定包含未知资源形式参数")
        formal_resources = {r.name: r.type for r in declaration.resources}
        resource_arguments = []
        for resource in target.resources:
            if resource.name not in explicit and resource.name in formal_resources:
                if formal_resources[resource.name] != resource.type:
                    raise ValidationError("声明和实现的资源类型不符")
                resource_arguments.append(resource.name)
            else:
                actual = explicit.get(resource.name, f"{slot}__{resource.name}")
                from pyqecclang.infrastructure.validation import name

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

    requirements, active = {}, set()

    def need(key):
        if key in active:
            raise ValidationError("绑定引入了循环调用")
        if key in requirements:
            return requirements[key]
        active.add(key)
        result = set()
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
    local_names, additions = {}, {}
    for key, module in modules.items():
        occupied = {r.name for r in module.registers + module.locals} | {
            r.name for r in module.resources
        }
        existing = {r.name: r.type for r in module.resources}
        names, extra = {}, []
        for resource in requirements[key]:
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

    def rewrite(nodes, owner):
        if nodes is None:
            return None
        result = []
        for node in nodes:
            if isinstance(node, Call):
                values = tuple(
                    local_names[owner][r[1:]] if r.startswith("@") else r for r in node.resources
                )
                values += tuple(local_names[owner][key] for key, _ in additions[node.module])
                result.append(replace(node, resources=values))
            elif isinstance(node, (Repeat, Control, Adjoint)):
                result.append(replace(node, body=rewrite(node.body, owner)))
            else:
                result.append(node)
        return tuple(result)

    linked = tuple(
        replace(
            module,
            resources=module.resources + tuple(r for _, r in additions[key]),
            body=rewrite(module.body, key),
        )
        for key, module in sorted(modules.items())
    )
    return validate(Program(program.entry, linked, VERSION))
