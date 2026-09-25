"Gap analysis of open oracles, partial binding and QRAM captured resource promotion."

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
    """Describe an unimplemented open oracle declaration and its call sites.

    Attributes:
        name: module name of the open declaration, i.e. the slot name used at binding time.
        paradigm: ``oracle_paradigm`` paradigm designated by the declaration.
        path: sequence of module names along the shortest call path from the entry module to the declaration.
        registers: formal register signature tuple of the declaration.
        attributes: attribute entry tuple of the declaration, each a key-value pair.
    """
    name: str
    paradigm: str
    path: tuple[str, ...]
    registers: tuple
    attributes: tuple


@dataclass(frozen=True)
class Binding:
    """Bind an open declaration slot to an implementing operation, optionally carrying a QRAM resource mapping.

    Attributes:
        operation: ``Operation`` providing the implementation; its module name must differ from the declaration slot name.
        resources: mapping from the implementation's formal resource parameters to entry logical resource names; resources not
            explicitly covered are auto-captured as ``slot name__resource name`` and promoted to entry resources.
    """
    operation: Operation
    resources: dict[str, str] | None = None


@dataclass(frozen=True)
class BindingIssue:
    """Structured reason for a binding failure; makes no claim of proving the oracle's mathematical semantics."""

    code: str
    path: tuple[str, ...]
    expected: object
    actual: object
    message: str


class BindingError(ValidationError):
    """Binding diagnostic that remains compatible with ValidationError."""

    def __init__(self, issues: tuple[BindingIssue, ...]) -> None:
        self.issues = issues
        super().__init__("; ".join(f"{i.code} {' -> '.join(i.path)}: {i.message}" for i in issues))


@dataclass(frozen=True)
class BindingReport:
    """Binding manifest, program fingerprints and remaining dependencies; stores only exchangeable data."""

    source_digest: str
    result_digest: str | None
    bindings: tuple[tuple[str, str], ...]
    resource_mappings: tuple[tuple[str, tuple[tuple[str, str], ...]], ...]
    remaining: tuple[OracleRequirement, ...]
    resources: tuple[Resource, ...]
    issues: tuple[BindingIssue, ...] = ()

    @property
    def ok(self) -> bool:
        """Whether this binding passed the structural checks; does not require every slot to be closed."""
        return not self.issues

    def to_dict(self) -> dict[str, object]:
        """Emit a JSON-friendly report, saved independently of the executable RIR."""
        return {"ok": self.ok, **asdict(self)}


@dataclass(frozen=True)
class BindingResult:
    """Program and report of one binding, avoiding a second linking pass after checking."""

    program: Program | None
    report: BindingReport

    def require(self) -> Program:
        """Obtain the program when binding succeeded, otherwise raise the structured issues."""
        if self.report.issues:
            raise BindingError(self.report.issues)
        if self.program is None:
            raise ValidationError("binding result lacks a program")
        return self.program


def bind_with_report(
    program: Program | Operation, bindings: dict[str, Binding | Operation]
) -> BindingResult:
    """Perform one pure-functional binding, returning a directly usable program plus a traceable report.

    Invalid input programs are still rejected by validate; binding errors of candidate
    implementations are recorded in the report. Data contents are not part of the RIR
    fingerprint, and run records must identify the QRAM data snapshot separately.
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
    """Yield all module call nodes in an instruction body, including nested structure blocks, in order of appearance.

    Args:
        nodes: instruction body to scan; ``None`` or an empty tuple is treated as no content.

    Returns:
        Iterator[Call]: ``Call`` nodes lazily yielded in order of appearance, including inside nested blocks.
    """
    for node in nodes or ():
        if isinstance(node, Call):
            yield node
        elif isinstance(node, (Repeat, Control, Adjoint)):
            yield from calls(node.body)


def stores(nodes: tuple[Instruction, ...] | None) -> Iterator[Store]:
    """List all Store side effects in an instruction body, including nested structure blocks.

    Args:
        nodes: instruction body to scan; ``None`` or an empty tuple is treated as no content.

    Returns:
        Iterator[Store]: ``Store`` nodes lazily yielded in order of appearance, including inside nested blocks.
    """
    for node in nodes or ():
        if isinstance(node, Store):
            yield node
        elif isinstance(node, (Repeat, Control, Adjoint)):
            yield from stores(node.body)


def uses_store(program: Program) -> bool:
    """Whether a QRAM random write exists in any entry-reachable module.

    Args:
        program: RIR program to check, traversed from the entry along the call graph.

    Returns:
        bool: True when any entry-reachable module body contains a ``Store`` node.
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
    """List structurally reachable unimplemented declarations from the entry, with the first shortest call path for each.

    Args:
        program: RIR program to analyze; structural validation runs first.

    Returns:
        tuple[OracleRequirement, ...]: list of open declaration requirements sorted
        by declaration name, each recording the first shortest call path from the
        entry.
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
                    # Validation guarantees the oracle_paradigm attribute of an open declaration is a str.
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
    """Derive the transformation capabilities of all modules in one traversal.

    Args:
        program: RIR program to analyze, traversing all its modules.

    Returns:
        dict[str, dict[str, bool]]: mapping from module names to capability
        dicts; modules containing ``Store`` side effects have all capabilities
        False.
    """
    cache: dict[str, dict[str, bool]]
    modules, cache = program.module_map, {}

    def infer(name: str) -> dict[str, bool]:
        """Derive one module's capability dict by conservatively conjoining its declaration with callee capabilities."""
        if name in cache:
            return cache[name]
        module = modules[name]
        attrs = dict(module.attributes)
        # Validation guarantees supports_adjoint/supports_controlled attributes are booleans whenever present.
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
    """Return the transformation capability dict of a given module, the entry by default.

    Capability values are derived by conservatively conjoining the module
    declaration with the capabilities of called modules; see
    ``capability_table``.

    Args:
        program: ``Program`` to analyze.
        key: target module name; the entry module is used when omitted.

    Returns:
        dict: boolean dict keyed by ``supports_adjoint`` and ``supports_controlled``.
    """
    return capability_table(program)[key or program.entry]


def bind(program: Program | Operation, bindings: dict[str, Binding | Operation]) -> Program:
    """Bind declared slots; newly added QRAM resources are promoted explicitly along the module graph, and other slots may stay open.

    Args:
        program: ``Program`` or ``Operation`` containing open declaration slots.
        bindings: mapping from declaration slot names to ``Binding`` or a bare ``Operation``.

    Returns:
        Program: revalidated program after binding and promoting captured resources; untouched slots stay open.
    """
    if isinstance(program, Operation):
        program = program.program()
    validate(program)
    modules = program.module_map.copy()
    global_types = {r.name: r.type for r in program.main.resources}
    captures = {}

    def reject(code: str, slot: str, expected: object, actual: object, message: str) -> NoReturn:
        """The diagnostic records the call path from the entry to the slot."""
        paths = {item.name: item.path for item in unresolved(program)}
        raise BindingError((BindingIssue(code, paths.get(slot, (program.entry, slot)), expected, actual, message),))

    def add(module: Module) -> None:
        """Merge an implementation module into the module table; raises on equal names with differing definitions."""
        existing = modules.get(module.name)
        if existing is not None and existing != module:
            raise ValidationError(f"module name conflict among bound implementations: {module.name}")
        modules[module.name] = module

    for slot, item in bindings.items():
        if slot not in modules or modules[slot].body is not None:
            reject("BIND_TARGET", slot, "open declaration", slot, f"binding target is not an open declaration: {slot}")
        declaration = modules[slot]
        binding = item if isinstance(item, Binding) else Binding(item)
        implementation = binding.operation
        implementation.program()
        target = implementation.module
        if target.name == slot:
            reject("BIND_NAME", slot, "distinct implementation name", target.name, "the implementation must use a module name different from the declaration slot")
        if tuple(r.type for r in target.registers) != tuple(r.type for r in declaration.registers):
            reject("BIND_SIGNATURE", slot, [asdict(r) for r in declaration.registers], [asdict(r) for r in target.registers], f"register type or width mismatch for {slot}; regenerate the algorithm for shape changes")
        expected, offered = dict(declaration.attributes), dict(target.attributes)
        for key in ("be_alpha", "fixed_width", "fixed_fraction", "fixed_signed", "rounding"):
            if key in expected and offered.get(key) != expected[key]:
                reject("BIND_ATTRIBUTE", slot, {key: expected[key]}, {key: offered.get(key)}, f"{key} mismatch for {slot}; regenerate the algorithm with the new constants")
        for key in ("zero_input", "clean_work"):
            # Historic annotations left at their default are still bindable; explicitly contradictory declarations cannot be overridden by a wrapper module.
            if expected.get(key) is True and offered.get(key) is False:
                reject("BIND_PROMISE", slot, {key: True}, {key: False}, f"conflicting {key} declarations for {slot}")
        role = offered.get("oracle_paradigm")
        if (
            role
            and role != expected["oracle_paradigm"]
            and expected["oracle_paradigm"] != "unitary"
        ):
            reject("BIND_ROLE", slot, expected["oracle_paradigm"], role, f"oracle paradigm mismatch for {slot}")
        provided_caps = capabilities(implementation.program())
        for capability in ("supports_adjoint", "supports_controlled"):
            if expected.get(capability, True) and not provided_caps[capability]:
                reject("BIND_CAPABILITY", slot, capability, False, f"{slot} lacks the required capability {capability}")
        for module in (*implementation.dependencies, target):
            add(module)
        explicit = binding.resources or {}
        if set(explicit) - {r.name for r in target.resources}:
            raise ValidationError("binding contains an unknown resource formal parameter")
        formal_resources = {r.name: r.type for r in declaration.resources}
        resource_arguments = []
        # resource first carries Resource formal objects then captured resource name strings; annotated as a union type.
        resource: Resource | str
        for resource in target.resources:
            if resource.name not in explicit and resource.name in formal_resources:
                if formal_resources[resource.name] != resource.type:
                    raise ValidationError("resource type mismatch between declaration and implementation")
                resource_arguments.append(resource.name)
            else:
                actual = explicit.get(resource.name, f"{slot}__{resource.name}")
                from oracq.infrastructure.validation import name

                name(actual)
                previous = global_types.get(actual)
                if previous is not None and previous != resource.type:
                    raise ValidationError(f"captured resource type conflict: {actual}")
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
        """Collect all captured resource names needed by a module and its callee chain, and detect cyclic calls."""
        if key in active:
            raise ValidationError("binding introduces a cyclic call")
        if key in requirements:
            return requirements[key]
        active.add(key)
        result: set[str] = set()
        for call in calls(modules[key].body):
            result.update(r[1:] for r in call.resources if r.startswith("@"))
            if call.module not in modules:
                raise ValidationError(f"unknown module after binding: {call.module}")
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
                    raise ValidationError("captured resource type conflict")
                names[resource] = local
                continue
            if key == program.entry:
                local = resource
                if local in existing:
                    if existing[local] != captures[resource]:
                        raise ValidationError("entry resource type conflict")
                    names[resource] = local
                    continue
                if local in occupied:
                    raise ValidationError("entry captured resource conflicts with a register name")
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
        """Rewrite captured resource bindings of calls in the body and append promoted resources to callee arguments."""
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
                # When the incoming node.body is a concrete tuple, rewrite must return a tuple (None only comes from a None input).
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
