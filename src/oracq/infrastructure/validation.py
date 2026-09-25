"Closed operation set, type, aliasing, resource and acyclic call checks for RIR."

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
"""Set of register kinds that are allowed."""
UNARY = {"h", "x", "y", "z", "s", "t", "rx", "ry", "rz", "phase"}
"""Set of allowed unary quantum primitive operation names."""
ROTATIONS = {"rx", "ry", "rz", "phase"}
"""Set of rotation primitive names that carry an angle argument; a subset of ``UNARY``."""
BINARY = {"xor", "swap"}
"""Set of allowed binary quantum primitive operation names."""
IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")
"""Regex pattern for valid identifiers: starts with a letter or underscore, followed by letters, digits or underscores."""


def require(condition: object, message: str) -> None:
    """Assert that a structural check condition holds, otherwise raise ``ValidationError``.

    Args:
        condition: check result interpreted with boolean semantics.
        message: error description written into the exception on failure.

    Raises:
        ValidationError: ``condition`` is false.
    """
    if not condition:
        raise ValidationError(message)


def integer(value: object, lo: int, hi: int, message: str) -> None:
    """Assert that ``value`` is an integer within the closed interval ``lo..hi``.

    Args:
        value: value to check.
        lo: minimum allowed value, inclusive.
        hi: maximum allowed value, inclusive.
        message: error description written into the exception on failure.

    Raises:
        ValidationError: ``value`` is not an ``int`` or lies outside the interval.
    """
    require(type(value) is int and lo <= value <= hi, message)


def name(value: object) -> None:
    """Assert that ``value`` is a valid identifier string matching the ``IDENTIFIER`` pattern.

    Args:
        value: identifier to check, such as a register name, resource name or module name.

    Raises:
        ValidationError: ``value`` is not a string or contains illegal characters.
    """
    require(
        isinstance(value, str) and IDENTIFIER.fullmatch(value) is not None, f"invalid identifier: {value!r}"
    )


def reg_type(dtype: RegType) -> None:
    """Assert that ``dtype`` is a register type with a known kind and width 0..64.

    Args:
        dtype: storage type object to check.

    Raises:
        ValidationError: ``dtype`` is not a ``RegType``, its kind is unknown, or its width is out of range.
    """
    require(isinstance(dtype, RegType) and dtype.kind in KINDS, "unknown register type")
    integer(dtype.width, 0, 64, "width of every register or view must be 0..64")


def locations(ref: Ref) -> tuple[tuple[str, int], ...]:
    """Expand all qubit coordinates covered by a view reference.

    Args:
        ref: register view to expand.

    Returns:
        tuple: tuple of ``(register name, bit index)`` pairs in view segment order.
    """
    return tuple((s.register, bit) for s in ref.parts for bit in range(s.start, s.start + s.width))


def capture_map(module: Module) -> dict[str, str]:
    """Parse the linker resource-origin attribute; ensures the referenced mapping exists and is one-to-one."""
    raw = dict(module.attributes).get("binding_captures", "{}")
    require(isinstance(raw, str), "binding_captures must be a JSON string")
    try:
        result = json.loads(cast(str, raw))
    except ValueError as exc:
        raise ValidationError("binding_captures is not valid JSON") from exc
    require(isinstance(result, dict), "binding_captures must be an object")
    existing = {r.name for r in module.resources}
    for logical, local in result.items():
        name(logical)
        name(local)
        require(local in existing, "binding_captures references a nonexistent resource")
    require(len(set(result.values())) == len(result), "binding_captures has duplicate local resource names")
    return cast(dict[str, str], result)


def _validate(program: Program) -> Program:
    """Run all cross-node structural checks on the program and return it unchanged on success."""
    require(isinstance(program, Program) and type(program.modules) is tuple, "an immutable Program is required")
    require(program.version == VERSION, f"unsupported RIR version: {program.version}")
    name(program.entry)
    modules = program.module_map
    require(len(modules) == len(program.modules), "duplicate module name")
    require(program.entry in modules, "entry module does not exist")
    graph: dict[str, set[str]] = {key: set() for key in modules}

    for module in program.modules:
        name(module.name)
        require(
            type(module.registers) is tuple and type(module.resources) is tuple,
            "module signature must be immutable",
        )
        require(type(module.attributes) is tuple, "module attributes must be immutable")
        require(
            all(type(pair) is tuple and len(pair) == 2 for pair in module.attributes),
            "attribute entries must be pairs",
        )
        attributes = dict(module.attributes)
        require(len(attributes) == len(module.attributes), "duplicate module attribute key")
        for key, value in attributes.items():
            name(key)
            require(type(value) in (str, int, float, bool), "attributes must be scalars")
            if type(value) is float:
                require(math.isfinite(value), "attribute floats must be finite")
        for capability in ("supports_adjoint", "supports_controlled"):
            if capability in attributes:
                require(type(attributes[capability]) is bool, "capability fields must be booleans")
        require(type(module.locals) is tuple, "local register list must be immutable")
        if module.body is None:
            require(not module.locals, "open declarations cannot define a private workspace")
        all_registers = module.registers + module.locals
        regs = {reg.name: reg.type for reg in all_registers}
        resources = {res.name: res.type for res in module.resources}
        require(len(regs) == len(all_registers), "duplicate register name")
        require(len(resources) == len(module.resources), "duplicate resource name")
        require(not regs.keys() & resources.keys(), "register and resource names conflict")
        for reg in all_registers:
            name(reg.name)
            reg_type(reg.type)
        require(sum(reg.type.width for reg in module.registers) > 0, "module must have a nonempty quantum interface")
        for resource in module.resources:
            name(resource.name)
            integer(resource.type.address_width, 1, 64, "QRAM address width must be 1..64")
            integer(resource.type.data_width, 1, 64, "QRAM data width must be 1..64")
        if "binding_captures" in attributes:
            capture_map(module)

        def check_ref(ref: Ref, regs: dict[str, RegType] = regs) -> set[tuple[str, int]]:
            """Validate a single view for legality and non-overlap, returning the set of qubits it covers."""
            require(isinstance(ref, Ref), "a register view is required")
            reg_type(ref.type)
            require(type(ref.parts) is tuple, "view must be immutable")
            require(sum(s.width for s in ref.parts) == ref.width, "view width does not match its segments")
            for span in ref.parts:
                require(span.register in regs, f"unknown register: {span.register}")
                width = regs[span.register].width
                integer(span.start, 0, width, "view start out of range")
                integer(span.width, 0, width - span.start, "view end out of range")
            locs = locations(ref)
            require(len(set(locs)) == len(locs), "overlapping qubits within a view")
            return set(locs)

        def distinct(refs: tuple[Ref, ...], protected: frozenset[tuple[str, int]]) -> None:
            """Assert that operands do not overlap each other and do not modify protected control qubits."""
            used: set[tuple[str, int]] = set()
            for ref in refs:
                current = check_ref(ref)
                require(not current & used, "operands alias or overlap")
                require(not current & protected, "operands modify a protected control register")
                used |= current

        def body(
            nodes: tuple[Instruction, ...],
            protected: frozenset[tuple[str, int]] = frozenset(),
            depth: int = 0,
            resources: dict[str, QRAM] = resources,
            module: Module = module,
            unitary: bool = True,
        ) -> None:
            """Recursively validate an instruction body: primitive arity, aliasing, control protection and call matching."""
            require(depth < 128, "nesting depth exceeds 127")
            require(type(nodes) is tuple, "instruction body must be immutable")
            for node in nodes:
                if isinstance(node, Primitive):
                    require(type(node.operands) is tuple, "primitive operands must be immutable")
                    require(node.op in UNARY | BINARY | {"add_const", "gphase"}, "unknown primitive")
                    expected = 2 if node.op in BINARY else 0 if node.op == "gphase" else 1
                    require(len(node.operands) == expected, "primitive operand count mismatch")
                    distinct(node.operands, protected)
                    if node.op in BINARY:
                        require(
                            node.operands[0].width == node.operands[1].width, "binary operation width mismatch"
                        )
                    if node.op in ROTATIONS | {"gphase"}:
                        require(
                            # short-circuiting of and guarantees angle is already of a finite real type when isfinite runs.
                            type(node.angle) in (float, int) and math.isfinite(cast(float, node.angle)),
                            "rotation angle must be a finite real number",
                        )
                    else:
                        require(node.angle is None, "this primitive does not accept an angle")
                    if node.op == "add_const":
                        require(node.operands[0].type.kind == "uint", "add_const requires uint")
                        integer(
                            node.value,
                            0,
                            (1 << node.operands[0].width) - 1,
                            "addition constant exceeds the register range",
                        )
                    else:
                        require(node.value is None, "this primitive does not accept an integer argument")
                elif isinstance(node, Load):
                    require(node.resource in resources, "QRAM resource is not declared")
                    spec = resources[node.resource]
                    distinct((node.address, node.data), protected)
                    require(node.address.width == spec.address_width, "QRAM address width mismatch")
                    require(node.data.width == spec.data_width, "QRAM data width mismatch")
                elif isinstance(node, Store):
                    require(unitary, "Store is a non-unitary side effect and cannot appear inside a Control or Adjoint body")
                    require(node.resource in resources, "QRAM resource is not declared")
                    spec = resources[node.resource]
                    distinct((node.address, node.data), protected)
                    require(node.address.width == spec.address_width, "QRAM address width mismatch")
                    require(node.data.width == spec.data_width, "QRAM data width mismatch")
                elif isinstance(node, Call):
                    require(
                        type(node.arguments) is tuple and type(node.resources) is tuple,
                        "call arguments must be immutable",
                    )
                    require(node.module in modules, f"unknown module: {node.module}")
                    target = modules[node.module]
                    require(len(node.arguments) == len(target.registers), "module quantum argument count mismatch")
                    require(len(node.resources) == len(target.resources), "module resource argument count mismatch")
                    distinct(node.arguments, protected)
                    for actual, formal in zip(node.arguments, target.registers, strict=True):
                        require(actual.type == formal.type, f"module argument type mismatch: {formal.name}")
                    # actual/formal are bound to Ref/Register in the loop above and carry str/Resource in this one.
                    for actual, formal in zip(  # type: ignore[assignment]
                        node.resources, target.resources, strict=True
                    ):
                        require(
                            actual in resources and resources[cast(str, actual)] == formal.type,
                            "module QRAM argument type mismatch",
                        )
                    graph[module.name].add(node.module)
                elif isinstance(node, Repeat):
                    integer(node.count, 0, 2**63 - 1, "repeat count must be 0..2^63-1")
                    body(node.body, protected, depth + 1, unitary=unitary)
                elif isinstance(node, Control):
                    locs = check_ref(node.register)
                    require(bool(locs), "control register cannot be empty")
                    integer(node.value, 0, (1 << node.register.width) - 1, "control value out of range")
                    require(not locs & protected, "nested control registers overlap")
                    body(node.body, protected | locs, depth + 1, unitary=False)
                elif isinstance(node, Adjoint):
                    body(node.body, protected, depth + 1, unitary=False)
                else:
                    raise ValidationError(f"unknown instruction type: {type(node).__name__}")

        if module.body is None:
            require(
                isinstance(attributes.get("oracle_paradigm"), str),
                "open declarations must specify oracle_paradigm",
            )
        else:
            body(module.body)

    visited: set[str]
    active: set[str]
    visited, active = set(), set()

    def visit(key: str) -> None:
        """Detect recursion and overly deep module calls by walking the call graph depth-first."""
        require(key not in active, "module call graph contains recursion")
        require(len(active) < 128, "module call depth exceeds 127")
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
        """Check that modules called in a controlled or adjoint context provide the corresponding capabilities."""
        for node in nodes or ():
            if isinstance(node, Call):
                if controlled:
                    require(
                        inferred[node.module]["supports_controlled"],
                        f"{node.module} lacks supports_controlled",
                    )
                if inverse:
                    require(
                        inferred[node.module]["supports_adjoint"],
                        f"{node.module} lacks supports_adjoint",
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
    """Run full structural validation on an RIR program and return it unchanged.

    Covers cross-node rules such as module signatures and attributes, register
    and view rules, primitive arity and aliasing constraints, control register
    protection, QRAM resource matching and an acyclic call graph, and checks
    that modules called in a controlled or adjoint context provide the
    corresponding capabilities; the rules stay consistent with
    docs/reference/rir.md.

    Args:
        program: ``Program`` to validate.
        require_closed: when true, additionally require the program to contain no unbound open oracle declarations.

    Returns:
        Program: the same ``program`` object that passed validation.

    Raises:
        ValidationError: the program violates a structural or semantic constraint; or a standard
            exception indicating that the input is not a valid RIR object was caught.
    """
    try:
        result = _validate(program)
        if require_closed:
            from oracq.infrastructure.linking import unresolved

            missing = unresolved(result)
            if missing:
                raise ValidationError(
                    "unbound oracle: "
                    + "; ".join(f"{r.name} {r.paradigm}; {' -> '.join(r.path)}" for r in missing)
                )
        return result
    except ValidationError:
        raise
    except (AttributeError, TypeError, KeyError, ValueError, OverflowError, RecursionError) as exc:
        raise ValidationError(f"invalid RIR object: {exc}") from exc
