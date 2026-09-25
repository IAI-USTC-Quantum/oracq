"On-demand traversal of the module graph, plus a small-scale reference executor based on tuples of register integers."

from __future__ import annotations

import cmath
import math
from collections.abc import Callable, Collection, Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

from oracq.infrastructure.ir import (
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
    Span,
    Store,
    ValidationError,
)
from oracq.infrastructure.validation import locations, validate

if TYPE_CHECKING:
    from oracq.infrastructure.native import NativeSite


def check_memory(
    program: Program,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None,
) -> dict[str, dict[int, int]]:
    """Validate and normalize the QRAM binding data of the entry.

    Args:
        program: closed RIR program.
        memory: data provided by resource name; each item is a word sequence or a sparse ``address -> word`` dict, and ``None`` is treated as all zeros.

    Returns:
        dict: dictionary from resource names to sparse cell maps; cells with value zero are dropped, and missing cells are interpreted as zero.

    Raises:
        ValidationError: the data is not a mapping by resource name, the resource set does not match the entry declaration, an address is out of range, or a word is not an unsigned integer within the bit width.
    """
    memory = {} if memory is None else memory
    if not isinstance(memory, Mapping):
        raise ValidationError("QRAM data must be provided as a mapping keyed by resource name")
    specs = {r.name: r.type for r in program.main.resources}
    if set(memory) != set(specs):
        raise ValidationError("data must be provided for every QRAM of the entry with no extra resources")
    result = {}
    for name, spec in specs.items():
        raw = memory[name]
        try:
            items = raw.items() if isinstance(raw, Mapping) else enumerate(raw)
        except TypeError as exc:
            raise ValidationError("each QRAM must provide a word sequence or a sparse dict") from exc
        cells = {}
        for address, value in items:
            if type(address) is not int or not 0 <= address < 1 << spec.address_width:
                raise ValidationError("QRAM address out of range")
            if type(value) is not int or not 0 <= value < 1 << spec.data_width:
                raise ValidationError("QRAM words must be unsigned integers within the bit width")
            if value:
                cells[address] = value
        result[name] = cells
    return result


def _counter(
    program: Program, limit: int, native_modules: Collection[str] = frozenset()
) -> Callable[[tuple[Instruction, ...]], int]:
    """Build a recursive step-counting function memoized by module, reused by ``events`` and budget checks."""
    modules, cached = program.module_map, {}

    def count(nodes: tuple[Instruction, ...]) -> int:
        """Accumulate steps node by node under the weighting rules, clamping once the total reaches ``limit + 1``."""
        total = 0
        for node in nodes:
            if isinstance(node, Primitive):
                cost = max(4, sum(r.width for r in node.operands) ** 2)
            elif isinstance(node, (Load, Store)):
                cost = 1
            elif isinstance(node, Call):
                if node.module in native_modules:
                    cost = 1
                    total = min(limit + 1, total + cost)
                    continue
                if node.module not in cached:
                    cached[node.module] = count(
                        cast("tuple[Instruction, ...]", modules[node.module].body)
                    )
                cost = 1 + cached[node.module]
            elif isinstance(node, Repeat):
                cost = node.count * (1 + count(node.body)) if node.body else 0
            elif isinstance(node, Control):
                cost = 1 + (1 + 2 * node.register.width) * count(node.body)
            else:
                cost = 1 + count(node.body)
            total = min(limit + 1, total + cost)
        return total

    return count


def expanded_steps(
    program: Program,
    limit: int = 1_000_000,
    native_modules: Collection[str] = frozenset(),
) -> int:
    """Estimate the execution step count of the fully expanded program, for budget checks before execution.

    Primitives are weighted by total operand bit width and Load/Store count as 1 step; module calls recursively include the callee body, and Repeat multiplies by the count.
    Once the count reaches the limit it is clamped at ``limit + 1`` and stops growing.

    Args:
        program: RIR program.
        limit: clamping limit of the count.
        native_modules: set of native module names counted as a single step; returns 1 directly when the entry itself is a native module.

    Returns:
        int: saturated estimate of the expanded step count.
    """
    if program.entry in native_modules:
        return 1
    return _counter(program, limit, native_modules)(
        cast("tuple[Instruction, ...]", program.main.body)
    )


@dataclass(frozen=True)
class LocalEnter:
    """Event of a local register entering scope; the register joins the execution state at value zero.

    Attributes:
        name: unique synthesized name assigned to the local register.
        type: type of the local register.
    """

    name: str
    type: RegType


@dataclass(frozen=True)
class LocalExit:
    """Event of a local register leaving scope; at exit all branches are required to be restored to zero.

    Attributes:
        name: name of the released local register.
        width: bit width of the local register.
    """

    name: str
    width: int


def remap(ref: Ref, mapping: Mapping[str, Ref]) -> Ref:
    """Rewrite a reference by a register-name substitution table, used at module calls to bind formal parameters to argument views.

    Args:
        ref: ``Ref`` to rewrite.
        mapping: mapping from register names to ``Ref``; each span of the original reference is sliced from the corresponding range of the replacement and concatenated in order.

    Returns:
        Ref: rewritten new reference, with the same type as the original.
    """
    parts: list[Span] = []
    for span in ref.parts:
        parts.extend(mapping[span.register][span.start : span.start + span.width].parts)
    return Ref(tuple(parts), ref.type)


def events(
    program: Program,
    *,
    max_steps: int = 1_000_000,
    native_modules: Collection[str] = frozenset(),
) -> Iterator[
    tuple[
        Primitive | Load | Store | LocalEnter | LocalExit | NativeSite,
        tuple[tuple[Ref, int], ...],
        bool,
    ]
]:
    """Expand calls layer by layer into an iterator while leaving the original RIR unchanged.

    Args:
        program: RIR program to expand for execution; must be closed when there are no native modules.
        max_steps: budget limit on the expanded instruction step count.
        native_modules: set of module names treated as natively executed and no longer inlined.

    Returns:
        Iterator: lazily yields ``(event, control chain, inverse flag)`` triples; the event is a
        primitive, a QRAM read/write, a local register entering or leaving, or a native call
        site, the control chain is a tuple of ``(Ref, int)`` pairs, and the inverse flag
        indicates an ``Adjoint`` context.
    """
    from oracq.infrastructure.ir import Span

    program = validate(program, require_closed=not bool(native_modules))
    counter = _counter(program, max_steps, native_modules)
    if expanded_steps(program, max_steps, native_modules) > max_steps:
        raise ValidationError("execution exceeds the expansion budget")
    modules = program.module_map

    def walk(
        nodes: tuple[Instruction, ...],
        mapping: Mapping[str, Ref],
        resources: Mapping[str, str],
        controls: tuple[tuple[Ref, int], ...] = (),
        inverse: bool = False,
    ) -> Iterator[
        tuple[
            Primitive | Load | Store | LocalEnter | LocalExit | NativeSite,
            tuple[tuple[Ref, int], ...],
            bool,
        ]
    ]:
        """Rewrite instruction nodes in order and lazily yield event triples, with the control chain and inverse flag passed along the recursion."""
        for node in reversed(nodes) if inverse else nodes:
            if isinstance(node, Call):
                target = modules[node.module]
                actuals = {
                    r.name: remap(a, mapping)
                    for r, a in zip(target.registers, node.arguments, strict=True)
                }
                bound = {
                    r.name: resources[a]
                    for r, a in zip(target.resources, node.resources, strict=True)
                }
                yield from enter_module(target, actuals, bound, controls, inverse)
            elif isinstance(node, Repeat):
                if counter(node.body):
                    for _ in range(node.count):
                        yield from walk(node.body, mapping, resources, controls, inverse)
            elif isinstance(node, Control):
                yield from walk(
                    node.body,
                    mapping,
                    resources,
                    controls + ((remap(node.register, mapping), node.value),),
                    inverse,
                )
            elif isinstance(node, Adjoint):
                yield from walk(node.body, mapping, resources, controls, not inverse)
            elif isinstance(node, Primitive):
                yield (
                    Primitive(
                        node.op,
                        tuple(remap(r, mapping) for r in node.operands),
                        node.angle,
                        node.value,
                    ),
                    controls,
                    inverse,
                )
            elif isinstance(node, Load):
                yield (
                    Load(
                        resources[node.resource],
                        remap(node.address, mapping),
                        remap(node.data, mapping),
                    ),
                    controls,
                    inverse,
                )
            elif isinstance(node, Store):
                yield (
                    Store(
                        resources[node.resource],
                        remap(node.address, mapping),
                        remap(node.data, mapping),
                    ),
                    controls,
                    inverse,
                )

    local_counter = 0

    def enter_module(
        module: Module,
        mapping: Mapping[str, Ref],
        resources: Mapping[str, str],
        controls: tuple[tuple[Ref, int], ...] = (),
        inverse: bool = False,
    ) -> Iterator[
        tuple[
            Primitive | Load | Store | LocalEnter | LocalExit | NativeSite,
            tuple[tuple[Ref, int], ...],
            bool,
        ]
    ]:
        """Enter a single module: bind the formal parameter mapping, synthesize local register names on demand and forward the body event stream."""
        nonlocal local_counter
        if module.name in native_modules:
            from oracq.infrastructure.native import NativeSite

            yield (
                NativeSite(
                    module,
                    tuple(mapping[r.name] for r in module.registers),
                    tuple(resources[r.name] for r in module.resources),
                ),
                controls,
                inverse,
            )
            return
        if module.body is None:
            raise ValidationError(f"unimplemented module: {module.name}")
        mapping = dict(mapping)
        allocated = []
        for register in module.locals:
            local_counter += 1
            name = f"__local_{local_counter}"
            mapping[register.name] = Ref((Span(name, 0, register.type.width),), register.type)
            if register.type.width:
                allocated.append((name, register.type.width))
                yield LocalEnter(name, register.type), (), False
        yield from walk(module.body, mapping, resources, controls, inverse)
        for name, width in reversed(allocated):
            yield LocalExit(name, width), (), False

    roots = {r.name: Ref((Span(r.name, 0, r.type.width),), r.type) for r in program.main.registers}
    resources = {r.name: r.name for r in program.main.resources}
    yield from enter_module(program.main, roots, resources)


def gate_matrix(
    op: str, angle: float | None = None, inverse: bool = False
) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
    """Return the standard 2x2 unitary matrix of a single-qubit gate.

    Args:
        op: gate name; one of ``h``, ``x``, ``y``, ``z``, ``s``, ``t``, ``phase``, ``rx``, ``ry``, ``rz``.
        angle: angle in radians for rotation and phase gates; ``s`` and ``t`` use fixed angles and ignore this argument.
        inverse: when true, return the conjugate transpose, i.e. the inverse matrix of the gate.

    Returns:
        tuple: 2x2 complex matrix given as a tuple of rows.

    Raises:
        ValidationError: the gate name is not in the supported list.
    """
    if op == "h":
        a = 1 / math.sqrt(2)
        matrix: tuple[tuple[complex, complex], tuple[complex, complex]] = ((a, a), (a, -a))
    elif op == "x":
        matrix = ((0, 1), (1, 0))
    elif op == "y":
        matrix = ((0, -1j), (1j, 0))
    elif op == "z":
        matrix = ((1, 0), (0, -1))
    elif op in {"s", "t", "phase"}:
        theta = {"s": math.pi / 2, "t": math.pi / 4}.get(op, angle)
        matrix = ((1, 0), (0, cmath.exp(1j * cast(float, theta))))
    elif op == "rz":
        matrix = (
            (cmath.exp(-0.5j * cast(float, angle)), 0),
            (0, cmath.exp(0.5j * cast(float, angle))),
        )
    elif op in {"rx", "ry"}:
        c, s = math.cos(cast(float, angle) / 2), math.sin(cast(float, angle) / 2)
        matrix = ((c, -1j * s), (-1j * s, c)) if op == "rx" else ((c, -s), (s, c))
    else:
        raise ValidationError(f"unknown single-qubit matrix: {op}")
    if inverse:
        return cast(
            "tuple[tuple[complex, complex], tuple[complex, complex]]",
            tuple(tuple(complex(matrix[j][i]).conjugate() for j in range(2)) for i in range(2)),
        )
    return matrix


@dataclass(frozen=True)
class RegisterState:
    """Final state of a reference execution: entry register layout and sparse amplitudes.

    Attributes:
        registers: register tuple of the entry module, fixing the order of value tuples.
        amplitudes: mapping from register value tuples to complex amplitudes.
    """

    registers: tuple
    amplitudes: dict[tuple[int, ...], complex]

    def statevector(self, max_qubits: int = 20) -> list[complex]:
        """Pack the sparse state into a dense state vector, with each register occupying an index bit segment LSB-first in declaration order.

        Args:
            max_qubits: maximum total bit count allowed.

        Returns:
            list: list of complex amplitudes of length ``2**total bits``.

        Raises:
            ValidationError: the total register bit count exceeds ``max_qubits``.
        """
        width = sum(r.type.width for r in self.registers)
        if width > max_qubits:
            raise ValidationError("dense state conversion exceeds the bit budget")
        vector = [0j] * (1 << width)
        for values, amplitude in self.amplitudes.items():
            index, offset = 0, 0
            for reg, value in zip(self.registers, values, strict=True):
                index |= value << offset
                offset += reg.type.width
            vector[index] = amplitude
        return vector


def simulate(
    program: Program,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None = None,
    *,
    initial: Mapping[str, int] | None = None,
    max_steps: int = 1_000_000,
    max_states: int = 65536,
) -> RegisterState:
    """Reference-execute a closed program on sparse register states, depending only on the standard library.

    The initial state is a single basis vector with each register at its given integer value (all zero by
    default), and the state is kept as a sparse dictionary from tuples of register integer values to complex
    amplitudes. The instruction stream is lazily expanded by ``events``: single-qubit gates act bit by bit on
    the first operand, while ``xor``, ``swap``, ``add_const`` and ``gphase`` update register values or
    amplitudes per RIR semantics; ``Control`` conditions restrict the affected branches by logical
    conjunction, ``Adjoint`` inverts in reverse order, and ``Repeat`` repeats by the count. A QRAM ``Load``
    is an XOR read, naturally supporting superposed addresses; ``Store`` is a random write requiring the
    address and data to be in a definite basis state. Local registers are zeroed on entry, and all branches
    are checked to be restored on exit.

    Args:
        program: closed RIR program to execute.
        memory: QRAM data provided by resource name, in the same format as ``check_memory``.
        initial: mapping from register names to initial integer values, all zero by default.
        max_steps: expansion step budget; exceeding it raises an error.
        max_states: sparse basis vector budget; exceeding it raises an error.

    Returns:
        RegisterState: entry register layout and final sparse amplitudes; components with magnitude at most 1e-15 are truncated.

    Raises:
        ValidationError: program or QRAM data validation failed, an initial register is unknown or out of range, a Store is not in a definite basis state,
            a local register is not restored to zero, or the step or sparse state count budget is exceeded.
    """
    program = validate(program, require_closed=True)
    memories = check_memory(program, memory)
    registers = program.main.registers
    index = {r.name: i for i, r in enumerate(registers)}
    initial = initial or {}
    if set(initial) - set(index):
        raise ValidationError("initial state contains an unknown register")
    values: list[int] | tuple[int, ...] = []
    for reg in registers:
        value = initial.get(reg.name, 0)
        if type(value) is not int or not 0 <= value < (1 << reg.type.width):
            raise ValidationError("initial register value out of range")
        cast("list[int]", values).append(value)
    state = {tuple(values): 1 + 0j}

    def read(ref: Ref, values: tuple[int, ...]) -> int:
        """Read an unsigned integer out of the value tuple by concatenating view segments from low to high bits."""
        result, offset = 0, 0
        for span in ref.parts:
            result |= (
                (values[index[span.register]] >> span.start) & ((1 << span.width) - 1)
            ) << offset
            offset += span.width
        return result

    def write(ref: Ref, values: tuple[int, ...], value: int) -> tuple[int, ...]:
        """Write an integer value into the value tuple through per-segment masks and return the new value tuple."""
        result, offset = list(values), 0
        for span in ref.parts:
            pos = index[span.register]
            mask = ((1 << span.width) - 1) << span.start
            result[pos] = (result[pos] & ~mask) | (
                ((value >> offset) & ((1 << span.width) - 1)) << span.start
            )
            offset += span.width
        return tuple(result)

    for node, controls, inverse in cast(
        "Iterator[tuple[Primitive | Load | Store | LocalEnter | LocalExit, tuple[tuple[Ref, int], ...], bool]]",
        events(program, max_steps=max_steps),
    ):
        if isinstance(node, LocalEnter):
            index[node.name] = len(next(iter(state)))
            state = {values + (0,): amplitude for values, amplitude in state.items()}
            continue
        if isinstance(node, LocalExit):
            position = index.pop(node.name)
            if any(values[position] != 0 for values in state):
                raise ValidationError(f"local register not restored to zero: {node.name}")
            state = {
                values[:position] + values[position + 1 :]: amplitude
                for values, amplitude in state.items()
            }
            continue

        def active(
            values: tuple[int, ...], controls: tuple[tuple[Ref, int], ...] = controls
        ) -> bool:
            """Decide whether a basis vector satisfies the comparison values of all control conditions."""
            return all(read(ref, values) == expected for ref, expected in controls)

        if isinstance(node, Store):
            observed = {
                (read(node.address, values), read(node.data, values)) for values in state
            }
            if len(observed) != 1:
                raise ValidationError("Store requires the address and data registers to be in a definite basis state at execution time")
            address, value = observed.pop()
            if value:
                memories[node.resource][address] = value
            else:
                memories[node.resource].pop(address, None)
            continue

        if isinstance(node, Primitive) and node.op not in {"xor", "swap", "add_const", "gphase"}:
            matrix = gate_matrix(node.op, node.angle, inverse)
            for register, bit in locations(node.operands[0]):
                result: dict[tuple[int, ...], complex] = {}
                pos = index[register]
                for values, amplitude in state.items():
                    if not active(values):
                        result[values] = result.get(values, 0j) + amplitude
                        continue
                    incoming = (values[pos] >> bit) & 1
                    for outgoing in range(2):
                        target: list[int] | tuple[int, ...] = list(values)
                        cast("list[int]", target)[pos] = (
                            (target[pos] & ~(1 << bit)) | (outgoing << bit)
                        )
                        target = tuple(target)
                        result[target] = (
                            result.get(target, 0j) + matrix[outgoing][incoming] * amplitude
                        )
                state = {k: v for k, v in result.items() if abs(v) > 1e-15}
                if len(state) > max_states:
                    raise ValidationError("reference simulation exceeds the sparse state count budget")
            continue
        result = {}
        for values, amplitude in state.items():
            target = values
            if active(values):
                if isinstance(node, Load):
                    value = memories[node.resource].get(read(node.address, values), 0)
                    target = write(node.data, values, read(node.data, values) ^ value)
                elif node.op == "gphase":
                    amplitude *= cmath.exp(1j * cast(float, node.angle) * (-1 if inverse else 1))
                elif node.op == "xor":
                    source, dest = node.operands
                    target = write(dest, values, read(source, values) ^ read(dest, values))
                elif node.op == "swap":
                    a, b = node.operands
                    av, bv = read(a, values), read(b, values)
                    target = write(b, write(a, values, bv), av)
                elif node.op == "add_const":
                    ref = node.operands[0]
                    value = read(ref, values) + cast(int, node.value) * (-1 if inverse else 1)
                    target = write(ref, values, value & ((1 << ref.width) - 1))
            result[target] = result.get(target, 0j) + amplitude
        state = {k: v for k, v in result.items() if abs(v) > 1e-15}
    return RegisterState(registers, state)
