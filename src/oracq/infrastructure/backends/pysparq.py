"Optional PySparQ execution adapter. RIR root registers correspond to native integer registers."

from __future__ import annotations

import cmath
from collections.abc import Iterator, Mapping, Sequence
from threading import Lock
from typing import cast

from oracq.infrastructure.execution import (
    LocalEnter,
    LocalExit,
    RegisterState,
    check_memory,
    events,
    gate_matrix,
)
from oracq.infrastructure.ir import Load, Primitive, Program, Ref, ValidationError
from oracq.infrastructure.native import NativeContext, NativeRegistry, NativeSite
from oracq.infrastructure.validation import locations, validate

_LOCK = Lock()


def run_pysparq(
    program: Program,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None = None,
    *,
    max_steps: int = 1_000_000,
    max_states: int = 65536,
    native_registry: NativeRegistry | None = None,
    report: dict[str, int | list[str] | str] | None = None,
) -> RegisterState:
    """Execute a program event by event through the PySparQ sparse register simulator.

    RIR root registers map to PySparQ native named registers, with gates
    and QRAM queries applied one by one on the sparse state. An
    interpreter with pysparq installed is required; a missing dependency
    raises ValidationError. The whole execution holds a module-level mutex,
    requires the PySparQ global register table to be empty before taking
    over, and cleans up the registers created by this adapter on normal or
    exceptional exit.

    Args:
        program: The program to execute; it must be closed when ``native_registry`` is not provided.
        memory: QRAM data keyed by resource name, either a full word sequence or a sparse address-to-word dictionary; all zeros when omitted, and when provided it must cover exactly all resources of the entry.
        max_steps: Expansion budget; execution is refused once the estimated step count exceeds it.
        max_states: Sparse basis-state budget; execution aborts once it is exceeded after a gate or native operator application.
        native_registry: Module-level native implementation registry; when provided, open oracles can be closed by native implementations.
        report: The passed-in dictionary is filled in place with execution statistics, keyed ``native_calls``,
            ``gate_events``, ``native_labels`` and ``correctness`` (the last initialized to ``pending``).

    Returns:
        RegisterState: Sparse complex amplitudes keyed by tuples of the entry registers' integer values, with components below 1e-15 dropped.

    Raises:
        ValidationError: The program is illegal, contains a runtime Store, a native implementation is missing or mismatches the module description,
        the QRAM materialization exceeds 2^20 entries, the expansion or sparse-state budget is exceeded, a local workspace is not restored to zero,
        the PySparQ global register table is non-empty, or pysparq is not installed.
    """
    from oracq.infrastructure.linking import uses_store

    program = validate(program, require_closed=native_registry is None)
    if uses_store(program):
        raise ValidationError("the PySparQ adapter does not yet support runtime QRAM writes via Store")
    native_modules = frozenset() if native_registry is None else native_registry.matching(program)
    if native_registry is not None and native_registry.missing(program):
        raise ValidationError("PySparQ is missing implementations: " + ", ".join(native_registry.missing(program)))
    report = {} if report is None else report
    report.update(native_calls=0, gate_events=0, native_labels=[], correctness="pending")
    memories = check_memory(program, memory)
    if any(r.type.address_width > 20 for r in program.main.resources):
        raise ValidationError("the current adapter refuses to materialize QRAMs with more than 2^20 entries")
    # Finish the expansion budget check before taking over the global register table.
    from oracq.infrastructure.execution import expanded_steps

    if expanded_steps(program, max_steps, native_modules) > max_steps:
        raise ValidationError("PySparQ execution exceeds the expansion budget")
    try:
        import pysparq as ps
    except ImportError as exc:
        raise ValidationError("PySparQ execution requires an environment with pysparq installed") from exc

    with _LOCK:
        if ps.System.get_activated_register_size():
            raise ValidationError("the PySparQ global register table is not empty; finish the existing simulation first")
        try:
            state = ps.SparseState()
            names = {reg.name: f"pyqec_{i}" for i, reg in enumerate(program.main.registers)}
            kinds = {
                "bits": ps.StateStorageType.General,
                "uint": ps.StateStorageType.UnsignedInteger,
                "sint": ps.StateStorageType.SignedInteger,
                "rational": ps.StateStorageType.Rational,
            }
            for reg in program.main.registers:
                if reg.type.width:
                    ps.AddRegister(names[reg.name], kinds[reg.type.kind], reg.type.width)(state)
            qrams: dict[str, object] = {}
            for res in program.main.resources:
                data = [0] * (1 << res.type.address_width)
                for address, value in memories[res.name].items():
                    data[address] = value
                qrams[res.name] = ps.QRAMCircuit_qutrit(
                    res.type.address_width, res.type.data_width, data
                )

            def apply(operator: object, controls: Sequence[tuple[str, int]] = ()) -> None:
                """Apply an operator to the current sparse state and check the sparse-state budget."""
                if controls:
                    # pysparq native operators are dynamic backend objects with no type stubs available.
                    operator.conditioned_by_bit(list(controls))  # type: ignore[attr-defined]
                operator(state)  # type: ignore[operator]
                if len(state.basis_states) > max_states:
                    raise ValidationError("PySparQ execution exceeds the sparse-state budget")

            def actual(ref: Ref) -> list[tuple[str, int]]:
                """Resolve a view into the list of PySparQ ``(register name, bit)`` coordinates."""
                return [(names[name], bit) for name, bit in locations(ref)]

            # Store was already rejected above, so only the other five node kinds actually appear in the event stream.
            for node, conditions, inverse in cast(
                "Iterator[tuple[Primitive | Load | NativeSite | LocalEnter | LocalExit, tuple[tuple[Ref, int], ...], bool]]",
                events(program, max_steps=max_steps, native_modules=native_modules),
            ):
                if isinstance(node, LocalEnter):
                    names[node.name] = "pyqec_" + node.name
                    ps.AddRegister(names[node.name], kinds[node.type.kind], node.type.width)(state)
                    continue
                if isinstance(node, LocalExit):
                    native_name = names[node.name]
                    rid = ps.System.get_id(native_name)
                    if any(
                        int(basis.get(rid).value) & ((1 << node.width) - 1)
                        for basis in state.basis_states
                    ):
                        raise ValidationError(f"the PySparQ local workspace was not restored to zero: {node.name}")
                    ps.RemoveRegister(native_name)(state)
                    del names[node.name]
                    continue
                controls, zeros = [], []
                for ref, expected in conditions:
                    for bit_index, pair in enumerate(actual(ref)):
                        controls.append(pair)
                        if not ((expected >> bit_index) & 1):
                            zeros.append(pair)
                for pair in zeros:
                    ps.Xgate_Bool(*pair)(state)
                if isinstance(node, NativeSite):
                    # NativeSite events are only produced when a native registry is provided, i.e. native_modules is non-empty.
                    entry = cast(NativeRegistry, native_registry).entries[node.module.name]
                    context = NativeContext(ps, node, names, qrams, memories)
                    untouched = (
                        ps.split_systems(
                            state,
                            [],
                            [],
                            [(ps.System.get_id(name), bit) for name, bit in controls],
                            [],
                        )
                        if controls
                        else None
                    )
                    try:
                        if state.size():
                            operator = entry.factory(context)
                            (operator.dag if inverse else operator)(state)
                            if len(state.basis_states) > max_states:
                                raise ValidationError("the native operator exceeds the sparse-state budget")
                    finally:
                        if untouched is not None:
                            ps.combine_systems(state, untouched)
                    cast("dict[str, int]", report)["native_calls"] += 1
                    if entry.label not in cast("list[str]", report["native_labels"]):
                        cast("list[str]", report["native_labels"]).append(entry.label)
                elif isinstance(node, Load):
                    # The view is copied into temporary integer registers via reversible XOR; fully uncomputed after loading.
                    a: str | tuple[str, int]
                    d: str
                    a, d = "pyqec_tmp_address", "pyqec_tmp_data"
                    ps.AddRegister(a, ps.StateStorageType.General, node.address.width)(state)
                    ps.AddRegister(d, ps.StateStorageType.General, node.data.width)(state)
                    for i, pair in enumerate(actual(node.address)):
                        apply(ps.Xgate_Bool(a, i), (pair,))
                    query = ps.QRAMLoad(qrams[node.resource], a, d)
                    query(state)
                    for i, pair in enumerate(actual(node.data)):
                        apply(ps.Xgate_Bool(*pair), tuple(controls) + ((d, i),))
                    query(state)
                    for i, pair in reversed(list(enumerate(actual(node.address)))):
                        apply(ps.Xgate_Bool(a, i), (pair,))
                    ps.RemoveRegister(d)(state)
                    ps.RemoveRegister(a)(state)
                elif node.op == "gphase":
                    theta = cast(float, node.angle) * (-1 if inverse else 1)
                    if controls:
                        apply(ps.Phase_Bool(*controls[-1], theta), controls[:-1])
                    else:
                        root = next(r for r in program.main.registers if r.type.width)
                        phase = cmath.exp(1j * theta)
                        apply(ps.Rot_Bool(names[root.name], 0, [phase, 0j, 0j, phase]))
                elif node.op in {"xor", "swap"}:
                    for a, b in zip(
                        actual(node.operands[0]), actual(node.operands[1]), strict=True
                    ):
                        apply(ps.Xgate_Bool(*b), tuple(controls) + (a,))
                        if node.op == "swap":
                            apply(ps.Xgate_Bool(*a), tuple(controls) + (b,))
                            apply(ps.Xgate_Bool(*b), tuple(controls) + (a,))
                elif node.op == "add_const":
                    bits = actual(node.operands[0])
                    value = (-cast(int, node.value) if inverse else cast(int, node.value)) % (
                        1 << len(bits)
                    )
                    for offset in range(len(bits)):
                        if (value >> offset) & 1:
                            for i in reversed(range(offset + 1, len(bits))):
                                apply(
                                    ps.Xgate_Bool(*bits[i]), tuple(controls) + tuple(bits[offset:i])
                                )
                            apply(ps.Xgate_Bool(*bits[offset]), controls)
                else:
                    matrix = gate_matrix(node.op, node.angle, inverse)
                    for pair in actual(node.operands[0]):
                        apply(ps.Rot_Bool(*pair, [v for row in matrix for v in row]), controls)
                if not isinstance(node, NativeSite):
                    cast("dict[str, int]", report)["gate_events"] += 1
                for pair in reversed(zeros):
                    ps.Xgate_Bool(*pair)(state)
            result: dict[tuple[int, ...], complex] = {}
            for basis in state.basis_states:
                key = tuple(
                    (
                        int(basis.get(ps.System.get_id(names[r.name])).value)
                        & ((1 << r.type.width) - 1)
                    )
                    if r.type.width
                    else 0
                    for r in program.main.registers
                )
                result[key] = result.get(key, 0j) + complex(basis.amplitude)
            return RegisterState(
                program.main.registers, {k: v for k, v in result.items() if abs(v) > 1e-15}
            )
        finally:
            ps.System.clear()


def run_pysparq_rir(
    program: Program,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None = None,
    *,
    max_steps: int = 1_000_000,
    max_states: int = 65536,
) -> RegisterState:
    """Execute through the PySparQ native RIR interpreter; an implementation independent of run_pysparq, used for cross-validation.

    Args:
        program: The closed RIR program to execute; it must not contain runtime QRAM writes via Store.
        memory: Binding from resource names to QRAM contents; sequences give data words by index, mappings by address.
        max_steps: The interpreter's instruction-step budget cap for expanded execution.
        max_states: The interpreter's budget cap on the number of maintained basis states.

    Returns:
        RegisterState: The final state over the entry's public register space, as sparse complex
        amplitudes keyed by tuples of the registers' integer values.
    """
    from oracq.infrastructure.linking import uses_store
    from oracq.infrastructure.serialization import dumps

    program = validate(program, require_closed=True)
    if uses_store(program):
        raise ValidationError("the PySparQ RIR interpreter does not yet support runtime QRAM writes via Store")
    memories = check_memory(program, memory)
    try:
        import pysparq as ps
    except ImportError as exc:
        raise ValidationError("PySparQ RIR execution requires an environment with pysparq installed") from exc
    result = ps.run_rir(dumps(program), memories, max_steps=max_steps, max_states=max_states)
    return RegisterState(program.main.registers, dict(result.amplitudes))
