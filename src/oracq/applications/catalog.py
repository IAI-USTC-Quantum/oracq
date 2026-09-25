"Reproducibly generated catalog of paradigm cases: open IR, binding plans and a finite set of concrete instances."

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, cast

from oracq.algorithms.basics.oracle_algorithms import deutsch_jozsa
from oracq.algorithms.common.estimation import phase_estimation
from oracq.algorithms.common.search import grover
from oracq.algorithms.common.transforms import (
    oblivious_amplification,
    qsvt_sequence,
    qubitization_walk,
)
from oracq.algorithms.input_model.block_encoding import (
    direct_sum,
    kronecker_sum,
    lcu,
    matrix_pauli_encoding,
    pad_signal,
    tensor,
)
from oracq.algorithms.input_model.operators import BlockEncoding, identity, pauli_x
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    XorDatabase,
    abstract_block_encoding,
    abstract_database,
    abstract_sparse_access,
    abstract_state_prep,
    declare,
    diagonal_block_encoding,
    gate_database,
    gate_state_prep,
    phase_marks,
    qram_database,
    qram_state_angles,
    qram_state_prep,
    sparse_entry,
    sparse_location_gate,
    sparse_location_qram,
    uniform_state,
)
from oracq.algorithms.input_model.sparse import (
    batch_lookup,
    reversible_lookup,
    sparse_block_encoding,
    word_rotation,
)
from oracq.algorithms.qlss.qlss import CostaConfig, costa_qlss, make_costa_qlss
from oracq.algorithms.qode.legacy import make_lchs_qode
from oracq.algorithms.qode.ode import make_euler_history_qode
from oracq.algorithms.qpde.pde import make_qpde
from oracq.applications.legacy import (
    qfvm_inputs,
    qfvm_step,
    qham_initial_vector,
    qham_lift_m1,
    qham_m1,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Program
from oracq.infrastructure.linking import Binding, bind
from oracq.infrastructure.readout import ReadoutAction

if TYPE_CHECKING:
    from oracq.infrastructure.backends import OriginIRArtifact


@dataclass(frozen=True)
class Case:
    """One paradigm reference case: an open program, implementation bindings and host memory data.

    Attributes:
        name: Case name, one of the entries of ``CASES``.
        program: The case's ``Program``; abstract slots await filling by ``bindings``.
        bindings: Mapping of slot names to ``Binding`` or ``Operation``; unlisted slots stay open.
        memory: Host-side memory tables of QRAM cases; keys are the names referenced by binding resource maps.
        source: Tuple of provenance identifiers (language specification, spec tests or reference workload entries).
        notes: Tuple of qualifying notes on the case's applicability boundary.
        readout: Tuple of terminal ``ReadoutAction`` readout actions executed in order after export.
    """

    name: str
    program: Program
    bindings: dict[str, Binding | Operation] = field(default_factory=dict)
    memory: dict[str, list[int] | dict[int, int]] = field(default_factory=dict)
    source: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    readout: tuple[ReadoutAction, ...] = ()

    def closed(self) -> Program:
        """Bind the open slots of ``program`` per ``bindings`` and return the closed ``Program``.

        Returns:
            Program: The closed program with each open slot replaced by its bound implementation.
        """
        return bind(self.program, self.bindings)

    def artifact(self) -> OriginIRArtifact:
        """Export the OriginIR artifact of the closed program.

        Returns:
            OriginIRArtifact: The OriginIR export artifact of the ``closed()`` result.

        The OriginIR exporter is imported lazily inside this method, so core
        code gains no backend dependency."""
        from oracq.infrastructure.backends import export_originir

        return export_originir(self.closed())


CASES = (
    "bell",
    "ghz",
    "python_generators",
    "trotter_hamsim",
    "dj_gate",
    "dj_qram",
    "grover_gate",
    "grover_qram",
    "stateprep_gate",
    "stateprep_qram",
    "sparse_gate",
    "sparse_qram",
    "costa_gate",
    "costa_qram",
    "costa_sparse_qram",
    "qfvm_gate",
    "qfvm_qram",
    "qham_qode",
    "qham_qpde",
    "qpe",
    "qsvt",
    "oaa",
    "lchs",
    "heat_qode",
    "poisson_qlss",
    "carleman_step",
    "schrodingerisation",
    "be_algebra",
    "arithmetic",
    "batch_qram",
    "banked_qram",
    "register_views",
    "measure_reset",
)
"""All case names accepted by ``build_case``.

Naming convention: most cases are prefixed with an application family (e.g.
``dj_``, ``grover_``, ``qham_``); the ``_gate`` suffix denotes a fixed
gate-level implementation, the ``_qram`` suffix denotes an implementation that
binds QRAM databases and carries host memory tables, and the
``_qode``/``_qpde`` suffixes distinguish solver wrapper layers; the remaining
unsuffixed names are demo cases of primitives or single algorithms.
"""


def _sparse_case(name: str) -> Case:
    """Build a sparse access case; with the ``costa`` prefix, additionally layer on the Costa QLSS solver assembly."""
    access = abstract_sparse_access("Sparse", 1, 2, 2)
    amplitude = declare(
        "SparseAmplitude", {"value": Bits(2), "amplitude": Bits(1)}, paradigm="reversible_function"
    )
    a = sparse_block_encoding(access, amplitude)
    binding: dict[str, Binding | Operation] = {"SparseAmplitude": word_rotation(2)}
    memory: dict[str, list[int] | dict[int, int]] = {}
    if name.endswith("qram"):
        binding["Sparse_position"] = Binding(
            sparse_location_qram(1), {"forward": "positions", "inverse": "inverse_positions"}
        )
        binding["Sparse_entry"] = Binding(
            sparse_entry(qram_database(2, 2), 1), {"db__table": "entries"}
        )
        memory = {
            "positions": [0, 1, 1, 0],
            "inverse_positions": [0, 1, 1, 0],
            "entries": [0, 1, 1, 0],
        }
    else:
        binding["Sparse_position"] = sparse_location_gate(1, [[0, 1], [1, 0]])
        binding["Sparse_entry"] = sparse_entry(gate_database(2, 2, [0, 1, 1, 0]), 1)
    if name.startswith("costa"):
        rhs = abstract_state_prep("Rhs", 1)
        binding["Rhs"] = uniform_state(1).operation
        operation = costa_qlss(a, rhs, CostaConfig(steps=1)).operation
    else:
        operation = a.operation
    return Case(
        name,
        operation.program(),
        binding,
        memory,
        ("language-spec-v2 §E.4; reference-workloads 4/5; CKS §1.1",),
        ("T†SWAP T sparse candidate construction; validation of the matrix and normalization is deferred to the next stage.",),
    )


def _qfvm_case(name: str) -> Case:
    """Build a QFVM paradigm case; the ``qram`` suffix binds the database slots to host memory tables."""
    inputs = qfvm_inputs()
    # QLSSProtocol.__call__ is loosely annotated on the algorithms side (*args -> SolveResult | StateOracle);
    # narrowed here to the solver callback as used by qfvm_step.
    state = qfvm_step(
        inputs, cast("Callable[..., StateOracle]", make_costa_qlss(CostaConfig(steps=1)))
    )
    bindings: dict[str, Binding | Operation] = {}
    memory: dict[str, list[int] | dict[int, int]] = {}
    database_slots = [
        (inputs.position, "position", [0, 1, 1, 0]),
        (inputs.reverse_slot, "reverse_slot", [0, 0, 1, 1]),
        (inputs.column_position, "column_position", [0, 1, 1, 0]),
        (inputs.column_reverse_slot, "column_reverse_slot", [0, 0, 1, 1]),
        (inputs.flow, "flow", [0, 1]),
        (inputs.boundary, "boundary", [1, 0]),
    ]
    for slot, label, data in database_slots:
        impl = (
            qram_database(slot.address_width, slot.data_width)
            if name.endswith("qram")
            else gate_database(slot.address_width, slot.data_width, data)
        )
        if name.endswith("qram"):
            bindings[slot.operation.module.name] = Binding(impl.operation, {"table": label})
            memory[label] = data
        else:
            bindings[slot.operation.module.name] = impl.operation
    ins = {"row": 1, "slot": 1, "column": 1, "flow": 1, "boundary": 1}
    outs = {"value": 2, "status": 1}
    table = {i: (((i & 1) + ((i >> 3) & 1) + 2 * ((i >> 4) & 1)) % 4) | 4 for i in range(32)}
    physical = reversible_lookup(
        ins, outs, qram_database(5, 3) if name.endswith("qram") else gate_database(5, 3, table)
    )
    bindings[inputs.physical_entry.module.name] = (
        Binding(physical, {"db__table": "physical_entries"}) if name.endswith("qram") else physical
    )
    if name.endswith("qram"):
        memory["physical_entries"] = table
    bindings[inputs.amplitude.module.name] = word_rotation(2)
    bindings[inputs.residual.operation.module.name] = uniform_state(1).operation
    return Case(
        name,
        state.operation.program(),
        bindings,
        memory,
        ("reference-workloads 1: QFVM sparse inputs, T_L†SWAP T_R, QLSS/filter",),
        (
            "The physical entry is bound to an explicit toy lookup function and does not pretend to be a full Roe implementation.",
            "This case validates bidirectional access, resource binding and filtered solver assembly; fluid and quantum correctness await the next stage.",
        ),
    )


def _qham_case(name: str) -> Case:
    """Build the QHAM m=1 lifted-system case; the ``qpde`` suffix instead solves through the QPDE wrapper."""
    linear_impl = matrix_pauli_encoding([[-0.2, 0.1], [0.1, -0.2]])
    fold_impl = matrix_pauli_encoding([[0, 0.5, 0, 0], [0, 0, -0.5, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    linear = abstract_block_encoding("QhamLinear", 1, linear_impl.signal_qubits, linear_impl.alpha)
    fold = abstract_block_encoding("QhamFold", 2, fold_impl.signal_qubits, fold_impl.alpha)
    lifted = qham_lift_m1(linear, fold)
    initial = abstract_state_prep("QhamInitial", lifted.width)
    # QLSSProtocol.__call__ is loosely annotated on the algorithms side; narrowed per the generator function protocol.
    qode = make_euler_history_qode(
        cast(
            "Callable[[BlockEncoding, StatePreparation], StateOracle]",
            make_costa_qlss(CostaConfig(steps=1)),
        ),
        steps=1,
    )
    solver = make_qpde(qode) if name.endswith("qpde") else qode
    state = qham_m1(linear, fold, initial, solver, via_pde=name.endswith("qpde"))
    return Case(
        name,
        state.operation.program(),
        {
            "QhamLinear": linear_impl.operation,
            "QhamFold": fold_impl.operation,
            "QhamInitial": gate_state_prep(qham_initial_vector([0.3, 0.1])).operation,
        },
        {},
        ("reference-workloads 6; generator-protocols §8; QHAM m=1 lifted system",),
        (
            "m=1 with two spatial points; the lifted dimension is 10, padded to 16.",
            "Generated through the QODE history system and Costa/filter; HAM convergence and the conditioned solution state are not validated at this stage.",
        ),
    )


def build_case(name: str) -> Case:
    """Build one paradigm reference case by name.

    All cases are generated from fixed constants and can be constructed
    reproducibly; QRAM cases also provide host memory tables.

    Args:
        name: Case name; must be one of the entries of ``CASES``.

    Returns:
        Case: The corresponding case, with implementation bindings, memory data and readout actions where applicable.

    Raises:
        KeyError: The case name is not in ``CASES``.
    """
    if name not in CASES:
        raise KeyError(name)
    if name in {"bell", "ghz"}:
        width = 2 if name == "bell" else 3
        b = Builder(name, {"q": Bits(width)})
        b.h(b["q"][0])
        for bit in range(1, width):
            b.xor(b["q"][0], b["q"][bit])
        return Case(name, b.finish().program(), source=("spec-tests 00-primitives",))
    if name == "python_generators":

        def count(n: int) -> int:
            """Return the n-th Fibonacci number by recursion, as a generation-time constant."""
            return 1 if n < 2 else count(n - 1) + count(n - 2)

        width = count(3)

        def generate(depth: int) -> Operation:
            """Recursively build a demo module with ``depth`` levels of nested ``Repeat`` calls."""
            b = Builder(f"recursive_{depth}", {"q": Bits(width)}, attributes={"const_alpha": 2.5})
            if depth:
                with b.repeat(2):
                    b.call(generate(depth - 1), q=b["q"])
            else:
                for bit, angle in enumerate(tuple(0.1 * (i + 1) for i in range(width))):
                    b.rz(b["q"][bit], angle)
            return b.finish()

        return Case(name, generate(2).program(), source=("spec-tests 02-generics-const",))
    if name == "trotter_hamsim":
        from oracq.algorithms.common.hamiltonian import trotter_hamsim

        return Case(
            name,
            trotter_hamsim([(0.5, "X"), (0.7, "Z")], 0.2).program(),
            source=("spec-tests 06-protocols/hamsim-trotter",),
        )
    if name == "carleman_step":
        implementation: BlockEncoding | XorDatabase = matrix_pauli_encoding(
            [[0, 0, 0, 0], [0, -0.2, 0.1, 0], [0, 0, -0.4, 0.2], [0, 0, 0, -0.6]]
        )
        generator = abstract_block_encoding(
            "CarlemanGenerator",
            2,
            cast(BlockEncoding, implementation).signal_qubits,
            cast(BlockEncoding, implementation).alpha,
        )
        initial = abstract_state_prep("CarlemanInitial", 2)
        solver = make_euler_history_qode(
            cast(
                "Callable[[BlockEncoding, StatePreparation], StateOracle]",
                make_costa_qlss(CostaConfig(steps=1)),
            ),
            steps=1,
        )
        return Case(
            name,
            solver(generator, initial, 0.1).operation.program(),
            {
                "CarlemanGenerator": implementation.operation,
                "CarlemanInitial": gate_state_prep([1, 0.3, 0.09, 0.027]).operation,
            },
            {},
            ("spec-tests 07-scientific/carleman-step",),
            ("Third-order truncated lift of a quadratic scalar model; only the assembly is accepted, the truncation error is not certified.",),
        )
    if name == "schrodingerisation":
        from oracq.algorithms.qode.legacy import make_schrodingerisation_qode

        hamiltonian = abstract_block_encoding("LiftedHamiltonian", 2, 1, 1.0)
        initial = abstract_state_prep("SchrodingerInitial", 1)
        solver = make_schrodingerisation_qode(
            lambda g: hamiltonian, lambda h, t: qsvt_sequence(h, [t, -t, t])
        )
        state = solver(identity(1), initial, 0.1)
        return Case(
            name,
            state.operation.program(),
            {
                "LiftedHamiltonian": pad_signal(pauli_x(2), 1).operation,
                "SchrodingerInitial": uniform_state(1).operation,
            },
            {},
            ("spec-tests 07-scientific/schrodingerization; QODE alternative",),
            ("The Hamiltonian lift is an explicit oracle boundary; the demo binding only checks the paradigm and does not claim to reproduce any specific PDE.",),
        )
    if name == "banked_qram":
        from oracq.algorithms.input_model.oracles import banked_database

        slot: Operation | XorDatabase | StatePreparation = banked_database(2, 96, abstract=True)
        impl: Operation | StatePreparation | BlockEncoding = banked_database(2, 96)
        return Case(
            name,
            cast(Operation, slot).program(),
            {cast(Operation, slot).module.name: Binding(cast(Operation, impl), {"bank0": "low", "bank1": "high"})},
            {"low": {0: 5}, "high": {0: 7}},
            ("reference-workloads 1 packed flow words; register-level storage",),
            ("The 96-bit logical data is split into two banks of 64/32 bits; this case accepts the description only and is not executed on the dense simulator.",),
        )
    if name == "register_views":
        from oracq.infrastructure.ir import fuse

        b = Builder("views_demo", {"input": Bits(4), "output": Bits(4)})
        b.h(b["input"][:2])
        b.xor(fuse(b["input"][2:], b["input"][:2]), b["output"])
        return Case(name, b.finish().program(), source=("spec-tests 01-registers / 00-primitives",))
    if name == "measure_reset":
        b = Builder("readout_demo", {"q": Bits(2)})
        b.h(b["q"][0])
        b.xor(b["q"][0], b["q"][1])
        return Case(
            name,
            b.finish().program(),
            source=("spec-tests 00-primitives/measure-reset",),
            readout=(ReadoutAction("measure", "q"), ReadoutAction("reset", "q")),
        )
    if name.startswith("qfvm_"):
        return _qfvm_case(name)
    if name.startswith("qham_"):
        return _qham_case(name)
    if name.startswith("sparse_") or name == "costa_sparse_qram":
        return _sparse_case(name)
    if name.startswith("dj_"):
        slot = abstract_database("BooleanFunction", 2, 1)
        implementation = (
            qram_database(2, 1) if name.endswith("qram") else gate_database(2, 1, [0, 1, 1, 0])
        )
        return Case(
            name,
            deutsch_jozsa(slot).program(),
            {
                "BooleanFunction": Binding(implementation.operation, {"table": "truth"})
                if name.endswith("qram")
                else implementation.operation
            },
            {"truth": [0, 1, 1, 0]} if name.endswith("qram") else {},
            ("Deutsch–Jozsa user-added workload",),
        )
    if name.startswith("grover_"):
        from oracq.algorithms.common.search import phase_from_database

        work = 1 if name.endswith("qram") else 0
        slot = declare(
            "Mark",
            {"target": Bits(2), **({"work": Bits(1)} if work else {})},
            paradigm="phase_oracle",
        )
        if work:
            impl = phase_from_database(qram_database(2, 1))
            binding: Binding | Operation = Binding(impl, {"db__table": "marks"})
            memory: dict[str, list[int] | dict[int, int]] = {"marks": [0, 0, 0, 1]}
        else:
            binding = phase_marks(2, [3])
            memory = {}
        return Case(
            name,
            grover(slot, 2).operation.program(),
            {"Mark": binding},
            memory,
            ("language-spec-v2 §15.1; spec-tests 03-oracle/grover2",),
        )
    if name.startswith("stateprep_"):
        amplitudes = (1, 2, 3, 4)
        if name.endswith("qram"):
            impl = qram_state_prep(2, 3)
            slot = abstract_state_prep("Prepare", 2, impl.work_width)
            binding = Binding(impl.operation, {"angles": "angles"})
            memory = {"angles": qram_state_angles(amplitudes, 3)}
        else:
            slot = abstract_state_prep("Prepare", 2)
            binding = gate_state_prep(amplitudes).operation
            memory = {}
        return Case(
            name,
            slot.operation.program(),
            {"Prepare": binding},
            memory,
            ("spec-tests 03-oracle/state-prep-isometry; reference-workloads 4",),
        )
    if name in {"costa_gate", "costa_qram", "poisson_qlss"}:
        qram = name.endswith("qram")
        a = abstract_block_encoding("A", 1, 3, 1.0)
        bp = abstract_state_prep("B", 1, 4 if qram else 0)
        a_impl = diagonal_block_encoding(
            qram_database(1, 2) if qram else gate_database(1, 2, [0, 1])
        )
        b_impl = qram_state_prep(1, 3) if qram else gate_state_prep([1, 2])
        if name == "poisson_qlss":
            a_impl = matrix_pauli_encoding([[2, -1], [-1, 2]])
            a = abstract_block_encoding("A", 1, a_impl.signal_qubits, a_impl.alpha)
            b_impl = gate_state_prep([1, 0])
        bindings: dict[str, Binding | Operation] = {
            "A": Binding(a_impl.operation, {"db__table": "matrix_angles"})
            if qram
            else a_impl.operation,
            "B": Binding(b_impl.operation, {"angles": "rhs_angles"}) if qram else b_impl.operation,
        }
        memory = (
            {"matrix_angles": [0, 1], "rhs_angles": qram_state_angles([1, 2], 3)} if qram else {}
        )
        return Case(
            name,
            costa_qlss(a, bp).operation.program(),
            bindings,
            memory,
            ("reference-workloads 4/5; language-spec-v2 §E.10",),
            ("general walk with actual unary Dolph–Chebyshev LCU filtering; numerical certification is deferred.",),
        )
    if name in {"qpe", "qsvt", "oaa", "lchs", "heat_qode"}:
        impl = (
            matrix_pauli_encoding([[-0.2, 0.1], [0.1, -0.2]])
            if name == "heat_qode"
            else pad_signal(pauli_x(1), 1)
        )
        a = abstract_block_encoding("Operator", 1, impl.signal_qubits, impl.alpha)
        bindings = {"Operator": impl.operation}
        if name == "qpe":
            op = phase_estimation(qubitization_walk(a), precision=2)
        elif name == "qsvt":
            op = qsvt_sequence(a, [0.1, -0.2, 0.3])
        elif name == "oaa":
            op = oblivious_amplification(a)
        else:
            initial = abstract_state_prep("Initial", 1)
            bindings["Initial"] = uniform_state(1).operation
            if name == "heat_qode":
                solver = make_euler_history_qode(
                    cast(
                        "Callable[[BlockEncoding, StatePreparation], StateOracle]",
                        make_costa_qlss(CostaConfig(steps=1)),
                    ),
                    steps=1,
                )
            else:
                solver = make_lchs_qode(
                    lambda g, t: qsvt_sequence(g, [t, -t, t]), [0.5, 1.0], [0.4, 0.6]
                )
            op = solver(a, initial, 0.1).operation
        return Case(
            name,
            op.program(),
            bindings,
            {},
            ("language-spec-v2 §E.7–E.11; spec-tests 05-qsvt-qpe/07-scientific",),
            ("This case validates the calling paradigm; binding a given phase or generator does not constitute an accuracy guarantee for PDEs or the exponential function.",),
        )
    if name == "be_algebra":
        a = abstract_block_encoding("A", 1, 0, 1.0)
        # The short variable name b is reused across branches (other branches use Builder); the proper fix is a rename, but this task forbids renames.
        b = abstract_block_encoding("B", 1, 1, 2.0)  # type: ignore[assignment]
        combined = lcu(
            [
                (1, direct_sum(a, cast(BlockEncoding, b))),
                (0.5, tensor(a, cast(BlockEncoding, b))),
                (-0.25, kronecker_sum(a)),
            ]
        )
        return Case(
            name,
            combined.operation.program(),
            {"A": identity(1).operation, "B": _one_signal_identity(2.0)},
            {},
            ("reference-workloads 2; operation tuple/array heterogeneous LCU",),
        )
    if name == "arithmetic":
        slot = declare(
            "Multiply",
            {"left": Bits(2), "right": Bits(2), "output": Bits(4)},
            paradigm="reversible_function",
            attributes={"fraction_bits": 1},
        )
        table = {i: (i & 3) * (i >> 2) for i in range(16)}
        impl = reversible_lookup({"left": 2, "right": 2}, {"output": 4}, gate_database(4, 4, table))
        b = Builder("arithmetic_demo", {"a": Bits(2), "b": Bits(2), "out": Bits(4)})
        b.h(b["a"])
        b.call(slot, left=b["a"], right=b["b"], output=b["out"])
        return Case(
            name,
            b.finish().program(),
            {"Multiply": impl},
            {},
            ("reference-workloads 3; reversible arithmetic and fixed-point interpretation",),
        )
    if name == "batch_qram":
        slot = abstract_database("BatchData", 4, 8)
        return Case(
            name,
            batch_lookup(slot, 8).program(),
            {"BatchData": Binding(qram_database(4, 8).operation, {"table": "batch"})},
            {"batch": {0: 7, 3: 12}},
            ("language-spec-v2 §E.13; register arrays / paged QRAM",),
        )
    raise AssertionError(name)


def _one_signal_identity(alpha: float) -> Operation:
    """Return an identity block encoding implementation scaled by ``alpha`` and padded with one signal bit."""
    from oracq.algorithms.input_model.operators import scale

    return pad_signal(scale(alpha, identity(1)), 1).operation
