"""QPDE/QODE tutorial: different input paradigms, batched binding, and replaceable generating functions.

Run: PYTHONPATH=src python examples/ode_input_models.py
Optional --native-parse parses the exports with real uniqc; --native-bindings cross-checks small QRAM bindings against a real backend.
These checks do not certify solution accuracy.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from functools import partial
from pathlib import Path
from typing import cast

from oracq import (
    Binding,
    BlockEncoding,
    FixedFormat,
    Operation,
    bind,
    dump_qram_yaml,
    dumps,
    export_originir,
    export_toffoli_u3_cz,
    identity,
    load_qram_yaml,
    loads,
    run_pysparq,
    scale,
    unresolved,
    zero,
)
from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    abstract_block_encoding,
    abstract_database,
    abstract_sparse_access,
    abstract_state_prep,
    diagonal_block_encoding,
    gate_database,
    gate_state_prep,
    qram_database,
    qram_state_angles,
    qram_state_prep,
    sparse_entry,
    sparse_location_gate,
    sparse_location_qram,
)
from oracq.algorithms.input_model.sparse import real_symmetric_sparse_encoding
from oracq.algorithms.qnlss.carleman import PolynomialODE, carleman_qode
from oracq.algorithms.qode.lchs import QuadraturePlan, lchs_qode
from oracq.algorithms.qode.ode import linear_qode
from oracq.algorithms.qode.ode_models import HermitianParts, LinearODE
from oracq.algorithms.qode.schrodingerization import SchrodingerPlan
from oracq.algorithms.qpde.pde import DiscretePDE, PDEInput, make_qpde, qpde_solver
from oracq.applications.qham import (
    Discretization,
    Field,
    Grid,
    PolynomialPDE,
    QHAMBindings,
    structured_fd_bindings,
)
from oracq.applications.qham.stencils import derivative_encoding


def save_case(
    root: Path,
    name: str,
    state: StateOracle,
    bindings: Mapping[str, Operation | Binding] | None = None,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None = None,
    *,
    notes: Mapping[str, object] | None = None,
    native_parse: bool = False,
) -> dict[str, object]:
    """Leave open/partially-bound/closed artifacts side by side; oracle bodies are not inlined."""
    folder = root / name
    folder.mkdir(parents=True, exist_ok=True)
    opened = state.operation.program()
    assert loads(dumps(opened)) == opened
    (folder / "open.rir.yaml").write_text(dumps(opened), encoding="utf-8")
    closed = opened
    for index, (slot, implementation) in enumerate((bindings or {}).items()):
        closed = bind(closed, {slot: implementation})
        if index == 0:
            (folder / "partial.rir.yaml").write_text(dumps(closed), encoding="utf-8")
    assert not unresolved(closed), unresolved(closed)
    assert loads(dumps(closed)) == closed
    (folder / "closed.rir.yaml").write_text(dumps(closed), encoding="utf-8")
    (folder / "modular.originir").write_text(export_originir(closed).text, encoding="utf-8")
    basis = export_toffoli_u3_cz(closed).text
    (folder / "toffoli_u3_cz.originir").write_text(basis, encoding="utf-8")
    (folder / "memory.qram.yaml").write_text(dump_qram_yaml(closed, memory), encoding="utf-8")
    if native_parse:
        from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser

        OriginIR_BaseParser().parse(basis)
    record = {
        "case": name,
        "open_slots": [r.name for r in unresolved(opened)],
        "modules": len(closed.modules),
        "resources": [r.name for r in closed.main.resources],
        "target_width": state.width,
        "signal_width": state.signal_qubits,
        "native_parse": "passed" if native_parse else "not_requested",
        "notes": notes or {},
        "correctness": "solver accuracy and success probability pending",
    }
    (folder / "report.json").write_text(
        json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        name,
        "modules",
        len(closed.modules),
        "open slots",
        len(cast(list[str], record["open_slots"])),
        flush=True,
    )
    return record


def verify_bindings(root: Path) -> None:
    """Compare the full complex amplitudes of the two bindings of the angle-database cases with real PySparQ."""
    records = []
    for stem in ("given_xor",):
        states = []
        for implementation in ("gate", "qram"):
            directory = root / f"{stem}_{implementation}_lchs"
            program = loads((directory / "closed.rir.yaml").read_text(encoding="utf-8"))
            memory = load_qram_yaml(directory / "memory.qram.yaml")
            states.append(run_pysparq(program, memory).amplitudes)
        keys = states[0].keys() | states[1].keys()
        error = max(abs(states[0].get(k, 0) - states[1].get(k, 0)) for k in keys)
        assert error < 1e-10, (stem, error)
        records.append(
            {
                "case": stem,
                "gate_vs_qram_max_amplitude_difference": error,
                "backend": "real pysparq",
                "scope": "same finite circuit construction, not ODE convergence",
            }
        )
        print(stem, "gate/QRAM amplitude difference", error, flush=True)
    (root / "binding-validation.json").write_text(json.dumps(records, indent=2), encoding="utf-8")


def polynomial_from_bindings(bindings: QHAMBindings) -> PolynomialODE:
    """Tutorial host adapter: merge the PDE multi-linear ports by degree and hand them to Carleman."""
    grouped = defaultdict(list)
    for _, port in bindings.ports:
        grouped[port.arity].append((1, port.encoding))
    return PolynomialODE(
        bindings.state_width,
        tuple((order, lcu(terms)) for order, terms in sorted(grouped.items())),
        bindings.initial,
        bindings.initial_norm,
    )


def shifted_solver(
    qode: Callable[[BlockEncoding, StatePreparation, float], StateOracle],
    recovery: dict[str, float],
) -> Callable[[BlockEncoding, StatePreparation, float], StateOracle]:
    """Apply a shift to the lifted linear system; the norm-recovery record goes into the host report."""

    def solve(generator: BlockEncoding, initial: StatePreparation, time: float) -> StateOracle:
        """Hand the generator shifted to ``G - mu*I`` to the underlying solver protocol for evolution."""
        mu = generator.alpha
        shifted = lcu([(1, generator), (-mu, identity(generator.width))])
        recovery.update(growth_shift=mu, log_amplitude_rescale=mu * time)
        return qode(shifted, initial, time)

    return solve


def main() -> None:
    """Build the cases for each input paradigm and write the open, partially-bound, and closed artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path, default=Path("out/ode-input-models"))
    parser.add_argument("--native-parse", action="store_true", help="parse all exports with real uniqc")
    parser.add_argument("--native-bindings", action="store_true", help="cross-check the angle-database bindings against PySparQ")
    args = parser.parse_args()
    write = partial(save_case, args.output, native_parse=args.native_parse)
    records = []
    time = 0.05
    hamiltonian_function = partial(taylor_hamiltonian, degree=1)
    quadrature = QuadraturePlan.cauchy(cutoff=1, spacing=1.0)
    lchs = linear_qode("lchs", plan=quadrature, hamiltonian_function=hamiltonian_function)
    schrodinger = linear_qode(
        "schrodingerization",
        plan=SchrodingerPlan(auxiliary_width=2, period=8.0, selected_index=1),
        hamiltonian_function=hamiltonian_function,
    )

    # 1. Directly provide the BE of G; G=-I. work=3 is reserved for the later QRAM initial-state implementation.
    initial = abstract_state_prep("Initial", 1, work_width=3)
    initial_gate = gate_state_prep([1, 1], work_width=3)
    generator = abstract_block_encoding("Generator", 1, 0, 1.0)
    state = lchs(generator, initial, time)
    records.append(
        write(
            "given_be_lchs",
            state,
            {
                "Generator": scale(-1, identity(1)).operation,
                "Initial": initial_gate.operation,
            },
        )
    )

    # 2. When A=L+iH is already given, use the parts entry directly instead of merging first and splitting again.
    dissipation = abstract_block_encoding("Dissipation", 1, 0, 1.0)
    hamiltonian = abstract_block_encoding("Hamiltonian", 1, 1, 1.0)
    state = lchs_qode(
        LinearODE(HermitianParts(dissipation, hamiltonian), initial),
        time,
        plan=quadrature,
        hamiltonian_function=hamiltonian_function,
    )
    records.append(
        write(
            "given_parts_lchs",
            state,
            {
                "Dissipation": identity(1).operation,
                "Hamiltonian": zero(1).operation,
                "Initial": initial_gate.operation,
            },
        )
    )

    # 3. Same open graph, bound to gate / QRAM databases and state preparations respectively.
    # The data is an integer encoding of rotation angles; A=diag(1, cos(pi/4)), G=-A.
    angles = abstract_database("DiagonalAngles", 1, 2)
    a = diagonal_block_encoding(angles, alpha=1.0)
    state = lchs(scale(-1, a), initial, time)
    records.append(
        write(
            "given_xor_gate_lchs",
            state,
            {
                "DiagonalAngles": gate_database(1, 2, [0, 1]).operation,
                "Initial": initial_gate.operation,
            },
        )
    )
    records.append(
        write(
            "given_xor_qram_lchs",
            state,
            {
                "DiagonalAngles": Binding(
                    qram_database(1, 2).operation, {"table": "diagonal_angles"}
                ),
                "Initial": Binding(qram_state_prep(1, 2).operation, {"angles": "initial_angles"}),
            },
            memory={
                "diagonal_angles": [0, 1],
                "initial_angles": qram_state_angles([1, 1], 2),
            },
        )
    )

    # 4. A=[[1,-.5],[-.5,1]]: real symmetric, non-negative diagonal, and A>=0.
    fmt = FixedFormat(4, 1)
    access = abstract_sparse_access("SparseA", 1, fmt.width, sparsity=2)
    a = real_symmetric_sparse_encoding(access, fmt, 1.0, diagonal_nonnegative=True)
    state = lchs(scale(-1, a), initial, time)
    values = [fmt.encode(v) for v in (1, -0.5, -0.5, 1)]
    records.append(
        write(
            "given_sparse_gate_lchs",
            state,
            {
                "SparseA_position": sparse_location_gate(1, [(0, 1), (0, 1)]),
                "SparseA_entry": sparse_entry(gate_database(2, fmt.width, values), 1),
                "Initial": initial_gate.operation,
            },
            notes={"sparsity": 2, "amax": 1.0, "alpha_A": a.alpha},
        )
    )
    records.append(
        write(
            "given_sparse_qram_lchs",
            state,
            {
                "SparseA_position": Binding(
                    sparse_location_qram(1),
                    {
                        "forward": "positions",
                        "inverse": "inverse_positions",
                    },
                ),
                "SparseA_entry": Binding(
                    sparse_entry(qram_database(2, fmt.width), 1),
                    {
                        "db__table": "entries",
                    },
                ),
                "Initial": initial_gate.operation,
            },
            memory={
                "positions": [0, 0, 1, 1],
                "inverse_positions": [0, 0, 1, 1],
                "entries": values,
            },
        )
    )

    # 5. Periodic four-point heat equation; the structured shift BE serves directly as the spatial discretization result.
    grid = Grid(("x",), (4,), (1.0,), boundary="periodic")
    heat = DiscretePDE(
        scale(0.1, derivative_encoding(grid, (("x", 2),))),
        gate_state_prep([1, 0, 0, 0]),
        label="periodic_heat",
    )
    # QODEProtocol.__call__ has the same shape as the three-argument solver callback; declare solver by the callback type
    # so the return value of shifted_solver below can share the same loop variable as the protocol instances.
    for method, solver in (
        ("lchs", cast("Callable[[BlockEncoding, StatePreparation, float], StateOracle]", lchs)),
        ("schrodingerization", schrodinger),
    ):
        records.append(write("heat_" + method, make_qpde(solver)(heat, time)))

    # 6. Ordinary PDE -> F_p ports -> Carleman -> replaceable linear solver.
    u = Field("u")
    pde = PolynomialPDE.from_equations({"u": 0.1 * u.d("x", 2) - u * u.d("x")})
    space = Discretization(pde, grid)
    ports = structured_fd_bindings(space, space.encode_fields({"u": [0.1, 0.2, 0, -0.1]}))
    concrete = polynomial_from_bindings(ports)
    coefficient_slots, bindings = [], {}
    for order, coefficient in concrete.coefficients:
        name = "BurgersF" + str(order)
        slot = abstract_block_encoding(
            name, coefficient.width, coefficient.signal_qubits, coefficient.alpha
        )
        coefficient_slots.append((order, slot))
        bindings[name] = coefficient.operation
    initial_slot = abstract_state_prep(
        "BurgersInitial", concrete.width, concrete.initial.work_width
    )
    bindings["BurgersInitial"] = concrete.initial.operation
    problem = PolynomialODE(
        concrete.width, tuple(coefficient_slots), initial_slot, concrete.initial_norm
    )
    recovery: dict[str, float] = {}
    for method, solver in (
        ("schrodingerization", schrodinger),
        ("shifted_lchs", shifted_solver(lchs, recovery)),
    ):
        nonlinear_solver = partial(carleman_qode, linear_solver=solver, cutoff=2)
        state = qpde_solver(nonlinear_solver)(PDEInput(problem, "burgers"), time)
        records.append(
            write(
                "burgers_carleman_" + method,
                state,
                bindings,
                notes={
                    "initial_norm": concrete.initial_norm,
                    "cutoff": 2,
                    **(recovery if method == "shifted_lchs" else {}),
                },
            )
        )

    # 7. Same QPDE, independently replacing the configuration of the internal Hamiltonian-function protocol.
    degree_two = linear_qode(
        "lchs", plan=quadrature, hamiltonian_function=partial(taylor_hamiltonian, degree=2)
    )
    records.append(write("heat_lchs_taylor2", make_qpde(degree_two)(heat, time)))
    (args.output / "index.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    if args.native_bindings:
        verify_bindings(args.output)


if __name__ == "__main__":
    main()
