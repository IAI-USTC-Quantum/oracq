"""Input-model variant tutorial: four solver families × structured / spectral-embedding / QRAM data paths.

Run: PYTHONPATH=src python examples/input_models.py
The artifacts record structural evidence only (open slots, resources, closure and export); solution accuracy is a later validation item.
Each family shows different input models on the same problem instance: structured gate ports (baseline),
spectral embedding (exact Pauli expansion), and variants whose data path goes through open angle databases/QRAM resources.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path
from typing import cast

from oracq import (
    Binding,
    FixedFormat,
    Operation,
    bind,
    dump_qram_yaml,
    dumps,
    export_originir,
    export_toffoli_u3_cz,
    loads,
    scale,
    unresolved,
)
from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.common.prepare_select import lcu_prepare_select, qram_prepare
from oracq.algorithms.input_model.block_encoding import lcu, matrix_pauli_encoding
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    abstract_block_encoding,
    abstract_sparse_access,
    abstract_state_prep,
    gate_database,
    gate_state_prep,
    qram_database,
    qram_state_angles,
    qram_state_prep,
    sparse_entry,
    sparse_location_qram,
)
from oracq.algorithms.input_model.sparse import real_symmetric_sparse_encoding
from oracq.algorithms.qnlss.carleman import PolynomialODE, carleman_qode
from oracq.algorithms.qode.lchs import QuadraturePlan
from oracq.algorithms.qode.ode import linear_qode
from oracq.algorithms.qpde.pde import DiscretePDE, make_qpde
from oracq.applications.qham import (
    Discretization,
    Field,
    Grid,
    Known,
    PolynomialPDE,
    QHAMBindings,
    QHAMPlan,
    gate_bindings,
    qham_input_model,
    qram_coefficient_encoding,
    qram_coefficient_memory,
    structured_fd_bindings,
)


def save_case(
    root: Path,
    name: str,
    state: StateOracle,
    bindings: Mapping[str, Operation | Binding] | None = None,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None = None,
    *,
    notes: Mapping[str, object] | None = None,
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
    record = {
        "case": name,
        "open_slots": [r.name for r in unresolved(opened)],
        "modules": len(closed.modules),
        "resources": [r.name for r in closed.main.resources],
        "target_width": state.width,
        "signal_width": state.signal_qubits,
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


def reopen_coefficients(
    problem: PolynomialODE, label: str
) -> tuple[PolynomialODE, dict[str, Operation | Binding]]:
    """Re-declare the concrete coefficient block encodings as open slots; returns (PolynomialODE, bindings)."""
    coefficient_slots, bindings = [], {}
    for order, coefficient in problem.coefficients:
        name = label + "F" + str(order)
        slot = abstract_block_encoding(
            name, coefficient.width, coefficient.signal_qubits, coefficient.alpha
        )
        coefficient_slots.append((order, slot))
        bindings[name] = coefficient.operation
    initial_slot = abstract_state_prep(
        label + "Initial", problem.width, problem.initial.work_width
    )
    bindings[label + "Initial"] = problem.initial.operation
    return (
        PolynomialODE(problem.width, tuple(coefficient_slots), initial_slot, problem.initial_norm),
        # Dict invariance: the actual value type is Operation; pass it out under the declared union type.
        cast("dict[str, Operation | Binding]", bindings),
    )


def main() -> None:
    """Assemble input-model variants for the four solver families and write the comparison artifacts and index to disk."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path, default=Path("out/input-models"))
    args = parser.parse_args()
    write = partial(save_case, args.output)
    records = []
    time = 0.05
    kernel = partial(taylor_hamiltonian, degree=1)
    quadrature = QuadraturePlan.cauchy(cutoff=1, spacing=1.0)
    lchs = linear_qode("lchs", plan=quadrature, hamiltonian_function=kernel)
    schrodinger = linear_qode("schrodingerization", hamiltonian_function=kernel)
    cbmd = linear_qode("cbmd", hamiltonian_function=kernel)

    # ---- QHAM family: forced Burgers (order 2, eta=-0.4), three input models ----------
    u = Field("u")
    pde = PolynomialPDE.from_equations(
        {"u": 0.1 * u.d("x", 2) - u * u.d("x") + Known("forcing")},
        label="forced_burgers",
    )
    plan = QHAMPlan(pde, order=2)
    grid = Grid(("x",), (4,), (1.0,), boundary="periodic")
    space = Discretization(pde, grid, {"forcing": [0.05, 0, -0.05, 0]})

    # 1. Structured gate ports + gate initial state: baseline (same input as examples/general_qham.py).
    initial = space.encode_fields({"u": [0.1, 0.2, 0, -0.1]})
    stencil_bindings = structured_fd_bindings(space, initial)
    model = qham_input_model(plan, stencil_bindings, eta=-0.4)
    records.append(write("qham_stencil_gate", model.solve(schrodinger, 0.01)))

    # 2. Spectral-embedding ports: exact Pauli-LCU expansion of the base port matrices + gate initial state.
    spectral_bindings = gate_bindings(space, initial)
    model = qham_input_model(plan, spectral_bindings, eta=-0.4)
    records.append(write("qham_spectral_ports", model.solve(schrodinger, 0.01)))

    # 3+4. Coefficient data stays in open angle-database slots and the initial state is declared as a slot: the same open
    # program is bound to a gate database and to a QRAM database respectively (the runtime angle tables are provided separately).
    # QRAM angle tables currently accept only non-negative amplitudes, so this pair of cases uses a non-negative profile.
    profile = [0.1, 0.2, 0.15, 0.05]
    qinitial = space.encode_fields({"u": profile})
    encoder = partial(qram_coefficient_encoding, angle_width=8)
    open_bindings = structured_fd_bindings(space, qinitial, coefficient_encoder=encoder)
    open_bindings = QHAMBindings(
        open_bindings.state_width,
        open_bindings.ports,
        abstract_state_prep("QhamInitial", open_bindings.state_width, 10),
        open_bindings.initial_norm,
    )
    model = qham_input_model(plan, open_bindings, eta=-0.4)
    state = model.solve(schrodinger, 0.01)
    coeff_slots = [
        r.name for r in unresolved(state.operation.program()) if r.name != "QhamInitial"
    ]
    assert len(coeff_slots) == 1, coeff_slots
    (db_slot,) = coeff_slots
    (forcing_monomial,) = (
        term.monomial for port in space.pde.ports for term in port.terms if term.monomial.known
    )
    angle_words = qram_coefficient_memory(space, forcing_monomial, angle_width=8)
    word_table = [angle_words.get(address, 0) for address in range(4)]
    records.append(
        write(
            "qham_coeff_gate",
            state,
            {
                db_slot: gate_database(2, 8, word_table).operation,
                "QhamInitial": gate_state_prep(profile, work_width=10).operation,
            },
        )
    )
    qram_memory: dict[str, Sequence[int] | Mapping[int, int]] = {
        "coeff_angles": word_table,
        "initial_angles": qram_state_angles(profile, 8),
    }
    qram_dict = {
        db_slot: Binding(qram_database(2, 8).operation, {"table": "coeff_angles"}),
        "QhamInitial": Binding(qram_state_prep(2, 8).operation, {"angles": "initial_angles"}),
    }
    records.append(write("qham_coeff_qram", state, qram_dict, memory=qram_memory))

    # ---- Carleman family: inviscid Burgers, cutoff=2, three input models -------------
    pde_c = PolynomialPDE.from_equations({"u": 0.1 * u.d("x", 2) - u * u.d("x")})
    space_c = Discretization(pde_c, grid)
    initial_c = space_c.encode_fields({"u": profile})

    # 5. Structured gate ports (slots re-abstracted, same pattern as case 6 in ode_input_models).
    concrete = polynomial_from_bindings(structured_fd_bindings(space_c, initial_c))
    problem, rebind = reopen_coefficients(concrete, "Burgers")
    records.append(
        write(
            "carleman_stencil_gate",
            carleman_qode(problem, time, schrodinger, cutoff=2),
            rebind,
        )
    )

    # 6. Spectral-embedding ports (Pauli-LCU base port matrices).
    concrete = polynomial_from_bindings(gate_bindings(space_c, initial_c))
    problem, rebind = reopen_coefficients(concrete, "Spectral")
    records.append(
        write(
            "carleman_spectral_ports",
            carleman_qode(problem, time, schrodinger, cutoff=2),
            rebind,
        )
    )

    # 7. Spatially varying viscosity: coefficients through an open angle database + QRAM initial state.
    pde_nu = PolynomialPDE.from_equations(
        {"u": 0.1 * Known("nu") * u.d("x", 2) - u * u.d("x")}, label="burgers_nu"
    )
    space_nu = Discretization(pde_nu, grid, {"nu": [1.0, 0.8, 0.6, 1.2]})
    concrete = polynomial_from_bindings(
        structured_fd_bindings(space_nu, initial_c, coefficient_encoder=encoder)
    )
    problem, rebind = reopen_coefficients(concrete, "QramBurgers")
    # The initial-state slot is declared with the register layout of the QRAM preparation (the gate implementation has work width 0).
    del rebind["QramBurgersInitial"]
    qinit_slot = abstract_state_prep("QramBurgersInitial", concrete.width, 10)
    problem = PolynomialODE(
        concrete.width, problem.coefficients, qinit_slot, concrete.initial_norm
    )
    extra = [
        r.name
        for _, coefficient in concrete.coefficients
        for r in unresolved(coefficient.operation.program())
    ]
    assert len(extra) == 1, extra
    (nu_slot,) = extra
    (nu_monomial,) = (
        term.monomial for port in space_nu.pde.ports for term in port.terms if term.monomial.known
    )
    nu_words = qram_coefficient_memory(space_nu, nu_monomial, angle_width=8)
    rebind[nu_slot] = Binding(qram_database(2, 8).operation, {"table": "nu_angles"})
    rebind["QramBurgersInitial"] = Binding(
        qram_state_prep(2, 8).operation, {"angles": "initial_angles"}
    )
    state = carleman_qode(problem, time, schrodinger, cutoff=2)
    records.append(
        write(
            "carleman_coeff_qram",
            state,
            rebind,
            memory={
                "nu_angles": [nu_words.get(address, 0) for address in range(4)],
                "initial_angles": qram_state_angles(profile, 8),
            },
        )
    )

    # ---- LCHS / CBMD: spectral embedding of a given generator (same matrix as the given_sparse case) ------
    a_matrix = [[1, -0.5], [-0.5, 1]]
    spectral_be = matrix_pauli_encoding(a_matrix)
    initial_slot = abstract_state_prep("SpectralInitial", 1, work_width=3)
    initial_gate = gate_state_prep([1, 1], work_width=3)
    gen_slot = abstract_block_encoding(
        "SpectralGenerator", 1, spectral_be.signal_qubits, spectral_be.alpha
    )
    for family, solver in (("lchs", lchs), ("cbmd", cbmd)):
        state = solver(scale(-1, gen_slot), initial_slot, time)
        records.append(
            write(
                family + "_given_spectral",
                state,
                {
                    "SpectralGenerator": scale(-1, spectral_be).operation,
                    "SpectralInitial": initial_gate.operation,
                },
            )
        )

    # Spectral embedding + QRAM PREPARE: Pauli coefficients loaded from a QRAM angle table. A = I - 0.5 X.
    terms = [(1.0, "I"), (-0.5, "X")]
    qprep = qram_prepare([1.0, -0.5], angle_width=8)
    qbe = lcu_prepare_select(terms, prepare=qprep.preparation)
    qgen_slot = abstract_block_encoding(
        "QramSpectralGenerator", 1, qbe.signal_qubits, qbe.alpha
    )
    qgen_op = scale(-1, qbe).operation
    (angles_formal,) = [r.name for r in qgen_op.program().main.resources]
    for family, solver in (("lchs", lchs), ("cbmd", cbmd)):
        state = solver(scale(-1, qgen_slot), initial_slot, time)
        records.append(
            write(
                family + "_spectral_qram_prepare",
                state,
                {
                    "QramSpectralGenerator": Binding(
                        qgen_op, {angles_formal: "spectral_angles"}
                    ),
                    "SpectralInitial": initial_gate.operation,
                },
                memory={"spectral_angles": qprep.memory["angles"]},
            )
        )

    # ---- LCHS / CBMD: QRAM grid discretization ---------------------------------------
    # Sparse access of the four-point periodic-ring Laplacian (location + entry tables from QRAM) + QRAM initial state.
    fmt = FixedFormat(4, 1)
    access = abstract_sparse_access("RingStencil", 2, fmt.width, sparsity=3)
    laplacian = real_symmetric_sparse_encoding(access, fmt, 2.0, diagonal_nonnegative=True)
    heat = DiscretePDE(
        scale(-0.1, laplacian),
        abstract_state_prep("GridInitial", 2, 10),
        label="ring_heat",
    )
    ring = [
        [2 if row == col else (-1 if (row - col) % 4 in (1, 3) else 0) for col in range(4)]
        for row in range(4)
    ]
    slots = {col: [row for row in range(4) if ring[row][col]] for col in range(4)}
    forward, inverse = [0] * 16, [0] * 16
    for col in range(4):
        for slot, row in enumerate(slots[col]):
            forward[col + slot * 4] = row
            inverse[col + row * 4] = slot
    entries = [fmt.encode(ring[row][col]) for col in range(4) for row in range(4)]
    grid_bindings = {
        "RingStencil_position": Binding(
            sparse_location_qram(2),
            {"forward": "positions", "inverse": "inverse_positions"},
        ),
        "RingStencil_entry": Binding(
            sparse_entry(qram_database(4, fmt.width), 2),
            {"db__table": "entries"},
        ),
        "GridInitial": Binding(qram_state_prep(2, 8).operation, {"angles": "grid_angles"}),
    }
    grid_memory: dict[str, Sequence[int] | Mapping[int, int]] = {
        "positions": forward,
        "inverse_positions": inverse,
        "entries": entries,
        "grid_angles": qram_state_angles([1, 0, 0, 0], 8),
    }
    for family, solver in (("lchs", lchs), ("cbmd", cbmd)):
        state = make_qpde(solver)(heat, time)
        records.append(
            write(family + "_qram_grid", state, grid_bindings, memory=grid_memory)
        )

    # ---- CBMD: QHAM open-coefficient input through a dissipative shift (input model reused across families) ------------
    model = qham_input_model(plan, open_bindings, eta=-0.4)
    state = model.dissipative_shift().solve(cbmd, 0.01)
    records.append(
        write("cbmd_qham_coeff_qram", state, qram_dict, memory=qram_memory)
    )

    (args.output / "index.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
