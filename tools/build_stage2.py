"""重建第二阶段开放/闭合 RIR 与模块化后端描述，产物一律放 out/。"""

from __future__ import annotations

import argparse
import json
from functools import partial
from pathlib import Path

from pyqecclang import bind, dumps, identity, scale
from pyqecclang.algorithms.arithmetic import FixedFormat, fixed_arithmetic
from pyqecclang.algorithms.block_encoding import matrix_pauli_encoding, pad_signal
from pyqecclang.algorithms.carleman import PolynomialODE, carleman_qode
from pyqecclang.algorithms.cbmd import ContourPlan
from pyqecclang.algorithms.hamiltonian import taylor_hamiltonian
from pyqecclang.algorithms.lchs import QuadraturePlan
from pyqecclang.algorithms.ode import linear_qode
from pyqecclang.algorithms.oracles import (
    abstract_block_encoding,
    abstract_state_prep,
    basis_state,
    gate_state_prep,
)
from pyqecclang.algorithms.pde import PDEInput, make_qpde, qpde_solver
from pyqecclang.algorithms.qlss import CostaConfig, SpectralPromise, make_costa_qlss
from pyqecclang.algorithms.schrodingerization import SchrodingerPlan
from pyqecclang.applications.flow_data import RoeFlowData
from pyqecclang.applications.legacy import qham_initial_vector, qham_m1
from pyqecclang.applications.qfvm import (
    bind_qfvm,
    qfvm_memories,
    roe_entry,
    roe_qfvm_block_encoding,
    roe_qfvm_inputs,
    roe_qfvm_step,
)
from pyqecclang.applications.roe import roe_face
from pyqecclang.infrastructure.backends.basis import export_toffoli_u3_cz
from pyqecclang.infrastructure.backends.originir import export_originir
from pyqecclang.infrastructure.layout import workspace_table
from pyqecclang.infrastructure.linking import unresolved


def save_case(root, name, opened, closed=None, memory=None):
    path = root / name
    path.mkdir(parents=True, exist_ok=True)
    (path / "open.rir.json").write_text(dumps(opened))
    requirements = [
        {"name": r.name, "paradigm": r.paradigm, "path": r.path} for r in unresolved(opened)
    ]
    report = {
        "name": name,
        "modules": len(opened.modules),
        "requirements": requirements,
        "correctness": "pending",
        "evidence": "construction and serialization",
        "private_workspace": workspace_table(opened)[opened.entry],
    }
    if closed is not None:
        (path / "closed.rir.json").write_text(dumps(closed))
        normal = export_originir(closed)
        strict = export_toffoli_u3_cz(closed)
        (path / "modular.originir").write_text(normal.text)
        (path / "toffoli_u3_cz.originir").write_text(strict.text)
        report.update(
            closed_modules=len(closed.modules),
            definitions=strict.text.count("DEF "),
            strict_bytes=len(strict.text.encode()),
            gate_alphabet=["TOFFOLI", "U3", "CZ"],
            qram_resources=list(strict.resources),
            strict_workspace_qubits=len(strict.workspace_qubits),
            backend_evidence="text emitted; native/parse evidence recorded separately",
        )
    if memory is not None:
        (path / "memory.json").write_text(json.dumps(memory, indent=2) + "\n")
    (path / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(
        name,
        report.get("closed_modules", len(opened.modules)),
        report.get("strict_bytes", 0),
        flush=True,
    )
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("out/stage2"))
    args = parser.parse_args()
    root = args.out
    reports = []
    fmt = FixedFormat(4, 1)
    for kind in (
        "add",
        "sub",
        "neg",
        "abs",
        "mul",
        "div",
        "reciprocal",
        "sqrt",
        "lt",
        "eq",
        "select",
        "and",
        "or",
        "xor",
    ):
        p = fixed_arithmetic(kind, fmt).program()
        reports.append(save_case(root, "arithmetic_" + kind, p, p))
    flow_format = FixedFormat(6, 2)
    p = roe_face(fmt=flow_format).program()
    reports.append(save_case(root, "roe_face", p, p))
    inputs = roe_qfvm_inputs(fmt=flow_format, angle_width=6)
    flow = RoeFlowData(
        [(1, 0, 2), (1.25, 0.125, 2.5), (1, 0, 2), (0.75, -0.125, 1.5)],
        fmt=flow_format,
        angle_width=inputs.angle_width,
    )
    memory = qfvm_memories(inputs, flow, amax=4.0)
    for name, op in (
        ("qfvm_roe_entry", roe_entry(inputs)),
        ("qfvm_roe_be", roe_qfvm_block_encoding(inputs, amax=4).operation),
        (
            "qfvm_costa_filter",
            roe_qfvm_step(
                inputs,
                make_costa_qlss(CostaConfig(steps=1, filter_degree=2)),
                spectrum=SpectralPromise(10, 0.5),
                rhs_norm=1,
                amax=4,
            ).operation,
        ),
    ):
        opened = op.program()
        closed = bind_qfvm(opened, inputs)
        selected = {r.name: memory[r.name] for r in closed.main.resources}
        reports.append(save_case(root, name, opened, closed, selected))
    patch = flow.update({1: (1.5, 0.5, 3.0)})
    (root / "qfvm_patch.json").write_text(
        json.dumps(
            {
                "version": patch.version,
                "changes": patch.changes,
                "faces": patch.recomputed_faces,
                "cells": patch.recomputed_cells,
                "native_update": "rebuild changed banks only",
            },
            indent=2,
        )
        + "\n"
    )

    concrete_a = pad_signal(scale(-0.5, identity(1)), 1)
    a = abstract_block_encoding("DifferentialA", 1, 1, concrete_a.alpha)
    initial = abstract_state_prep("DifferentialInitial", 1)
    bindings = {
        a.operation.module.name: concrete_a.operation,
        initial.operation.module.name: basis_state(1).operation,
    }
    hf = partial(taylor_hamiltonian, degree=1)
    method_options = {
        "lchs": {"plan": QuadraturePlan.cauchy(cutoff=1)},
        "cbmd": {"plan": ContourPlan(cutoff=1)},
        "schrodingerization": {"plan": SchrodingerPlan(auxiliary_width=2)},
    }
    for method, options in method_options.items():
        qode = linear_qode(method, hamiltonian_function=hf, **options)
        state = qode(a, initial, 0.1)
        p = state.operation.program()
        reports.append(save_case(root, method + "_qode", p, bind(p, bindings)))
        # PDE 入口只负责空间离散化到开放算子；同一 oracle 图可来自非矩阵输入。
        from pyqecclang.algorithms.pde import DiscretePDE

        state = make_qpde(qode)(DiscretePDE(a, initial, "heat_equation_open_space"), 0.1)
        p = state.operation.program()
        reports.append(save_case(root, method + "_qpde", p, bind(p, bindings)))

    actual_f2 = matrix_pauli_encoding(
        [[0, 0.125, 0, 0], [0, 0, 0.125, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    )
    f2 = abstract_block_encoding("DifferentialF2", 2, actual_f2.signal_qubits, actual_f2.alpha)
    polynomial = PolynomialODE(1, ((1, a), (2, f2)), initial)
    carleman = partial(
        carleman_qode,
        linear_solver=linear_qode("cbmd", hamiltonian_function=hf, plan=ContourPlan(cutoff=0)),
        cutoff=2,
    )
    more = {**bindings, f2.operation.module.name: actual_f2.operation}
    for name, state in (
        ("carleman_qode", carleman(polynomial, 0.1)),
        (
            "carleman_qpde",
            qpde_solver(carleman)(PDEInput(polynomial, "burgers_polynomial_input"), 0.1),
        ),
    ):
        p = state.operation.program()
        reports.append(save_case(root, name, p, bind(p, more)))

    qham_prep = gate_state_prep(qham_initial_vector([1.0, 0.0]))
    for method in ("cbmd", "schrodingerization"):
        qode = linear_qode(
            method,
            hamiltonian_function=hf,
            **({"plan": ContourPlan(cutoff=0)} if method == "cbmd" else {}),
        )
        state = qham_m1(a, f2, qham_prep, make_qpde(qode), via_pde=True)
        p = state.operation.program()
        available = {r.name for r in unresolved(p)}
        reports.append(
            save_case(
                root,
                "qham_" + method + "_qpde",
                p,
                bind(p, {k: v for k, v in more.items() if k in available}),
            )
        )
    (root / "index.json").write_text(json.dumps(reports, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
