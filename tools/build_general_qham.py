"""自动生成一般 QHAM 的推导、数学见证、开放/闭合量子描述。"""

import json
import math
import random
from functools import partial
from pathlib import Path

from pyqecclang import bind, dumps, export_originir, export_toffoli_u3_cz, unresolved
from pyqecclang.algorithms.cbmd import ContourPlan
from pyqecclang.algorithms.hamiltonian import taylor_hamiltonian
from pyqecclang.algorithms.ode import linear_qode
from pyqecclang.algorithms.schrodingerization import SchrodingerPlan
from pyqecclang.applications.qham import (
    Discretization,
    Grid,
    QHAMBindings,
    QHAMPlan,
    qham_input_model,
    structured_fd_bindings,
    taylor_qode,
)
from pyqecclang.applications.qham.examples import example_pde
from pyqecclang.applications.qham.report import export_derivation


def main():
    root = Path("out/qham-general")
    root.mkdir(exist_ok=True)
    cases = (
        ("burgers", 2, (4,)),
        ("kdv", 3, (4,)),
        ("reaction", 2, (2,)),
        ("coupled", 2, (2,)),
        ("vector_burgers_2d", 2, (4, 4)),
    )
    reports = []
    for name, order, shape in cases:
        pde = example_pde(name)
        grid = Grid(pde.axes, shape, (1.0,) * len(shape))
        known = (
            {"f": [0.03 * math.cos(2 * math.pi * i / grid.size) for i in range(grid.size)]}
            if any(p.name == "F" for p in pde.ports)
            else {}
        )
        disc = Discretization(pde, grid, known)
        initial = disc.encode_fields(
            {
                field: [
                    0.1 * math.cos(2 * math.pi * i / grid.size + component * 0.3)
                    for i in range(grid.size)
                ]
                for component, field in enumerate(pde.fields)
            }
        )
        plan = QHAMPlan(pde, order)
        path = root / (name + "_m" + str(order))
        export_derivation(plan, path, state_width=disc.width, eta=-0.4)
        rng = random.Random(17)
        values = [[rng.uniform(-0.1, 0.1) for _ in range(disc.dimension)] for _ in range(order + 1)]
        lifted = disc.lift(plan, values)
        direct = disc.chain_rule(plan, values, -0.4)
        action = disc.linear_action(plan, -0.4, lifted)
        residual = max(abs(a - b) for a, b in zip(action, direct, strict=True))
        numerical = structured_fd_bindings(disc, initial)
        specs = {
            key: (port.encoding.alpha, port.encoding.signal_qubits) for key, port in numerical.ports
        }
        declared = QHAMBindings.declare(
            plan, disc.width, specs, initial_norm=numerical.initial_norm, prefix="Input_" + name
        )
        model = qham_input_model(plan, declared, eta=-0.4)
        (path / "generator-open.rir.json").write_text(dumps(model.generator.operation.program()))
        (path / "initial-open.rir.json").write_text(dumps(model.initial.operation.program()))
        (path / "bindings.json").write_text(
            json.dumps(
                {
                    key: {
                        "alpha": port.encoding.alpha,
                        "signal": port.encoding.signal_qubits,
                        "arity": port.arity,
                    }
                    for key, port in numerical.ports
                },
                indent=2,
            )
            + "\n"
        )
        solution = model.solve(partial(taylor_qode, degree=1), 0.01)
        mapping = {
            dict(declared.ports)[key].encoding.operation.module.name: port.encoding.operation
            for key, port in numerical.ports
        }
        mapping[declared.initial.operation.module.name] = numerical.initial.operation
        opened = solution.operation.program()
        closed = bind(opened, mapping)
        (path / "solution-open.rir.json").write_text(dumps(opened))
        (path / "solution-closed.rir.json").write_text(dumps(closed))
        (path / "solution.originir").write_text(export_originir(closed).text)
        strict = export_toffoli_u3_cz(closed)
        (path / "toffoli_u3_cz.originir").write_text(strict.text)
        report = {
            "case": name,
            "order": order,
            "degree": pde.degree,
            "fields": pde.fields,
            "axes": pde.axes,
            "state_dimension": disc.dimension,
            "raw_dimension": len(lifted),
            "blocks": plan.block_count,
            "chain_rule_residual": residual,
            "generator_alpha": model.generator.alpha,
            "log_initial_norm": model.log_initial_norm,
            "rir_modules": len(closed.modules),
            "open_oracles": [r.name for r in unresolved(opened)],
            "strict_bytes": len(strict.text.encode()),
            "port_implementation": "structured shifts/component projections/contractions; no dense matrix",
            "global_matrix_materialized": False,
            "solver_accuracy": "pending",
        }
        if name == "burgers":
            alternatives = {}
            hf = partial(taylor_hamiltonian, degree=1)
            for method in ("schrodingerization", "cbmd"):
                effective = model.dissipative_shift() if method == "cbmd" else model
                options = (
                    {"plan": ContourPlan(cutoff=0)}
                    if method == "cbmd"
                    else {"plan": SchrodingerPlan(auxiliary_width=2, selected_index=1)}
                )
                result = effective.solve(
                    linear_qode(method, hamiltonian_function=hf, **options), 0.01
                )
                program = bind(result.operation.program(), mapping)
                (path / (method + ".rir.json")).write_text(dumps(program))
                (path / (method + ".originir")).write_text(export_originir(program).text)
                alternatives[method] = {
                    "modules": len(program.modules),
                    "growth_shift": effective.growth_shift,
                    "growth_rescale_log": effective.growth_shift * 0.01,
                    "accuracy": "pending",
                }
            report["qode_alternatives"] = alternatives
        (path / "verification.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=2) + "\n"
        )
        reports.append(report)
        print(
            name,
            {
                "dimension": len(lifted),
                "blocks": plan.block_count,
                "residual": residual,
                "modules": len(closed.modules),
            },
            flush=True,
        )
    (root / "index.json").write_text(json.dumps(reports, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
