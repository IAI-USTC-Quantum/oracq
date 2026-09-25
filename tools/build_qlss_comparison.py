"""Reviewable artifacts for switching the same QFVM sparse problem between CKS and Costa."""

import hashlib
import json
from pathlib import Path

from oracq import (
    FixedFormat,
    dump_qram_yaml,
    dumps,
    estimate_resources,
    export_originir,
    export_toffoli_u3_cz,
    unresolved,
)
from oracq.algorithms.qlss.qlss import (
    CKSConfig,
    CostaConfig,
    SpectralPromise,
    make_cks_qlss,
    make_costa_qlss,
)
from oracq.applications.flow_data import RoeFlowData
from oracq.applications.qfvm import bind_qfvm, qfvm_memories, roe_qfvm_inputs, roe_qfvm_problem


def main():
    root = Path("out/qlss-audit")
    root.mkdir(parents=True, exist_ok=True)
    inputs = roe_qfvm_inputs(fmt=FixedFormat(6, 2), angle_width=6)
    flow = RoeFlowData(
        [(1, 0, 2), (1.25, 0.25, 2.5), (1, 0, 2), (0.75, -0.25, 1.5)],
        fmt=inputs.fmt,
        angle_width=inputs.angle_width,
    )
    problem = roe_qfvm_problem(
        inputs, spectrum=SpectralPromise(10, 0.5), rhs_norm=flow.rhs_norm, amax=4
    )
    memory = qfvm_memories(inputs, flow, amax=4)
    reports = []
    for protocol in (make_cks_qlss(CKSConfig(order=2)), make_costa_qlss(CostaConfig(steps=1))):
        result = protocol(problem)
        path = root / protocol.name
        path.mkdir(exist_ok=True)
        opened = result.operation.program()
        closed = bind_qfvm(opened, inputs)
        (path / "open.rir.yaml").write_text(dumps(opened))
        (path / "closed.rir.yaml").write_text(dumps(closed))
        (path / "modular.originir").write_text(export_originir(closed).text)
        (path / "toffoli_u3_cz.originir").write_text(export_toffoli_u3_cz(closed).text)
        probe = bind_qfvm(result.norm_probe.operation.program(), inputs)
        (path / "norm-probe.rir.yaml").write_text(dumps(probe))
        memory_text = dump_qram_yaml(closed, {r.name: memory[r.name] for r in closed.main.resources})
        (path / "memory.qram.yaml").write_text(memory_text, encoding="utf-8")
        reports.append(
            {
                "solver": protocol.name,
                "open_cost": estimate_resources(opened, require_closed=False).to_dict(),
                "closed_cost": estimate_resources(closed).to_dict(),
                "program_sha256": hashlib.sha256(dumps(closed).encode()).hexdigest(),
                "memory_sha256": hashlib.sha256(memory_text.encode()).hexdigest(),
                "scope": "input adaptation, binding, and resource analysis; this script does not certify solution accuracy or success probability",
                "input_model": protocol.input_model,
                "adapter_trace": result.adapter_trace,
                "alpha": result.input_alpha,
                "encoded_inverse_bound": result.encoded_inverse_bound,
                "physical_width": result.state.width,
                "signal_width": result.state.signal_qubits,
                "rhs_norm": flow.rhs_norm,
                "closed_modules": len(closed.modules),
                "open_requirements": [r.name for r in unresolved(opened)],
                "kernel_status": result.kernel_status,
                "norm_recovery": "rhs_norm / (alpha * sqrt(p_joint / p_solver))",
            }
        )
        cost = estimate_resources(closed)
        print(protocol.name, "qubits", cost.qubits, "Toffoli", cost.toffoli, "QRAM", cost.qram_total, flush=True)
    (root / "comparison.json").write_text(json.dumps(reports, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
