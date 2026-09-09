"""同一 QFVM 稀疏问题切换 CKS / Costa 的可审阅产物。"""

import json
from pathlib import Path

from pyqecclang import FixedFormat, dumps, export_originir, export_toffoli_u3_cz, unresolved
from pyqecclang.algorithms.cks import CKSConfig, make_cks_qlss
from pyqecclang.algorithms.costa import CostaConfig, make_costa_qlss
from pyqecclang.flow_data import RoeFlowData
from pyqecclang.qfvm import bind_qfvm, qfvm_memories, roe_qfvm_inputs, roe_qfvm_problem
from pyqecclang.qlss import SpectralPromise


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
        (path / "open.rir.json").write_text(dumps(opened))
        (path / "closed.rir.json").write_text(dumps(closed))
        (path / "modular.originir").write_text(export_originir(closed).text)
        (path / "toffoli_u3_cz.originir").write_text(export_toffoli_u3_cz(closed).text)
        probe = bind_qfvm(result.norm_probe.operation.program(), inputs)
        (path / "norm-probe.rir.json").write_text(dumps(probe))
        (path / "memory.json").write_text(
            json.dumps({r.name: memory[r.name] for r in closed.main.resources}, indent=2) + "\n"
        )
        reports.append(
            {
                "solver": protocol.name,
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
        print(protocol.name, reports[-1], flush=True)
    (root / "comparison.json").write_text(json.dumps(reports, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
