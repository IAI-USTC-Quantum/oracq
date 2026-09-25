"""Record of real-backend consumption of the pure-function compilation artifacts."""

import json
from pathlib import Path

from oracq import FixedFormat, arithmetic_native_registry, run_pysparq
from oracq.applications.flow_data import RoeFlowData
from oracq.applications.qfvm import bind_qfvm, qfvm_memories, roe_entry, roe_qfvm_inputs


def main():
    from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser

    root = Path("out/math-functions")
    records = []
    for name in ("pressure", "guarded_reciprocal", "polynomial"):
        parser = OriginIR_BaseParser()
        parser.parse((root / name / "toffoli_u3_cz.originir").read_text())
        records.append(
            {
                "case": name,
                "evidence": "real OriginIR_BaseParser",
                "qubits": parser.n_qubit,
                "status": "passed",
            }
        )
        print(name, "parsed", flush=True)
    inputs = roe_qfvm_inputs(fmt=FixedFormat(6, 2), angle_width=6)
    flow = RoeFlowData(
        [(1, 0, 2), (1.25, 0.125, 2.5), (1, 0, 2), (0.75, -0.125, 1.5)],
        fmt=inputs.fmt,
        angle_width=6,
    )
    program = bind_qfvm(roe_entry(inputs).program(), inputs)
    memory = qfvm_memories(inputs, flow)
    memory = {r.name: memory[r.name] for r in program.main.resources}
    report = {}
    state = run_pysparq(
        program, memory, native_registry=arithmetic_native_registry(program), report=report
    )
    records.append(
        {
            "case": "qfvm_entry_compiled_roe",
            "evidence": "real PySparQ, QRAM and compiled arithmetic",
            "execution": report,
            "states": len(state.amplitudes),
            "status": "passed",
            "numeric_correctness": "pending",
        }
    )
    (root / "backend-validation.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2) + "\n"
    )
    print("qfvm_entry", report, flush=True)


if __name__ == "__main__":
    main()
