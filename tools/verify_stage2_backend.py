"""对阶段产物使用真实 UnifiedQuantum 解析器；不运行大规模状态向量。"""

import json
from pathlib import Path

from pyqecclang import Program, loads
from pyqecclang.backends.basis import export_toffoli_u3_cz


def main():
    from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser

    root = Path("out/stage2")
    records = []
    names = (
        "arithmetic_add",
        "arithmetic_mul",
        "arithmetic_div",
        "arithmetic_sqrt",
        "lchs_qode",
        "cbmd_qode",
        "schrodingerization_qode",
        "carleman_qode",
        "qham_cbmd_qpde",
        "qham_schrodingerization_qpde",
    )
    for name in names:
        text = (root / name / "toffoli_u3_cz.originir").read_text()
        parser = OriginIR_BaseParser()
        parser.parse(text)
        records.append(
            {
                "case": name,
                "backend": "real OriginIR_BaseParser",
                "status": "passed",
                "qubits": parser.n_qubit,
                "correctness": "pending",
            }
        )
        print(name, "parsed", parser.n_qubit, flush=True)
    closed = loads((root / "qfvm_roe_be" / "closed.rir.json").read_text())
    structural = next(
        m
        for m in closed.modules
        if dict(m.attributes).get("algorithm") == "cks_real_symmetric_isometry"
    )
    program = Program(structural.name, closed.modules)
    artifact = export_toffoli_u3_cz(program)
    parser = OriginIR_BaseParser()
    parser.parse(artifact.text)
    (root / "qfvm_sparse_T.originir").write_text(artifact.text)
    records.append(
        {
            "case": "qfvm_sparse_T",
            "backend": "real OriginIR_BaseParser",
            "status": "passed",
            "qubits": parser.n_qubit,
            "qram": True,
            "correctness": "pending",
        }
    )
    (root / "backend-validation.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
