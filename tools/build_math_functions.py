"""Rebuild the pure-function compilation cases; export the quantum implementations and MIR to out/ as well."""

import json
from pathlib import Path

from oracq import (
    FixedFormat,
    MathConfig,
    compile_function,
    dumps,
    export_originir,
    export_toffoli_u3_cz,
)
from oracq.applications.qfvm import bind_qfvm, roe_qfvm_block_encoding, roe_qfvm_inputs
from oracq.applications.roe_formulas import frozen_roe_face
from oracq.infrastructure.layout import workspace_table
from oracq.infrastructure.mathfunc import Index


def main():
    root = Path("out/math-functions")
    root.mkdir(parents=True, exist_ok=True)
    source = Path("examples/math_functions.py").read_text()
    records = []
    for name in ("pressure", "roe_speed", "phase_response", "guarded_reciprocal", "polynomial"):
        compiled = compile_function(
            source, entry=name, fmt=FixedFormat(8, 3), config=MathConfig(degree=3)
        )
        path = root / name
        path.mkdir(exist_ok=True)
        (path / "function.mir.json").write_text(compiled.math_ir.dumps())
        (path / "function.rir.yaml").write_text(dumps(compiled.program()))
        (path / "function.originir").write_text(export_originir(compiled.program()).text)
        strict = export_toffoli_u3_cz(compiled.program())
        (path / "toffoli_u3_cz.originir").write_text(strict.text)
        record = {
            "function": name,
            "math_modules": len(compiled.math_ir.functions),
            "rir_modules": len(compiled.program().modules),
            "workspace": workspace_table(compiled.program())[compiled.program().entry],
            "strict_bytes": len(strict.text.encode()),
            "correctness": "pending",
        }
        records.append(record)
        print(name, record["rir_modules"], record["strict_bytes"], flush=True)
    compiled = compile_function(
        frozen_roe_face,
        fmt=FixedFormat(6, 2),
        inputs={
            **{key: "real" for key in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r")},
            "row": Index(2),
            "col": Index(2),
        },
        constants={"gamma": 1.4, "entropy_delta": 0.125},
        output_names=("left", "right"),
    )
    path = root / "roe_face"
    path.mkdir(exist_ok=True)
    (path / "function.mir.json").write_text(compiled.math_ir.dumps())
    (path / "function.rir.yaml").write_text(dumps(compiled.program()))
    (path / "toffoli_u3_cz.originir").write_text(export_toffoli_u3_cz(compiled.program()).text)
    records.append(
        {
            "function": "roe_face",
            "math_modules": len(compiled.math_ir.functions),
            "rir_modules": len(compiled.program().modules),
            "correctness": "pending",
        }
    )
    inputs = roe_qfvm_inputs(fmt=FixedFormat(6, 2))
    be = roe_qfvm_block_encoding(inputs)
    path = root / "qfvm"
    path.mkdir(exist_ok=True)
    (path / "open.rir.yaml").write_text(dumps(be.operation.program()))
    closed = bind_qfvm(be.operation.program(), inputs)
    (path / "closed.rir.yaml").write_text(dumps(closed))
    (path / "toffoli_u3_cz.originir").write_text(export_toffoli_u3_cz(closed).text)
    records.append(
        {
            "function": "qfvm_with_compiled_roe",
            "rir_modules": len(closed.modules),
            "correctness": "pending",
        }
    )
    (root / "index.json").write_text(json.dumps(records, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
