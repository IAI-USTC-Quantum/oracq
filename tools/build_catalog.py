"""生成设计案例面板与开放/绑定/后端产物。"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from pyqecclang import bind, dumps, unresolved
from pyqecclang.linking import Binding
from pyqecclang.workloads import CASES, build_case


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("out/catalog"))
    parser.add_argument("--case", action="append", choices=CASES)
    parser.add_argument("--native-parse", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.native_parse:
        from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser
    rows = []
    for name in args.case or CASES:
        case = build_case(name)
        directory = args.out / name
        directory.mkdir(parents=True, exist_ok=True)
        closed = case.closed()
        remaining = unresolved(closed)
        if remaining:
            raise RuntimeError(f"{name} 的具体实例仍有未绑定槽：{remaining}")
        (directory / "open.rir.json").write_text(dumps(case.program))
        (directory / "closed.rir.json").write_text(dumps(closed))
        partial_keys = sorted(case.bindings)[: len(case.bindings) // 2]
        partial = bind(case.program, {key: case.bindings[key] for key in partial_keys})
        (directory / "partial.rir.json").write_text(dumps(partial))
        artifact = case.artifact()
        (directory / "program.originir").write_text(artifact.text)
        (directory / "memory.json").write_text(json.dumps(case.memory, indent=2) + "\n")
        manifest = {}
        for slot, value in sorted(case.bindings.items()):
            item = value if isinstance(value, Binding) else Binding(value)
            path = "bindings/" + slot + ".rir.json"
            (directory / "bindings").mkdir(exist_ok=True)
            (directory / path).write_text(dumps(item.operation.program()))
            manifest[slot] = {"program": path, "resources": item.resources or {}}
        (directory / "readout.json").write_text(
            json.dumps([asdict(a) for a in case.readout], indent=2) + "\n"
        )
        (directory / "bindings.json").write_text(json.dumps(manifest, indent=2) + "\n")
        flat_operations = None
        if args.native_parse:
            native = OriginIR_BaseParser()
            native.parse(artifact.text)
            assert native.n_qubit == sum(r.type.width for r in closed.main.registers)
            flat_operations = len(native.program_body)
            if case.readout:
                from uniqc.circuit_builder.classical_program import parse_originir_ext_dynamic

                from pyqecclang.readout import export_with_readout

                execution = export_with_readout(closed, case.readout)
                parse_originir_ext_dynamic(execution.text)
                (directory / "execution-with-readout.originir").write_text(execution.text)
        row = {
            "case": name,
            "status": "paradigm_complete",
            "correctness": "not_assessed",
            "unresolved": [asdict(item) for item in unresolved(case.program)],
            "partial_unresolved": [r.name for r in unresolved(partial)],
            "closed_modules": len(closed.modules),
            "qubits": sum(r.type.width for r in closed.main.registers),
            "qram_resources": [r.name for r in closed.main.resources],
            "originir_bytes": len(artifact.text.encode()),
            "originir_definitions": artifact.text.count("DEF "),
            "native_parsed": args.native_parse,
            "readout_bridge": "explicit_flat_dynamic" if case.readout else None,
            "native_operations": flat_operations,
            "source": case.source,
            "notes": case.notes,
            "readout": [asdict(action) for action in case.readout],
        }
        (directory / "report.json").write_text(json.dumps(row, ensure_ascii=False, indent=2) + "\n")
        rows.append(row)
        print(
            f"{name}: slots={len(row['unresolved'])}, modules={row['closed_modules']}, "
            f"OriginIR={row['originir_bytes']} bytes, native={args.native_parse}",
            flush=True,
        )
    (args.out / "index.json").write_text(
        json.dumps({"phase": "paradigm", "cases": rows}, ensure_ascii=False, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
