"""Compare three oracle implementations from one open program; --native enables cross-checking against the real backend."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path

from oracq import bind_with_report, dump_qram_yaml, dumps, estimate_resources, loads, simulate
from oracq.applications.oracle_study import oracle_study, oracle_study_reference


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", action="store_true")
    parser.add_argument("--widths", type=int, nargs="+", default=[2, 3, 4])
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("-o", "--output", type=Path, default=Path("out/research-workflow"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    for width in args.widths:
        opened, variants = oracle_study(width, args.repetitions)
        text = dumps(opened)
        (args.output / f"width-{width}.open.rir.yaml").write_text(text)
        opened = loads(text)
        expected = oracle_study_reference(width, args.repetitions)
        for name, (binding, memory) in variants.items():
            linked = bind_with_report(opened, {"AngleWord": binding})
            program = linked.require()
            actual = dict(simulate(program, memory).amplitudes)
            errors = {"reference": max(abs(actual.get(k, 0) - expected.get(k, 0)) for k in actual.keys() | expected.keys())}
            if args.native:
                from oracq import run_originir, run_pysparq, run_pysparq_rir

                for label, runner in (("adapter_pysparq", run_pysparq), ("rir_pysparq", run_pysparq_rir)):
                    state = dict(runner(program, memory).amplitudes)
                    errors[label] = max(abs(state.get(k, 0) - expected.get(k, 0)) for k in state.keys() | expected.keys())
                vector = run_originir(program, memory)
                expected_vector = {target + (signal << width): a for (target, signal), a in expected.items()}
                errors["originir"] = max(abs(complex(a) - expected_vector.get(i, 0)) for i, a in enumerate(vector))
            if max(errors.values()) > 1e-9:
                raise AssertionError((width, name, errors))
            prefix = args.output / f"width-{width}.{name}"
            prefix.with_suffix(f".{name}.rir.yaml").write_text(dumps(program))
            memory_text = dump_qram_yaml(program, memory)
            prefix.with_suffix(f".{name}.memory.qram.yaml").write_text(memory_text, encoding="utf-8")
            record = {
                "width": width, "repetitions": args.repetitions, "implementation": name,
                "binding": linked.report.to_dict(),
                "memory_sha256": hashlib.sha256(memory_text.encode()).hexdigest(),
                "open_cost": estimate_resources(opened, require_closed=False).to_dict(),
                "closed_cost": estimate_resources(program).to_dict(),
                "amplitude_errors": errors, "tolerance": 1e-9,
                "success_probability": sum(abs(a) ** 2 for (target, signal), a in actual.items() if signal == 0),
            }
            records.append(record)
            cost = estimate_resources(program)
            print(width, name, "Toffoli", cost.toffoli, "QRAM", cost.qram_total, errors, flush=True)
    report = {
        "python": platform.python_version(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "backend_versions": {
            name: importlib.metadata.version(name)
            for name in ("pysparq", "unified-quantum", "uniqc-cppsimulator")
        } if args.native else {},
        "reference": "cos/sin of repetitions*pi*((address+1) mod 2**width)/2**width",
        "scope": "all target basis branches from a uniform input; public word and private work return to zero",
        "native_requested": args.native, "cases": records,
    }
    (args.output / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
