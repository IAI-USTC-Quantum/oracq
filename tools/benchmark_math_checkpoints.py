"""Measure clean-region workspace and gate tradeoffs at fixed numeric precision.

Example: python tools/benchmark_math_checkpoints.py --roe --output out/math-checkpoints.json
The retained SSA circuit (size 0) is the reference. No precision, formula,
or resource counting rules change between rows.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from oracq import FixedFormat, estimate_resources
from oracq.applications.roe_formulas import frozen_roe_face
from oracq.infrastructure.mathfunc import Index, MathConfig, compile_function


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width", type=int, default=6)
    parser.add_argument("--fraction", type=int, default=2)
    parser.add_argument("--sizes", nargs="+", type=int, default=[0, 2, 4, 8, 16, 32])
    parser.add_argument("--roe", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("out/math-checkpoints.json"))
    args = parser.parse_args()
    source = "def f(x):\n a=x+1\n b=a*a\n c=b+1\n d=c*c\n e=d+1\n return e*e"
    options = {}
    if args.roe:
        source = frozen_roe_face
        options = {
            "inputs": {**{k: "real" for k in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r")},
                       "row": Index(2), "col": Index(2)},
            "constants": {"gamma": 1.4, "entropy_delta": 0.5},
            "output_names": ("left", "right"),
        }
    rows = []
    for size in args.sizes:
        started = time.perf_counter()
        result = compile_function(source, fmt=FixedFormat(args.width, args.fraction),
                                  config=MathConfig(checkpoint_size=size), **options)
        generated = time.perf_counter()
        estimate = estimate_resources(result.program())
        row = {"checkpoint_size": size, "qubits": estimate.qubits,
               "mcx_ancilla": estimate.mcx_ancilla,
               "lowered_qubits": estimate.qubits + estimate.mcx_ancilla,
               "toffoli": estimate.toffoli, "clifford": estimate.clifford,
               "modules": len(result.program().modules),
               "generate_seconds": generated - started,
               "estimate_seconds": time.perf_counter() - generated}
        rows.append(row)
        print(json.dumps(row), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"schema": "oracq.math-checkpoints.v1", "roe": args.roe,
        "width": args.width, "fraction": args.fraction, "rows": rows}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
