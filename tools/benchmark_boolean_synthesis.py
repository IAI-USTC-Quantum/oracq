"""Compare exact SSA and affine Boolean circuits at unchanged numeric precision.

Run ``python tools/benchmark_boolean_synthesis.py --widths 3 6 8 --output out/boolean-synthesis.json``.
Reports real emitted-circuit counts and total lowered qubits. Generation/analysis
wall times are measurements, not quantum runtime estimates. No hardware backend,
network access, changed approximation tolerance, or idealized gate counts are used.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from oracq import estimate_resources
from oracq.algorithms.common.arithmetic import BooleanNetwork, FixedFormat, fixed_arithmetic


def measure(network: BooleanNetwork, *, optimize: bool) -> dict:
    """Compile and count the selected implementation of the same Boolean graph."""
    started = time.perf_counter()
    operation = network.operation(optimize=optimize)
    generated = time.perf_counter()
    estimate = estimate_resources(operation.program())
    counted = time.perf_counter()
    return {
        "qubits": estimate.qubits,
        "mcx_ancilla": estimate.mcx_ancilla,
        "lowered_qubits": estimate.qubits + estimate.mcx_ancilla,
        "toffoli": estimate.toffoli,
        "clifford": estimate.clifford,
        "gate_total": estimate.gate_total,
        "generate_seconds": generated - started,
        "estimate_seconds": counted - generated,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--widths", nargs="+", type=int, default=[3, 6, 8])
    parser.add_argument("--output", type=Path, default=Path("out/boolean-synthesis.json"))
    args = parser.parse_args()
    rows = []
    for width in args.widths:
        fmt = FixedFormat(width, max(0, width // 3))
        for kind in ("add", "mul", "div", "reciprocal", "sqrt"):
            arithmetic = fixed_arithmetic(kind, fmt)
            network = BooleanNetwork.from_payload(dict(arithmetic.module.attributes)["arithmetic_network"])
            before = measure(network, optimize=False)
            after = measure(network, optimize=True)
            rows.append({"kind": kind, "width": width, "fraction": fmt.fraction, "ssa": before, "affine": after})
            print(f"{kind} w={width}: qubits {before['lowered_qubits']} -> {after['lowered_qubits']}; "
                  f"Toffoli {before['toffoli']} -> {after['toffoli']}; Clifford {before['clifford']} -> {after['clifford']}", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"schema": "oracq.boolean-synthesis-benchmark.v1", "rows": rows}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
