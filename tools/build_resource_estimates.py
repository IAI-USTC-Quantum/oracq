"""Toffoli+Clifford+T+QRAM resource scaling experiments for representative algorithms.

For each algorithm family, performs compositional resource estimation at
several sizes n (Repeat symbols are multiplied symbolically without
exponential expansion), fits the scaling exponent of each resource metric
log-log and compares it against the theoretical expectation; artifacts are
written to out/resource-estimates/. QRAM query counts are listed separately
as a first-class metric independent of gate-level cost (the gate-synthesized
QROM control group shows the two implementation costs of the same data
source).

Run: PYTHONPATH=src python tools/build_resource_estimates.py
"""

import json
import math
from pathlib import Path

from oracq import (
    QRAM,
    Bits,
    Builder,
    FixedFormat,
    UInt,
    estimate_resources,
)
from oracq.algorithms.common.arithmetic import fixed_arithmetic
from oracq.algorithms.common.estimation import phase_estimation
from oracq.algorithms.common.fourier import qft
from oracq.algorithms.common.search import grover
from oracq.algorithms.input_model.data_loading import qrom_lookup
from oracq.algorithms.input_model.oracles import (
    diagonal_block_encoding,
    gate_database,
    gate_state_prep,
    qram_database,
)
from oracq.applications.roe_formulas import frozen_roe_face
from oracq.infrastructure.mathfunc import Index, MathConfig, compile_function


def _prog(obj):
    "Operations use .program() directly; OracleView uses .operation.program()."
    return obj.program() if hasattr(obj, "program") else obj.operation.program()

OUT = Path("out/resource-estimates")


def _slope(xs, ys):
    "log2-log2 least-squares slope and coefficient of determination R²; non-positive samples are skipped."
    pairs = [(x, y) for x, y in zip(xs, ys, strict=True) if x > 0 and y > 0]
    if len(pairs) < 2:
        return None
    lx = [math.log2(x) for x, _ in pairs]
    ly = [math.log2(y) for _, y in pairs]
    mx, my = sum(lx) / len(lx), sum(ly) / len(ly)
    var = sum((x - mx) ** 2 for x in lx)
    if not var:
        return None
    slope = sum((x - mx) * (y - my) for x, y in zip(lx, ly, strict=True)) / var
    residual = sum((y - (my + slope * (x - mx))) ** 2 for x, y in zip(lx, ly, strict=True))
    total = sum((y - my) ** 2 for y in ly)
    r2 = 1 - residual / total if total else 1.0
    return {"slope": slope, "r2": r2}


def _run_group(name, sizes, build, expectation, note=""):
    records = []
    for n in sizes:
        program = build(n)
        estimate = estimate_resources(program)
        records.append({"n": n, **estimate.to_dict()})
        print(
            f"  {name} n={n}: qubits={estimate.qubits} toffoli={estimate.toffoli} "
            f"t_exact={estimate.t_exact} rot={len(estimate.rotations)} "
            f"qram={estimate.qram_total}",
            flush=True,
        )
    metrics = ("qubits", "toffoli", "clifford", "t_exact", "rotations", "t_total", "qram_total", "gate_total")
    slopes = {
        metric: _slope([r["n"] for r in records], [r[metric] for r in records])
        for metric in metrics
    }
    payload = {
        "group": name,
        "expectation": expectation,
        "note": note,
        "records": records,
        "slopes": slopes,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"[artifact] {OUT / f'{name}.json'} slopes={slopes}", flush=True)
    return payload


def row_add_const():
    def build(n):
        b = Builder(f"add_const_{n}", {"w": UInt(n)})
        b.add_const(b["w"], (1 << n) - 1)
        return b.finish().program()

    return _run_group(
        "add_const",
        [4, 8, 16, 32, 64],
        build,
        "ladder emission for constant addition: Toffoli ~ Σ(2d−3) = O(n^3) (current state of this compiler; literature optimum O(n), see gidney2018halving)",
    )


def row_fixed_mul():
    return _run_group(
        "fixed_mul",
        [4, 6, 8, 12, 16, 24, 32],
        lambda n: _prog(fixed_arithmetic("mul", FixedFormat(n, n // 2))),
        "AND network of the shift-and-add multiplier: Toffoli ≈ #AND = O(n^2)",
    )


def row_qft():
    return _run_group(
        "qft",
        [3, 4, 6, 8, 10, 12],
        lambda n: _prog(qft(n)),
        "n(n−1)/2 controlled rotations: rotation primitives O(n^2), T estimate = rotations × ceil(3 log2(1/eps))",
    )


def row_qpe_add_const():
    def build(p):
        unitary_builder = Builder("qpe_u", {"w": UInt(4)})
        unitary_builder.add_const(unitary_builder["w"], 1)
        return _prog(phase_estimation(unitary_builder.finish(), precision=p))

    return _run_group(
        "qpe_add_const",
        [3, 4, 6, 8, 10, 12],
        build,
        "controlled U^{2^k} via symbolic Repeat: total gate count ~ 2^p × cost(add_const)",
        note="shows the estimator counting without expanding the 2^k repetitions",
    )


def row_grover():
    def build(n):
        marked = (1 << (n - 1)) | 1
        oracle = Builder(f"mark_{n}", {"target": Bits(n)})
        with oracle.control(oracle["target"], marked):
            oracle.global_phase(math.pi)
        iterations = max(1, int(math.pi / 4 * math.sqrt(1 << n)))
        return _prog(grover(oracle.finish(), n, iterations=iterations))

    return _run_group(
        "grover",
        [4, 6, 8, 10, 12, 14],
        build,
        "⌊π/4·2^{n/2}⌋ iterations × per-iteration (marking + diffusion) cost: total gate count ~ 2^{n/2}·poly(n)",
    )


def row_state_prep():
    def build(n):
        state = 0x9E3779B97F4A7C15
        amplitudes = []
        for _ in range(1 << n):
            state = (state * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
            amplitudes.append(complex((state % 1000) / 1000, ((state >> 20) % 1000) / 1000))
        norm = math.sqrt(sum(abs(a) ** 2 for a in amplitudes))
        amplitudes = [a / norm for a in amplitudes]
        return _prog(gate_state_prep(amplitudes))

    return _run_group(
        "state_prep_dense",
        [2, 3, 4, 5, 6, 7, 8, 9, 10],
        build,
        "rotation cascade for dense state preparation: rotation primitives O(2^n)",
    )


def row_qram_vs_qrom():
    def qram_build(n):
        b = Builder(f"qram_{n}", {"a": UInt(n), "d": UInt(8)}, {"mem": QRAM(n, 8)})
        b.qram("mem", b["a"], b["d"])
        return b.finish().program()

    payload_qram = _run_group(
        "qram_lookup",
        [4, 8, 12, 16],
        qram_build,
        "QRAM as a first-class resource: each Load = 1 query, gate-level cost 0 (query count is size-independent)",
    )

    def qrom_build(n):
        table = [(i * 37 + 11) % 256 for i in range(1 << n)]
        return _prog(qrom_lookup(table))

    payload_qrom = _run_group(
        "qrom_lookup_gate",
        [4, 6, 8, 10, 12],
        qrom_build,
        "gate-level QROM synthesis of the same data source: Toffoli O(2^n), contrasted with the single query of qram_lookup",
    )
    return payload_qram, payload_qrom


def row_diagonal_be():
    def gate_build(n):
        table = [((i * 73 + 5) % 7) + 1 for i in range(1 << n)]
        return _prog(diagonal_block_encoding(gate_database(n, 3, table)))

    payload_gate = _run_group(
        "diagonal_be_gate",
        [2, 3, 4, 5, 6, 7, 8],
        gate_build,
        "diagonal block encoding (gate reuses rotations): rotation primitives O(2^n)",
    )

    def qram_build(n):
        return _prog(diagonal_block_encoding(qram_database(n, 3)))

    payload_qram = _run_group(
        "diagonal_be_qram",
        [2, 4, 6, 8, 10, 12],
        qram_build,
        "diagonal block encoding (QRAM angle database): QRAM queries O(1), gate-level cost weakly size-dependent",
    )
    return payload_gate, payload_qram


def row_roe_face():
    def build(w):
        compiled = compile_function(
            frozen_roe_face,
            fmt=FixedFormat(w, w // 2),
            inputs={
                **{key: "real" for key in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r")},
                "row": Index(2),
                "col": Index(2),
            },
            constants={"gamma": 1.4, "entropy_delta": 0.125},
            output_names=("left", "right"),
        )
        return compiled.program()

    return _run_group(
        "roe_face",
        [6, 8, 10, 12, 16],
        build,
        "fixed-point arithmetic network for the Roe flux face: gate count grows polynomially with word length w (application-level arithmetic scaling)",
    )


def row_math_polynomial():
    def build(w):
        compiled = compile_function(
            "def f(x):\n return x*x*x+0.5*x+0.25",
            fmt=FixedFormat(w, w // 2),
            config=MathConfig(degree=3),
        )
        return compiled.program()

    return _run_group(
        "math_polynomial",
        [4, 6, 8, 10, 12, 16],
        build,
        "cubic polynomial via compile_function: gate count grows polynomially with word length w",
    )


def main():
    print("== resource estimate scaling experiments ==", flush=True)
    row_add_const()
    row_fixed_mul()
    row_qft()
    row_qpe_add_const()
    row_grover()
    row_state_prep()
    row_qram_vs_qrom()
    row_diagonal_be()
    row_roe_face()
    row_math_polynomial()
    print("done", flush=True)


if __name__ == "__main__":
    main()
