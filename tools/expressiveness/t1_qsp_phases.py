"""oracq side of T1: QSP phase sequence for matrix inversion (κ=8, ε=1e-2).

The specification and pass thresholds are in ~/projects/oracq-dev/benchmarks/t1/SPEC.md.
Standalone run:

    cd ~/projects/qcfd-dev/oracq && \
    PYTHONPATH=src ~/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python \
    tools/expressiveness/t1_qsp_phases.py

The target polynomial follows the library documentation convention (constructed
internally by qsvt_matrix_inversion in algorithms/qsvt.py):
c·J_b(x) = c·(1−(1−x²)^b)/x, b = ceil(ln(1/ε)/−ln(1−1/κ²)), scaling c = 1/(3·sup J_b).
"""

import math
from importlib.metadata import version

import oracq
from oracq.algorithms.common.qsvt import qsp_phases, qsp_response
from oracq.infrastructure.ir import ValidationError

KAPPA, EPS, MAX_DEGREE, THRESHOLD = 8, 1e-2, 40, 1e-2


def inversion_poly(b):
    """Real coefficients of c·J_b in ascending powers (constant term first) and the scaling c."""
    f = [0.0] * (2 * b)
    for m in range(1, b + 1):
        f[2 * m - 1] = ((-1.0) ** (m + 1)) * math.comb(b, m)
    grid = [math.cos(math.pi * j / 2048) for j in range(2049)]
    sup = max(abs(sum(c * x**i for i, c in enumerate(f))) for x in grid)
    scale = 1.0 / (3.0 * sup)
    return [v * scale for v in f], scale


def synthesize(b):
    f, scale = inversion_poly(b)
    return qsp_phases(tuple(f), imag=(0.0, math.sqrt(1.0 - f[-1] ** 2))), f, scale


def main():
    print(f"oracq {version('oracq')} ({oracq.__file__})")
    b_spec = math.ceil(math.log(1.0 / EPS) / -math.log(1.0 - 1.0 / KAPPA**2))
    try:
        synthesize(b_spec)
        b = b_spec
        print(f"spec b={b_spec} (degree {2 * b_spec - 1}) synthesized")
    except ValidationError as err:
        print(f"spec b={b_spec} (degree {2 * b_spec - 1}) rejected: {err}")
        b = 0
        while True:
            try:
                synthesize(b + 1)
                b += 1
            except ValidationError:
                break
        print(f"max feasible b={b} (degree {2 * b - 1})")
    phases, f, scale = synthesize(b)
    grid = [math.cos(math.pi * j / 512) for j in range(513)]

    def target(x):
        return complex(sum(c * x**i for i, c in enumerate(f)),
                       math.sqrt(1.0 - f[-1] ** 2) * x)

    residual = max(abs(qsp_response(x, phases) - target(x)) for x in grid)
    approx = max(abs(target(x).real - scale / x) for x in grid if abs(x) >= 1.0 / KAPPA)
    approx_hi = max(abs(target(x).real - scale / x) for x in grid if abs(x) >= 0.5)
    print(f"degree={2 * b - 1} phases={len(phases)}")
    print(f"roundtrip_residual={residual:.3e} (threshold {THRESHOLD})")
    print(f"approx_err_vs_c_over_x={approx:.3e} on [1/{KAPPA},1]; {approx_hi:.3e} on [0.5,1]")
    print("status:", "ok" if residual <= THRESHOLD else "failed")


if __name__ == "__main__":
    main()
