"""oracq side of T1: QSP phase sequence for matrix inversion (κ=8, ε=1e-2).

Pass thresholds: the round-trip residual against the target polynomial and the
approximation error against c/x on [1/κ, 1], both at most 1e-2 (the fixture
recorded in the paper's discussion section).

Standalone run (the pyqsp route needs the optional ``pyqsp`` extra):

    cd oracq && uv run --extra pyqsp python tools/expressiveness/t1_qsp_phases.py

The target polynomial follows the library documentation convention (constructed
internally by qsvt_matrix_inversion in algorithms/common/qsvt.py):
c·J_b(x) = c·(1−(1−x²)^b)/x, b = ceil(ln(1/ε)/−ln(1−1/κ²)), scaling c = 1/(3·sup J_b).

Two routes are measured:

- bundled: qsp_phases' root-finding + layer-stripping route, degree-guarded at 40;
  the spec degree 585 is rejected and the tool reports the maximum feasible b.
- pyqsp: the replaceable-synthesizer route through
  oracq.algorithms.common.qsp_pyqsp, with the target constructed and validated in
  the Chebyshev basis. The library's monomial-coefficient pipeline loses
  reliability around degree 40 independently of the synthesis route, so the
  high-degree target is never round-tripped through monomial coefficients.
"""

import importlib.util
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


def synthesize_pyqsp(b):
    """Chebyshev-native route: c·J_b built stably in the Chebyshev basis, phases from the
    pyqsp adapter; returns (phases, scale), with the scale computed from the closed form."""
    import numpy as np
    import numpy.polynomial.chebyshev as cheb

    from oracq.algorithms.common.qsp_pyqsp import chebyshev_phases_via_pyqsp

    d = 2 * b - 1

    def target(x):
        x = np.asarray(x, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(np.abs(x) < 1e-8, b * x, (1.0 - (1.0 - x * x) ** b) / x)

    grid = np.cos(np.pi * np.arange(4097) / 4096)
    scale = 1.0 / (3.0 * float(np.max(np.abs(target(grid)))))
    interp = cheb.Chebyshev.interpolate(lambda xs: scale * target(xs), 4 * d)
    coef = np.array(interp.coef[: d + 1])
    coef[0::2] = 0.0  # odd parity exactly
    return chebyshev_phases_via_pyqsp(tuple(coef)), scale


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
    print(f"[bundled] degree={2 * b - 1} phases={len(phases)}")
    print(f"[bundled] roundtrip_residual={residual:.3e} (threshold {THRESHOLD})")
    print(f"[bundled] approx_err_vs_c_over_x={approx:.3e} on [1/{KAPPA},1]; {approx_hi:.3e} on [0.5,1]")
    print("[bundled] status:", "ok" if residual <= THRESHOLD else "failed")

    if importlib.util.find_spec("pyqsp") is None:
        print("[pyqsp] skipped: the optional pyqsp extra is not installed")
        return
    phases_p, scale_p = synthesize_pyqsp(b_spec)

    def target_re(x):
        return scale_p * (b_spec * x if abs(x) < 1e-8 else (1.0 - (1.0 - x * x) ** b_spec) / x)

    # Free completion: the adapter pins the real part only, so the residual is on Re p.
    residual_p = max(abs(qsp_response(x, phases_p).real - target_re(x)) for x in grid)
    approx_p = max(abs(qsp_response(x, phases_p).real - scale_p / x) for x in grid
                   if abs(x) >= 1.0 / KAPPA)
    approx_p_hi = max(abs(qsp_response(x, phases_p).real - scale_p / x) for x in grid
                      if abs(x) >= 0.5)
    print(f"[pyqsp] degree={2 * b_spec - 1} phases={len(phases_p)}")
    print(f"[pyqsp] roundtrip_residual_re={residual_p:.3e} (threshold {THRESHOLD})")
    print(f"[pyqsp] approx_err_vs_c_over_x={approx_p:.3e} on [1/{KAPPA},1];"
          f" {approx_p_hi:.3e} on [0.5,1]")
    print("[pyqsp] status:", "ok" if residual_p <= THRESHOLD and approx_p <= THRESHOLD else "failed")


if __name__ == "__main__":
    main()
