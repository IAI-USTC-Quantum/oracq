"""pyqsp-backed phase synthesizer: convention mapping, contract conformance, and the kappa=8 case.

The suite is skipped unless the optional ``pyqsp`` extra is installed; nothing in the
core depends on it.
"""

import importlib.util
import math
import unittest

from oracq import ValidationError, simulate
from oracq.algorithms.common.qsvt import qsp_phases, qsp_response, qsvt_matrix_inversion
from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding

pyqsp_available = importlib.util.find_spec("pyqsp") is not None


def chebyshev_t(n):
    """Monomial coefficients in ascending powers of T_n."""
    if n == 0:
        return (1.0,)
    if n == 1:
        return (0.0, 1.0)
    a, b = (1.0,), (0.0, 1.0)
    for _ in range(2, n + 1):
        nb = tuple(2 * v for v in (0.0,) + b)
        nb = tuple(nb[i] - (a[i] if i < len(a) else 0.0) for i in range(len(nb)))
        a, b = b, nb
    return b


def peval(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))


@unittest.skipUnless(pyqsp_available, "the optional pyqsp extra is not installed")
class PyqspSynthesizerTests(unittest.TestCase):
    def setUp(self):
        from oracq.algorithms.common.qsp_pyqsp import PyqspSynthesizer

        self.synth = PyqspSynthesizer()

    def test_convention_mapping_matches_target_several_parities(self):
        # The Wx-to-reflection mapping is checked by qsp_phases itself on a 512-point
        # grid; here the response is compared against an independent evaluation of the target.
        for n in (2, 3, 5, 8, 11, 16):
            f = tuple(0.5 * v for v in chebyshev_t(n))
            phases = qsp_phases(f, synthesizer=self.synth)
            self.assertEqual(len(phases), n + 1)
            for i in range(41):
                x = -0.99 + 1.98 * i / 40
                self.assertAlmostEqual(
                    qsp_response(x, phases).real, peval(f, x), delta=1e-8
                )

    def test_pinned_completion_mode_rejected(self):
        with self.assertRaisesRegex(ValidationError, "free-completion"):
            qsp_phases((0.0, 0.5), imag=(0.0, math.sqrt(0.75)), synthesizer=self.synth)

    def test_library_consumer_beyond_bundled_guard(self):
        # kappa=3 at error=0.05 requires degree 51: rejected by the bundled route's degree
        # guard, synthesized by the adapter, and the block encoding verifies against the
        # true inverse at the requested precision.
        be = matrix_pauli_encoding([[0.4, -0.2], [-0.2, 0.4]])  # eigenvalues 0.2, 0.6
        with self.assertRaisesRegex(ValidationError, "synthesis limit"):
            qsvt_matrix_inversion(be, 3.0, error=0.05)
        inv = qsvt_matrix_inversion(be, 3.0, error=0.05, synthesizer=self.synth)
        attrs = dict(inv.operation.module.attributes)
        self.assertEqual(attrs["qsp_degree"], 51)
        scale = attrs["inverse_scale"]
        inverse = [[10.0 / 3.0, 5.0 / 3.0], [5.0 / 3.0, 10.0 / 3.0]]
        for col in range(2):
            state = simulate(inv.operation.program(), initial={"target": col})
            for row in range(2):
                amp = state.amplitudes.get((row, 0), 0j)
                self.assertAlmostEqual(amp, scale * inverse[row][col], delta=0.02)

    def test_kappa8_coefficient_pipeline_boundary_is_explicit(self):
        # kappa=8 at error=1e-2 needs the degree-585 polynomial: the adapter lifts the
        # synthesis-degree guard, but the monomial-coefficient pipeline of the library-level
        # entry loses reliability around degree 40, so the failure is caught explicitly by
        # the round-trip guard rather than producing a silently wrong scaling.
        be = matrix_pauli_encoding([[0.6, -0.2], [-0.2, 0.6]])
        with self.assertRaises(ValidationError):
            qsvt_matrix_inversion(be, 8.0, error=1e-2, synthesizer=self.synth)

    def test_kappa8_phase_level_target_reached(self):
        # The retained kappa=8, eps=1e-2 inversion target, reached through the adapter with
        # the target constructed and validated in the Chebyshev basis: degree 585, round-trip
        # residual ~1e-13, and approximation error < 1e-2 against c/x on |x| >= 1/8.
        import numpy as np
        import numpy.polynomial.chebyshev as cheb

        from oracq.algorithms.common.qsp_pyqsp import chebyshev_phases_via_pyqsp

        kappa, eps = 8, 1e-2
        b = max(1, math.ceil(math.log(1.0 / eps) / -math.log(1.0 - 1.0 / kappa**2)))
        d = 2 * b - 1
        self.assertEqual(d, 585)

        def target_closed_form(x):
            return b * x if abs(x) < 1e-8 else (1.0 - (1.0 - x * x) ** b) / x

        grid = np.cos(np.pi * np.arange(4097) / 4096)
        sup = float(max(abs(target_closed_form(float(x))) for x in grid))
        c_scale = 1.0 / (3.0 * sup)

        with np.errstate(divide="ignore", invalid="ignore"):
            interp = cheb.Chebyshev.interpolate(
                lambda xs: c_scale * np.where(
                    np.abs(xs) < 1e-8, b * xs, (1.0 - (1.0 - xs**2) ** b) / xs
                ),
                4 * d,
            )
        coef = np.array(interp.coef[: d + 1])
        coef[0::2] = 0.0  # odd parity exactly
        phases = chebyshev_phases_via_pyqsp(tuple(coef))
        self.assertEqual(len(phases), d + 1)

        check = np.cos(np.pi * np.arange(2049) / 2048)
        residual = max(
            abs(qsp_response(float(x), phases).real - c_scale * target_closed_form(float(x)))
            for x in check
        )
        self.assertLessEqual(residual, eps)
        band = [float(x) for x in check if abs(x) >= 1.0 / kappa]
        approx = max(
            abs(qsp_response(x, phases).real - c_scale / x) for x in band
        )
        self.assertLessEqual(approx, eps)


if __name__ == "__main__":
    unittest.main()
