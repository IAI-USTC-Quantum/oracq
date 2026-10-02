"""Numerical witnesses for the QSVT standard transform library: phase synthesis round trips, pinned conventions, inversion, filtering, simulation, and fixed-point search."""

import cmath
import math
import random
import unittest

from oracq import ValidationError, dumps, simulate
from oracq.algorithms.common.qsvt import (
    eigenstate_filter,
    fixed_point_search,
    fixed_point_search_phases,
    qsp_phases,
    qsp_response,
    qsvt_hamiltonian_simulation,
    qsvt_matrix_inversion,
)
from oracq.algorithms.common.transforms import qsvt_sequence
from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding


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


class PhaseSynthesisTests(unittest.TestCase):
    def test_roundtrip_chebyshev(self):
        for n in range(1, 7):
            phases = qsp_phases(chebyshev_t(n))
            self.assertEqual(len(phases), n + 1)
            for i in range(51):
                x = -0.99 + 1.98 * i / 50
                self.assertAlmostEqual(
                    qsp_response(x, phases), peval(chebyshev_t(n), x), delta=1e-8
                )

    def test_roundtrip_with_explicit_imaginary_part(self):
        # P(x) = (0.5 + i·√0.75)·x: |P|² = x², Q = 1, exactly realizable.
        phases = qsp_phases((0.0, 0.5), imag=(0.0, math.sqrt(0.75)))
        for i in range(21):
            x = -1.0 + 0.1 * i
            self.assertAlmostEqual(
                qsp_response(x, phases), complex(0.5, math.sqrt(0.75)) * x, delta=1e-9
            )

    def test_convention_matches_qsvt_sequence(self):
        # Under random phases, the zero-signal block of the qsvt_sequence circuit must equal qsp_response pointwise.
        rng = random.Random(7)
        phases = tuple(rng.uniform(-math.pi, math.pi) for _ in range(6))
        matrix = [[0.85, 0.0], [0.0, -0.55]]
        be = matrix_pauli_encoding(matrix)
        op = qsvt_sequence(be, phases)
        for col in range(2):
            x = matrix[col][col] / be.alpha
            state = simulate(op.program(), initial={"target": col})
            amp = state.amplitudes.get((col, 0), 0j)
            self.assertAlmostEqual(amp, qsp_response(x, phases), delta=1e-10)

    def test_negated_phases_conjugate_polynomial(self):
        # Property relied on by the real-part-extraction LCU: −Φ realizes P̄.
        rng = random.Random(11)
        phases = tuple(rng.uniform(-math.pi, math.pi) for _ in range(5))
        for i in range(21):
            x = -1.0 + 0.1 * i
            self.assertAlmostEqual(
                qsp_response(x, tuple(-p for p in phases)),
                qsp_response(x, phases).conjugate(),
                delta=1e-12,
            )

    def test_input_validation(self):
        with self.assertRaises(ValidationError):  # parity mismatch
            qsp_phases((0.0, 0.5, 0.5))
        with self.assertRaises(ValidationError):  # upper bound violation
            qsp_phases((0.0, 2.0))
        with self.assertRaises(ValidationError):  # endpoints not saturated and no imaginary part given
            qsp_phases((0.0, 0.5))
        with self.assertRaises(ValidationError):  # non-real coefficients
            qsp_phases((0.0, 1.0 + 1.0j))
        with self.assertRaises(ValidationError):  # imaginary-part parity mismatch
            qsp_phases((0.0, 0.5), imag=(0.5,))


class PhaseSynthesizerTests(unittest.TestCase):
    """Replaceable phase synthesizers: injection, the convention guard, and the degree-guard bypass."""

    @staticmethod
    def toy_chebyshev(coeffs, imag=None):
        # Reflection-convention phases realizing T_n exactly, mapped from the all-zero
        # sequence of the Wx (rotation) convention: reverse, subtract pi/4 from the two
        # end phases and pi/2 from the middle ones, then add n*pi/2 to the first-applied
        # phase so that p = e^{i n pi/2} i^{-n} T_n = T_n.
        n = len(tuple(coeffs)) - 1
        phases = [-math.pi / 2] * (n + 1)
        phases[0] = -math.pi / 4 + n * math.pi / 2
        phases[n] = -math.pi / 4
        return tuple(phases)

    @staticmethod
    def bundled_style(coeffs, imag=None):
        """A free-completion synthesizer rebuilt from the bundled route's own scaffolding."""
        if imag is not None:
            return qsp_phases(coeffs, imag)
        f = tuple(coeffs)
        v = 0.0
        for c in reversed(f):  # Horner at x = 1, matching the library's evaluation order
            v += c
        sat = math.sqrt(max(0.0, 1.0 - v * v))
        h = (0.0, sat) if (len(f) - 1) % 2 else (0.0, 0.0, sat)
        return qsp_phases(f, imag=h)

    def test_injected_synthesizer_realizes_chebyshev(self):
        for n in (3, 5):
            phases = qsp_phases(chebyshev_t(n), synthesizer=self.toy_chebyshev)
            self.assertEqual(len(phases), n + 1)
            for i in range(21):
                x = -0.99 + 1.98 * i / 20
                self.assertAlmostEqual(
                    qsp_response(x, phases).real, peval(chebyshev_t(n), x), delta=1e-12
                )

    def test_degree_guard_bypassed_with_synthesizer(self):
        # The degree limit guards the bundled numerics only: with a synthesizer the guard
        # does not fire and the call reaches the synthesizer (a well-conditioned target is
        # used so the coefficient-domain input checks still pass at degree 41).
        calls = []

        def recorder(coeffs, imag=None):
            calls.append(len(tuple(coeffs)) - 1)
            raise ValidationError("sentinel: synthesizer reached")

        with self.assertRaisesRegex(ValidationError, "sentinel"):
            qsp_phases((0.0,) * 41 + (0.5,), synthesizer=recorder)
        self.assertEqual(calls, [41])
        with self.assertRaisesRegex(ValidationError, "synthesis limit"):
            qsp_phases((0.0,) * 41 + (0.5,))  # the bundled route stays degree-guarded

    def test_roundtrip_guard_rejects_wrong_phases(self):
        with self.assertRaisesRegex(ValidationError, "round-trip"):
            qsp_phases(chebyshev_t(3), synthesizer=lambda coeffs, imag=None: (0.1,) * 4)

    def test_synthesizer_output_shape_checked(self):
        with self.assertRaisesRegex(ValidationError, "expected 4"):
            qsp_phases(chebyshev_t(3), synthesizer=lambda coeffs, imag=None: (0.0, 0.0))
        with self.assertRaisesRegex(ValidationError, "non-real or non-finite"):
            qsp_phases(
                chebyshev_t(3),
                synthesizer=lambda coeffs, imag=None: (0.0, float("nan"), 0.0, 0.0),
            )

    def test_pinned_completion_mode_is_exact(self):
        # With imag given, the synthesizer must honor p = f + i·h exactly, not just Re p = f.
        target, completion = (0.0, 0.5), (0.0, math.sqrt(0.75))

        def shifted(coeffs, imag=None):
            phases = qsp_phases(coeffs, imag)
            return (phases[0] + 0.3,) + phases[1:]

        with self.assertRaisesRegex(ValidationError, "round-trip"):
            qsp_phases(target, imag=completion, synthesizer=shifted)
        phases = qsp_phases(
            target, imag=completion, synthesizer=lambda c, imag=None: qsp_phases(c, imag)
        )
        x = 0.37
        self.assertAlmostEqual(
            qsp_response(x, phases), complex(0.5, math.sqrt(0.75)) * x, delta=1e-9
        )

    def test_upstream_threading_matches_default_bitwise(self):
        # A synthesizer reproducing the bundled scaffolding yields byte-identical programs.
        be = matrix_pauli_encoding([[0.6, -0.2], [-0.2, 0.6]])
        default = qsvt_matrix_inversion(be, 2.0, error=0.15)
        custom = qsvt_matrix_inversion(be, 2.0, error=0.15, synthesizer=self.bundled_style)
        self.assertEqual(
            dumps(default.operation.program()), dumps(custom.operation.program())
        )
        be2 = matrix_pauli_encoding([[0.05, 0.0], [0.0, 0.5]])
        default_f = eigenstate_filter(be2, 0.2, 8, center=0.1)
        custom_f = eigenstate_filter(be2, 0.2, 8, center=0.1, synthesizer=self.bundled_style)
        self.assertEqual(
            dumps(default_f.operation.program()), dumps(custom_f.operation.program())
        )

    def test_fixed_point_search_with_synthesizer(self):
        default = fixed_point_search_phases(0.4, 5)
        custom = fixed_point_search_phases(
            0.4, 5, synthesizer=lambda coeffs, imag=None: qsp_phases(coeffs, imag)
        )
        self.assertEqual(len(custom), 6)
        for i in range(1, 41):
            x = i / 41
            self.assertAlmostEqual(
                abs(qsp_response(x, custom)) ** 2,
                abs(qsp_response(x, default)) ** 2,
                delta=1e-6,
            )

    def test_kappa8_inversion_boundary_is_explicit(self):
        # The retained boundary: kappa=8 at error=1e-2 requires the degree-585 polynomial,
        # which the bundled route rejects by its degree guard (and the monomial-coefficient
        # pipeline itself loses reliability around degree 40, independently of synthesis).
        be = matrix_pauli_encoding([[0.6, -0.2], [-0.2, 0.6]])
        with self.assertRaisesRegex(ValidationError, "synthesis limit"):
            qsvt_matrix_inversion(be, 8.0, error=1e-2)


class TransformWitnessTests(unittest.TestCase):
    def test_matrix_inversion_block(self):
        matrix = [[0.6, -0.2], [-0.2, 0.6]]  # eigenvalues 0.4, 0.8: κ = 2
        be = matrix_pauli_encoding(matrix)
        inv = qsvt_matrix_inversion(be, 2.0, error=0.15)
        attrs = dict(inv.operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "qsvt_matrix_inversion")
        self.assertEqual(attrs["be_alpha"], 1.0)
        scale = attrs["inverse_scale"]
        inverse = [[1.875, 0.625], [0.625, 1.875]]
        for col in range(2):
            state = simulate(inv.operation.program(), initial={"target": col})
            for row in range(2):
                amp = state.amplitudes.get((row, 0), 0j)
                expected = scale * inverse[row][col]
                self.assertAlmostEqual(amp, expected, delta=0.35 * abs(expected))

    def test_eigenstate_filter_isolates_eigenvalue(self):
        # Eigenvalues 0.05 and 0.5 (α = 0.5); the filter center shifts to 0.1 (spectral variable).
        be = matrix_pauli_encoding([[0.05, 0.0], [0.0, 0.5]])
        flt = eigenstate_filter(be, 0.2, 8, center=0.1)
        attrs = dict(flt.operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "eigenstate_filter")
        self.assertLess(attrs["suppression"], 0.08)
        # Shifted eigenvalues: (0.05−0.1)/0.6 ≈ −0.083 (passband) and (0.5−0.1)/0.6 ≈ 0.667 (stopband).
        state = simulate(flt.operation.program(), initial={"target": 0})
        self.assertGreater(abs(state.amplitudes.get((0, 0), 0j)), 0.7)
        state = simulate(flt.operation.program(), initial={"target": 1})
        self.assertLess(abs(state.amplitudes.get((1, 0), 0j)), 0.1)

    def test_hamiltonian_simulation_block(self):
        matrix = [[0.5, 0.0], [0.0, -0.25]]
        be = matrix_pauli_encoding(matrix)
        t = 0.7
        hs = qsvt_hamiltonian_simulation(be, t, error=0.01)
        attrs = dict(hs.operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "qsvt_hamiltonian_simulation")
        scale = attrs["sim_scale"]
        for col in range(2):
            x = matrix[col][col] / be.alpha
            state = simulate(hs.operation.program(), initial={"target": col})
            amp = state.amplitudes.get((col, 0), 0j)
            self.assertAlmostEqual(amp, cmath.exp(1j * t * x) / scale, delta=5e-4)

    def test_fixed_point_search_phase_properties(self):
        delta = 0.4
        thresholds = []
        for degree in (3, 5, 7):
            phases = fixed_point_search_phases(delta, degree)
            self.assertEqual(len(phases), degree + 1)
            c = math.cosh(math.acosh(1.0 / delta) / degree)
            theta = math.sqrt(1.0 - 1.0 / c**2)
            thresholds.append(theta)
            for i in range(1, 401):
                x = i / 400
                p2 = abs(qsp_response(x, phases)) ** 2
                self.assertLessEqual(p2, 1.0 + 1e-9)
                # matches the YLC closed-form success probability
                v = c * math.sqrt(1.0 - x * x)
                t_l = math.cos(degree * math.acos(v)) if v <= 1 else math.cosh(
                    degree * math.acosh(v)
                )
                self.assertAlmostEqual(p2, 1.0 - delta**2 * t_l**2, delta=1e-5)
                if x >= theta:
                    self.assertGreaterEqual(p2, 1.0 - delta**2 - 1e-9)
        # Fixed-point property: the threshold decreases monotonically with the query count and tends to 0
        self.assertGreater(thresholds[0], thresholds[1])
        self.assertGreater(thresholds[1], thresholds[2])

    def test_fixed_point_search_circuit_amplifies(self):
        be = matrix_pauli_encoding([[0.6, 0.0], [0.0, 0.1]])  # α = 0.6: x = 1.0 and 1/6
        delta, degree = 0.4, 5
        fp = fixed_point_search(be, delta, degree)
        attrs = dict(fp.operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "fixed_point_search")
        self.assertAlmostEqual(attrs["threshold"], math.sqrt(1.0 - 1.0 / 1.04955**2), places=3)
        state = simulate(fp.operation.program(), initial={"target": 0})
        self.assertGreaterEqual(abs(state.amplitudes.get((0, 0), 0j)) ** 2, 1.0 - delta**2 - 1e-9)

    def test_transform_input_validation(self):
        be = matrix_pauli_encoding([[0.5, 0.0], [0.0, 0.25]])
        with self.assertRaises(ValidationError):
            qsvt_matrix_inversion(be, 0.5)
        with self.assertRaises(ValidationError):
            qsvt_matrix_inversion(be, 2.0, error=1.0)
        with self.assertRaises(ValidationError):  # required degree exceeds the synthesis cap
            qsvt_matrix_inversion(be, 10.0, error=1e-9)
        with self.assertRaises(ValidationError):
            eigenstate_filter(be, 1.5, 4)
        with self.assertRaises(ValidationError):
            eigenstate_filter(be, 0.2, 0)
        with self.assertRaises(ValidationError):
            qsvt_hamiltonian_simulation(be, 0.0)
        with self.assertRaises(ValidationError):
            fixed_point_search_phases(1.5, 5)
        with self.assertRaises(ValidationError):  # degree must be odd (even given)
            fixed_point_search_phases(0.3, 4)


if __name__ == "__main__":
    unittest.main()
