"""Numerical witnesses for QPCA and density matrix exponentiation (the LMR protocol)."""

import math
import unittest

from oracq import ValidationError, simulate
from oracq.algorithms.input_model.density import (
    gate_purification,
    partial_trace,
    trace_distance,
)
from oracq.algorithms.input_model.oracles import gate_state_prep
from oracq.algorithms.qml.qpca import (
    density_matrix_exponentiation,
    eigenvalue_from_phase,
    qpca,
)

SQRT2 = math.sqrt(2)
PLUS = [1 / SQRT2, 1 / SQRT2]
MINUS = [1 / SQRT2, -1 / SQRT2]


def dense_system(state, copies_width, copies):
    vector = [0j] * (1 << (1 + copies_width * copies))
    for key, amplitude in state.amplitudes.items():
        vector[key[0] + (key[1] << 1)] += amplitude
    return partial_trace(vector, 1, copies_width * copies)


def exact_evolved_pure(time):
    """Exact result of e^{-it·|+⟩⟨+|} acting on ``|0⟩⟨0|``."""
    import cmath

    a0 = (1 + cmath.exp(-1j * time)) / 2
    a1 = -(1 - cmath.exp(-1j * time)) / 2
    return (
        (a0 * a0.conjugate(), a0 * a1.conjugate()),
        (a1 * a0.conjugate(), a1 * a1.conjugate()),
    )


def phase_mode(state):
    distribution = {}
    for key, amplitude in state.amplitudes.items():
        distribution[key[-1]] = distribution.get(key[-1], 0.0) + abs(amplitude) ** 2
    return max(distribution, key=distribution.get), distribution


def circular_eigenvalue(distribution, precision, step):
    """Circular-mean decoding of the phase distribution: a λ estimate robust to a broadened peak."""
    import cmath

    z = sum(
        p * cmath.exp(2j * math.pi * v / (1 << precision))
        for v, p in distribution.items()
    )
    return -2 * math.pi * (cmath.phase(z) / (2 * math.pi)) / step


class DensityMatrixExponentiationTests(unittest.TestCase):
    def test_small_step_first_order_accurate(self):
        # Single-step error O(Δt²): at Δt = 0.05 the trace distance stays within the Δt² scale (coefficient about 0.53).
        step = 0.05
        prep = gate_state_prep(PLUS)
        operation = density_matrix_exponentiation(prep, time=step, copies=1)
        reduced = dense_system(simulate(operation.program()), 1, 1)
        self.assertLess(trace_distance(reduced, exact_evolved_pure(step)), 0.6 * step * step)

    def test_error_halves_with_copies(self):
        # LMR first-order scaling: with the total time t fixed, doubling the copies roughly halves the error.
        time = 0.4
        prep = gate_state_prep(PLUS)
        distances = []
        for copies in (1, 2, 4):
            operation = density_matrix_exponentiation(prep, time=time, copies=copies)
            reduced = dense_system(simulate(operation.program()), 1, copies)
            distances.append(trace_distance(reduced, exact_evolved_pure(time)))
        self.assertLess(distances[1], distances[0] * 0.6)
        self.assertLess(distances[2], distances[1] * 0.6)


class QpcaTests(unittest.TestCase):
    def test_pure_state_eigenvalues(self):
        # ρ = |+⟩⟨+|: eigenvalues 1 and 0 are read out with the system input in |+⟩ and |−⟩, respectively.
        step = math.pi / 4  # λ=1 gives φ = 7/8; deterministic readout at precision=3.
        prep = gate_state_prep(PLUS)
        operation = qpca(prep, precision=3, step_time=step, system=gate_state_prep(PLUS))
        attrs = dict(operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "qpca")
        self.assertEqual(attrs["copies"], 7)
        # |+⟩ is the copy state itself: each partial swap acts on a SWAP-symmetric eigenstate, exact for any Δt.
        mode, distribution = phase_mode(simulate(operation.program()))
        self.assertAlmostEqual(distribution[mode], 1.0, places=9)
        self.assertAlmostEqual(eigenvalue_from_phase(mode, 3, step), 1.0, places=12)

        # |−⟩ is a generic eigenstate: at large Δt the LMR error appears as peak broadening, but the
        # circular mean of the distribution is still centered at φ = 0 (λ = 0).
        operation = qpca(prep, precision=3, step_time=step, system=gate_state_prep(MINUS))
        _, distribution = phase_mode(simulate(operation.program()))
        self.assertAlmostEqual(circular_eigenvalue(distribution, 3, step), 0.0, delta=0.15)

    def test_mixed_state_eigenvalue_via_purification(self):
        # ρ = diag(0.75, 0.25) adapted via purification; the system input |0⟩ reads out λ = 0.75.
        step = 2 * math.pi / 6  # λ=0.75 gives φ = 7/8.
        purification = gate_purification(((0.75 + 0j, 0j), (0j, 0.25 + 0j)))
        prep = purification.as_state_preparation()
        operation = qpca(prep, precision=3, step_time=step, swap_width=1)
        mode, distribution = phase_mode(simulate(operation.program()))
        # The large-Δt LMR error only broadens the peak without shifting it: the mode decodes to exactly 0.75.
        self.assertGreater(distribution[mode], 0.4)
        self.assertAlmostEqual(eigenvalue_from_phase(mode, 3, step), 0.75, places=12)
        self.assertAlmostEqual(circular_eigenvalue(distribution, 3, step), 0.75, delta=0.15)

    def test_invalid_inputs_fail_at_generation(self):
        prep = gate_state_prep(PLUS)
        cases = [
            lambda: density_matrix_exponentiation(prep, time=0.0, copies=1),
            lambda: density_matrix_exponentiation(prep, time=1.0, copies=0),
            lambda: density_matrix_exponentiation(object(), time=1.0, copies=1),
            lambda: qpca(prep, precision=0, step_time=1.0),
            lambda: qpca(prep, precision=7, step_time=1.0),
            lambda: qpca(prep, precision=3, step_time=0.0),
            lambda: qpca(prep, precision=3, step_time=1.0, swap_width=2),
            lambda: qpca(
                prep, precision=3, step_time=1.0, system=gate_state_prep([1, 0, 0, 0])
            ),
            lambda: eigenvalue_from_phase(0, 3, 0.0),
            lambda: eigenvalue_from_phase(8, 3, 1.0),
        ]
        for case in cases:
            with self.subTest(case=case):
                with self.assertRaises(ValidationError):
                    case()


if __name__ == "__main__":
    unittest.main()
