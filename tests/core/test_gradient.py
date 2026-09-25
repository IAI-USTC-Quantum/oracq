"""Numerical witnesses for Jordan gradient estimation: exact linear readout, perturbation convergence, and binding consistency."""

import math
import unittest

from oracq import ValidationError, bind, simulate, unresolved
from oracq.algorithms.optimization.gradient import (
    abstract_phase_oracle,
    function_phase_oracle,
    gate_phase_oracle,
    gradient_estimation,
    gradient_from_readout,
)


def linear_oracle(coefficients, grid_bits, *, name=None):
    """Explicit phase-table oracle for the linear function f(x) = Σ_i c_i x_i (Jordan scaling)."""
    dimension = len(coefficients)
    grid_points = 1 << grid_bits
    angles = []
    for value in range(1 << (dimension * grid_bits)):
        phase = 0.0
        for i, coefficient in enumerate(coefficients):
            chunk = (value >> (i * grid_bits)) & (grid_points - 1)
            phase += coefficient * chunk / grid_points
        angles.append(2 * math.pi * grid_points * phase)
    return gate_phase_oracle(
        dimension * grid_bits, angles, phase_scale=grid_points, name=name
    )


def readout_distribution(state):
    """Readout distribution when target is the only root register."""
    return {key[0]: abs(amplitude) ** 2 for key, amplitude in state.amplitudes.items()}


def mode(distribution):
    return max(distribution, key=distribution.get)


class GradientTests(unittest.TestCase):
    def test_linear_function_exact_two_dimensions(self):
        # Gradient components take exactly representable grid values; the readout register lands deterministically on the encoded value.
        grid_bits = 4
        oracle = linear_oracle((3 / 16, -2 / 16), grid_bits)
        operation = gradient_estimation(oracle, dimension=2, grid_bits=grid_bits)
        attrs = dict(operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "jordan_gradient")
        self.assertEqual(attrs["oracle_queries"], 1)
        measured = readout_distribution(simulate(operation.program()))
        self.assertEqual(len(measured), 1)
        self.assertAlmostEqual(measured[mode(measured)], 1.0, places=12)
        self.assertEqual(
            gradient_from_readout(mode(measured), dimension=2, grid_bits=grid_bits),
            (3 / 16, -2 / 16),
        )

    def test_perturbed_linear_concentrates_with_grid_bits(self):
        # f(x) = a·x + x²/N² (N = 2**m): the phase perturbation is o(1/N) and the probability of a correct readout rises monotonically.
        a_numerator, a_denominator = 3, 8
        probabilities = []
        for grid_bits in (3, 4, 5):
            grid_points = 1 << grid_bits
            angles = [
                2
                * math.pi
                * grid_points
                * (
                    a_numerator / a_denominator * value / grid_points
                    + (value / grid_points) ** 2 / grid_points**2
                )
                for value in range(grid_points)
            ]
            oracle = gate_phase_oracle(grid_bits, angles, phase_scale=grid_points)
            operation = gradient_estimation(oracle, dimension=1, grid_bits=grid_bits)
            measured = readout_distribution(simulate(operation.program()))
            exact = a_numerator * grid_points // a_denominator
            self.assertEqual(mode(measured), exact)
            probabilities.append(measured[exact])
        self.assertLess(probabilities[0], probabilities[1])
        self.assertLess(probabilities[1], probabilities[2])
        # Failure decay rate: under an o(1/N) phase perturbation the leakage outside the peak
        # converges roughly quadratically as the grid refines (q decays by about 1/4).
        # Measured ratios q1/q0 ≈ 0.292 and q2/q1 ≈ 0.268; use the headroom upper bound 0.34.
        failures = [1 - p for p in probabilities]
        self.assertLess(failures[1], 0.34 * failures[0])
        self.assertLess(failures[2], 0.34 * failures[1])

    def test_abstract_oracle_binds_to_gate_implementation(self):
        grid_bits = 3
        oracle = abstract_phase_oracle("JordanPhase", 2 * grid_bits, phase_scale=1 << grid_bits)
        operation = gradient_estimation(oracle, dimension=2, grid_bits=grid_bits)
        program = operation.program()
        self.assertEqual([r.name for r in unresolved(program)], ["JordanPhase"])
        bound = bind(
            program,
            {"JordanPhase": linear_oracle((1 / 8, 1 / 4), grid_bits).operation},
        )
        self.assertFalse(unresolved(bound))
        measured = readout_distribution(simulate(bound))
        self.assertEqual(
            gradient_from_readout(mode(measured), dimension=2, grid_bits=grid_bits),
            (1 / 8, 1 / 4),
        )

    def test_function_phase_oracle_from_mathfunc(self):
        # mathfunc route: f(x) = 0.25·x is exact in the fixed-point format; the readout recovers 1/4.
        # Gradient components are only resolvable mod 1 (the phase e^{2πi·N·a·j} is identical for
        # a and a±1), so the witness takes components with |a| < 1/2.
        grid_bits = 2
        oracle = function_phase_oracle(
            "def f(x0):\n    return 0.25 * x0\n",
            dimension=1,
            grid_bits=grid_bits,
        )
        self.assertEqual(oracle.phase_scale, 1 << grid_bits)
        operation = gradient_estimation(oracle, dimension=1, grid_bits=grid_bits)
        measured = readout_distribution(simulate(operation.program()))
        self.assertEqual(
            gradient_from_readout(mode(measured), dimension=1, grid_bits=grid_bits),
            (0.25,),
        )

    def test_invalid_inputs_fail_at_generation(self):
        grid_bits = 3
        oracle = linear_oracle((1 / 8,), grid_bits)
        cases = [
            lambda: gradient_estimation(oracle, dimension=2, grid_bits=grid_bits),
            lambda: gradient_estimation(oracle, dimension=1, grid_bits=grid_bits + 1),
            lambda: gradient_estimation(oracle, dimension=0, grid_bits=grid_bits),
            lambda: gradient_estimation(oracle, dimension=1, grid_bits=0),
            lambda: gate_phase_oracle(2, [0.0, 0.0], phase_scale=4),
            lambda: gate_phase_oracle(2, [0.0] * 4, phase_scale=0),
            lambda: abstract_phase_oracle("Bad", 65, phase_scale=1),
            lambda: function_phase_oracle(
                "def f(x0):\n    return x0\n", dimension=1, grid_bits=0
            ),
            lambda: gradient_from_readout(-1, dimension=1, grid_bits=2),
        ]
        for case in cases:
            with self.subTest(case=case):
                with self.assertRaises(ValidationError):
                    case()


if __name__ == "__main__":
    unittest.main()
