"""首批算法的独立数学见证，覆盖相位、位序、概率和经典后处理。"""

import cmath
import importlib
import math
import unittest

from oracq import Bits, Builder, ValidationError, simulate
from oracq.algorithms.basics.number_theory import (
    factors_from_phase,
    modular_multiply,
    order_finding,
)
from oracq.algorithms.basics.oracle_algorithms import (
    affine_boolean_oracle,
    bernstein_vazirani,
    simon_nullspace,
    simon_sample,
)
from oracq.algorithms.common.estimation import (
    amplitude_estimation,
    amplitude_from_phase,
    hadamard_test,
    swap_test,
)
from oracq.algorithms.common.fourier import fourier_add, qft
from oracq.algorithms.common.search import amplify_success
from oracq.algorithms.common.walks import cycle_walk
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    basis_state,
    gate_database,
    uniform_state,
)
from oracq.algorithms.optimization.variational import (
    hardware_efficient_ansatz,
    qaoa_maxcut,
    vqe_measurements,
)
from oracq.algorithms.qec.error_correction import repetition_encode, repetition_recover
from oracq.infrastructure.ir import fuse


def distribution(state, index):
    result = {}
    for key, amplitude in state.amplitudes.items():
        result[key[index]] = result.get(key[index], 0) + abs(amplitude) ** 2
    return result


class AlgorithmExpansionTests(unittest.TestCase):
    def test_legacy_imports_reference_canonical_objects(self):
        from oracq.algorithms.qlss.qlss import make_costa_qlss

        self.assertIs(
            importlib.import_module("oracq.algorithms.costa").make_costa_qlss, make_costa_qlss
        )
        self.assertIs(importlib.import_module("oracq.ir").Bits, Bits)
        self.assertIs(importlib.import_module("oracq.algorithms.elementary").qft, qft)

    def test_qft_matches_positive_fourier_matrix(self):
        size = 8
        operation = qft(3)
        for column in range(size):
            state = simulate(operation.program(), initial={"target": column})
            for row in range(size):
                expected = cmath.exp(2j * math.pi * column * row / size) / math.sqrt(size)
                self.assertAlmostEqual(state.amplitudes.get((row,), 0), expected, places=10)

    def test_fourier_add_all_basis_inputs(self):
        operation = fourier_add(3)
        for a in range(8):
            for b in range(8):
                state = simulate(operation.program(), initial={"a": a, "b": b})
                self.assertAlmostEqual(state.amplitudes.get((a, (a + b) % 8), 0), 1, places=10)

    def test_bernstein_vazirani_recovers_secret_with_affine_bias(self):
        for secret in range(8):
            for bias in (0, 1):
                state = simulate(
                    bernstein_vazirani(affine_boolean_oracle(3, secret, bias=bias)).program()
                )
                self.assertAlmostEqual(distribution(state, 0).get(secret, 0), 1, places=10)

    def test_simon_sampling_and_rank_aware_postprocess(self):
        state = simulate(simon_sample(gate_database(2, 1, [0, 1, 1, 0])).program())
        p = distribution(state, 0)
        self.assertAlmostEqual(p[0], 0.5)
        self.assertAlmostEqual(p[3], 0.5)
        self.assertEqual(simon_nullspace([3], 2), (3,))
        self.assertEqual(len(simon_nullspace([], 3)), 3)
        self.assertEqual(simon_nullspace([1, 2, 4], 3), ())

    def test_amplification_of_known_quarter_probability(self):
        b = Builder("QuarterSuccess", {"target": Bits(1), "signal": Bits(1)})
        b.ry(b["signal"], 2 * math.pi / 3)
        result = amplify_success(StateOracle(b.finish()))
        state = simulate(result.operation.program())
        self.assertAlmostEqual(distribution(state, 1).get(0, 0), 1, places=10)

    def test_hadamard_test_real_and_imaginary(self):
        b = Builder("PhaseWitness", {"target": Bits(1)})
        b.gate("phase", b["target"], 0.6)
        unitary = b.finish()
        for component, expected in (("real", math.cos(0.6)), ("imag", math.sin(0.6))):
            state = simulate(
                hadamard_test(unitary, basis_state(1, 1), component=component).program()
            )
            p = distribution(state, 2)
            self.assertAlmostEqual(p.get(0, 0) - p.get(1, 0), expected, places=10)

    def test_swap_test_overlap(self):
        for first, second, expected in ((0, 0, 1), (0, 1, 0.5)):
            state = simulate(swap_test(basis_state(1, first), basis_state(1, second)).program())
            self.assertAlmostEqual(distribution(state, 4).get(0, 0), expected, places=10)

    def test_amplitude_estimation_half_probability(self):
        operation = amplitude_estimation(uniform_state(1), (1,), precision=3)
        p = distribution(simulate(operation.program()), 2)
        self.assertAlmostEqual(p.get(2, 0) + p.get(6, 0), 1, places=10)
        for phase in (2, 6):
            self.assertAlmostEqual(amplitude_from_phase(phase, 3), 0.5)

    def test_qaoa_single_edge_optimal_layer(self):
        operation = qaoa_maxcut(2, ((0, 1, 1.0),), (math.pi / 2,), (math.pi / 8,))
        p = distribution(simulate(operation.program()), 0)
        self.assertAlmostEqual(p.get(1, 0) + p.get(2, 0), 1, places=10)

    def test_vqe_measurement_of_y_eigenstate(self):
        b = Builder("PlusI", {"target": Bits(1)})
        b.h(b["target"])
        b.gate("phase", b["target"], math.pi / 2)
        coefficient, operation = vqe_measurements(b.finish(), ((2, "Y"),))[0]
        self.assertEqual(coefficient, 2)
        self.assertAlmostEqual(distribution(simulate(operation.program()), 0).get(0, 0), 1)

    def test_ansatz_zero_angles_and_walk_one_step(self):
        operation = hardware_efficient_ansatz(2, (((0, 0), (0, 0)),))
        self.assertAlmostEqual(simulate(operation.program()).amplitudes[(0,)], 1)
        state = simulate(cycle_walk(2).program())
        for key in ((1, 0), (3, 1)):
            self.assertAlmostEqual(state.amplitudes[key], 1 / math.sqrt(2))

    def test_modular_multiplication_total_permutation_and_order(self):
        op = modular_multiply(2, 5)
        for x in range(8):
            expected = 2 * x % 5 if x < 5 else x
            self.assertAlmostEqual(
                simulate(op.program(), initial={"target": x}).amplitudes[(expected,)], 1
            )
        p = distribution(simulate(order_finding(2, 3, precision=2).program()), 1)
        self.assertAlmostEqual(p[0], 0.5)
        self.assertAlmostEqual(p[2], 0.5)
        self.assertEqual(factors_from_phase(1, 2, 2, 15), (3, 5))
        self.assertIsNone(factors_from_phase(0, 2, 2, 15))

    def test_repetition_codes_preserve_arbitrary_logical_amplitudes(self):
        for error in ("bit", "phase"):
            for position in range(3):
                b = Builder(
                    "CodeWitness_" + error + str(position), {"target": Bits(1), "syndrome": Bits(2)}
                )
                b.ry(b["target"], 0.73)
                b.rz(b["target"], 0.29)
                b.call(repetition_encode(error=error), target=b["target"], syndrome=b["syndrome"])
                b.gate("x" if error == "bit" else "z", fuse(b["target"], b["syndrome"])[position])
                b.call(repetition_recover(error=error), target=b["target"], syndrome=b["syndrome"])
                state = simulate(b.finish().program())
                syndromes = {key[1] for key in state.amplitudes}
                self.assertEqual(len(syndromes), 1)
                syndrome = next(iter(syndromes))
                self.assertAlmostEqual(
                    state.amplitudes.get((0, syndrome), 0),
                    cmath.exp(-0.145j) * math.cos(0.365),
                    places=10,
                )
                self.assertAlmostEqual(
                    state.amplitudes.get((1, syndrome), 0),
                    cmath.exp(0.145j) * math.sin(0.365),
                    places=10,
                )

    def test_bad_inputs_fail_at_generation(self):
        cases = [
            lambda: modular_multiply(2, 4),
            lambda: qft(True),
            lambda: qaoa_maxcut(2, ((0, 0, 1),), (1,), (1,)),
            lambda: hardware_efficient_ansatz(2, (((0, 0),),)),
            lambda: repetition_encode(error="unknown"),
        ]
        for case in cases:
            with self.assertRaises(ValidationError):
                case()

    def test_measurement_algorithms_require_clean_preparation(self):
        from oracq.algorithms.input_model.oracles import StatePreparation, annotate

        dirty = StatePreparation(annotate(basis_state(1, work_width=1).operation,
                                         "state_prep_isometry", clean_work=False))
        with self.assertRaisesRegex(ValidationError, "clean_work"):
            swap_test(dirty, basis_state(1))
