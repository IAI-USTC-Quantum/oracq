"""Structural witnesses and small-instance numerical cross-checks for DQI (decoded quantum interferometry).

The numerical part verifies the conclusion of the paper (arXiv:2408.08292): measuring the
syndrome register yields assignment x with probability proportional to ``K_l(u(x))**2``,
where u(x) is the number of unsatisfied constraints and ``K_l`` is a Krawtchouk polynomial.
"""

import math
import unittest
from math import comb

from oracq import ValidationError, bind, simulate, unresolved
from oracq.algorithms.optimization.dqi import (
    XorSatInstance,
    abstract_decoder,
    bruteforce_decoder,
    dicke_state,
    dqi,
    table_decoder,
)


def krawtchouk(m, weight, j):
    return sum(
        ((-1) ** t) * comb(j, t) * comb(m - j, weight - t) for t in range(weight + 1)
    )


def expected_distribution(instance, weight):
    m, n = instance.num_constraints, instance.num_variables
    weights = {}
    for x in range(1 << n):
        unsatisfied = m - instance.satisfied_count(x)
        weights[x] = krawtchouk(m, weight, unsatisfied) ** 2
    total = sum(weights.values())
    return {x: w / total for x, w in weights.items()}


def distribution(state, index):
    result = {}
    for key, amplitude in state.amplitudes.items():
        result[key[index]] = result.get(key[index], 0) + abs(amplitude) ** 2
    return result


# An n=3 instance using all 7 nonzero rows, with the right-hand side planted from assignment x*=0b101;
# rows are pairwise distinct and nonzero, so weight-1 errors are all uniquely decodable.
PLANTED = XorSatInstance(
    ((0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)),
    (1, 0, 1, 1, 0, 1, 0),
    3,
)

# An m=n=5 instance with B the identity; every weight is uniquely decodable.
IDENTITY = XorSatInstance(tuple((j,) for j in range(5)), (1, 0, 1, 1, 0), 5)


class DqiTests(unittest.TestCase):
    def check_distribution(self, instance, weight):
        operation = dqi(instance, bruteforce_decoder(instance, max_weight=weight), weight=weight)
        state = simulate(operation.program())
        # All errors of weight ≤ weight are uniquely decoded; the error register is deterministically uncomputed.
        self.assertAlmostEqual(distribution(state, 0).get(0, 0), 1, places=12)
        measured = distribution(state, 1)
        expected = expected_distribution(instance, weight)
        self.assertEqual(set(measured), set(range(1 << instance.num_variables)))
        for x, probability in expected.items():
            self.assertAlmostEqual(measured[x], probability, places=10)
        satisfied = sum(
            x * instance.satisfied_count(assignment) for assignment, x in measured.items()
        )
        return satisfied

    def test_planted_instance_beats_random_guessing(self):
        satisfied = self.check_distribution(PLANTED, 1)
        # Random guessing has expected satisfied count m/2 = 3.5; the formula gives 6.5.
        self.assertAlmostEqual(satisfied, 6.5, places=10)
        self.assertGreater(satisfied, PLANTED.num_constraints / 2)

    def test_identity_instance_with_weight_two(self):
        satisfied = self.check_distribution(IDENTITY, 2)
        expected = expected_distribution(IDENTITY, 2)
        reference = sum(p * IDENTITY.satisfied_count(x) for x, p in expected.items())
        self.assertAlmostEqual(satisfied, reference, places=10)
        # The identity decoder gives no advantage: the expectation returns exactly to the random baseline m/2 (floating-point jitter allowed).
        self.assertGreaterEqual(satisfied, IDENTITY.num_constraints / 2 - 1e-9)

    def test_dicke_state_weight_and_uniformity(self):
        m, weight = 5, 2
        state = simulate(dicke_state(m, weight).operation.program())
        self.assertEqual(len(state.amplitudes), comb(m, weight))
        for key, amplitude in state.amplitudes.items():
            self.assertEqual(key[1], 0)
            self.assertEqual(key[0].bit_count(), weight)
            self.assertAlmostEqual(amplitude, 1 / math.sqrt(comb(m, weight)), places=12)

    def test_abstract_decoder_binds_to_witness(self):
        n, m = PLANTED.num_variables, PLANTED.num_constraints
        decoder = abstract_decoder("DqiSyndromeDecoder", n, m)
        operation = dqi(PLANTED, decoder, weight=1)
        attrs = dict(operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "dqi")
        self.assertEqual(attrs["field"], "GF(2)")
        self.assertEqual(attrs["num_constraints"], m)
        self.assertEqual(attrs["num_variables"], n)
        self.assertEqual(attrs["dicke_weight"], 1)
        program = operation.program()
        self.assertEqual([r.name for r in unresolved(program)], ["DqiSyndromeDecoder"])
        bound = bind(program, {"DqiSyndromeDecoder": bruteforce_decoder(PLANTED, max_weight=1).operation})
        self.assertFalse(unresolved(bound))
        measured = distribution(simulate(bound), 1)
        for x, probability in expected_distribution(PLANTED, 1).items():
            self.assertAlmostEqual(measured[x], probability, places=10)

    def test_invalid_inputs_fail_at_generation(self):
        n, m = PLANTED.num_variables, PLANTED.num_constraints
        cases = [
            lambda: XorSatInstance(((0, 3),), (0,), 2),
            lambda: XorSatInstance(((0, 0),), (0,), 2),
            lambda: XorSatInstance(((0,),), (0, 1), 1),
            lambda: XorSatInstance(((0,),), (2,), 1),
            lambda: XorSatInstance((), (), 2),
            lambda: dicke_state(4, 5),
            lambda: table_decoder(2, 2, {4: 1}),
            lambda: table_decoder(2, 2, {1: 4}),
            lambda: dqi(PLANTED, bruteforce_decoder(PLANTED), weight=m + 1),
            lambda: dqi(PLANTED, table_decoder(n + 1, m, {}), weight=1),
            lambda: dqi(PLANTED, object(), weight=1),
            lambda: dqi(PLANTED, bruteforce_decoder(PLANTED), weight=-1),
        ]
        for case in cases:
            with self.assertRaises(ValidationError):
                case()


if __name__ == "__main__":
    unittest.main()
