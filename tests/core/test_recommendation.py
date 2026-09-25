"""Numerical validation of the KP quantum recommendation system: exact rank-one reference and threshold filtering."""

import math
import unittest

from oracq import ValidationError, dumps, estimate_resources, loads, simulate
from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.input_model.qdata import QMatrix
from oracq.algorithms.qml.recommendation import (
    KPRecommendationConfig,
    kp_recommendation,
    sigma_from_phase,
)

FMT = FixedFormat(8, 4)
AW = 12
CONFIG = KPRecommendationConfig(precision=4, sigma=0.45)


def normalized(vector):
    norm = math.sqrt(sum(x * x for x in vector))
    return [x / norm for x in vector]


class SigmaDecodingTests(unittest.TestCase):
    def test_mirror_phases_share_singular_value(self):
        for t in range(16):
            self.assertAlmostEqual(
                sigma_from_phase(t, 4, 1.0), sigma_from_phase((16 - t) % 16, 4, 1.0)
            )

    def test_endpoints(self):
        self.assertAlmostEqual(sigma_from_phase(0, 4, 2.0), 2.0)
        self.assertAlmostEqual(sigma_from_phase(8, 4, 2.0), 0.0)
        with self.assertRaises(ValidationError):
            sigma_from_phase(16, 4, 1.0)


class RecommendationTests(unittest.TestCase):
    def test_two_by_two_rank_one_is_exact(self):
        matrix = QMatrix([[0.5, 0.25], [0.5, 0.25]], fmt=FMT, angle_width=AW)
        result = kp_recommendation(matrix, 0, KPRecommendationConfig(precision=4, sigma=0.5))
        success, distribution = result.readout(
            simulate(result.operation.program(), result.memories())
        )
        v1 = normalized((0.5, 0.25))
        self.assertAlmostEqual(success, 1.0, delta=2e-3)
        for item, probability in sorted(distribution.items()):
            self.assertAlmostEqual(probability, v1[item] ** 2, delta=2e-3)

    def test_four_by_four_rank_one_for_every_user(self):
        base = (0.5, 0.25, 0.25, 0.0)
        matrix = QMatrix(
            [[scale * value for value in base] for scale in (1.0, 0.5, 0.5, 1.0)],
            fmt=FMT,
            angle_width=AW,
        )
        v1 = normalized(base)
        for user in range(4):
            with self.subTest(user=user):
                result = kp_recommendation(matrix, user, CONFIG)
                success, distribution = result.readout(
                    simulate(result.operation.program(), result.memories())
                )
                self.assertAlmostEqual(success, 1.0, delta=2e-3)
                for item, probability in sorted(distribution.items()):
                    self.assertAlmostEqual(probability, v1[item] ** 2, delta=3e-3)

    def test_threshold_filters_small_singular_values(self):
        base = (0.5, 0.25, 0.25, 0.0)
        scales = (1.0, 0.5, 0.5, 1.0)
        eps = 0.0625
        matrix = QMatrix(
            [
                [scales[i] * base[j] + (eps if i == j else 0.0) for j in range(4)]
                for i in range(4)
            ],
            fmt=FMT,
            angle_width=AW,
        )
        result = kp_recommendation(matrix, 1, KPRecommendationConfig(precision=5, sigma=0.45))
        success, distribution = result.readout(
            simulate(result.operation.program(), result.memories())
        )
        v1 = normalized(base)
        self.assertGreater(success, 0.9)
        self.assertLess(success, 1.0)
        for item, probability in sorted(distribution.items()):
            self.assertAlmostEqual(probability, v1[item] ** 2, delta=0.02)

    def test_validation(self):
        matrix = QMatrix([[0.5, 0.25], [0.5, 0.25]], fmt=FMT, angle_width=AW)
        with self.assertRaises(ValidationError):
            kp_recommendation(matrix, 0, KPRecommendationConfig(precision=0))
        with self.assertRaises(ValidationError):
            kp_recommendation(matrix, 0, KPRecommendationConfig(precision=13))
        with self.assertRaises(ValidationError):
            kp_recommendation(matrix, 0, KPRecommendationConfig(sigma=2.0 * matrix.frobenius))
        with self.assertRaises(ValidationError):
            kp_recommendation([[1, 2], [3, 4]], 0)

    def test_program_closes_and_estimates(self):
        matrix = QMatrix([[0.5, 0.25], [0.5, 0.25]], fmt=FMT, angle_width=AW)
        result = kp_recommendation(matrix, 0, KPRecommendationConfig(precision=3, sigma=0.5))
        program = result.operation.program()
        self.assertEqual(loads(dumps(program)), program)
        estimate = estimate_resources(program)
        self.assertGreater(estimate.qram_total, 0)


if __name__ == "__main__":
    unittest.main()
