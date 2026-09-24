"""DQI（解码量子干涉优化）的结构见证与小实例数值对拍。

数值部分验证论文（arXiv:2408.08292）的结论：测量 syndrome 寄存器得到赋值 x
的概率正比于 ``K_l(u(x))**2``，其中 u(x) 是未满足约束数、``K_l`` 是 Krawtchouk
多项式。
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


# 全部 7 个非零行的 n=3 实例，右端项由赋值 x*=0b101 植入；
# 行两两不同且非零，因此权重 1 的错误都可唯一译码。
PLANTED = XorSatInstance(
    ((0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)),
    (1, 0, 1, 1, 0, 1, 0),
    3,
)

# B 取单位阵的 m=n=5 实例，任意权重都可唯一译码。
IDENTITY = XorSatInstance(tuple((j,) for j in range(5)), (1, 0, 1, 1, 0), 5)


class DqiTests(unittest.TestCase):
    def check_distribution(self, instance, weight):
        operation = dqi(instance, bruteforce_decoder(instance, max_weight=weight), weight=weight)
        state = simulate(operation.program())
        # 所有权重不超过 weight 的错误都被唯一译码，error 寄存器确定性复净。
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
        # 随机猜测的期望满足数是 m/2 = 3.5；公式给出 6.5。
        self.assertAlmostEqual(satisfied, 6.5, places=10)
        self.assertGreater(satisfied, PLANTED.num_constraints / 2)

    def test_identity_instance_with_weight_two(self):
        satisfied = self.check_distribution(IDENTITY, 2)
        expected = expected_distribution(IDENTITY, 2)
        reference = sum(p * IDENTITY.satisfied_count(x) for x, p in expected.items())
        self.assertAlmostEqual(satisfied, reference, places=10)
        # identity 译码器不提供优势，期望恰好回到随机基线 m/2（允许浮点抖动）。
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
