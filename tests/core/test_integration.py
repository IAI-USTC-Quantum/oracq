"""Heinrich 量子求和/积分的数值见证：比较器构造、QAE 读出、三层绑定一致性与收敛率。"""

import unittest

from pyqecclang import ValidationError, bind, simulate, unresolved
from pyqecclang.algorithms.integration import (
    heinrich_rate,
    integral_from_phase,
    mean_from_phase,
    quantum_integral,
    quantum_sum,
    sum_preparation,
    table_loader,
)
from pyqecclang.algorithms.oracles import abstract_database


def mode(distribution):
    return max(distribution, key=distribution.get)


def phase_distribution(state):
    """phase 是最后一个公开寄存器。"""
    result = {}
    for key, amplitude in state.amplitudes.items():
        result[key[-1]] = result.get(key[-1], 0.0) + abs(amplitude) ** 2
    return result


class SumPreparationTests(unittest.TestCase):
    def test_flag_probability_matches_mean(self):
        # 好状态概率恰为 E[v]/2**w（线性构造的核心恒等式）。
        values = (0, 1, 2, 3, 4, 5, 6, 7)
        loader = table_loader(values)
        prep = sum_preparation(loader.database)
        state = simulate(prep.program())
        n, w = 3, 3
        flag_bit = 1 << (n + w)
        probability = sum(
            abs(amp) ** 2
            for key, amp in state.amplitudes.items()
            if key[0] & flag_bit
        )
        self.assertAlmostEqual(probability, sum(values) / len(values) / (1 << w), places=12)

    def test_constant_table_exact(self):
        loader = table_loader([5] * 8, data_width=3)
        prep = sum_preparation(loader.database)
        state = simulate(prep.program())
        flag_bit = 1 << 6
        probability = sum(
            abs(amp) ** 2 for key, amp in state.amplitudes.items() if key[0] & flag_bit
        )
        self.assertAlmostEqual(probability, 5 / 8, places=12)

    def test_qram_binding_matches_gate(self):
        values = (3, 1, 4, 1, 5, 9, 2, 6)
        gate = sum_preparation(table_loader(values).database)
        qram_bundle = table_loader(values, backend="qram")
        qram = sum_preparation(qram_bundle.database)
        flag_bit = 1 << 6
        gate_probability = sum(
            abs(a) ** 2
            for k, a in simulate(gate.program()).amplitudes.items()
            if k[0] & flag_bit
        )
        qram_probability = sum(
            abs(a) ** 2
            for k, a in simulate(qram.program(), memory=qram_bundle.memory).amplitudes.items()
            if k[0] & flag_bit
        )
        self.assertAlmostEqual(gate_probability, qram_probability, places=12)

    def test_abstract_loader_binds(self):
        values = (1, 2, 3, 4)
        abstract = abstract_database("SumLoader", 2, 3)
        prep = sum_preparation(abstract)
        self.assertEqual([r.name for r in unresolved(prep.program())], ["SumLoader"])
        bound = bind(
            prep.program(), {"SumLoader": table_loader(values).database.operation}
        )
        self.assertFalse(unresolved(bound))
        flag_bit = 1 << 5
        probability = sum(
            abs(a) ** 2 for k, a in simulate(bound).amplitudes.items() if k[0] & flag_bit
        )
        self.assertAlmostEqual(probability, sum(values) / 4 / 8, places=12)


class QuantumSumTests(unittest.TestCase):
    def test_mean_on_qae_grid_is_exact(self):
        # E[v]/2**w = 1/2 落在 QAE 栅格上，读出确定。
        loader = table_loader([2] * 8, data_width=2)
        operation = quantum_sum(loader.database, precision=4)
        attrs = dict(operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "quantum_sum")
        measured = phase_distribution(simulate(operation.program()))
        for value, probability in measured.items():
            if probability > 1e-9:
                self.assertAlmostEqual(
                    mean_from_phase(value, 4, 2), 2.0, places=9
                )

    def test_ramp_mean_within_qae_resolution(self):
        values = tuple(range(8))
        loader = table_loader(values)
        precision = 5
        operation = quantum_sum(loader.database, precision=precision)
        measured = phase_distribution(simulate(operation.program()))
        estimate = mean_from_phase(mode(measured), precision, 3)
        self.assertAlmostEqual(estimate, 3.5, delta=0.5)

    def test_quantum_integral_trapezoid_scale(self):
        # f(x) = x 在 [0,1] 上 8 个中点网格：积分真值 0.5。
        grid = 8
        data_width = 4
        values = tuple(round((i + 0.5) / grid * ((1 << data_width) - 1)) for i in range(grid))
        loader = table_loader(values, data_width=data_width)
        precision = 4
        operation = quantum_integral(loader.database, precision=precision, interval=1.0)
        attrs = dict(operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "quantum_integral")
        measured = phase_distribution(simulate(operation.program()))
        estimate = integral_from_phase(mode(measured), precision, data_width)
        exact = sum(values) / grid / ((1 << data_width) - 1)
        self.assertAlmostEqual(estimate, exact, delta=0.06)
        self.assertAlmostEqual(estimate, 0.5, delta=0.09)


class RateTests(unittest.TestCase):
    def test_heinrich_rate_known_values(self):
        rates = heinrich_rate(1, 1)
        self.assertEqual(rates["deterministic"], 1.0)
        self.assertEqual(rates["randomized"], 1.5)
        self.assertEqual(rates["quantum"], 2.0)
        rates = heinrich_rate(2, 4)
        self.assertEqual(
            rates["quantum"] - rates["randomized"],
            rates["randomized"] - rates["deterministic"],
        )

    def test_invalid_inputs_fail_at_generation(self):
        loader = table_loader([1, 2, 3])
        cases = [
            lambda: table_loader([]),
            lambda: table_loader([1, -1]),
            lambda: table_loader([1.5, 2]),
            lambda: table_loader([8], data_width=3),
            lambda: table_loader([1], backend="memory"),
            lambda: quantum_sum(loader.database, precision=0),
            lambda: quantum_integral(loader.database, interval=0.0),
            lambda: heinrich_rate(0, 1),
            lambda: heinrich_rate(1, 0),
            lambda: mean_from_phase(0, 3, 0),
            lambda: sum_preparation(object()),
        ]
        for case in cases:
            with self.subTest(case=case):
                with self.assertRaises(ValidationError):
                    case()


if __name__ == "__main__":
    unittest.main()
