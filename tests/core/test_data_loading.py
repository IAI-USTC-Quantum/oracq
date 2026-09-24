"""Select-Swap QROM 数据加载的数学见证：逐地址对拍、资源公式与复净检查。"""

import unittest

from oracq import Bits, Builder, ValidationError, simulate
from oracq.algorithms.input_model.data_loading import (
    QromCost,
    qrom_cost,
    qrom_lookup,
    select_swap_qrom,
)
from oracq.algorithms.input_model.oracles import gate_database

TABLE16 = (3, 0, 5, 2, 7, 1, 6, 4, 0, 2, 1, 7, 5, 3, 6, 4)


def probe(operation, address, data, address_bits, data_bits):
    """在基态初值 ``|address, data>`` 上执行数据库调用并返回输出寄存器对。"""
    b = Builder(
        f"Probe_{operation.module.name}_{address}_{data}",
        {"address": Bits(address_bits), "data": Bits(data_bits)},
    )
    b.call(operation, address=b["address"], data=b["data"])
    state = simulate(b.finish().program(), initial={"address": address, "data": data})
    (key, amplitude), = state.amplitudes.items()
    assert abs(amplitude - 1) < 1e-12
    return key


class DataLoadingTests(unittest.TestCase):
    def test_select_swap_matches_gate_database_per_address(self):
        baseline = gate_database(4, 3, TABLE16).operation
        expected = {
            address: probe(baseline, address, 0, 4, 3) for address in range(16)
        }
        for partitions in (1, 2, 4, 8, 16):
            database = select_swap_qrom(TABLE16, partitions=partitions)
            self.assertEqual(database.address_width, 4)
            self.assertEqual(database.data_width, 3)
            actual = {
                address: probe(database.operation, address, 0, 4, 3)
                for address in range(16)
            }
            self.assertEqual(actual, expected)
            self.assertEqual(
                [actual[a] for a in range(16)], [(a, TABLE16[a]) for a in range(16)]
            )

    def test_xor_semantics_with_nonzero_data_register(self):
        database = select_swap_qrom(TABLE16, partitions=4).operation
        for address in range(16):
            out_address, out_data = probe(database, address, 5, 4, 3)
            self.assertEqual(out_address, address)
            self.assertEqual(out_data, TABLE16[address] ^ 5)

    def test_superposition_query_restores_clean_work(self):
        database = select_swap_qrom(TABLE16, partitions=4).operation
        b = Builder("SuperpositionProbe", {"address": Bits(4), "data": Bits(3)})
        b.h(b["address"])
        b.call(database, address=b["address"], data=b["data"])
        state = simulate(b.finish().program())
        # 地址叠加上的查询保持地址不变并逐点 XOR 表字；窗口局部寄存器由模拟器
        # 在 LocalExit 处强制复净，此处地址边缘分布均匀即说明 compute/uncompute 干净。
        self.assertEqual(len(state.amplitudes), 16)
        for address in range(16):
            amplitude = state.amplitudes.get((address, TABLE16[address]), 0)
            self.assertAlmostEqual(abs(amplitude) ** 2, 1 / 16, places=12)

    def test_result_invariant_across_partitions_and_qrom_lookup_baseline(self):
        expected = [probe(qrom_lookup(TABLE16).operation, a, 0, 4, 3)[1] for a in range(16)]
        self.assertEqual(expected, list(TABLE16))
        for partitions in (1, 2, 4, 8, 16):
            database = select_swap_qrom(TABLE16, partitions=partitions).operation
            actual = [probe(database, a, 0, 4, 3)[1] for a in range(16)]
            self.assertEqual(actual, expected)

    def test_cost_model_matches_formulas_and_tradeoff_curve(self):
        cost = qrom_cost(16, 3, 4)
        self.assertIsInstance(cost, QromCost)
        self.assertEqual(cost.address_bits, 4)
        self.assertEqual(cost.select_toffoli, 3)
        self.assertEqual(cost.swap_toffoli, 9)
        self.assertEqual(cost.t_count, 4 * (3 + 9))
        self.assertEqual(cost.round_trip_t_count, 2 * cost.t_count)
        self.assertEqual(cost.work_qubits, 12)
        self.assertEqual(cost.fanout_qubits, 6)
        self.assertEqual(cost.ancilla_qubits, 18)
        self.assertEqual(cost.t_depth, 3 + 3)
        curve = {p: qrom_cost(16, 3, p).t_count for p in (1, 2, 4, 8, 16)}
        self.assertEqual(curve[1], 4 * 15)
        self.assertEqual(curve[16], 4 * 45)
        self.assertLess(curve[4], curve[1])
        self.assertLess(curve[4], curve[16])
        self.assertEqual(
            curve, {p: 4 * ((16 // p) - 1 + 3 * (p - 1)) for p in curve}
        )

    def test_generated_operations_carry_cost_attributes(self):
        database = select_swap_qrom(TABLE16, partitions=4)
        attrs = dict(database.operation.module.attributes)
        cost = qrom_cost(16, 3, 4)
        self.assertEqual(attrs["implementation"], "select_swap")
        self.assertEqual(attrs["qrom_partitions"], 4)
        self.assertEqual(attrs["select_toffoli"], cost.select_toffoli)
        self.assertEqual(attrs["swap_toffoli"], cost.swap_toffoli)
        self.assertEqual(attrs["t_count"], cost.t_count)
        self.assertEqual(attrs["work_qubits"], cost.work_qubits)
        self.assertEqual(attrs["dirty_fanout_qubits"], cost.fanout_qubits)
        baseline = qrom_lookup(TABLE16)
        baseline_attrs = dict(baseline.operation.module.attributes)
        self.assertEqual(baseline_attrs["implementation"], "qrom_unary_iteration")
        self.assertEqual(baseline_attrs["qrom_partitions"], 1)
        self.assertEqual(baseline_attrs["t_count"], qrom_cost(16, 3, 1).t_count)

    def test_cost_report_is_structured(self):
        report = qrom_cost(1024, 32, 8).to_dict()
        self.assertEqual(report["n_addresses"], 1024)
        self.assertEqual(report["partitions"], 8)
        self.assertEqual(report["t_count"], 4 * (127 + 32 * 7))
        self.assertEqual(report["ancilla_qubits"], 8 * 32 + 7 * 3)

    def test_invalid_inputs_fail_at_generation(self):
        cases = [
            lambda: select_swap_qrom([1, 0], partitions=3),
            lambda: select_swap_qrom([1, 0], partitions=0),
            lambda: select_swap_qrom([1, 0], partitions=4),
            lambda: select_swap_qrom([], partitions=1),
            lambda: select_swap_qrom([1, -2], partitions=1),
            lambda: select_swap_qrom([1, 4], partitions=1, data_bits=2),
            lambda: select_swap_qrom({-1: 3}, partitions=1),
            lambda: select_swap_qrom([1.5, 0], partitions=1),
            lambda: qrom_lookup([]),
            lambda: qrom_lookup([1, 2], data_bits=1),
            lambda: qrom_cost(16, 3, 3),
            lambda: qrom_cost(16, 3, 32),
            lambda: qrom_cost(0, 3, 1),
            lambda: qrom_cost(16, 0, 1),
        ]
        for case in cases:
            with self.assertRaises(ValidationError):
                case()


if __name__ == "__main__":
    unittest.main()
