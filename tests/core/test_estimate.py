"资源估计与 strict 网表的手算计数及交叉验证。"

import math
import unittest
from collections import Counter

from pyqecclang import (
    QRAM,
    Bits,
    Builder,
    UInt,
    estimate_resources,
    export_strict,
)


def _strict_atom_counts(artifact):
    "展开 strict 网表（单 DEF 平铺程序）逐行计数。"
    counts, qram = Counter(), Counter()
    for line in artifact.text.splitlines():
        keyword = line.split(" ", 1)[0]
        if keyword.startswith("ram_"):
            qram[keyword] += 1
        elif keyword in {"TOFFOLI", "CZ", "H", "S", "SDG", "T", "TDG", "X", "Y", "Z", "RZ", "RY"}:
            counts[keyword.lower()] += 1
    return counts, qram


class HandCountTests(unittest.TestCase):
    def test_single_atoms(self):
        b = Builder("atoms", {"q": Bits(4)})
        b.x(b["q"][0])
        b.gate("t", b["q"][1])
        b.h(b["q"][2])
        b.rz(b["q"][3], 0.37)
        estimate = estimate_resources(b.finish().program())
        self.assertEqual(estimate.atoms, Counter({"x": 3, "t": 1, "h": 1}))
        # rz(0.37) 全局相位锚定额外贡献 2 个相位原子与 2 个 X
        self.assertEqual(len(estimate.rotations), 3)
        self.assertEqual(estimate.rotations[0][0], "rz")
        self.assertEqual(estimate.qubits, 4)
        self.assertEqual(estimate.qram_total, 0)
        self.assertEqual(
            estimate.t_total(1e-10), 1 + 3 * estimate.synthesis_t_per_rotation(1e-10)
        )

    def test_add_const_staircase(self):
        b = Builder("add", {"w": UInt(3)})
        b.add_const(b["w"], 1)
        estimate = estimate_resources(b.finish().program())
        # 阶梯：i=2 需要两位进位控制（Toffoli），i=1 一个 CNOT，末位一个 X
        self.assertEqual(estimate.atoms, Counter({"toffoli": 1, "h": 2, "cz": 1, "x": 1}))
        b = Builder("add7", {"w": UInt(3)})
        b.add_const(b["w"], 7)
        estimate = estimate_resources(b.finish().program())
        self.assertEqual(estimate.atoms, Counter({"toffoli": 1, "h": 4, "cz": 2, "x": 3}))

    def test_controlled_gate_uses_controlled_u3_recipe(self):
        b = Builder("ch", {"q": Bits(2)})
        with b.control(b["q"][1], 1):
            b.h(b["q"][0])
        estimate = estimate_resources(b.finish().program())
        self.assertEqual(estimate.toffoli, 0)
        self.assertEqual(len(estimate.rotations), 0)
        self.assertEqual(estimate.atoms["cz"], 2)
        self.assertEqual(estimate.t_exact, 2)  # ry(±π/4) 各贡献一个 T/TDG
        self.assertEqual(estimate.mcx_ancilla, 0)

    def test_multi_control_toffoli_ladder(self):
        b = Builder("ccx", {"c": Bits(3), "t": Bits(1)})
        with b.control(b["c"], 7):
            b.x(b["t"])
        estimate = estimate_resources(b.finish().program())
        self.assertEqual(estimate.toffoli, 2 * 3 - 3)
        self.assertEqual(estimate.mcx_ancilla, 2)

    def test_zero_control_flip_cost(self):
        b = Builder("zc", {"c": Bits(2), "t": Bits(1)})
        with b.control(b["c"], 1):  # 01：一个零位，发射行组两侧各翻一次
            b.x(b["t"])
        estimate = estimate_resources(b.finish().program())
        self.assertEqual(estimate.atoms, Counter({"x": 2, "toffoli": 1}))

    def test_symbolic_repeat_multiplies(self):
        b = Builder("rep", {"w": UInt(2)})
        with b.repeat(2**40):
            b.add_const(b["w"], 1)
        estimate = estimate_resources(b.finish().program())
        base = Counter({"h": 2, "cz": 1, "x": 1})
        self.assertEqual(estimate.atoms, Counter({k: v * 2**40 for k, v in base.items()}))
        self.assertEqual(estimate.gate_total, 4 * 2**40)

    def test_qram_query_count(self):
        b = Builder("qram", {"a": UInt(2), "d": UInt(3)}, {"mem": QRAM(2, 3)})
        b.qram("mem", b["a"], b["d"])
        estimate = estimate_resources(b.finish().program())
        self.assertEqual(estimate.qram_queries, Counter({"mem": 1}))
        self.assertEqual(estimate.gate_total, 0)
        b2 = Builder("qram2", {"a": UInt(2), "d": UInt(3)}, {"mem": QRAM(2, 3)})
        with b2.repeat(17):
            b2.qram("mem", b2["a"], b2["d"])
        estimate = estimate_resources(b2.finish().program())
        self.assertEqual(estimate.qram_total, 17)

    def test_gphase_classification(self):
        b = Builder("gp", {"q": Bits(1)})
        b.global_phase(math.pi / 4)
        estimate = estimate_resources(b.finish().program())
        self.assertEqual(estimate.t_exact, 2)
        self.assertEqual(estimate.atoms["x"], 2)
        self.assertEqual(len(estimate.rotations), 0)

    def test_open_module_rejected(self):
        from pyqecclang import ValidationError
        from pyqecclang.algorithms.input_model.oracles import abstract_database

        oracle = abstract_database("open_db", 2, 2)
        open_main = Builder("main", {"q": Bits(4)})
        open_main.call(oracle.operation, address=open_main["q"][:2], data=open_main["q"][2:])
        with self.assertRaises(ValidationError):
            estimate_resources(open_main.finish().program())

    def test_operation_method_and_to_dict(self):
        b = Builder("meth", {"w": UInt(2)})
        b.h(b["w"])
        op = b.finish()
        estimate = op.estimate()
        self.assertEqual(estimate.atoms["h"], 2)
        payload = estimate.to_dict()
        self.assertEqual(payload["qubits"], 2)
        self.assertEqual(payload["clifford"], 2)
        self.assertGreaterEqual(payload["t_total"], payload["t_exact"])


class StrictNetlistCrossTests(unittest.TestCase):
    def compare(self, program):
        estimate = estimate_resources(program)
        artifact = export_strict(program)
        counts, qram = _strict_atom_counts(artifact)
        strict_rotations = counts.pop("rz", 0) + counts.pop("ry", 0)
        estimate_atoms = {k: v for k, v in estimate.atoms.items() if v}
        self.assertEqual(estimate_atoms, {k: v for k, v in counts.items() if v})
        self.assertEqual(len(estimate.rotations), strict_rotations)
        self.assertEqual(estimate.qram_total, sum(qram.values()))

    def test_flat_gates_and_rotations(self):
        b = Builder("flat", {"q": Bits(3)})
        b.h(b["q"][0])
        b.gate("t", b["q"][1])
        b.ry(b["q"][2], 0.61)
        b.gate("phase", b["q"][0], math.pi / 2)
        self.compare(b.finish().program())

    def test_flat_add_const_controlled(self):
        b = Builder("flatadd", {"w": UInt(3), "c": Bits(2)})
        with b.control(b["c"], 3):
            b.add_const(b["w"], 5)
        program = b.finish().program()
        self.compare(program)
        estimate = estimate_resources(program)
        qinit = next(
            int(line.split()[1])
            for line in export_strict(program).text.splitlines()
            if line.startswith("QINIT ")
        )
        self.assertEqual(estimate.qubits + estimate.mcx_ancilla, qinit)

    def test_flat_xor_swap(self):
        b = Builder("flatxs", {"a": UInt(2), "b": UInt(2), "c": Bits(1)})
        with b.control(b["c"], 1):
            b.xor(b["a"], b["b"])
            b.swap(b["a"], b["b"])
        self.compare(b.finish().program())

    def test_flat_qram_and_gphase(self):
        b = Builder("flatq", {"a": UInt(2), "d": UInt(2)}, {"mem": QRAM(2, 2)})
        b.h(b["a"])
        b.qram("mem", b["a"], b["d"])
        b.global_phase(0.23)
        self.compare(b.finish().program())

    def test_controlled_rotation_counts(self):
        b = Builder("cr", {"c": Bits(1), "q": Bits(1)})
        with b.control(b["c"], 1):
            b.rz(b["q"], math.pi / 3)
        program = b.finish().program()
        estimate = estimate_resources(program)
        self.compare(program)
        # RZ 的 (θ,φ,λ)=(0,0,λ) 参数化下受控形式只需两个 rz 旋转原子
        self.assertEqual(len(estimate.rotations), 2)
        self.assertEqual(estimate.atoms, Counter({"h": 4, "cz": 2}))


if __name__ == "__main__":
    unittest.main()
