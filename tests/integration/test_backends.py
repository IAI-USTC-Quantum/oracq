"""真实 UnifiedQuantum 与 PySparQ 对拍；用具有这两个依赖的解释器运行。"""

import cmath
import math
import unittest

from pyqecclang import (
    QRAM,
    Bits,
    Builder,
    UInt,
    export_originir,
    identity,
    linear_combination,
    pauli_x,
    run_originir,
    run_pysparq,
    run_pysparq_rir,
    simulate,
)


class BackendTests(unittest.TestCase):
    def test_native_register_word_64_and_total_width_over_64(self):
        b = Builder("wide", {"word": UInt(64), "flag": Bits(2)})
        b.x(b["word"])
        b.add_const(b["word"], 1)
        b.x(b["word"][63])
        b.x(b["flag"][1])
        p = b.finish().program()
        self.assertEqual(run_pysparq(p).amplitudes, simulate(p).amplitudes)

    def test_pysparq_preserves_existing_registry(self):
        import pysparq as ps

        from pyqecclang import ValidationError

        ps.System.add_register("existing", ps.StateStorageType.General, 1)
        try:
            with self.assertRaisesRegex(ValidationError, "非空"):
                run_pysparq(identity(1).operation.program())
            self.assertEqual(ps.System.get_activated_register_size(), 1)
        finally:
            ps.System.clear()

    def test_pysparq_budget_cleans_its_registry(self):
        import pysparq as ps

        from pyqecclang import ValidationError

        b = Builder("budget", {"q": Bits(5)})
        b.h(b["q"])
        with self.assertRaisesRegex(ValidationError, "数量预算"):
            run_pysparq(b.finish().program(), max_states=4)
        self.assertEqual(ps.System.get_activated_register_size(), 0)

    def test_signed_and_rational_storage(self):
        from pyqecclang import Rational, SInt

        b = Builder("storage", {"s": SInt(3), "r": Rational(2)})
        b.x(b["s"])
        b.ry(b["r"], 0.41)
        self.compare(b.finish().program())

    def compare(self, program, memory=None):
        expected = simulate(program, memory).statevector()
        native = run_pysparq(program, memory).statevector()
        origin = list(run_originir(program, memory))
        self.assertEqual(len(expected), len(native))
        self.assertEqual(len(expected), len(origin))
        for index, (a, b, c) in enumerate(zip(expected, native, origin, strict=True)):
            with self.subTest(index=index):
                self.assertAlmostEqual(a, b, places=11)
                self.assertAlmostEqual(a, complex(c), places=11)

    def test_nested_def_qram_arbitrary_data_and_lsb(self):
        leaf = Builder("load", {"a": UInt(2), "d": UInt(3)}, {"table": QRAM(2, 3)})
        leaf.qram("table", leaf["a"], leaf["d"])
        op = leaf.finish()
        wrapper = Builder("wrap", {"a": UInt(2), "d": UInt(3)}, {"table": QRAM(2, 3)})
        wrapper.call(op, a=wrapper["a"], d=wrapper["d"], resources={"table": "table"})
        main = Builder("main", {"a": UInt(2), "d": UInt(3)}, {"data": QRAM(2, 3)})
        main.h(main["a"])
        main.x(main["d"][0])
        main.call(wrapper.finish(), a=main["a"], d=main["d"], resources={"table": "data"})
        p = main.finish().program()
        self.assertEqual(export_originir(p).text.count("DEF "), 3)
        self.compare(p, {"data": [1, 2, 4, 7]})

    def test_qram_views_and_zero_control(self):
        b = Builder("main", {"word": Bits(7)}, {"mem": QRAM(2, 3)})
        b.h(b["word"][:2])
        b.h(b["word"][6])
        b.x(b["word"][3])
        with b.control(b["word"][6], 0):
            b.qram("mem", b["word"][:2], b["word"][3:6])
        self.compare(b.finish().program(), {"mem": [3, 1, 7, 5]})

    def test_repeat_adjoint_arithmetic(self):
        step = Builder("step", {"w": UInt(3)})
        step.add_const(step["w"], 3)
        step.rz(step["w"][1], 0.37)
        op = step.finish()
        b = Builder("main", {"w": UInt(3), "ctrl": Bits(1)})
        b.h(b["w"])
        b.h(b["ctrl"])
        with b.control(b["ctrl"], 1):
            with b.repeat(5):
                b.call(op, w=b["w"])
            with b.adjoint():
                with b.repeat(3):
                    b.call(op, w=b["w"])
        self.compare(b.finish().program())

    def test_all_one_bit_gates_and_global_phase(self):
        b = Builder("main", {"q": Bits(3)})
        for op in ["h", "x", "y", "z", "s", "t", "rx", "ry", "rz", "phase"]:
            angle = 0.23 if op in {"rx", "ry", "rz", "phase"} else None
            b.gate(op, b["q"][0], angle)
            with b.adjoint():
                b.gate(op, b["q"][1], angle)
        b.h(b["q"][2])
        with b.control(b["q"][2], 0):
            b.global_phase(-0.48)
        b.global_phase(0.12)
        self.compare(b.finish().program())

    def test_complex_lcu(self):
        be = linear_combination(1j, identity(1), -2, pauli_x(1))
        self.compare(be.operation.program())

    def test_qram_resource_specialization(self):
        leaf = Builder("load", {"a": UInt(1), "d": UInt(2)}, {"m": QRAM(1, 2)})
        leaf.qram("m", leaf["a"], leaf["d"])
        op = leaf.finish()
        b = Builder(
            "main", {"a": UInt(1), "d": UInt(2)}, {"first": QRAM(1, 2), "second": QRAM(1, 2)}
        )
        b.h(b["a"])
        for resource in ["first", "second"]:
            b.call(op, a=b["a"], d=b["d"], resources={"m": resource})
        self.compare(b.finish().program(), {"first": [1, 3], "second": [2, 0]})

    def test_origin_global_phase_all_control_bits(self):
        b = Builder("main", {"q": Bits(2)})
        b.h(b["q"])
        with b.control(b["q"], 2):
            b.global_phase(math.pi / 4)
        p = b.finish().program()
        expected = list(run_originir(p))
        self.assertAlmostEqual(expected[2], 0.5 * cmath.exp(1j * math.pi / 4))
        self.compare(p)

    def test_pysparq_native_rir_cross_validation(self):
        """PySparQ 原生 RIR 解释器与既有三条路径独立实现对拍。"""
        leaf = Builder("load", {"a": UInt(2), "d": UInt(3)}, {"table": QRAM(2, 3)})
        leaf.qram("table", leaf["a"], leaf["d"])
        op = leaf.finish()
        step = Builder("step", {"w": UInt(3)})
        step.add_const(step["w"], 3)
        step.rz(step["w"][1], 0.37)
        step_op = step.finish()
        b = Builder("main", {"a": UInt(2), "d": UInt(3), "w": UInt(3), "ctrl": Bits(1)},
                    {"data": QRAM(2, 3)})
        b.h(b["a"])
        b.h(b["ctrl"])
        b.call(op, a=b["a"], d=b["d"], resources={"table": "data"})
        with b.control(b["ctrl"], 1):
            with b.repeat(5):
                b.call(step_op, w=b["w"])
            with b.adjoint():
                with b.repeat(3):
                    b.call(step_op, w=b["w"])
        p = b.finish().program()
        memory = {"data": [1, 2, 4, 7]}
        expected = simulate(p, memory).amplitudes
        native = run_pysparq_rir(p, memory).amplitudes
        self.assertEqual(set(expected), set(native))
        for key in expected:
            with self.subTest(key=key):
                self.assertAlmostEqual(expected[key], native[key], places=11)


if __name__ == "__main__":
    unittest.main()
