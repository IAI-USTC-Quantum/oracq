"""qmem 指针式寻址的真实后端对拍；用具有 uniqc 和 pysparq 的解释器运行。"""

import unittest

from pyqecclang import (
    QRAM,
    Builder,
    QMem,
    UInt,
    ValidationError,
    run_originir,
    run_pysparq,
    run_pysparq_rir,
    simulate,
)


class QmemBackendTests(unittest.TestCase):
    def compare(self, program, memory=None):
        expected = simulate(program, memory).statevector()
        native = run_pysparq(program, memory).statevector()
        origin = list(run_originir(program, memory))
        rir = run_pysparq_rir(program, memory).statevector()
        for index, (a, b, c, d) in enumerate(zip(expected, native, origin, rir, strict=True)):
            with self.subTest(index=index):
                self.assertAlmostEqual(a, b, places=11)
                self.assertAlmostEqual(a, complex(c), places=11)
                self.assertAlmostEqual(a, d, places=11)

    def test_constant_offset_pointer_load(self):
        b = Builder("const_off", {"out": UInt(4)}, {"rom": QRAM(4, 4)})
        mem = QMem(b, "rom")
        (mem.ptr() + 5).load(b["out"])
        self.compare(b.finish().program(), {"rom": list(range(16))})

    def test_quantum_pointer_and_offset(self):
        b = Builder("qptr", {"idx": UInt(2), "out": UInt(4)}, {"rom": QRAM(2, 4)})
        mem = QMem(b, "rom")
        b.h(b["idx"])
        (mem.ptr(b["idx"]) + 1).load(b["out"])
        self.compare(b.finish().program(), {"rom": [3, 5, 7, 9]})

    def test_pointer_chain_with_ripple_adder(self):
        b = Builder("chain", {"i": UInt(2), "j": UInt(2), "out": UInt(4)}, {"rom": QRAM(2, 4)})
        mem = QMem(b, "rom")
        b.h(b["i"])
        b.h(b["j"])
        (mem.ptr() + b["i"] + b["j"]).load(b["out"])
        self.compare(b.finish().program(), {"rom": [3, 5, 7, 9]})

    def test_two_dimensional_grid(self):
        b = Builder("grid", {"row": UInt(2), "col": UInt(2), "out": UInt(4)}, {"rom": QRAM(4, 4)})
        grid = QMem(b, "rom", shape=(4, 4))
        b.h(b["row"])
        b.h(b["col"])
        grid[b["row"], b["col"]].load(b["out"])
        self.compare(b.finish().program(), {"rom": list(range(16))})

    def test_slice_view(self):
        b = Builder("slice", {"col": UInt(2), "out": UInt(4)}, {"rom": QRAM(4, 4)})
        grid = QMem(b, "rom", shape=(4, 4))
        b.h(b["col"])
        grid[2:][1, b["col"]].load(b["out"])
        self.compare(b.finish().program(), {"rom": list(range(16))})

    def test_store_is_rejected_by_real_backends(self):
        b = Builder("wr", {"a": UInt(2), "v": UInt(4)}, {"ram": QRAM(2, 4)})
        QMem(b, "ram")[b["a"]].store(b["v"])
        program = b.finish().program()
        # 存储单元按经典单元建模：文本执行器（UnifiedQuantum、PySparQ）暂不接受运行期写，
        # 必须在执行入口明确报错；参考模拟器仍然可以执行。
        simulate(program, {"ram": [0, 0, 0, 0]})
        with self.assertRaisesRegex(ValidationError, "UnifiedQuantum 执行暂不支持"):
            run_originir(program, {"ram": [0, 0, 0, 0]})
        with self.assertRaisesRegex(ValidationError, "PySparQ 适配器暂不支持"):
            run_pysparq(program, {"ram": [0, 0, 0, 0]})
        with self.assertRaisesRegex(ValidationError, "PySparQ RIR 解释器暂不支持"):
            run_pysparq_rir(program, {"ram": [0, 0, 0, 0]})


if __name__ == "__main__":
    unittest.main()
