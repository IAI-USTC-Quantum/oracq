"""Compare full complex amplitudes on both real PySparQ adapters."""

import sys
import unittest
from pathlib import Path

from oracq import Bits, Builder, compile_qram, run_pysparq, run_pysparq_rir, simulate

sys.path.insert(0, str(Path(__file__).parents[1] / "core"))
from test_qram_compile import rotation_program


class QRAMCompilerNativeTests(unittest.TestCase):
    def compare(self, program):
        compiled, memory, _ = compile_qram(program)
        expected = simulate(program).amplitudes
        for execute in (run_pysparq, run_pysparq_rir):
            actual = execute(compiled, memory).amplitudes
            self.assertLess(
                max(
                    abs(expected.get(k, 0) - actual.get(k, 0))
                    for k in expected.keys() | actual.keys()
                ),
                1e-10,
            )

    def test_entangled_variable_rotation_and_negative_coherent_controls(self):
        for width in [2, 3]:
            self.compare(rotation_program(width, 2, controlled=True))
            self.compare(rotation_program(width, 2, negative=True))

    def test_modular_arithmetic_and_nonzero_xor_targets(self):
        b = Builder("native_boolean", {"word": Bits(4), "target": Bits(1)})
        b.h(b["word"])
        b.h(b["target"])
        with b.control(b["word"], 9):
            b.x(b["target"])
        b.add_const(b["word"].reinterpret("uint"), 7)
        self.compare(b.finish().program())

    def test_multi_controlled_h_and_non_dyadic_rotation(self):
        b = Builder("native_controlled_h", {"controls": Bits(3), "target": Bits(1)})
        b.h(b["controls"])
        with b.control(b["controls"], 2):
            b.h(b["target"])
            b.rz(b["target"], 0.217)
        self.compare(b.finish().program())
