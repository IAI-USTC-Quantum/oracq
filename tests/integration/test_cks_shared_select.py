"""Execute shared CKS SELECT through the actual PySparQ native gate adapter."""

import math
import unittest

from oracq import Bits, Builder, run_pysparq, simulate
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import annotate
from oracq.algorithms.qlss.qlss import CKSConfig, cks_inverse_block


class NativeSharedSelectTests(unittest.TestCase):
    def test_complex_coherent_controls_and_padded_selectors(self):
        for terms in (3, 4, 7):
            matrix = Builder("native_involution", {"target": Bits(1), "signal": Bits(1)})
            matrix.ry(matrix["signal"], 0.71)
            with matrix.control(matrix["signal"], 0):
                matrix.z(matrix["target"])
            with matrix.control(matrix["signal"]):
                matrix.x(matrix["target"])
            matrix.ry(matrix["signal"], -0.71)
            a = BlockEncoding(annotate(matrix.finish(), "block_encoding", be_alpha=1.0,
                                       self_adjoint_extension=True))
            operation = cks_inverse_block(a, CKSConfig(order=8, terms=terms)).operation
            b = Builder("coherent_native_select", {r.name: r.type for r in operation.module.registers} |
                        {"control": Bits(1)})
            for name in ("target", "signal", "control"):
                b.h(b[name])
            b.rz(b["target"], math.pi/2)
            with b.control(b["control"], 0):
                b.call(operation, target=b["target"], signal=b["signal"])
            program = b.finish().program()
            expected, actual = simulate(program), run_pysparq(program)
            with self.subTest(terms=terms):
                self.assertEqual(expected.registers, actual.registers)
                for key in expected.amplitudes.keys() | actual.amplitudes.keys():
                    self.assertAlmostEqual(expected.amplitudes.get(key, 0),
                                           actual.amplitudes.get(key, 0), places=10)


if __name__ == "__main__":
    unittest.main()
