"""quantikz 导出的结构与文本断言：结构保持 RIR，不展开调用。"""

import math
import unittest

from pyqecclang import QRAM, Bits, Builder, UInt, bind, quantikz
from pyqecclang.algorithms.basics.oracle_algorithms import affine_boolean_oracle, bernstein_vazirani
from pyqecclang.algorithms.input_model.oracles import abstract_database


class QuantikzTests(unittest.TestCase):
    def test_bound_bv_keeps_oracle_box(self):
        given = abstract_database("BooleanFunction", 3, 1)
        closed = bind(
            bernstein_vazirani(given).program(),
            {"BooleanFunction": affine_boolean_oracle(3, secret=5, bias=1).operation},
        )
        text = quantikz(closed)
        self.assertIn("\\begin{quantikz}", text)
        self.assertIn("\\gate{H}", text)
        self.assertIn("\\mathtt{BooleanFunction}", text)
        self.assertIn("\\gate[wires=4]", text)
        self.assertNotIn("dashed", text.split("\\mathtt{BooleanFunction}")[0])

    def test_open_call_is_dashed(self):
        given = abstract_database("BooleanFunction", 1, 1)
        opened = bernstein_vazirani(given).program()
        text = quantikz(opened)
        line = next(line for line in text.splitlines() if "BooleanFunction" in line)
        self.assertIn("dashed", line)

    def test_repeat_adjoint_control_and_xor(self):
        b = Builder("demo", {"q": Bits(2)})
        with b.repeat(3):
            b.h(b["q"][0])
            b.xor(b["q"][0], b["q"][1])
        with b.control(b["q"][0]):
            b.x(b["q"][1])
        with b.control(b["q"][1], 0):
            b.x(b["q"][0])
        with b.adjoint():
            b.h(b["q"][0])
        text = quantikz(b.finish())
        self.assertIn("\\times 3", text)
        self.assertIn("\\ctrl{1}", text)
        self.assertIn("\\octrl{", text)
        self.assertIn("\\targ{}", text)
        self.assertIn("\\dagger", text)
        self.assertIn("remember picture", text)
        self.assertIn("fit=", text)

    def test_qram_load_and_power_of_two_repeat(self):
        b = Builder("rep", {"q": Bits(1), "a": UInt(2), "d": UInt(3)}, {"table": QRAM(2, 3)})
        with b.repeat(1 << 40):
            b.x(b["q"])
        b.qram("table", b["a"], b["d"])
        text = quantikz(b.finish())
        self.assertIn("2^{40}", text)
        self.assertIn("\\mathrm{QRAM}_{\\mathtt{table}}", text)
        self.assertEqual(text.count("\\times"), 1)

    def test_angles_and_add_const_and_phase(self):
        b = Builder("rot", {"q": Bits(1), "acc": UInt(3)})
        b.ry(b["q"], math.pi / 4)
        b.rz(b["q"], -math.pi / 2)
        b.add_const(b["acc"], 5)
        b.global_phase(math.pi)
        text = quantikz(b.finish())
        self.assertIn("R_y(\\frac{\\pi}{4})", text)
        self.assertIn("R_z(-\\frac{\\pi}{2})", text)
        self.assertIn("+5", text)
        self.assertIn("e^{i\\pi}", text)

    def test_method_form_matches_function(self):
        given = abstract_database("BooleanFunction", 1, 1)
        operation = bernstein_vazirani(given)
        self.assertEqual(
            operation.quantikz(name="bvcheck"), quantikz(operation, name="bvcheck")
        )


if __name__ == "__main__":
    unittest.main()
