"""第二阶段结构测试；不宣称数值算法已获验证。"""

import unittest
from dataclasses import replace

from pyqecclang import Bits, Builder, dumps, loads
from pyqecclang.arithmetic import BooleanNetwork, FixedFormat, fixed_arithmetic
from pyqecclang.backends.basis import export_toffoli_u3_cz
from pyqecclang.execution import simulate
from pyqecclang.ir import ValidationError
from pyqecclang.layout import workspace_table


class Stage2StructureTests(unittest.TestCase):
    def test_private_workspace_roundtrip_and_lifecycle(self):
        child = Builder("private_child", {"a": Bits(1)})
        tmp = child.local("scratch", Bits(1))
        child.xor(child["a"], tmp)
        child.xor(child["a"], tmp)
        child = child.finish()
        root = Builder("private_root", {"a": Bits(1)})
        root.x(root["a"])
        with root.repeat(3):
            root.call(child, a=root["a"])
        program = root.finish().program()
        self.assertEqual(loads(dumps(program)), program)
        self.assertEqual(workspace_table(program)[program.entry], 1)
        self.assertEqual(simulate(program).amplitudes, {(1,): 1 + 0j})
        with self.assertRaises(ValidationError):
            dumps(replace(program, version="0.2"))

    def test_private_workspace_does_not_replace_public_interface(self):
        b = Builder("private_only", {})
        b.local("scratch", Bits(1))
        with self.assertRaisesRegex(ValidationError, "非空量子接口"):
            b.finish()

    def test_dirty_workspace_rejected_at_return(self):
        b = Builder("dirty", {"a": Bits(1)})
        b.x(b.local("scratch", Bits(1)))
        with self.assertRaises(ValidationError):
            simulate(b.finish().program())

    def test_arithmetic_construction_and_basis(self):
        fmt = FixedFormat(4, 1)
        for kind in (
            "add",
            "sub",
            "neg",
            "abs",
            "mul",
            "div",
            "reciprocal",
            "sqrt",
            "lt",
            "eq",
            "select",
            "and",
            "or",
            "xor",
        ):
            with self.subTest(kind=kind):
                op = fixed_arithmetic(kind, fmt)
                self.assertTrue(op.module.locals)
                self.assertEqual(loads(dumps(op.program())), op.program())
                text = export_toffoli_u3_cz(op.program()).text
                for line in text.splitlines():
                    if line.startswith(("DEF ", "m_", "QINIT ", "CREG ")) or line in {
                        "ENDDEF",
                        "DAGGER",
                        "ENDDAGGER",
                    }:
                        continue
                    self.assertIn(line.split()[0], {"TOFFOLI", "U3", "CZ"})
                    self.assertNotIn("controlled_by", line)

    def test_gate_network_is_executable_small_example(self):
        op = fixed_arithmetic("add", FixedFormat(2, 0, False))
        b = Builder("add_example", {r.name: r.type for r in op.module.registers})
        b.x(b["a"][0])
        b.x(b["b"][1])
        b.call(op, **{r.name: b[r.name] for r in op.module.registers})
        result = simulate(b.finish().program(), max_steps=10000)
        self.assertEqual(result.amplitudes, {(1, 2, 3, 0): 1 + 0j})
        net = BooleanNetwork.from_payload(dict(op.module.attributes)["arithmetic_network"])
        self.assertEqual(net.evaluate(a=1, b=2), {"out": 3, "status": 0})
