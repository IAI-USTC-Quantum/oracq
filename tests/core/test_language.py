"""语言核心与 BE 组合的行为验证，全部测试只需要 Python 标准库。"""

import json
import math
import unittest

from pyqecclang import (
    QRAM,
    Bits,
    BlockEncoding,
    Builder,
    Call,
    Module,
    Operation,
    Primitive,
    Program,
    Rational,
    Ref,
    Register,
    SInt,
    Span,
    UInt,
    ValidationError,
    block_encoding,
    dumps,
    export_originir,
    fuse,
    identity,
    linear_combination,
    loads,
    pauli_x,
    product,
    scale,
    simulate,
    validate,
    zero,
)


def lookup():
    b = Builder("lookup", {"a": UInt(2), "d": UInt(3)}, {"mem": QRAM(2, 3)})
    b.qram("mem", b["a"], b["d"])
    return b.finish()


class LanguageTests(unittest.TestCase):
    def test_huge_empty_module_call_still_has_execution_cost(self):
        empty = Builder("empty", {"q": Bits(1)}).finish()
        b = Builder("main", {"q": Bits(1)})
        with b.repeat(1 << 40):
            b.call(empty, q=b["q"])
        program = b.finish().program()
        self.assertLess(len(export_originir(program).text), 50000)
        with self.assertRaisesRegex(ValidationError, "预算"):
            simulate(program)

    def test_mutable_instruction_operands_are_rejected(self):
        b = Builder("main", {"q": Bits(1)})
        b.emit(Primitive("x", [b["q"]]))
        with self.assertRaisesRegex(ValidationError, "不可变"):
            b.finish()

    def test_sparse_state_budget(self):
        b = Builder("main", {"q": Bits(8)})
        b.h(b["q"])
        with self.assertRaisesRegex(ValidationError, "数量预算"):
            simulate(b.finish().program(), max_states=4)

    def test_register_broadcast_is_one_ir_node(self):
        b = Builder("main", {"word": UInt(64)})
        b.x(b["word"])
        p = b.finish().program()
        self.assertEqual(len(p.main.body), 1)
        self.assertEqual(simulate(p).amplitudes, {((1 << 64) - 1,): 1 + 0j})

    def test_add_wraps_at_64_bits(self):
        b = Builder("main", {"word": UInt(64)})
        b.add_const(b["word"], 1)
        p = b.finish().program()
        self.assertEqual(simulate(p, initial={"word": (1 << 64) - 1}).amplitudes, {(0,): 1 + 0j})

    def test_more_than_64_total_qubits(self):
        b = Builder("main", {"a": UInt(64), "b": Bits(2)})
        b.x(b["a"][63])
        b.x(b["b"][1])
        self.assertEqual(simulate(b.finish().program()).amplitudes, {(1 << 63, 2): 1 + 0j})

    def test_register_width_and_kind_validation(self):
        for dtype in [Bits(65), UInt(-1), Bits(True), Bits(2.5)]:
            with self.subTest(dtype=dtype), self.assertRaises(ValidationError):
                Builder("bad", {"q": dtype}).finish()

    def test_all_storage_kinds_have_raw_bit_semantics(self):
        b = Builder("main", {"a": SInt(3), "b": Rational(3)})
        b.x(b["a"])
        b.x(b["b"])
        self.assertEqual(simulate(b.finish().program()).amplitudes, {(7, 7): 1 + 0j})

    def test_cross_register_view_and_slice(self):
        b = Builder("main", {"a": Bits(3), "b": Bits(3)})
        view = fuse(b["b"][1:], b["a"][:2])
        b.x(view[1:3])
        self.assertEqual(simulate(b.finish().program()).amplitudes, {(1, 4): 1 + 0j})

    def test_view_rejects_overlap(self):
        b = Builder("main", {"a": Bits(3)})
        b.x(fuse(b["a"][:2], b["a"][1:]))
        with self.assertRaisesRegex(ValidationError, "重叠"):
            b.finish()

    def test_view_rejects_invalid_slices(self):
        b = Builder("main", {"a": Bits(3)})
        for key in [-1, 3, slice(0, 4), slice(None, None, 2)]:
            with self.subTest(key=key), self.assertRaises(ValidationError):
                _ = b["a"][key]

    def test_binary_alias(self):
        b = Builder("main", {"a": Bits(3)})
        b.xor(b["a"][:2], b["a"][1:])
        with self.assertRaisesRegex(ValidationError, "重叠"):
            b.finish()

    def test_protected_control(self):
        b = Builder("main", {"q": Bits(3)})
        with b.control(b["q"][0]):
            b.x(b["q"])
        with self.assertRaisesRegex(ValidationError, "控制"):
            b.finish()

    def test_control_value_bounds(self):
        b = Builder("main", {"q": Bits(3)})
        with b.control(b["q"][:1], 2):
            b.x(b["q"][1])
        with self.assertRaises(ValidationError):
            b.finish()

    def test_adjoint_reverses_order_and_angles(self):
        b = Builder("main", {"q": Bits(1)})
        b.ry(b["q"], 0.7)
        b.rz(b["q"], -1.3)
        with b.adjoint():
            b.ry(b["q"], 0.7)
            b.rz(b["q"], -1.3)
        state = simulate(b.finish().program())
        self.assertAlmostEqual(state.amplitudes[(0,)], 1 + 0j)

    def test_repeat_preserves_graph_and_export_is_logarithmic(self):
        b = Builder("main", {"word": UInt(8)})
        with b.repeat(1 << 40):
            b.add_const(b["word"], 1)
        p = b.finish().program()
        original = dumps(p)
        text = export_originir(p).text
        self.assertEqual(dumps(p), original)
        self.assertEqual(len(p.main.body), 1)
        self.assertLess(text.count("DEF "), 45)
        self.assertLess(len(text), 50000)
        with self.assertRaisesRegex(ValidationError, "预算"):
            simulate(p)

    def test_empty_repeat_does_not_iterate(self):
        b = Builder("main", {"q": Bits(1)})
        with b.repeat((1 << 63) - 1):
            pass
        self.assertEqual(simulate(b.finish().program()).amplitudes, {(0,): 1 + 0j})

    def test_builder_block_rollback(self):
        b = Builder("main", {"q": Bits(1)})
        try:
            with b.repeat(3):
                b.x(b["q"])
                raise RuntimeError("abort")
        except RuntimeError:
            pass
        self.assertEqual(len(b.finish().module.body), 0)

    def test_builder_freezes_after_finish(self):
        b = Builder("main", {"q": Bits(1)})
        b.finish()
        with self.assertRaises(ValidationError):
            b.x(b["q"])

    def test_qram_nonzero_data_in_superposition(self):
        p = lookup().program()
        state = simulate(p, {"mem": [1, 2, 4, 7]}, initial={"a": 3, "d": 5})
        self.assertEqual(state.amplitudes, {(3, 2): 1 + 0j})

    def test_qram_double_load_is_identity(self):
        b = Builder("main", {"a": UInt(2), "d": UInt(3)}, {"table": QRAM(2, 3)})
        b.h(b["a"])
        b.x(b["d"])
        for _ in range(2):
            b.call(lookup(), a=b["a"], d=b["d"], resources={"mem": "table"})
        state = simulate(b.finish().program(), {"table": [1, 2, 4, 7]})
        for a in range(4):
            self.assertAlmostEqual(state.amplitudes[(a, 7)], 0.5)

    def test_qram_data_errors(self):
        p = lookup().program()
        for memory in [None, {}, {"mem": {4: 0}}, {"mem": {1: 8}}, {"mem": {True: 1}}]:
            with self.subTest(memory=memory), self.assertRaises(ValidationError):
                simulate(p, memory)

    def test_qram_shape_error(self):
        b = Builder("main", {"q": UInt(4)}, {"m": QRAM(2, 3)})
        b.qram("m", b["q"][:2], b["q"][2:])
        with self.assertRaises(ValidationError):
            b.finish()

    def test_call_shape_and_resource_checks(self):
        b = Builder("main", {"a": Bits(2), "d": UInt(3)}, {"table": QRAM(2, 3)})
        b.call(lookup(), a=b["a"], d=b["d"], resources={"mem": "table"})
        with self.assertRaisesRegex(ValidationError, "类型"):
            b.finish()

    def test_same_module_specialized_for_two_qrams(self):
        b = Builder(
            "main", {"a": UInt(2), "d": UInt(3)}, {"first": QRAM(2, 3), "second": QRAM(2, 3)}
        )
        for mem in ["first", "second"]:
            b.call(lookup(), a=b["a"], d=b["d"], resources={"mem": mem})
        p = b.finish().program()
        text = export_originir(p).text
        self.assertEqual(len(p.modules), 2)
        self.assertEqual(text.count("DEF "), 3)
        self.assertEqual(simulate(p, {"first": [1], "second": [3]}).amplitudes, {(0, 2): 1 + 0j})

    def test_call_cycle_rejected(self):
        ref = Ref((Span("q", 0, 1),), Bits(1))
        m = Module("a", (Register("q", Bits(1)),), (), (Call("a", (ref,)),))
        with self.assertRaisesRegex(ValidationError, "递归"):
            validate(Program("a", (m,)))

    def test_conflicting_modules_rejected(self):
        a = Builder("same", {"q": Bits(1)}).finish()
        b = Builder("same", {"q": Bits(1)})
        b.x(b["q"])
        c = Builder("main", {"q": Bits(1)})
        c.call(a, q=c["q"])
        with self.assertRaisesRegex(ValidationError, "冲突"):
            c.call(b.finish(), q=c["q"])

    def test_json_roundtrip_retains_modules(self):
        b = Builder("main", {"a": UInt(2), "d": UInt(3)}, {"table": QRAM(2, 3)})
        with b.repeat(123):
            b.call(lookup(), a=b["a"], d=b["d"], resources={"mem": "table"})
        p = b.finish().program()
        text = dumps(p)
        self.assertEqual(dumps(loads(text)), text)
        self.assertEqual(loads(text), p)

    def test_json_rejects_unknown_version_fields_and_duplicates(self):
        text = dumps(identity(1).operation.program())
        data = json.loads(text)
        for key, value in [("version", "9.0"), ("extra", 1), ("tag", "Unknown")]:
            mutated = dict(data)
            mutated[key] = value
            with self.subTest(key=key), self.assertRaises(ValidationError):
                loads(json.dumps(mutated))
        with self.assertRaisesRegex(ValidationError, "重复"):
            loads('{"tag":"Program","tag":"Program"}')

    def test_unknown_primitive_and_nan_rejected(self):
        for op, angle in [("unknown", None), ("ry", float("nan"))]:
            b = Builder("main", {"q": Bits(1)})
            b.emit(Primitive(op, (b["q"],), angle))
            with self.assertRaises(ValidationError):
                b.finish()

    def test_global_phase_is_preserved(self):
        b = Builder("main", {"q": Bits(1)})
        b.h(b["q"])
        with b.control(b["q"], 0):
            b.global_phase(math.pi / 3)
        state = simulate(b.finish().program())
        self.assertAlmostEqual(state.amplitudes[(1,)], 1 / math.sqrt(2))
        self.assertAlmostEqual(state.amplitudes[(0,)].imag, math.sqrt(3 / 8))


def matrix(be):
    columns = []
    for value in range(1 << be.width):
        state = simulate(be.operation.program(), initial={"target": value})
        columns.append([state.amplitudes.get((i, 0), 0j) * be.alpha for i in range(1 << be.width)])
    return [[columns[j][i] for j in range(len(columns))] for i in range(len(columns))]


class BlockEncodingTests(unittest.TestCase):
    def assertMatrix(self, actual, expected):
        for arow, erow in zip(actual, expected, strict=True):
            for a, e in zip(arow, erow, strict=True):
                self.assertAlmostEqual(a, e)

    def test_unequal_normalization_and_signed_sum(self):
        a, b = scale(2, identity(1)), scale(3, pauli_x(1))
        combined = linear_combination(-0.5, a, 2, b)
        self.assertEqual(combined.alpha, 7)
        self.assertMatrix(matrix(combined), [[-1, 6], [6, -1]])

    def test_complex_coefficients(self):
        combined = linear_combination(1j, identity(1), -0.5, pauli_x(1))
        self.assertMatrix(matrix(combined), [[1j, -0.5], [-0.5, 1j]])

    def test_product_order_and_independent_signal_spaces(self):
        b = Builder("h_be", {"target": Bits(1), "signal": Bits(1)})
        b.h(b["target"])
        b.ry(b["signal"], 0.8)
        a = block_encoding(b.finish(), 2)
        result = product(a, pauli_x(1))
        k = 2 * math.cos(0.4) / math.sqrt(2)
        self.assertMatrix(matrix(result), [[k, k], [-k, k]])

    def test_zero_and_zero_coefficient(self):
        self.assertMatrix(matrix(zero(1)), [[0, 0], [0, 0]])
        self.assertMatrix(
            matrix(linear_combination(0, identity(1), 0, pauli_x(1))), [[0, 0], [0, 0]]
        )

    def test_alpha_survives_ir_serialization(self):
        be = scale(4, identity(2))
        p = loads(dumps(be.operation.program()))
        restored = BlockEncoding(
            Operation(p.main, tuple(m for m in p.modules if m.name != p.entry))
        )
        self.assertEqual(restored.alpha, 4)

    def test_invalid_alpha(self):
        for alpha in [0, -1, True, float("inf")]:
            with self.subTest(alpha=alpha), self.assertRaises(ValidationError):
                block_encoding(identity(1).operation, alpha)


if __name__ == "__main__":
    unittest.main()
