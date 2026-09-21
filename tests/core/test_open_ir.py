"""开放 IR 和绑定验证。只验证范式、ABI、闭合与导出，不认证算法数值。"""

import unittest

from pyqecclang import (
    Binding,
    Bits,
    Builder,
    ValidationError,
    bind,
    dumps,
    export_originir,
    loads,
    unresolved,
)
from pyqecclang.algorithms.input_model.oracles import (
    abstract_database,
    abstract_state_prep,
    basis_state,
    qram_database,
)


class OpenIRTests(unittest.TestCase):
    def test_open_roundtrip_and_export_boundary(self):
        slot = abstract_database("Function", 2, 1)
        from pyqecclang.algorithms.basics.oracle_algorithms import deutsch_jozsa

        p = deutsch_jozsa(slot).program()
        restored = loads(dumps(p))
        self.assertEqual(unresolved(restored)[0].name, "Function")
        self.assertIsNone(restored.module_map["Function"].body)
        with self.assertRaisesRegex(ValidationError, "Function"):
            export_originir(restored)

    def test_partial_binding_and_lifted_resources(self):
        first, second = abstract_database("First", 1, 1), abstract_database("Second", 1, 1)
        b = Builder("Top", {"a": Bits(1), "d": Bits(1)})
        for slot in (first, second):
            b.call(slot.operation, address=b["a"], data=b["d"])
        p = b.finish().program()
        q = bind(p, {"First": Binding(qram_database(1, 1).operation, {"table": "one"})})
        self.assertEqual([r.name for r in unresolved(q)], ["Second"])
        r = bind(q, {"Second": Binding(qram_database(1, 1).operation, {"table": "two"})})
        self.assertFalse(unresolved(r))
        self.assertEqual({r.name for r in r.main.resources}, {"one", "two"})
        self.assertEqual(export_originir(r).text.count("QRAMDECL"), 2)

    def test_shared_capture_has_one_entry_resource(self):
        a, bslot = abstract_database("A", 1, 1), abstract_database("B", 1, 1)
        b = Builder("Top", {"a": Bits(1), "d": Bits(1)})
        for slot in (a, bslot):
            b.call(slot.operation, address=b["a"], data=b["d"])
        impl = qram_database(1, 1).operation
        q = bind(
            b.finish().program(),
            {"A": Binding(impl, {"table": "shared"}), "B": Binding(impl, {"table": "shared"})},
        )
        self.assertEqual(len(q.main.resources), 1)

    def test_shape_mismatch_does_not_mutate_source(self):
        slot = abstract_state_prep("B", 1)
        original = dumps(slot.operation.program())
        with self.assertRaisesRegex(ValidationError, "宽度"):
            bind(slot.operation, {"B": basis_state(2).operation})
        self.assertEqual(dumps(slot.operation.program()), original)

    def test_identity_is_not_an_unfinished_operation(self):
        b = Builder("Identity", {"q": Bits(1)}).finish()
        self.assertEqual(b.module.body, ())
        self.assertFalse(unresolved(b.program()))

    def test_isometry_requires_declared_adjoint_capability(self):
        slot = abstract_state_prep("Prep", 1, reversible=False)
        b = Builder("Use", {"q": Bits(1), "work": Bits(0)})
        with b.adjoint():
            b.call(slot.operation, target=b["q"], work=b["work"])
        with self.assertRaisesRegex(ValidationError, "supports_adjoint"):
            b.finish()

    def test_unused_open_declaration_does_not_block_entry(self):
        from dataclasses import replace

        slot = abstract_state_prep("Unused", 1)
        p = basis_state(1).operation.program()
        q = replace(p, modules=p.modules + (slot.operation.module,))
        self.assertFalse(unresolved(q))
        self.assertTrue(export_originir(q).text)
