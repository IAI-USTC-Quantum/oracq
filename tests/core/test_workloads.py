"""Paradigm acceptance: generation, open serialization, binding, and backend description; solving accuracy is not checked."""

import unittest

from oracq import ValidationError, bind, dumps, export_originir, loads, unresolved, validate
from oracq.algorithms.input_model.block_encoding import pad_signal
from oracq.algorithms.input_model.operators import identity, scale
from oracq.algorithms.input_model.oracles import abstract_block_encoding, banked_database
from oracq.algorithms.qlss.qlss import dolph_chebyshev_plan
from oracq.applications.catalog import CASES, build_case
from oracq.infrastructure.linking import calls


class WorkloadTests(unittest.TestCase):
    def test_every_catalog_case_has_a_closed_description(self):
        for name in CASES:
            with self.subTest(case=name):
                case = build_case(name)
                opened = loads(dumps(case.program))
                self.assertEqual(
                    [r.name for r in unresolved(opened)], [r.name for r in unresolved(case.program)]
                )
                closed = case.closed()
                validate(closed, require_closed=True)
                artifact = case.artifact()
                self.assertIn("DEF ", artifact.text)
                self.assertIn("QINIT ", artifact.text)
                self.assertFalse(unresolved(closed))

    def test_costa_keeps_filtering_structure_and_only_input_holes(self):
        case = build_case("costa_qram")
        self.assertEqual({r.name for r in unresolved(case.program)}, {"A", "B"})
        filters = [
            m
            for m in case.program.modules
            if dict(m.attributes).get("algorithm") == "coherent_lcu_filter"
        ]
        self.assertEqual(len(filters), 1)
        self.assertIsNotNone(filters[0].body)
        self.assertTrue(any("unary_prepare" in c.module for c in calls(filters[0].body)))
        plan = dolph_chebyshev_plan(4, 0.2)
        self.assertEqual((len(plan.weights), plan.stride, plan.offset), (5, 2, -4))

    def test_qfvm_can_stop_at_physics_without_losing_the_outer_algorithm(self):
        case = build_case("qfvm_qram")
        remaining = {k: v for k, v in case.bindings.items() if k != "QfvmPhysicalEntry"}
        partial = bind(case.program, remaining)
        holes = unresolved(partial)
        self.assertEqual([r.name for r in holes], ["QfvmPhysicalEntry"])
        self.assertGreater(len(partial.modules), 10)
        self.assertGreater(len(partial.main.resources), 0)
        with self.assertRaisesRegex(ValidationError, "QfvmPhysicalEntry"):
            export_originir(partial)

    def test_qham_depends_on_lifted_operator_components(self):
        for name in ("qham_qode", "qham_qpde"):
            with self.subTest(case=name):
                case = build_case(name)
                self.assertEqual(
                    {r.name for r in unresolved(case.program)},
                    {"QhamLinear", "QhamFold", "QhamInitial"},
                )
                self.assertEqual(case.program.main.registers[0].type.width, 1)
                self.assertTrue(
                    any(
                        dict(m.attributes).get("algorithm") == "costa_qlss"
                        for m in case.program.modules
                    )
                )

    def test_banked_words_do_not_violate_register_storage_width(self):
        op = banked_database(2, 96, abstract=True)
        self.assertEqual([r.type.width for r in op.module.registers], [2, 64, 32])
        closed = build_case("banked_qram").closed()
        self.assertEqual(len(closed.main.resources), 2)

    def test_alpha_change_requires_regeneration(self):
        slot = abstract_block_encoding("A", 1, 1, 1.0)
        implementation = pad_signal(scale(2, identity(1)), 1)
        with self.assertRaisesRegex(ValidationError, "be_alpha"):
            bind(slot.operation, {"A": implementation.operation})
