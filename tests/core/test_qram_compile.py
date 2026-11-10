"""Exact coherent QRAM assisted compilation, including workspace cleanup."""

import math
import unittest

from oracq import Bits, Builder, compile_qram, dumps, estimate_resources, loads, simulate


def rotation_program(width=3, count=3, controlled=False, negative=False):
    b = Builder(
        f"rotation_words_{width}", {"angle": Bits(width), "target": Bits(1), "enabled": Bits(1)}
    )
    b.h(b["angle"])
    b.h(b["target"])
    b.h(b["enabled"])

    def rotations():
        for _ in range(count):
            for bit in range(width):
                with b.control(b["angle"][bit]):
                    b.ry(
                        b["target"],
                        (-1 if negative else 1) * 2 * math.pi * (1 << bit) / (1 << width),
                    )
            b.rz(b["target"], 0.17)

    if controlled:
        with b.control(b["enabled"], 0):
            rotations()
    else:
        rotations()
    return b.finish().program()


class QRAMCompilationTests(unittest.TestCase):
    def compare(self, program, memory=None, **options):
        compiled, extra, report = compile_qram(program, **options)
        actual_memory = {**(memory or {}), **extra}
        before = simulate(program, memory).amplitudes
        after = simulate(compiled, actual_memory).amplitudes
        self.assertLess(
            max(abs(before.get(k, 0) - after.get(k, 0)) for k in before.keys() | after.keys()),
            1e-10,
        )
        self.assertEqual(loads(dumps(compiled, format="json")), compiled)
        return compiled, report

    def test_shared_gradient_restores_zero_with_entangled_angle_and_negative_control(self):
        for negative in [False, True]:
            with self.subTest(negative=negative):
                compiled, report = self.compare(
                    rotation_program(controlled=True, negative=negative)
                )
                self.assertEqual(report["phase_gradient_bits"], 4)
                self.assertEqual(report["rotation_word_groups"], 3)
                self.assertTrue(report["gradient_preparation_and_cleanup_in_circuit"])
                self.assertEqual(estimate_resources(compiled).toffoli, 0)

    def test_boolean_truth_table_preserves_nonzero_target_and_negative_controls(self):
        b = Builder("boolean", {"control": Bits(5), "target": Bits(2)})
        b.h(b["control"])
        b.h(b["target"])
        with b.control(b["control"], 9):
            b.x(b["target"])
        b.add_const(b["control"].reinterpret("uint"), 11)
        with b.control(b["target"], 2):
            b.swap(b["control"][:2], b["control"][2:4])
        compiled, report = self.compare(b.finish().program())
        self.assertGreater(report["boolean_xor_replacements"], 0)
        self.assertEqual(estimate_resources(compiled).toffoli, 0)
        self.assertEqual(report["truth_table_words"], 4)

    def test_repeat_and_adjoint_remain_compressed(self):
        b = Builder("compressed", {"angle": Bits(3), "target": Bits(1), "enabled": Bits(1)})
        op = rotation_program(3, 1)
        from oracq import Operation

        with b.repeat(2):
            b.call(
                Operation.from_program(op),
                angle=b["angle"],
                target=b["target"],
                enabled=b["enabled"],
            )
        with b.adjoint():
            b.call(
                Operation.from_program(op),
                angle=b["angle"],
                target=b["target"],
                enabled=b["enabled"],
            )
        compiled, _ = self.compare(b.finish().program())
        self.assertIn('"count": 2', dumps(compiled, format="json"))

    def test_non_dyadic_rotations_and_modes_are_preserved(self):
        for phase_gradient, boolean in [(True, False), (False, True), (False, False)]:
            with self.subTest(phase_gradient=phase_gradient, boolean=boolean):
                self.compare(rotation_program(2, 1), phase_gradient=phase_gradient, boolean=boolean)

    def test_idempotency_and_readonly_boundary(self):
        p, memory, _ = compile_qram(rotation_program(2, 1))
        same, extra, report = compile_qram(p)
        self.assertEqual(p, same)
        self.assertEqual(extra, {})
        self.assertFalse(report["applied"])
        self.assertEqual(memory, {"__qr_truth": {3: 1}})

    def test_readonly_namespace_and_controlled_nonclassical_gates(self):
        from oracq import QRAM, ValidationError

        b = Builder("write", {"a": Bits(2), "v": Bits(1)}, {"bank": QRAM(2, 1)})
        b.store("bank", b["a"], b["v"])
        with self.assertRaises(ValidationError):
            compile_qram(b.finish().program())
        b = Builder("__qr_conflict", {"target": Bits(1)})
        with self.assertRaises(ValidationError):
            compile_qram(b.finish().program())
        b = Builder("controlled_h", {"controls": Bits(3), "target": Bits(1)})
        b.h(b["controls"])
        with b.control(b["controls"], 2):
            b.h(b["target"])
            b.rz(b["target"], 0.217)
        self.compare(b.finish().program())

    def test_shared_gradient_aligns_different_angle_word_widths(self):
        from oracq import Operation

        b = Builder(
            "mixed_words",
            {
                "small": Bits(2),
                "large": Bits(3),
                "target": Bits(1),
                "enabled": Bits(1),
                "outer": Bits(1),
            },
        )
        b.h(b["outer"])
        for width, key in [(2, "small"), (3, "large")]:
            with b.control(b["outer"], 0):
                b.call(
                    Operation.from_program(rotation_program(width, 1, negative=(width == 2))),
                    angle=b[key],
                    target=b["target"],
                    enabled=b["enabled"],
                )
        _, report = self.compare(b.finish().program())
        self.assertEqual(report["phase_gradient_bits"], 4)
