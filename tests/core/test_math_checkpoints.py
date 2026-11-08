"""Checkpoint circuits preserve exact fixed-point values, flags, and clean workspace."""

import unittest

from oracq import (
    Bits,
    Builder,
    FixedFormat,
    ValidationError,
    dumps,
    estimate_resources,
    loads,
    simulate,
)
from oracq.infrastructure.mathfunc import Index, MathConfig, compile_function

SOURCE = """def f(x,y):
 a=x+y
 b=x-y
 c=a*b
 d=c+1
 return d if x<y else 1/(b+1)
"""


class MathCheckpointTests(unittest.TestCase):
    def test_all_small_words_with_nonzero_outputs_and_path_flags(self):
        reference = compile_function(SOURCE, fmt=FixedFormat(3, 1), config=MathConfig(checkpoint_size=0))
        checkpoint = compile_function(SOURCE, fmt=reference.fmt, config=MathConfig(checkpoint_size=2))
        self.assertTrue(any("math_checkpoint" in dict(m.attributes) for m in checkpoint.program().modules))
        for x in range(8):
            for y in range(8):
                initial = {"x": x, "y": y, "out": 5, "status": 3}
                with self.subTest(x=x, y=y):
                    expected = simulate(reference.program(), initial=initial)
                    actual = simulate(checkpoint.program(), initial=initial)
                    self.assertEqual(actual.amplitudes, expected.amplitudes)

    def test_superposition_negative_control_inverse_and_roundtrip(self):
        def build(size):
            compiled = compile_function(SOURCE, fmt=FixedFormat(3, 1), config=MathConfig(checkpoint_size=size))
            op = compiled.operation
            b = Builder("coherent_checkpoints", {**{r.name: r.type for r in op.module.registers}, "control": Bits(1)})
            b.h(b["x"][0])
            b.h(b["control"])
            b.x(b["y"][0])
            b.x(b["out"][1])
            b.x(b["status"][1])
            with b.control(b["control"], 0):
                b.call(op, **{r.name: b[r.name] for r in op.module.registers})
            return b, op

        reference, _ = build(0)
        checkpoint, operation = build(2)
        program = checkpoint.finish().program()
        self.assertEqual(simulate(program).amplitudes, simulate(reference.finish().program()).amplitudes)
        self.assertEqual(simulate(loads(dumps(program))).amplitudes, simulate(program).amplitudes)
        b = Builder("inverse_checkpoint", {r.name: r.type for r in operation.module.registers})
        b.h(b["x"][0])
        b.x(b["out"][0])
        arguments = {r.name: b[r.name] for r in operation.module.registers}
        b.call(operation, **arguments)
        with b.adjoint():
            b.call(operation, **arguments)
        state = simulate(b.finish().program())
        self.assertEqual(set(state.amplitudes), {(0, 0, 1, 0), (1, 0, 1, 0)})

    def test_index_complex_helpers_and_multiple_outputs(self):
        source = """def helper(z:complex):
 return z*z+1, z*z-1
def f(i,z:complex):
 a,b=helper(z)
 c=a+b
 return (a if i==0 else c), b
"""
        options = {"fmt": FixedFormat(3, 1), "inputs": {"i": Index(1), "z": "complex"}}
        before = compile_function(source, config=MathConfig(checkpoint_size=0), **options)
        after = compile_function(source, config=MathConfig(checkpoint_size=1), **options)
        for i, re, im in ((0, 0, 0), (1, 1, 2), (0, 6, 3)):
            initial = {"i": i, "z_real": re, "z_imag": im, "out_0_real": 3, "status": 1}
            self.assertEqual(simulate(before.program(), initial=initial).amplitudes,
                             simulate(after.program(), initial=initial).amplitudes)

    def test_long_function_reduces_actual_peak_space(self):
        source = "def f(x):\n a=x+1\n b=a*a\n c=b+1\n d=c*c\n e=d+1\n return e*e"
        programs = [compile_function(source, fmt=FixedFormat(6, 2), config=MathConfig(checkpoint_size=size)).program()
                    for size in (0, 2)]
        before, after = [estimate_resources(program) for program in programs]
        self.assertLess(after.qubits + after.mcx_ancilla, before.qubits + before.mcx_ancilla)
        self.assertGreater(after.toffoli, before.toffoli)
        for size in (-1, True, 1025):
            with self.assertRaises(ValidationError):
                MathConfig(checkpoint_size=size)


if __name__ == "__main__":
    unittest.main()
