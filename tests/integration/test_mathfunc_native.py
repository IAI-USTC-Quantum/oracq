"""自动纯函数的真实 PySparQ 执行与 OriginIR 消费。"""

import unittest

from pyqecclang import (
    Builder,
    FixedFormat,
    arithmetic_native_registry,
    export_toffoli_u3_cz,
    run_pysparq,
    simulate,
)
from pyqecclang.infrastructure.mathfunc import MathConfig, compile_function


class MathFunctionNativeTests(unittest.TestCase):
    def run_compiled(self, source, fmt, initial, **options):
        result = compile_function(source, fmt=fmt, **options)
        op = result.operation
        b = Builder("math_native_input", {r.name: r.type for r in op.module.registers})
        for name, value in initial.items():
            for bit in range(b[name].width):
                if (value >> bit) & 1:
                    b.x(b[name][bit])
        b.call(op, **{r.name: b[r.name] for r in op.module.registers})
        p = b.finish().program()
        report = {}
        actual = run_pysparq(p, native_registry=arithmetic_native_registry(p), report=report)
        return actual, report

    def test_complex_nonzero_output(self):
        state, report = self.run_compiled(
            "def f(z:complex):\n return z*z+1",
            FixedFormat(6, 2),
            {"z_real": 2, "z_imag": 4, "out_real": 7},
        )
        self.assertEqual(set(state.amplitudes), {(2, 4, 6, 4, 0)})
        self.assertGreater(report["native_calls"], 0)

    def test_guarded_singular_branch_on_superposition(self):
        result = compile_function(
            "def safe(x):\n if x==0:\n  return 0\n return 1/x", fmt=FixedFormat(6, 2)
        )
        b = Builder(
            "guarded_superposition", {r.name: r.type for r in result.operation.module.registers}
        )
        b.h(b["x"][0])
        b.call(result.operation, **{r.name: b[r.name] for r in result.operation.module.registers})
        p = b.finish().program()
        actual = run_pysparq(p, native_registry=arithmetic_native_registry(p))
        self.assertEqual(set(actual.amplitudes), {(0, 0, 0), (1, 16, 0)})

    def test_elementary_polynomial_native_smoke(self):
        fmt = FixedFormat(8, 3)
        state, report = self.run_compiled(
            "import math\ndef f(x):\n return math.exp(x)",
            fmt,
            {},
            config=MathConfig(degree=3, intervals=(("exp", -1.0, 1.0),)),
        )
        basis = next(iter(state.amplitudes))
        self.assertAlmostEqual(fmt.decode(basis[1]), 1.0, delta=0.3)
        self.assertEqual(basis[2], 0)
        self.assertGreater(report["native_calls"], 10)
        outside, _ = self.run_compiled(
            "import math\ndef f(x):\n return math.sin(x)",
            fmt,
            {"x": fmt.encode(4)},
            config=MathConfig(degree=2),
        )
        self.assertTrue(next(iter(outside.amplitudes))[-1] & 2)

    def test_actual_originir_function_description(self):
        from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser
        from uniqc.simulator import Simulator

        formula = compile_function("def copy_value(x):\n return x", fmt=FixedFormat(2, 0))
        b = Builder(
            "function_originir", {r.name: r.type for r in formula.operation.module.registers}
        )
        b.h(b["x"])
        b.x(b["out"][1])
        b.call(formula.operation, **{r.name: b[r.name] for r in formula.operation.module.registers})
        p = b.finish().program()
        text = export_toffoli_u3_cz(p).text
        vector = Simulator(least_qubit_remapping=False).simulate_statevector(text)
        expected = simulate(p).statevector()
        for value, ref in zip(vector[: len(expected)], expected, strict=True):
            self.assertAlmostEqual(complex(value), ref, places=11)
        self.assertLess(sum(abs(v) ** 2 for v in vector[len(expected) :]), 1e-20)
        polynomial = compile_function(
            "import math\ndef f(x):\n return math.sin(x)+math.exp(x)",
            fmt=FixedFormat(8, 3),
            config=MathConfig(degree=2),
        )
        artifact = export_toffoli_u3_cz(polynomial.program())
        parser = OriginIR_BaseParser()
        parser.parse(artifact.text)
        self.assertGreater(parser.n_qubit, 8)
