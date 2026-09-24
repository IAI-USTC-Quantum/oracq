"""纯数学函数编译契约、路径状态和模块边界。"""

import unittest

from oracq import Bits, Builder, FixedFormat, ValidationError, dumps, loads, simulate
from oracq.infrastructure.mathfunc import (
    Index,
    MathConfig,
    MathProgram,
    compile_function,
    lower_math_ir,
)


class MathFunctionTests(unittest.TestCase):
    def test_nonzero_output_xor_superposition_and_inverse(self):
        compiled = compile_function("def square_plus_one(x):\n return x*x+1", fmt=FixedFormat(4, 1))
        op = compiled.operation
        b = Builder("auto_function_xor", {r.name: r.type for r in op.module.registers})
        b.x(b["x"][0])
        b.h(b["x"][1])
        b.x(b["out"][:2])
        b.x(b["status"][0])
        b.call(op, **{r.name: b[r.name] for r in op.module.registers})
        result = simulate(b.finish().program())
        self.assertEqual(set(result.amplitudes), {(1, 1, 1), (3, 5, 1)})
        inverse = Builder("auto_function_inverse", {r.name: r.type for r in op.module.registers})
        inverse.x(inverse["x"][0])
        inverse.x(inverse["out"][0])
        for _ in range(2):
            inverse.call(op, **{r.name: inverse[r.name] for r in op.module.registers})
        self.assertEqual(simulate(inverse.finish().program()).amplitudes, {(1, 1, 0): 1 + 0j})

    def test_unused_branch_error_is_masked(self):
        op = compile_function(
            "def safe_inverse(x):\n return 0 if x == 0 else 1 / x", fmt=FixedFormat(4, 1)
        ).operation
        result = simulate(op.program())
        self.assertEqual(result.amplitudes, {(0, 0, 0): 1 + 0j})
        bad = compile_function("def inverse(x):\n return 1 / x", fmt=FixedFormat(4, 1)).operation
        self.assertTrue(next(iter(simulate(bad.program()).amplitudes))[-1] & 1)

    def test_source_is_not_executed_and_effects_rejected(self):
        for source in (
            "def f(x):\n print(x)\n return x",
            "def f(x):\n x.real = 2\n return x",
            "def f(x):\n while x > 0:\n  x=x-1\n return x",
            "def f(x):\n return f(x)",
            "def f(x):\n return open('must_not_create','w')",
        ):
            with self.subTest(source=source):
                with self.assertRaises(ValidationError):
                    compile_function(source)
        self.assertFalse(__import__("pathlib").Path("must_not_create").exists())

    def test_future_imports_and_module_docstring_are_ignored(self):
        result = compile_function(
            '"""模块说明。"""\n'
            "from __future__ import annotations\n"
            "import math\n"
            "def norm(x: float) -> float:\n"
            " return math.sqrt(x)",
            fmt=FixedFormat(8, 3),
        )
        parameters = result.math_ir.function_map[result.math_ir.entry].parameters
        self.assertEqual([(p.name, p.kind) for p in parameters], [("x", "real")])

    def test_helper_module_reuse_and_roundtrip(self):
        result = compile_function(
            "def helper(x):\n return x*x\ndef f(x):\n y=helper(x)\n return y+helper(x)"
        )
        self.assertEqual(len(result.math_ir.functions), 2)
        graph = result.math_ir.function_map[result.math_ir.entry]
        self.assertEqual(sum(n.op == "call" for n in graph.nodes), 1)
        self.assertEqual(MathProgram.loads(result.math_ir.dumps()), result.math_ir)
        rebuilt = lower_math_ir(MathProgram.loads(result.math_ir.dumps()), fmt=result.fmt)
        self.assertEqual(dumps(rebuilt.program()), dumps(result.program()))
        self.assertEqual(loads(dumps(result.program())), result.program())

    def test_static_loop_constants_and_index(self):
        result = compile_function(
            "def polynomial(x, n=3):\n y=0\n for i in range(n):\n  y=y+x**i\n return y"
        )
        self.assertEqual(
            [p.name for p in result.math_ir.function_map[result.math_ir.entry].parameters], ["x"]
        )
        index = compile_function(
            "def select(i,x):\n return x if i==2 else 0",
            inputs={"i": Index(2), "x": "real"},
            fmt=FixedFormat(6, 2),
        )
        self.assertEqual(index.operation.module.registers[0].type.width, 2)

    def test_all_elementary_families_construct(self):
        real = (
            "sqrt",
            "exp",
            "log",
            "log10",
            "sin",
            "cos",
            "tan",
            "asin",
            "acos",
            "atan",
            "sinh",
            "cosh",
            "tanh",
            "asinh",
            "acosh",
            "atanh",
        )
        for module in ("math", "cmath"):
            for name in real:
                with self.subTest(module=module, name=name):
                    source = f"import {module}\ndef f(x):\n return {module}.{name}(x)"
                    result = compile_function(
                        source,
                        inputs={"x": "complex" if module == "cmath" else "real"},
                        fmt=FixedFormat(8, 3),
                        config=MathConfig(degree=2),
                    )
                    self.assertFalse(any(m.body is None for m in result.program().modules))
        result = compile_function(
            "import cmath\ndef f(z:complex):\n r,p=cmath.polar(z)\n return cmath.rect(r,p),z.conjugate()"
        )
        self.assertEqual(len(result.output_layout), 2)

    def test_different_configs_can_coexist(self):
        src = "import math\ndef f(x):\n return math.sin(x)"
        a = compile_function(src, config=MathConfig(degree=2))
        b = compile_function(src, config=MathConfig(degree=3))
        self.assertNotEqual(a.operation.module.name, b.operation.module.name)
        program = Builder("both", {"x": Bits(12), "a": Bits(12), "b": Bits(12), "flags": Bits(4)})
        program.call(a.operation, x=program["x"], out=program["a"], status=program["flags"][:2])
        program.call(b.operation, x=program["x"], out=program["b"], status=program["flags"][2:])
        program.finish()

    def test_roe_is_a_compiled_pure_function(self):
        from oracq.applications.roe import roe_face

        op = roe_face(fmt=FixedFormat(4, 1))
        self.assertEqual(dict(op.module.attributes)["math_function"], "frozen_roe_face")
        labels = {dict(m.attributes).get("math_function") for m in op.program().modules}
        self.assertTrue(
            {"pick3", "conserved_to_primitive", "entropy_absolute", "euler_entry"} <= labels
        )
