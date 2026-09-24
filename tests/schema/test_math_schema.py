"""MIR JSON 形状与语义约束；MIR 文本格式保持 JSON。"""

import json
import unittest
from pathlib import Path

from jsonschema import Draft202012Validator

from oracq import ValidationError
from oracq.infrastructure.mathfunc import MathProgram, compile_function


class MathSchemaTests(unittest.TestCase):
    def test_generated_math_graph(self):
        root = Path(__file__).resolve().parents[2]
        validator = Draft202012Validator(
            json.loads((root / "docs/reference/schemas/math-ir.schema.json").read_text())
        )
        result = compile_function(
            "import cmath\ndef square(z):\n return z*z\ndef f(z:complex):\n return cmath.exp(square(z)),z.real>0"
        )
        validator.validate(json.loads(result.math_ir.dumps()))

    def test_semantic_types_and_unknown_fields_rejected(self):
        result = compile_function("def f(x):\n return x*x+1")
        data = json.loads(result.math_ir.dumps())
        data["functions"][0]["nodes"][-1]["kind"] = "bool"
        with self.assertRaises(ValidationError):
            MathProgram.loads(json.dumps(data))
        data = json.loads(result.math_ir.dumps())
        data["functions"][0]["extra"] = "ignored"
        with self.assertRaises(ValidationError):
            MathProgram.loads(json.dumps(data))

    def test_intrinsic_and_future_reference_rejected(self):
        result = compile_function("import math\ndef f(x):\n return math.sin(x)")
        data = json.loads(result.math_ir.dumps())
        data["functions"][0]["nodes"][-1]["data"] = ["unknown"]
        with self.assertRaises(ValidationError):
            MathProgram.loads(json.dumps(data))
        data = json.loads(result.math_ir.dumps())
        data["functions"][0]["nodes"][-1]["args"] = [999]
        with self.assertRaises(ValidationError):
            MathProgram.loads(json.dumps(data))
