"""Standalone JSON representations and constraints of PDE/QCL."""

import json
import unittest
from pathlib import Path

from jsonschema import Draft202012Validator

from oracq import ValidationError
from oracq.applications.qham import Field, Known, PolynomialPDE, QHAMPlan


class QhamSchemaTests(unittest.TestCase):
    def test_schemas_and_roundtrip(self):
        root = Path(__file__).resolve().parents[2]
        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": u.d("x", 2) - u * u.d("x") + Known("f")})
        plan = QHAMPlan(pde, 3)
        for name, value in (("pde", pde), ("qcl-plan", plan)):
            schema = json.loads((root / f"docs/reference/schemas/{name}.schema.json").read_text())
            Draft202012Validator.check_schema(schema)
            Draft202012Validator(schema).validate(json.loads(value.dumps()))

    def test_unknown_atoms_and_fields_rejected(self):
        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": u * u})
        data = json.loads(pde.dumps())
        data["terms"][0]["monomial"]["fields"][0]["name"] = "v"
        with self.assertRaises(ValidationError):
            PolynomialPDE.loads(json.dumps(data))
        data = json.loads(pde.dumps())
        data["terms"][0]["monomial"]["extra"] = 1
        with self.assertRaises(ValidationError):
            PolynomialPDE.loads(json.dumps(data))

    def test_negative_order_rejected(self):
        pde = PolynomialPDE.from_equations({"u": Field("u")})
        with self.assertRaises(ValidationError):
            QHAMPlan(pde, -1)
