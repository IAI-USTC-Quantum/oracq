"""Check specification serialization instances with a real JSON Schema validator; the schema itself stays JSON."""

import json
import unittest
from pathlib import Path

import yaml
from jsonschema import Draft202012Validator

from oracq import Bits, Builder, dumps, identity, linear_combination, pauli_x

ROOT = Path(__file__).resolve().parents[2]


class SchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema = json.loads((ROOT / "docs/reference/schemas/rir.schema.json").read_text())
        Draft202012Validator.check_schema(cls.schema)
        cls.validator = Draft202012Validator(cls.schema)

    def test_open_application_modules_are_schema_valid(self):
        from oracq.applications.catalog import build_case

        for name in ("costa_qram", "qfvm_qram", "qham_qode"):
            with self.subTest(case=name):
                self.validator.validate(yaml.safe_load(dumps(build_case(name).program)))

    def test_qram_fixture(self):
        from oracq import loads

        texts = {
            name: (ROOT / "examples" / name).read_text() for name in ("qram.rir.json", "qram.rir.yaml")
        }
        for name, text in texts.items():
            with self.subTest(fixture=name):
                self.validator.validate(yaml.safe_load(text))
        self.assertEqual(loads(texts["qram.rir.yaml"]), loads(texts["qram.rir.json"]))

    def test_all_structured_nodes_and_metadata(self):
        be = linear_combination(1j, identity(1), -2, pauli_x(1))
        b = Builder(
            "schema_example",
            {"target": Bits(1), "signal": Bits(be.signal_qubits), "control": Bits(1)},
        )
        with b.repeat(13):
            with b.adjoint():
                with b.control(b["control"], 0):
                    b.call(be.operation, target=b["target"], signal=b["signal"])
        self.validator.validate(yaml.safe_load(dumps(b.finish().program())))

    def test_store_node_is_schema_valid(self):
        from oracq import QRAM, Builder, QMem, UInt

        b = Builder("schema_store", {"a": UInt(2), "v": UInt(4)}, {"ram": QRAM(2, 4)})
        QMem(b, "ram")[b["a"]].store(b["v"])
        data = yaml.safe_load(dumps(b.finish().program()))
        self.validator.validate(data)
        node = data["modules"][0]["body"][-1]
        self.assertEqual(node["tag"], "Store")

    def test_rejects_invalid_version_and_unknown_fields(self):
        data = yaml.safe_load(dumps(identity(1).operation.program()))
        data["version"] = "1.0"
        self.assertTrue(list(self.validator.iter_errors(data)))
        data["version"] = "0.1"
        data["modules"][0]["extra"] = "wrong"
        self.assertTrue(list(self.validator.iter_errors(data)))

    def test_repeat_limit_is_exact(self):
        self.assertEqual(
            self.schema["$defs"]["Repeat"]["properties"]["count"]["maximum"], (1 << 63) - 1
        )

    def test_binding_capture_metadata_survives_schema_roundtrip(self):
        from oracq import bind, loads
        from oracq.applications.oracle_study import oracle_study

        opened, variants = oracle_study(2)
        closed = bind(opened, {"AngleWord": variants["qram_table"][0]})
        self.validator.validate(yaml.safe_load(dumps(closed)))
        self.assertEqual(loads(dumps(closed)), closed)
        self.assertIn("binding_captures", dict(closed.main.attributes))


class QramMemorySchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema = json.loads(
            (ROOT / "docs/reference/schemas/qram-memory.schema.json").read_text()
        )
        Draft202012Validator.check_schema(cls.schema)
        cls.validator = Draft202012Validator(cls.schema)

    def test_memory_fixture_is_schema_valid(self):
        data = yaml.safe_load((ROOT / "examples/memory.qram.yaml").read_text(encoding="utf-8"))
        self.validator.validate(data)

    def test_rejects_unknown_type_and_extra_field(self):
        data = yaml.safe_load((ROOT / "examples/memory.qram.yaml").read_text(encoding="utf-8"))
        data["qram_segments"][0]["type"] = "float"
        self.assertTrue(list(self.validator.iter_errors(data)))
        data["qram_segments"][0]["type"] = "uint"
        data["qram_segments"][0]["extra"] = 1
        self.assertTrue(list(self.validator.iter_errors(data)))

    def test_type_specific_data_constraints(self):
        def check(type_, words, invalid):
            data = {
                "qram_segments": [
                    {"name": "s", "address_length": 1, "word_length": 3, "type": type_, "data": words}
                ]
            }
            if invalid:
                self.assertTrue(list(self.validator.iter_errors(data)), (type_, words))
            else:
                self.validator.validate(data)

        check("uint", [0, 7], invalid=False)
        check("uint", [-1], invalid=True)
        check("uint", [1.5], invalid=True)
        check("sint", [-4, 3], invalid=False)
        check("sint", [0.5], invalid=True)
        check("fixedpoint", [0, 0.5, 0.999], invalid=False)
        check("fixedpoint", [1], invalid=True)
        check("fixedpoint", [-0.1], invalid=True)
