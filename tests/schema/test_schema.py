"""用真实 JSON Schema 验证器检查规范和序列化实例。"""

import json
import unittest
from pathlib import Path

from jsonschema import Draft202012Validator

from pyqecclang import Bits, Builder, dumps, identity, linear_combination, pauli_x

ROOT = Path(__file__).resolve().parents[2]


class SchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema = json.loads((ROOT / "docs/reference/schemas/rir.schema.json").read_text())
        Draft202012Validator.check_schema(cls.schema)
        cls.validator = Draft202012Validator(cls.schema)

    def test_open_application_modules_are_schema_valid(self):
        from pyqecclang.applications.catalog import build_case

        for name in ("costa_qram", "qfvm_qram", "qham_qode"):
            with self.subTest(case=name):
                self.validator.validate(json.loads(dumps(build_case(name).program)))

    def test_qram_fixture(self):
        data = json.loads((ROOT / "examples/qram.rir.json").read_text())
        self.validator.validate(data)

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
        self.validator.validate(json.loads(dumps(b.finish().program())))

    def test_store_node_is_schema_valid(self):
        from pyqecclang import QRAM, Builder, QMem, UInt

        b = Builder("schema_store", {"a": UInt(2), "v": UInt(4)}, {"ram": QRAM(2, 4)})
        QMem(b, "ram")[b["a"]].store(b["v"])
        data = json.loads(dumps(b.finish().program()))
        self.validator.validate(data)
        node = data["modules"][0]["body"][-1]
        self.assertEqual(node["tag"], "Store")

    def test_rejects_invalid_version_and_unknown_fields(self):
        data = json.loads(dumps(identity(1).operation.program()))
        data["version"] = "1.0"
        self.assertTrue(list(self.validator.iter_errors(data)))
        data["version"] = "0.1"
        data["modules"][0]["extra"] = "wrong"
        self.assertTrue(list(self.validator.iter_errors(data)))

    def test_repeat_limit_is_exact(self):
        self.assertEqual(
            self.schema["$defs"]["Repeat"]["properties"]["count"]["maximum"], (1 << 63) - 1
        )
