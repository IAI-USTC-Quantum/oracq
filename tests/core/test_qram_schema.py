"""Parsing, dump scenarios, and the error matrix for qram YAML memory definitions."""

import tempfile
import unittest
from pathlib import Path

from oracq import (
    QRAM,
    Builder,
    UInt,
    ValidationError,
    dump_qram_yaml,
    identity,
    load_qram_yaml,
    simulate,
)
from oracq.infrastructure.qram_schema import parse_qram_yaml

ROOT = Path(__file__).resolve().parents[2]

DOCUMENT = """\
qram_segments:
  - name: values
    address_length: 2
    word_length: 3
    type: uint
    data: [1, 2, 4, 7]
"""


def load_program():
    """Entry program containing a single QRAM(2,3) resource."""
    b = Builder("main", {"a": UInt(2), "d": UInt(3)}, {"values": QRAM(2, 3)})
    b.qram("values", b["a"], b["d"])
    return b.finish().program()


class QramYamlTests(unittest.TestCase):
    def test_example_file_loads(self):
        self.assertEqual(
            load_qram_yaml(ROOT / "examples" / "memory.qram.yaml"), {"values": [1, 2, 4, 7]}
        )

    def test_parse_returns_segments_with_source_data(self):
        document = (
            "qram_segments:\n"
            "  - {name: t, address_length: 1, word_length: 3, type: sint, data: [-1, 3]}\n"
            "  - {name: f, address_length: 1, word_length: 2, type: fixedpoint, data: [0, 0.5]}\n"
        )
        first, second = parse_qram_yaml(document)
        self.assertEqual((first.name, first.address_length, first.word_length, first.type), ("t", 1, 3, "sint"))
        self.assertEqual(first.data, (-1, 3))
        self.assertEqual((second.name, second.type), ("f", "fixedpoint"))
        self.assertEqual(second.data, (0, 0.5))

    def test_empty_segments_allowed(self):
        self.assertEqual(parse_qram_yaml("qram_segments: []\n"), ())

    def test_loaded_memory_compatible_with_simulate(self):
        program = load_program()
        loaded = load_qram_yaml(ROOT / "examples" / "memory.qram.yaml")
        inline = {"values": [1, 2, 4, 7]}
        self.assertEqual(
            simulate(program, loaded, initial={"a": 3, "d": 5}),
            simulate(program, inline, initial={"a": 3, "d": 5}),
        )

    def _load_document(self, document: str) -> dict[str, list[int]]:
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "memory.qram.yaml"
            path.write_text(document, encoding="utf-8")
            return load_qram_yaml(path)

    def test_sint_encodes_twos_complement(self):
        document = (
            "qram_segments:\n"
            "  - {name: s, address_length: 2, word_length: 3, type: sint, data: [-4, -1, 0, 3]}\n"
        )
        self.assertEqual(self._load_document(document), {"s": [4, 7, 0, 3]})

    def test_fixedpoint_encodes_by_truncation_and_pads(self):
        document = (
            "qram_segments:\n"
            "  - {name: f, address_length: 2, word_length: 3, type: fixedpoint, data: [0, 0.5, 0.999]}\n"
        )
        self.assertEqual(self._load_document(document), {"f": [0, 4, 7, 0]})

    def test_short_uint_data_is_padded(self):
        document = (
            "qram_segments:\n  - {name: p, address_length: 2, word_length: 3, type: uint, data: [5]}\n"
        )
        self.assertEqual(self._load_document(document), {"p": [5, 0, 0, 0]})

    def test_dump_round_trip_with_sparse_input(self):
        text = dump_qram_yaml(load_program(), {"values": {1: 2, 3: 7}})
        memory = {s.name: list(s.data) for s in parse_qram_yaml(text)}
        self.assertEqual(memory, {"values": [0, 2, 0, 7]})
        self.assertEqual({s.type for s in parse_qram_yaml(text)}, {"uint"})

    def test_dump_without_resources(self):
        program = identity(1).operation.program()
        self.assertEqual(parse_qram_yaml(dump_qram_yaml(program, None)), ())

    def test_dump_rejects_resource_mismatch(self):
        with self.assertRaises(ValidationError):
            dump_qram_yaml(load_program(), {})

    def test_parse_rejects_bad_documents(self):
        fields = "address_length: 1, word_length: 1"
        cases = {
            "top level missing qram_segments": "{}",
            "extra top-level key": "qram_segments: []\nextra: 1",
            "top level not a mapping": "- 1",
            "segment list not a list": "qram_segments: 3",
            "segment not a mapping": "qram_segments: [1]",
            "segment missing field": f"qram_segments: [{{name: a, {fields}, type: uint}}]",
            "segment extra field": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [0], extra: 1}}]"
            ),
            "empty segment name": f"qram_segments: [{{name: '', {fields}, type: uint, data: [0]}}]",
            "duplicate segment name": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [0]}}, "
                f"{{name: a, {fields}, type: uint, data: [0]}}]"
            ),
            "address_length is zero": (
                "qram_segments: [{name: a, address_length: 0, word_length: 1, "
                "type: uint, data: [0]}]"
            ),
            "word_length is 65": (
                "qram_segments: [{name: a, address_length: 1, word_length: 65, "
                "type: uint, data: [0]}]"
            ),
            "address_length is boolean": (
                "qram_segments: [{name: a, address_length: true, word_length: 1, "
                "type: uint, data: [0]}]"
            ),
            "word_length is string": (
                'qram_segments: [{name: a, address_length: 1, word_length: "1", '
                "type: uint, data: [0]}]"
            ),
            "unsupported data type": (
                f"qram_segments: [{{name: a, {fields}, type: float, data: [0]}}]"
            ),
            "data not a list": f"qram_segments: [{{name: a, {fields}, type: uint, data: 1}}]",
            "uint word overflow": (
                "qram_segments: [{name: a, address_length: 1, word_length: 3, "
                "type: uint, data: [8]}]"
            ),
            "uint word negative": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [-1]}}]"
            ),
            "uint word is boolean": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [true]}}]"
            ),
            "uint word is string": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [x]}}]"
            ),
            "data list too long": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [0, 0, 0]}}]"
            ),
            "sint below lower bound": (
                "qram_segments: [{name: a, address_length: 1, word_length: 3, "
                "type: sint, data: [-5]}]"
            ),
            "sint at upper bound": (
                "qram_segments: [{name: a, address_length: 1, word_length: 3, "
                "type: sint, data: [4]}]"
            ),
            "fixedpoint negative": (
                "qram_segments: [{name: a, address_length: 1, word_length: 3, "
                "type: fixedpoint, data: [-0.1]}]"
            ),
            "fixedpoint equals 1": (
                "qram_segments: [{name: a, address_length: 1, word_length: 3, "
                "type: fixedpoint, data: [1]}]"
            ),
            "fixedpoint is boolean": (
                f"qram_segments: [{{name: a, {fields}, type: fixedpoint, data: [true]}}]"
            ),
            "fixedpoint is string": (
                f"qram_segments: [{{name: a, {fields}, type: fixedpoint, data: [x]}}]"
            ),
            "invalid YAML": "qram_segments: [unclosed",
        }
        for label, document in cases.items():
            with self.subTest(label=label), self.assertRaises(ValidationError):
                parse_qram_yaml(document)


if __name__ == "__main__":
    unittest.main()
