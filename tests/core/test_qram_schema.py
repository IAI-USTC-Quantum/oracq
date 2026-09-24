"""qram YAML 内存定义的解析、写出现场与错误矩阵。"""

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
    """入口含单个 QRAM(2,3) 资源的装载程序。"""
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
            "顶层缺少 qram_segments": "{}",
            "顶层多余键": "qram_segments: []\nextra: 1",
            "顶层不是映射": "- 1",
            "段列表不是列表": "qram_segments: 3",
            "段不是映射": "qram_segments: [1]",
            "段缺少字段": f"qram_segments: [{{name: a, {fields}, type: uint}}]",
            "段多余字段": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [0], extra: 1}}]"
            ),
            "段名为空": f"qram_segments: [{{name: '', {fields}, type: uint, data: [0]}}]",
            "段名重复": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [0]}}, "
                f"{{name: a, {fields}, type: uint, data: [0]}}]"
            ),
            "address_length 为零": (
                "qram_segments: [{name: a, address_length: 0, word_length: 1, "
                "type: uint, data: [0]}]"
            ),
            "word_length 为 65": (
                "qram_segments: [{name: a, address_length: 1, word_length: 65, "
                "type: uint, data: [0]}]"
            ),
            "address_length 为布尔": (
                "qram_segments: [{name: a, address_length: true, word_length: 1, "
                "type: uint, data: [0]}]"
            ),
            "word_length 为字符串": (
                'qram_segments: [{name: a, address_length: 1, word_length: "1", '
                "type: uint, data: [0]}]"
            ),
            "数据类型不受支持": (
                f"qram_segments: [{{name: a, {fields}, type: float, data: [0]}}]"
            ),
            "data 不是数组": f"qram_segments: [{{name: a, {fields}, type: uint, data: 1}}]",
            "uint 字溢出": (
                "qram_segments: [{name: a, address_length: 1, word_length: 3, "
                "type: uint, data: [8]}]"
            ),
            "uint 字为负": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [-1]}}]"
            ),
            "uint 字为布尔": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [true]}}]"
            ),
            "uint 字为字符串": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [x]}}]"
            ),
            "数组超长": (
                f"qram_segments: [{{name: a, {fields}, type: uint, data: [0, 0, 0]}}]"
            ),
            "sint 低于下界": (
                "qram_segments: [{name: a, address_length: 1, word_length: 3, "
                "type: sint, data: [-5]}]"
            ),
            "sint 达到上界": (
                "qram_segments: [{name: a, address_length: 1, word_length: 3, "
                "type: sint, data: [4]}]"
            ),
            "fixedpoint 为负": (
                "qram_segments: [{name: a, address_length: 1, word_length: 3, "
                "type: fixedpoint, data: [-0.1]}]"
            ),
            "fixedpoint 达到 1": (
                "qram_segments: [{name: a, address_length: 1, word_length: 3, "
                "type: fixedpoint, data: [1]}]"
            ),
            "fixedpoint 为布尔": (
                f"qram_segments: [{{name: a, {fields}, type: fixedpoint, data: [true]}}]"
            ),
            "fixedpoint 为字符串": (
                f"qram_segments: [{{name: a, {fields}, type: fixedpoint, data: [x]}}]"
            ),
            "非法 YAML": "qram_segments: [unclosed",
        }
        for label, document in cases.items():
            with self.subTest(label=label), self.assertRaises(ValidationError):
                parse_qram_yaml(document)


if __name__ == "__main__":
    unittest.main()
