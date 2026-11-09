"""Sparse memory interchange must preserve addresses without visiting the address space."""

import json
import tempfile
import unittest
from pathlib import Path

from oracq import QRAM, Bits, Builder, ValidationError, dump_qram_yaml, load_qram_yaml, simulate
from oracq.infrastructure.qram_schema import parse_qram_yaml


class SparseMemoryTests(unittest.TestCase):
    def test_64_bit_addresses_roundtrip_without_densification(self):
        b = Builder("sparse_memory", {"a": Bits(64), "v": Bits(8)}, {"bank": QRAM(64, 8)})
        b.qram("bank", b["a"], b["v"])
        program = b.finish().program()
        expected = {"bank": {0: 7, (1 << 64) - 1: 255}}
        text = dump_qram_yaml(program, expected, sparse=True)
        self.assertLess(len(text), 500)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "memory.qram.yaml"
            path.write_text(text, encoding="utf-8")
            self.assertEqual(load_qram_yaml(path), expected)
        segment = parse_qram_yaml(text)[0]
        self.assertEqual(segment.addresses, (0, (1 << 64) - 1))

    def test_dense_legacy_and_sparse_loaded_execution_are_equivalent(self):
        b = Builder("interference", {"a": Bits(3), "v": Bits(3)}, {"bank": QRAM(3, 3)})
        b.h(b["a"][0])
        b.qram("bank", b["a"], b["v"])
        b.h(b["a"][0])
        program = b.finish().program()
        memory = {"bank": {0: 1, 1: 7, 4: 0}}
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / f"{sparse}.qram.yaml" for sparse in (False, True)]
            for sparse, path in zip((False, True), paths, strict=True):
                path.write_text(dump_qram_yaml(program, memory, sparse=sparse))
            dense, sparse = map(load_qram_yaml, paths)
        self.assertIsInstance(dense["bank"], list)
        self.assertIsInstance(sparse["bank"], dict)
        self.assertEqual(simulate(program, dense), simulate(program, sparse))

    def test_sparse_signed_and_fixedpoint_values_keep_the_existing_encoding(self):
        text = """qram_segments:
- {name: s, address_length: 32, word_length: 3, type: sint, encoding: sparse, data: [[4294967295, -1]]}
- {name: f, address_length: 32, word_length: 3, type: fixedpoint, encoding: sparse, data: [[4, 0.5]]}
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "memory.qram.yaml"
            path.write_text(text)
            self.assertEqual(load_qram_yaml(path), {"s": {4294967295: 7}, "f": {4: 4}})

    def test_invalid_pairs_ranges_duplicates_and_encoding_are_rejected(self):
        base = {"name": "s", "address_length": 3, "word_length": 3, "type": "uint", "encoding": "sparse"}
        for pairs in ([[1]], [[1, 2, 3]], [[True, 1]], [[-1, 1]], [[8, 1]], [[1, 8]],
                      [[1, -1]], [[1, True]], [[1, 1], [1, 2]]):
            with self.subTest(pairs=pairs), self.assertRaises(ValidationError):
                parse_qram_yaml(json.dumps({"qram_segments": [{**base, "data": pairs}]}))
        with self.assertRaises(ValidationError):
            parse_qram_yaml(json.dumps({"qram_segments": [{**base, "encoding": "other", "data": []}]}))


if __name__ == "__main__":
    unittest.main()
