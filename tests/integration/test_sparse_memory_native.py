"""Sparse interchange is checked on the real backend within its existing materialization limit."""

import tempfile
import unittest
from pathlib import Path

from oracq import (
    QRAM,
    Builder,
    UInt,
    ValidationError,
    dump_qram_yaml,
    load_qram_yaml,
    run_pysparq,
    simulate,
)


class SparseMemoryNativeTests(unittest.TestCase):
    def test_high_addresses_and_interference_match_inline_memory(self):
        b = Builder("high_address_memory", {"a": UInt(12), "v": UInt(3)},
                    {"bank": QRAM(12, 3)})
        for bit in range(1, 12):
            b.x(b["a"][bit])
        b.x(b["v"][0])
        b.x(b["v"][2])
        b.h(b["a"][0])
        b.qram("bank", b["a"], b["v"])
        b.h(b["a"][0])
        program = b.finish().program()
        memory = {"bank": {(1 << 12) - 2: 1, (1 << 12) - 1: 7}}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "memory.qram.yaml"
            path.write_text(dump_qram_yaml(program, memory, sparse=True))
            loaded = load_qram_yaml(path)
        expected = simulate(program, memory)
        report = {}
        actual = run_pysparq(program, loaded, report=report)
        for key in expected.amplitudes.keys() | actual.amplitudes.keys():
            self.assertAlmostEqual(abs(expected.amplitudes.get(key, 0) - actual.amplitudes.get(key, 0)), 0, places=11)
        self.assertEqual(report["native_calls"], 0)

    def test_backend_materialization_limit_is_retained(self):
        b = Builder("large_memory_limit", {"a": UInt(32), "v": UInt(3)},
                    {"bank": QRAM(32, 3)})
        b.qram("bank", b["a"], b["v"])
        with self.assertRaisesRegex(ValidationError, r"more than 2\^20"):
            run_pysparq(b.finish().program(), {"bank": {(1 << 32) - 1: 7}})


if __name__ == "__main__":
    unittest.main()
