"""Sparse interchange and wide sparse banks are checked on the real backend."""

import tempfile
import unittest
from pathlib import Path

from oracq import (
    QRAM,
    Builder,
    UInt,
    dump_qram_yaml,
    load_qram_yaml,
    run_pysparq,
    run_pysparq_rir,
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

    def test_wide_sparse_bank_executes_via_logical_word_interface(self):
        # Since the 64-bit logical width cap was lifted, banks with address
        # width over 20 no longer raise a materialization refusal: the dense
        # path stays structurally bounded at 2^20 entries (it is only chosen
        # for address widths <= 20) and wider banks route through the PySparQ
        # logical-word sparse interface instead of expanding the address space.
        b = Builder("large_memory_query", {"a": UInt(32), "v": UInt(3)},
                    {"bank": QRAM(32, 3)})
        b.x(b["a"][31])
        b.x(b["a"][29])
        b.qram("bank", b["a"], b["v"])
        program = b.finish().program()
        memory = {"bank": {(1 << 31) + (1 << 29): 7, 0: 5}}
        expected = simulate(program, memory).amplitudes
        for actual in (run_pysparq(program, memory).amplitudes,
                       run_pysparq_rir(program, memory).amplitudes):
            for key in expected.keys() | actual.keys():
                self.assertAlmostEqual(abs(expected.get(key, 0) - actual.get(key, 0)), 0, places=11)


if __name__ == "__main__":
    unittest.main()
