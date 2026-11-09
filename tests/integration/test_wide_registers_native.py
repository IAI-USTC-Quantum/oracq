"""Wide logical registers must execute on actual native sparse states."""

import unittest

from oracq import QRAM, Bits, Builder, UInt, run_pysparq, run_pysparq_rir, simulate


class WideNativeTests(unittest.TestCase):
    def compare(self, program, memory=None):
        expected = simulate(program, memory)
        for actual in (run_pysparq(program, memory), run_pysparq_rir(program, memory)):
            for key in expected.amplitudes.keys() | actual.amplitudes.keys():
                self.assertAlmostEqual(abs(expected.amplitudes.get(key, 0)-actual.amplitudes.get(key, 0)), 0)

    def test_high_word_entanglement_phase_and_global_phase(self):
        b = Builder("wide_coherence", {"q": Bits(129), "label": Bits(1)})
        b.h(b["q"][128])
        with b.control(b["q"][128]):
            b.x(b["label"])
        b.z(b["q"][128])
        b.global_phase(.37)
        self.compare(b.finish().program())

    def test_wide_qram_with_shared_zero_address_control_and_unquery(self):
        b = Builder("wide_query", {"address": Bits(65), "data": Bits(129)}, {"bank": QRAM(65, 129)})
        b.h(b["address"][64])
        with b.control(b["address"][64], 0):
            b.qram("bank", b["address"], b["data"])
        memory = {"bank": {0: (1 << 128)+3, 1 << 64: (1 << 65)+1}}
        self.compare(b.finish().program(), memory)

    def test_dense_qram_zero_control_shared_with_address(self):
        b = Builder("dense_query", {"address": Bits(4), "data": Bits(3)}, {"bank": QRAM(4, 3)})
        b.h(b["address"][3])
        with b.control(b["address"][3], 0):
            b.qram("bank", b["address"], b["data"])
        program = b.finish().program()
        memory = {"bank": {0: 5, 3: 7, 1 << 3: 1}}
        expected = simulate(program, memory).amplitudes
        # The gate-level adapter resolves the shared zero control against the
        # logical address word, so it matches the reference here.
        actual = run_pysparq(program, memory).amplitudes
        for key in expected.keys() | actual.keys():
            self.assertAlmostEqual(abs(expected.get(key, 0)-actual.get(key, 0)), 0)

    @unittest.expectedFailure
    def test_dense_rir_zero_control_shared_address_defect(self):
        # Reproduced against pysparq 0.3.0 and 0.4.0: the native run_rir
        # executor toggles zero-control bits around a query with X gates, so a
        # dense QRAM whose address register shares a bit with a zero control
        # reads the redirected cell (address 8 instead of 0). Tracked with the
        # SparQSim owner; when the upstream fix lands this expected failure
        # becomes an unexpected success and the marker must be removed.
        b = Builder("dense_query_defect", {"address": Bits(4), "data": Bits(3)}, {"bank": QRAM(4, 3)})
        b.h(b["address"][3])
        with b.control(b["address"][3], 0):
            b.qram("bank", b["address"], b["data"])
        program = b.finish().program()
        memory = {"bank": {0: 5, 3: 7, 1 << 3: 1}}
        expected = simulate(program, memory).amplitudes
        actual = run_pysparq_rir(program, memory).amplitudes
        for key in expected.keys() | actual.keys():
            self.assertAlmostEqual(abs(expected.get(key, 0)-actual.get(key, 0)), 0)

    def test_large_addition_and_clean_wide_workspace(self):
        b = Builder("wide_carries", {"q": UInt(130)})
        scratch = b.local("scratch", Bits(129))
        b.h(b["q"][63])
        with b.control(b["q"][63]):
            b.x(scratch[128])
        with b.control(b["q"][63]):
            b.x(scratch[128])
        b.add_const(b["q"], (1 << 96)+3)
        self.compare(b.finish().program())


if __name__ == "__main__":
    unittest.main()
