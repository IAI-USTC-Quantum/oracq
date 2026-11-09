"""A logical register is a qubit sequence, not one host uint64 word."""

import unittest

from oracq import (
    QRAM,
    Bits,
    Builder,
    UInt,
    ValidationError,
    dumps,
    estimate_resources,
    export_originir,
    loads,
    simulate,
)


class WideRegisterTests(unittest.TestCase):
    def test_coherent_high_bit_phase_survives_serialization(self):
        b = Builder("wide_phase", {"q": Bits(129), "label": Bits(1)})
        b.h(b["q"][128])
        with b.control(b["q"][128]):
            b.x(b["label"])
        b.z(b["q"][128])
        program = loads(dumps(b.finish().program(), format="json"))
        actual = simulate(program).amplitudes
        self.assertEqual(set(actual), {(0, 0), (1 << 128, 1)})
        self.assertAlmostEqual(actual[(0, 0)].real, 2**-.5)
        self.assertAlmostEqual(actual[(1 << 128, 1)].real, -2**-.5)
        self.assertEqual(estimate_resources(program).qubits, 130)

    def test_unsigned_addition_and_view_controls_cross_host_words(self):
        b = Builder("wide_add", {"q": UInt(130), "out": Bits(1)})
        amount = (1 << 96) + 3
        b.add_const(b["q"], amount)
        expected = ((1 << 129) - 1 + amount) % (1 << 130)
        with b.control(b["q"][64:], expected >> 64):
            b.x(b["out"])
        program = loads(dumps(b.finish().program(), format="json"))
        actual = simulate(program, initial={"q": (1 << 129) - 1}).amplitudes
        self.assertEqual(actual, {(expected, 1): 1 + 0j})

    def test_query_address_may_share_controls_with_zero_expectation(self):
        # A query reads its address, so sharing control qubits with it is legal;
        # the interpreter must keep reading the original address per branch.
        b = Builder("shared_address", {"address": Bits(4), "data": Bits(3)}, {"bank": QRAM(4, 3)})
        b.h(b["address"][3])
        with b.control(b["address"][3], 0):
            b.qram("bank", b["address"], b["data"])
        memory = {"bank": {0: 5, 1 << 3: 1}}
        # The |address=0> branch must read bank[0]=5, not the toggled bank[8]=1.
        actual = simulate(b.finish().program(), memory).amplitudes
        self.assertEqual(set(actual), {(0, 5), (1 << 3, 0)})
        self.assertAlmostEqual(actual[(0, 5)].real, 2**-.5)
        self.assertAlmostEqual(actual[(1 << 3, 0)].real, 2**-.5)
        with self.assertRaises(ValidationError):
            # The data word is still written and stays protected.
            b2 = Builder("bad", {"address": Bits(4), "data": Bits(3)}, {"bank": QRAM(4, 3)})
            with b2.control(b2["data"][0], 0):
                b2.qram("bank", b2["address"], b2["data"])
            b2.finish()

    def test_originir_rejects_zero_control_on_query_address(self):
        zero = Builder("zero", {"address": Bits(4), "data": Bits(3)}, {"bank": QRAM(4, 3)})
        with zero.control(zero["address"][3], 0):
            zero.qram("bank", zero["address"], zero["data"])
        with self.assertRaisesRegex(ValidationError, "zero control"):
            export_originir(zero.finish().program())
        one = Builder("one", {"address": Bits(4), "data": Bits(3)}, {"bank": QRAM(4, 3)})
        with one.control(one["address"][3], 1):
            one.qram("bank", one["address"], one["data"])
        self.assertIn("controlled_by", export_originir(one.finish().program()).text)


if __name__ == "__main__":
    unittest.main()
