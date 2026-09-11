"""witness 断言库的自证测试：正确程序通过、故意错误的构造抛 AssertionError。"""

import math
import unittest

import witness

from pyqecclang import Bits, Builder
from pyqecclang.algorithms.block_encoding import matrix_pauli_encoding
from pyqecclang.algorithms.prepare_select import abstract_prepare, gate_prepare

COEFFICIENTS = (0.6, -0.8, 0.3j, -0.5)


def bell_program():
    """两比特 Bell 电路：H 后 CNOT，是酉程序的最小正例。"""
    b = Builder("bell", {"q": Bits(2)})
    b.h(b["q"][0])
    b.xor(b["q"][0], b["q"][1])
    return b.finish().program()


class _FakeState:
    """simulate 返回值的替身：witness 只依赖 .amplitudes。"""

    def __init__(self, amplitudes):
        self.amplitudes = amplitudes


class UnitaryWitnessTests(unittest.TestCase):
    def test_unitary_program_passes(self):
        witness.assert_unitary(self, bell_program())

    def test_non_normalized_column_fails(self):
        # 原语只能拼出酉程序，非酉反例通过替换 witness.simulate 的返回构造。
        with self.assertRaises(AssertionError):
            original = witness.simulate
            witness.simulate = lambda program, initial=None: _FakeState({(0, 0): 0.5})
            try:
                witness.assert_unitary(self, bell_program())
            finally:
                witness.simulate = original

    def test_non_orthogonal_columns_fail(self):
        # 所有列相同：各自归一但两两内积为 1。
        with self.assertRaises(AssertionError):
            original = witness.simulate
            witness.simulate = lambda program, initial=None: _FakeState(
                {(0, 0): 1 / math.sqrt(2), (3, 0): 1 / math.sqrt(2)}
            )
            try:
                witness.assert_unitary(self, bell_program())
            finally:
                witness.simulate = original


class UncomputationWitnessTests(unittest.TestCase):
    def test_uncomputed_work_register_passes(self):
        b = Builder("clean", {"source": Bits(1), "work": Bits(1)})
        b.xor(b["source"], b["work"])
        b.xor(b["source"], b["work"])
        witness.assert_uncomputation(
            self, b.finish().program(), initial={"source": 1}, work_registers=["work"]
        )

    def test_cleaned_local_passes(self):
        b = Builder("clean_local", {"source": Bits(1)})
        scratch = b.local("scratch", Bits(1))
        b.xor(b["source"], scratch)
        b.xor(b["source"], scratch)
        witness.assert_uncomputation(self, b.finish().program(), initial={"source": 1})

    def test_dirty_work_register_fails(self):
        b = Builder("dirty_work", {"source": Bits(1), "work": Bits(1)})
        b.xor(b["source"], b["work"])
        with self.assertRaises(AssertionError):
            witness.assert_uncomputation(
                self, b.finish().program(), initial={"source": 1}, work_registers=["work"]
            )

    def test_dirty_local_fails(self):
        # simulate 在 LocalExit 处抛 ValidationError，原语须转成 AssertionError。
        b = Builder("dirty_local", {"source": Bits(1)})
        b.x(b.local("scratch", Bits(1)))
        with self.assertRaises(AssertionError):
            witness.assert_uncomputation(self, b.finish().program())


class BindInvariantWitnessTests(unittest.TestCase):
    def test_equivalent_candidates_pass(self):
        slot = abstract_prepare(COEFFICIENTS)
        witness.assert_bind_invariant(
            self,
            slot.operation.program(),
            {
                "gate_a": gate_prepare(COEFFICIENTS).operation,
                "gate_b": gate_prepare(COEFFICIENTS).operation,
            },
        )

    def test_different_candidate_fails(self):
        slot = abstract_prepare(COEFFICIENTS)
        with self.assertRaises(AssertionError):
            witness.assert_bind_invariant(
                self,
                slot.operation.program(),
                {
                    "faithful": gate_prepare(COEFFICIENTS).operation,
                    "wrong": gate_prepare((1.0, 1.0, 1.0, 1.0)).operation,
                },
            )


class BlockEqualsWitnessTests(unittest.TestCase):
    def test_correct_block_passes(self):
        be = matrix_pauli_encoding(((1.0, 0.0), (0.0, -1.0)))
        witness.assert_block_equals(self, be, ((1.0, 0.0), (0.0, -1.0)))

    def test_wrong_matrix_fails(self):
        be = matrix_pauli_encoding(((1.0, 0.0), (0.0, -1.0)))
        with self.assertRaises(AssertionError):
            witness.assert_block_equals(self, be, ((0.0, 1.0), (1.0, 0.0)))


if __name__ == "__main__":
    unittest.main()
