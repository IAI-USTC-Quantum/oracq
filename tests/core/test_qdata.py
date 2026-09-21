"""QVector/QMatrix 量子数据结构的语义测试：角树真值、局部更新与制备振幅。"""

import math
import unittest

from pyqecclang import QRAM, Builder, QMem, UInt, ValidationError, simulate
from pyqecclang.algorithms.common.arithmetic import FixedFormat
from pyqecclang.algorithms.input_model.oracles import qram_state_angles
from pyqecclang.algorithms.input_model.qdata import QMatrix, QVector

FMT = FixedFormat(8, 4)
AW = 10
TOL = 4e-3  # angle_width=10 的旋转角量化误差量级


class QVectorTests(unittest.TestCase):
    def test_angle_bank_matches_qram_state_angles(self):
        values = (1.0, 2.0, 0.0, 3.0)
        vector = QVector(values, fmt=FMT, angle_width=AW)
        self.assertEqual(vector.angles, qram_state_angles(values, AW))

    def test_preparation_amplitudes_match_classical(self):
        for values, signed in (
            ((1.0, 2.0, 0.0, 3.0), False),
            ((1.0, -2.0, 0.5, 0.0), True),
            ((0.5, 0.25, 0.125, 0.0625, 0.5, 0.5, 0.0, 1.0), False),
        ):
            with self.subTest(values=values, signed=signed):
                vector = QVector(values, fmt=FMT, angle_width=AW)
                state = simulate(
                    vector.preparation(signed=signed).operation.program(),
                    vector.snapshot(signed=signed),
                )
                expected = vector.amplitudes()
                for index, value in enumerate(expected):
                    keys = [key for key in state.amplitudes if key[0] == index]
                    if value == 0:
                        self.assertFalse(keys)
                    else:
                        self.assertEqual(len(keys), 1)
                        self.assertAlmostEqual(state.amplitudes[keys[0]], value, delta=TOL)

    def test_update_matches_full_rebuild(self):
        vector = QVector((1.0, 2.0, 0.0, 3.0), fmt=FMT, angle_width=AW)
        changed = vector.update(1, 0.5)
        rebuilt = QVector((1.0, 0.5, 0.0, 3.0), fmt=FMT, angle_width=AW)
        self.assertEqual(vector.angles, rebuilt.angles)
        self.assertEqual(vector.tree, rebuilt.tree)
        for address, word in changed.items():
            self.assertEqual(rebuilt.angles[address], word)

    def test_norm_and_validation(self):
        vector = QVector((3.0, 4.0), fmt=FMT, angle_width=AW)
        self.assertAlmostEqual(vector.norm, 5.0)
        with self.assertRaises(ValidationError):
            QVector((1.0, 2.0, 3.0), fmt=FMT, angle_width=AW)
        with self.assertRaises(ValidationError):
            QVector((1.0, 2.0), fmt=FMT, angle_width=0)


class QMatrixTests(unittest.TestCase):
    def setUp(self):
        self.matrix = QMatrix(
            [
                [0.5, 0.5, 0.25, 0.0],
                [0.5, 0.0, 0.25, 0.5],
                [0.25, 0.5, 0.0, 0.5],
                [0.0, 0.0, 0.5, 0.5],
            ],
            fmt=FMT,
            angle_width=AW,
        )

    def test_query_returns_entry_words(self):
        program = self.matrix.query().operation.program()
        for row in range(4):
            for col in range(4):
                with self.subTest(row=row, col=col):
                    state = simulate(
                        program,
                        {"entries": self.matrix.snapshot()["entries"]},
                        initial={"address": (row << 2) | col},
                    )
                    word = self.matrix.words[row][col]
                    self.assertEqual(state.amplitudes, {((row << 2) | col, word): 1 + 0j})
    def test_row_preparation_matches_classical_rows(self):
        program = self.matrix.row_preparation().program()
        for row in range(4):
            with self.subTest(row=row):
                state = simulate(
                    program,
                    {"row_angles": self.matrix.snapshot()["row_angles"]},
                    initial={"row": row},
                )
                expected = self.matrix.row_amplitudes(row)
                for key, amplitude in state.amplitudes.items():
                    self.assertAlmostEqual(amplitude, expected[key[1]], delta=TOL)

    def test_amplitude_preparation_matches_row_norms(self):
        state = simulate(
            self.matrix.amplitude_preparation().program(),
            {"root_angles": self.matrix.snapshot()["root_angles"]},
        )
        expected = self.matrix.user_amplitudes()
        for key, amplitude in state.amplitudes.items():
            self.assertAlmostEqual(amplitude, expected[key[0]], delta=TOL)

    def test_row_state_prep_matches_row(self):
        state = simulate(
            self.matrix.row_state_prep(2).operation.program(),
            {"row_angles": self.matrix.snapshot()["row_angles"]},
        )
        expected = self.matrix.row_amplitudes(2)
        for key, amplitude in state.amplitudes.items():
            self.assertAlmostEqual(amplitude, expected[key[0]], delta=TOL)

    def test_frobenius_and_updates(self):
        frobenius = math.sqrt(
            sum(FMT.decode(word) ** 2 for row in self.matrix.words for word in row)
        )
        self.assertAlmostEqual(self.matrix.frobenius, frobenius)
        changed = self.matrix.update(1, 2, 0.75)
        rebuilt = QMatrix(
            [
                [0.5, 0.5, 0.25, 0.0],
                [0.5, 0.0, 0.75, 0.5],
                [0.25, 0.5, 0.0, 0.5],
                [0.0, 0.0, 0.5, 0.5],
            ],
            fmt=FMT,
            angle_width=AW,
        )
        self.assertEqual(self.matrix.row_angles, rebuilt.row_angles)
        self.assertEqual(self.matrix.root_angles, rebuilt.root_angles)
        self.assertEqual(changed["entries"], {6: FMT.encode(0.75)})

    def test_entries_must_be_nonnegative_representable(self):
        with self.assertRaises(ValidationError):
            QMatrix([[1.0, -0.5], [0.0, 1.0]], fmt=FMT, angle_width=AW)
        with self.assertRaises(ValidationError):
            QMatrix([[99.0, 0.5], [0.0, 1.0]], fmt=FMT, angle_width=AW)


class TwoDimensionalBankLayoutTests(unittest.TestCase):
    """(row, node) 二维寻址与手工 flat 地址逐点一致。"""

    def test_row_angle_bank_flat_layout(self):
        matrix = QMatrix(
            [[0.5, 0.5, 0.25, 0.0], [0.5, 0.0, 0.25, 0.5], [0.25, 0.5, 0.0, 0.5], [0.0, 0.0, 0.5, 0.5]],
            fmt=FMT,
            angle_width=AW,
        )
        bank = matrix.snapshot()["row_angles"]
        for row in range(4):
            flat = qram_state_angles([FMT.decode(w) for w in matrix.words[row]], AW)
            for address, word in flat.items():
                self.assertEqual(bank[row * 4 + address], word)

    def test_qmem_layout_hits_expected_cells(self):
        b = Builder("layout", {"row": UInt(2), "col": UInt(2), "out": UInt(8)}, {"rom": QRAM(4, 8)})
        mem = QMem(b, "rom", shape=(4, 4))
        b.h(b["row"])
        mem[b["row"], 3].load(b["out"])
        state = simulate(
            b.finish().program(),
            {"rom": {(i << 2) | 3: 16 * i + 7 for i in range(4)}},
        )
        for key, amplitude in state.amplitudes.items():
            row = key[0]
            self.assertEqual(key[2], 16 * row + 7)
            self.assertAlmostEqual(abs(amplitude) ** 2, 0.25)


if __name__ == "__main__":
    unittest.main()
