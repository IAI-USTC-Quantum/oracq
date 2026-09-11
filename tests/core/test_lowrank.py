"""化学低秩分解（DF/THC）块编码的数值见证。"""

import math
import unittest

from witness import assert_block_equals, block_column

from pyqecclang import ValidationError, simulate, unresolved
from pyqecclang.algorithms.block_encoding import matrix_pauli_encoding
from pyqecclang.algorithms.lowrank import (
    DoubleFactorization,
    THCDecomposition,
    diagonalize_symmetric,
    double_factorized_encoding,
    thc_encoding,
)
from pyqecclang.algorithms.transforms import qubitization_walk

IDENTITY2 = ((1.0, 0.0), (0.0, 1.0))
HADAMARD = ((1 / math.sqrt(2), 1 / math.sqrt(2)), (1 / math.sqrt(2), -1 / math.sqrt(2)))


def matmul(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0])))
        for i in range(len(a))
    )


def transpose_conj(a):
    return tuple(tuple(a[j][i].conjugate() for j in range(len(a))) for i in range(len(a[0])))


def add(*matrices):
    return tuple(
        tuple(sum(m[i][j] for m in matrices) for j in range(len(matrices[0])))
        for i in range(len(matrices[0]))
    )


def scale_matrix(c, a):
    return tuple(tuple(c * v for v in row) for row in a)


def pauli_l1(matrix):
    """2×2 矩阵的 Pauli l1 范数（独立闭式，不经过 matrix_pauli_encoding）。

    M = cI·I + x·X + y·Y + z·Z，系数由 M = ((a, b), (c, d)) 的元素线性表出。
    """
    a, b = matrix[0]
    c, d = matrix[1]
    return abs((a + d) / 2) + abs((b + c) / 2) + abs((c - b) / 2) + abs((a - d) / 2)


class DiagonalizeSymmetricTests(unittest.TestCase):
    def test_reconstructs_factor(self):
        factor = ((1.2, 0.4), (0.4, 0.8))
        values, vectors = diagonalize_symmetric(factor)
        self.assertAlmostEqual(sum(values), 2.0, places=12)
        reconstructed = tuple(
            tuple(
                sum(values[k] * vectors[i][k] * vectors[j][k] for k in range(2))
                for j in range(2)
            )
            for i in range(2)
        )
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(reconstructed[i][j], factor[i][j], places=10)
        # 特征向量矩阵正交。
        for i in range(2):
            for j in range(2):
                dot = sum(vectors[k][i] * vectors[k][j] for k in range(2))
                self.assertAlmostEqual(dot, 1.0 if i == j else 0.0, places=10)

    def test_rejects_non_symmetric(self):
        with self.assertRaises(ValidationError):
            diagonalize_symmetric(((1.0, 0.4), (0.3, 0.8)))


class DoubleFactorizationTests(unittest.TestCase):
    def test_single_rank_block_equals_hamiltonian(self):
        hamiltonian = ((1.2, 0.4), (0.4, 0.8))
        df = DoubleFactorization.from_symmetric(0.0, (IDENTITY2,), (hamiltonian,))
        self.assertEqual(df.rank, 1)
        self.assertEqual(df.width, 1)
        be = double_factorized_encoding(df)
        attrs = dict(be.operation.module.attributes)
        self.assertEqual(attrs["be_form"], "double_factorization")
        self.assertAlmostEqual(be.alpha, sum(abs(v) for v in (1.4472135954999579, 0.5527864045000421)))
        assert_block_equals(self, be, hamiltonian)

    def test_scalar_and_rotated_ranks(self):
        g1 = ((0.9, 0.3), (0.3, 0.5))
        g2 = ((0.3, -0.1), (-0.1, 0.4))
        scalar = 0.2
        df = DoubleFactorization.from_symmetric(scalar, (IDENTITY2, HADAMARD), (g1, g2))
        rotated = matmul(matmul(HADAMARD, g2), HADAMARD)
        expected = add(
            scale_matrix(scalar, IDENTITY2),
            g1,
            rotated,
        )
        be = double_factorized_encoding(df)
        self.assertAlmostEqual(
            be.alpha,
            abs(scalar)
            + sum(abs(v) for v in diagonalize_symmetric(g1)[0])
            + sum(abs(v) for v in diagonalize_symmetric(g2)[0]),
            places=10,
        )
        assert_block_equals(self, be, expected)

    def test_matches_pauli_encoding_block(self):
        # 同一哈密顿量的 DF 编码与 Pauli LCU 编码的 (0,0) 块一致（alpha 各自不同）。
        hamiltonian = ((1.2, 0.4), (0.4, 0.8))
        df = DoubleFactorization.from_symmetric(0.0, (IDENTITY2,), (hamiltonian,))
        df_be = double_factorized_encoding(df)
        pauli_be = matrix_pauli_encoding(hamiltonian)
        for column in range(2):
            df_column = block_column(df_be, column)
            pauli_state = simulate(
                pauli_be.operation.program(), initial={"target": column}
            )
            for row in range(2):
                self.assertAlmostEqual(
                    df_column[row],
                    pauli_state.amplitudes.get((row, 0), 0) * pauli_be.alpha,
                    places=9,
                )

    def test_feeds_qubitization_walk(self):
        hamiltonian = ((1.2, 0.4), (0.4, 0.8))
        df = DoubleFactorization.from_symmetric(0.0, (IDENTITY2,), (hamiltonian,))
        be = double_factorized_encoding(df)
        walk = qubitization_walk(be)
        self.assertFalse(unresolved(walk.program()))
        names = [r.name for r in walk.module.registers]
        self.assertEqual(names, ["target", "signal"])

    def test_alpha_matches_closed_form_eigenvalues(self):
        # 独立闭式切断共用实现：2×2 对称矩阵 λ = (t ± √(t²−4d))/2，
        # t 为迹、d 为行列式；be.alpha 须等于 Σ|λ|（实测 2.0，两端独立计算）。
        hamiltonian = ((1.2, 0.4), (0.4, 0.8))
        df = DoubleFactorization.from_symmetric(0.0, (IDENTITY2,), (hamiltonian,))
        be = double_factorized_encoding(df)
        trace = hamiltonian[0][0] + hamiltonian[1][1]
        determinant = hamiltonian[0][0] * hamiltonian[1][1] - hamiltonian[0][1] ** 2
        gap = math.sqrt(trace * trace - 4 * determinant)
        expected = (abs(trace + gap) + abs(trace - gap)) / 2
        self.assertAlmostEqual(be.alpha, expected, places=10)

    def test_df_alpha_tighter_than_pauli(self):
        # 非对角项主导的正定 g：PSD 故 DF α = Σ|λ| = tr(g) = 1.4（实测），
        # Pauli LCU α = tr/2 + |g01| + |g00−g11|/2 = 1.5（实测）。理论条件：
        # 单比特 PSD 矩阵 Pauli α − DF α = |g01| + |Δ/2| − tr/2，非对角主导时为正。
        g = ((1.1, 0.4), (0.4, 0.3))
        df_be = double_factorized_encoding(
            DoubleFactorization.from_symmetric(0.0, (IDENTITY2,), (g,))
        )
        pauli_be = matrix_pauli_encoding(g)
        self.assertAlmostEqual(df_be.alpha, 1.4, places=10)
        self.assertLessEqual(df_be.alpha, pauli_be.alpha)


class ThcTests(unittest.TestCase):
    def test_thc_block_equals_hamiltonian(self):
        l0 = ((1.0, 0.2), (0.0, 1.0))
        l1 = ((0.5, 0.0), (0.3, 1.0))
        zeta = ((0.7, 0.1), (0.1, 0.4))
        thc = THCDecomposition(zeta, (l0, l1))
        self.assertEqual(thc.leaf_count, 2)
        leaves = (l0, l1)
        expected = add(
            *(
                scale_matrix(
                    zeta[mu][nu],
                    matmul(leaves[mu], transpose_conj(leaves[nu])),
                )
                for mu in range(2)
                for nu in range(2)
            )
        )
        be = thc_encoding(thc)
        attrs = dict(be.operation.module.attributes)
        self.assertEqual(attrs["be_form"], "thc")
        self.assertEqual(attrs["thc_leaves"], 2)
        assert_block_equals(self, be, expected)

    def test_thc_alpha_matches_hand_computed_bound(self):
        # thc_encoding 的 α 口径：Σ_{μν} |ζ_{μν}|·α_μ·α_ν，α_μ 为叶矩阵的
        # Pauli l1 上界。此处用 pauli_l1 闭式独立手算：α_0 = 1.2、α_1 = 1.3，
        # α = 0.7·1.44 + 0.1·1.56 + 0.1·1.56 + 0.4·1.69 = 1.996（实测一致）。
        l0 = ((1.0, 0.2), (0.0, 1.0))
        l1 = ((0.5, 0.0), (0.3, 1.0))
        zeta = ((0.7, 0.1), (0.1, 0.4))
        be = thc_encoding(THCDecomposition(zeta, (l0, l1)))
        alphas = (pauli_l1(l0), pauli_l1(l1))
        expected = sum(
            abs(zeta[mu][nu]) * alphas[mu] * alphas[nu]
            for mu in range(2)
            for nu in range(2)
        )
        self.assertAlmostEqual(be.alpha, expected, places=10)

    def test_invalid_inputs_fail_at_generation(self):
        bad_rotation = ((1.0, 1.0), (0.0, 1.0))
        cases = [
            lambda: DoubleFactorization(0.0, (), ()),
            lambda: DoubleFactorization(0.0, (IDENTITY2,), ((1.0,),)),
            lambda: DoubleFactorization(0.0, (bad_rotation,), ((1.0, 1.0),)),
            lambda: DoubleFactorization("x", (IDENTITY2,), ((1.0, 1.0),)),
            lambda: DoubleFactorization.from_symmetric(0.0, (IDENTITY2,), ()),
            lambda: THCDecomposition(((0.7, 0.2), (0.1, 0.4)), (IDENTITY2, IDENTITY2)),
            lambda: THCDecomposition(((0.7, 0.1),), (IDENTITY2, IDENTITY2)),
            lambda: thc_encoding(THCDecomposition(((0.0, 0.0), (0.0, 0.0)), (IDENTITY2, IDENTITY2))),
            lambda: diagonalize_symmetric(((1.0, 0.4), (0.3, 0.8))),
        ]
        for case in cases:
            with self.subTest(case=case):
                with self.assertRaises(ValidationError):
                    case()


if __name__ == "__main__":
    unittest.main()
