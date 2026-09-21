"""QSVT 标准变换库的数值见证：相位合成往返、约定钉死、求逆、过滤、模拟与定点搜索。"""

import cmath
import math
import random
import unittest

from pyqecclang import ValidationError, simulate
from pyqecclang.algorithms.common.qsvt import (
    eigenstate_filter,
    fixed_point_search,
    fixed_point_search_phases,
    qsp_phases,
    qsp_response,
    qsvt_hamiltonian_simulation,
    qsvt_matrix_inversion,
)
from pyqecclang.algorithms.common.transforms import qsvt_sequence
from pyqecclang.algorithms.input_model.block_encoding import matrix_pauli_encoding


def chebyshev_t(n):
    """T_n 的升幂单项式系数。"""
    if n == 0:
        return (1.0,)
    if n == 1:
        return (0.0, 1.0)
    a, b = (1.0,), (0.0, 1.0)
    for _ in range(2, n + 1):
        nb = tuple(2 * v for v in (0.0,) + b)
        nb = tuple(nb[i] - (a[i] if i < len(a) else 0.0) for i in range(len(nb)))
        a, b = b, nb
    return b


def peval(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))


class PhaseSynthesisTests(unittest.TestCase):
    def test_roundtrip_chebyshev(self):
        for n in range(1, 7):
            phases = qsp_phases(chebyshev_t(n))
            self.assertEqual(len(phases), n + 1)
            for i in range(51):
                x = -0.99 + 1.98 * i / 50
                self.assertAlmostEqual(
                    qsp_response(x, phases), peval(chebyshev_t(n), x), delta=1e-8
                )

    def test_roundtrip_with_explicit_imaginary_part(self):
        # P(x) = (0.5 + i·√0.75)·x：|P|² = x²，Q = 1，精确可实现。
        phases = qsp_phases((0.0, 0.5), imag=(0.0, math.sqrt(0.75)))
        for i in range(21):
            x = -1.0 + 0.1 * i
            self.assertAlmostEqual(
                qsp_response(x, phases), complex(0.5, math.sqrt(0.75)) * x, delta=1e-9
            )

    def test_convention_matches_qsvt_sequence(self):
        # 随机相位下，qsvt_sequence 电路的零信号块必须逐点等于 qsp_response。
        rng = random.Random(7)
        phases = tuple(rng.uniform(-math.pi, math.pi) for _ in range(6))
        matrix = [[0.85, 0.0], [0.0, -0.55]]
        be = matrix_pauli_encoding(matrix)
        op = qsvt_sequence(be, phases)
        for col in range(2):
            x = matrix[col][col] / be.alpha
            state = simulate(op.program(), initial={"target": col})
            amp = state.amplitudes.get((col, 0), 0j)
            self.assertAlmostEqual(amp, qsp_response(x, phases), delta=1e-10)

    def test_negated_phases_conjugate_polynomial(self):
        # 实部提取 LCU 所依赖的性质：−Φ 实现 P̄。
        rng = random.Random(11)
        phases = tuple(rng.uniform(-math.pi, math.pi) for _ in range(5))
        for i in range(21):
            x = -1.0 + 0.1 * i
            self.assertAlmostEqual(
                qsp_response(x, tuple(-p for p in phases)),
                qsp_response(x, phases).conjugate(),
                delta=1e-12,
            )

    def test_input_validation(self):
        with self.assertRaises(ValidationError):  # 奇偶性不符
            qsp_phases((0.0, 0.5, 0.5))
        with self.assertRaises(ValidationError):  # 上界违例
            qsp_phases((0.0, 2.0))
        with self.assertRaises(ValidationError):  # 端点未饱和且未给虚部
            qsp_phases((0.0, 0.5))
        with self.assertRaises(ValidationError):  # 非实系数
            qsp_phases((0.0, 1.0 + 1.0j))
        with self.assertRaises(ValidationError):  # 虚部奇偶性不符
            qsp_phases((0.0, 0.5), imag=(0.5,))


class TransformWitnessTests(unittest.TestCase):
    def test_matrix_inversion_block(self):
        matrix = [[0.6, -0.2], [-0.2, 0.6]]  # 本征值 0.4, 0.8：κ = 2
        be = matrix_pauli_encoding(matrix)
        inv = qsvt_matrix_inversion(be, 2.0, error=0.15)
        attrs = dict(inv.operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "qsvt_matrix_inversion")
        self.assertEqual(attrs["be_alpha"], 1.0)
        scale = attrs["inverse_scale"]
        inverse = [[1.875, 0.625], [0.625, 1.875]]
        for col in range(2):
            state = simulate(inv.operation.program(), initial={"target": col})
            for row in range(2):
                amp = state.amplitudes.get((row, 0), 0j)
                expected = scale * inverse[row][col]
                self.assertAlmostEqual(amp, expected, delta=0.35 * abs(expected))

    def test_eigenstate_filter_isolates_eigenvalue(self):
        # 本征值 0.05 与 0.5（α = 0.5），过滤中心平移到 0.1（谱变量）。
        be = matrix_pauli_encoding([[0.05, 0.0], [0.0, 0.5]])
        flt = eigenstate_filter(be, 0.2, 8, center=0.1)
        attrs = dict(flt.operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "eigenstate_filter")
        self.assertLess(attrs["suppression"], 0.08)
        # 平移后本征值：(0.05−0.1)/0.6 ≈ −0.083（通带内）与 (0.5−0.1)/0.6 ≈ 0.667（阻带）。
        state = simulate(flt.operation.program(), initial={"target": 0})
        self.assertGreater(abs(state.amplitudes.get((0, 0), 0j)), 0.7)
        state = simulate(flt.operation.program(), initial={"target": 1})
        self.assertLess(abs(state.amplitudes.get((1, 0), 0j)), 0.1)

    def test_hamiltonian_simulation_block(self):
        matrix = [[0.5, 0.0], [0.0, -0.25]]
        be = matrix_pauli_encoding(matrix)
        t = 0.7
        hs = qsvt_hamiltonian_simulation(be, t, error=0.01)
        attrs = dict(hs.operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "qsvt_hamiltonian_simulation")
        scale = attrs["sim_scale"]
        for col in range(2):
            x = matrix[col][col] / be.alpha
            state = simulate(hs.operation.program(), initial={"target": col})
            amp = state.amplitudes.get((col, 0), 0j)
            self.assertAlmostEqual(amp, cmath.exp(1j * t * x) / scale, delta=5e-4)

    def test_fixed_point_search_phase_properties(self):
        delta = 0.4
        thresholds = []
        for degree in (3, 5, 7):
            phases = fixed_point_search_phases(delta, degree)
            self.assertEqual(len(phases), degree + 1)
            c = math.cosh(math.acosh(1.0 / delta) / degree)
            theta = math.sqrt(1.0 - 1.0 / c**2)
            thresholds.append(theta)
            for i in range(1, 401):
                x = i / 400
                p2 = abs(qsp_response(x, phases)) ** 2
                self.assertLessEqual(p2, 1.0 + 1e-9)
                # 与 YLC 成功概率闭式一致
                v = c * math.sqrt(1.0 - x * x)
                t_l = math.cos(degree * math.acos(v)) if v <= 1 else math.cosh(
                    degree * math.acosh(v)
                )
                self.assertAlmostEqual(p2, 1.0 - delta**2 * t_l**2, delta=1e-5)
                if x >= theta:
                    self.assertGreaterEqual(p2, 1.0 - delta**2 - 1e-9)
        # 不动点性质：阈值随查询次数单调下降并趋于 0
        self.assertGreater(thresholds[0], thresholds[1])
        self.assertGreater(thresholds[1], thresholds[2])

    def test_fixed_point_search_circuit_amplifies(self):
        be = matrix_pauli_encoding([[0.6, 0.0], [0.0, 0.1]])  # α = 0.6：x = 1.0 与 1/6
        delta, degree = 0.4, 5
        fp = fixed_point_search(be, delta, degree)
        attrs = dict(fp.operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "fixed_point_search")
        self.assertAlmostEqual(attrs["threshold"], math.sqrt(1.0 - 1.0 / 1.04955**2), places=3)
        state = simulate(fp.operation.program(), initial={"target": 0})
        self.assertGreaterEqual(abs(state.amplitudes.get((0, 0), 0j)) ** 2, 1.0 - delta**2 - 1e-9)

    def test_transform_input_validation(self):
        be = matrix_pauli_encoding([[0.5, 0.0], [0.0, 0.25]])
        with self.assertRaises(ValidationError):
            qsvt_matrix_inversion(be, 0.5)
        with self.assertRaises(ValidationError):
            qsvt_matrix_inversion(be, 2.0, error=1.0)
        with self.assertRaises(ValidationError):  # 需要度数超过合成上限
            qsvt_matrix_inversion(be, 10.0, error=1e-9)
        with self.assertRaises(ValidationError):
            eigenstate_filter(be, 1.5, 4)
        with self.assertRaises(ValidationError):
            eigenstate_filter(be, 0.2, 0)
        with self.assertRaises(ValidationError):
            qsvt_hamiltonian_simulation(be, 0.0)
        with self.assertRaises(ValidationError):
            fixed_point_search_phases(1.5, 5)
        with self.assertRaises(ValidationError):  # 度数必须为偶数之外（奇数）
            fixed_point_search_phases(0.3, 4)


if __name__ == "__main__":
    unittest.main()
