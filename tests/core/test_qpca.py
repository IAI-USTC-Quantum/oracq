"""QPCA 与密度矩阵指数化的数值见证（LMR 协议）。"""

import math
import unittest

from pyqecclang import ValidationError, simulate
from pyqecclang.algorithms.density import gate_purification, partial_trace, trace_distance
from pyqecclang.algorithms.oracles import gate_state_prep
from pyqecclang.algorithms.qpca import (
    density_matrix_exponentiation,
    eigenvalue_from_phase,
    qpca,
)

SQRT2 = math.sqrt(2)
PLUS = [1 / SQRT2, 1 / SQRT2]
MINUS = [1 / SQRT2, -1 / SQRT2]


def dense_system(state, copies_width, copies):
    vector = [0j] * (1 << (1 + copies_width * copies))
    for key, amplitude in state.amplitudes.items():
        vector[key[0] + (key[1] << 1)] += amplitude
    return partial_trace(vector, 1, copies_width * copies)


def exact_evolved_pure(time):
    """e^{-it·|+⟩⟨+|} 作用在 ``|0⟩⟨0|`` 上的精确结果。"""
    import cmath

    a0 = (1 + cmath.exp(-1j * time)) / 2
    a1 = -(1 - cmath.exp(-1j * time)) / 2
    return (
        (a0 * a0.conjugate(), a0 * a1.conjugate()),
        (a1 * a0.conjugate(), a1 * a1.conjugate()),
    )


def phase_mode(state):
    distribution = {}
    for key, amplitude in state.amplitudes.items():
        distribution[key[-1]] = distribution.get(key[-1], 0.0) + abs(amplitude) ** 2
    return max(distribution, key=distribution.get), distribution


def circular_eigenvalue(distribution, precision, step):
    """相位分布的圆周均值解码：对展宽峰稳健的 λ 估计。"""
    import cmath

    z = sum(
        p * cmath.exp(2j * math.pi * v / (1 << precision))
        for v, p in distribution.items()
    )
    return -2 * math.pi * (cmath.phase(z) / (2 * math.pi)) / step


class DensityMatrixExponentiationTests(unittest.TestCase):
    def test_small_step_first_order_accurate(self):
        # 单步误差 O(Δt²)：Δt = 0.05 时迹距离在 Δt² 量级内（系数约 0.53）。
        step = 0.05
        prep = gate_state_prep(PLUS)
        operation = density_matrix_exponentiation(prep, time=step, copies=1)
        reduced = dense_system(simulate(operation.program()), 1, 1)
        self.assertLess(trace_distance(reduced, exact_evolved_pure(step)), 0.6 * step * step)

    def test_error_halves_with_copies(self):
        # LMR 一阶标度：固定总时间 t，拷贝数翻倍时误差近似减半。
        time = 0.4
        prep = gate_state_prep(PLUS)
        distances = []
        for copies in (1, 2, 4):
            operation = density_matrix_exponentiation(prep, time=time, copies=copies)
            reduced = dense_system(simulate(operation.program()), 1, copies)
            distances.append(trace_distance(reduced, exact_evolved_pure(time)))
        self.assertLess(distances[1], distances[0] * 0.6)
        self.assertLess(distances[2], distances[1] * 0.6)


class QpcaTests(unittest.TestCase):
    def test_pure_state_eigenvalues(self):
        # ρ = |+⟩⟨+|：本征值 1 与 0 分别由系统输入 |+⟩ 与 |−⟩ 读出。
        step = math.pi / 4  # λ=1 时 φ = 7/8，precision=3 下读出确定。
        prep = gate_state_prep(PLUS)
        operation = qpca(prep, precision=3, step_time=step, system=gate_state_prep(PLUS))
        attrs = dict(operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "qpca")
        self.assertEqual(attrs["copies"], 7)
        # |+⟩ 是拷贝态本身：每步部分交换作用在 SWAP 对称本征态上，对任意 Δt 都精确。
        mode, distribution = phase_mode(simulate(operation.program()))
        self.assertAlmostEqual(distribution[mode], 1.0, places=9)
        self.assertAlmostEqual(eigenvalue_from_phase(mode, 3, step), 1.0, places=12)

        # |−⟩ 是一般本征态：大 Δt 下 LMR 误差表现为峰展宽，但分布的圆周均值
        # 仍以 φ = 0（λ = 0）为中心。
        operation = qpca(prep, precision=3, step_time=step, system=gate_state_prep(MINUS))
        _, distribution = phase_mode(simulate(operation.program()))
        self.assertAlmostEqual(circular_eigenvalue(distribution, 3, step), 0.0, delta=0.15)

    def test_mixed_state_eigenvalue_via_purification(self):
        # ρ = diag(0.75, 0.25) 经纯化适配；系统输入 |0⟩ 读出 λ = 0.75。
        step = 2 * math.pi / 6  # λ=0.75 时 φ = 7/8。
        purification = gate_purification(((0.75 + 0j, 0j), (0j, 0.25 + 0j)))
        prep = purification.as_state_preparation()
        operation = qpca(prep, precision=3, step_time=step, swap_width=1)
        mode, distribution = phase_mode(simulate(operation.program()))
        # 大 Δt 的 LMR 误差只展宽峰而不移峰位：mode 解码恰为 0.75。
        self.assertGreater(distribution[mode], 0.4)
        self.assertAlmostEqual(eigenvalue_from_phase(mode, 3, step), 0.75, places=12)
        self.assertAlmostEqual(circular_eigenvalue(distribution, 3, step), 0.75, delta=0.15)

    def test_invalid_inputs_fail_at_generation(self):
        prep = gate_state_prep(PLUS)
        cases = [
            lambda: density_matrix_exponentiation(prep, time=0.0, copies=1),
            lambda: density_matrix_exponentiation(prep, time=1.0, copies=0),
            lambda: density_matrix_exponentiation(object(), time=1.0, copies=1),
            lambda: qpca(prep, precision=0, step_time=1.0),
            lambda: qpca(prep, precision=7, step_time=1.0),
            lambda: qpca(prep, precision=3, step_time=0.0),
            lambda: qpca(prep, precision=3, step_time=1.0, swap_width=2),
            lambda: qpca(
                prep, precision=3, step_time=1.0, system=gate_state_prep([1, 0, 0, 0])
            ),
            lambda: eigenvalue_from_phase(0, 3, 0.0),
            lambda: eigenvalue_from_phase(8, 3, 1.0),
        ]
        for case in cases:
            with self.subTest(case=case):
                with self.assertRaises(ValidationError):
                    case()


if __name__ == "__main__":
    unittest.main()
