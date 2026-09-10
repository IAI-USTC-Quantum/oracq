"""Jordan 梯度估计的数值见证：线性精确读出、扰动收敛与绑定一致性。"""

import math
import unittest

from pyqecclang import ValidationError, bind, simulate, unresolved
from pyqecclang.algorithms.gradient import (
    abstract_phase_oracle,
    function_phase_oracle,
    gate_phase_oracle,
    gradient_estimation,
    gradient_from_readout,
)


def linear_oracle(coefficients, grid_bits, *, name=None):
    """线性函数 f(x) = Σ_i c_i x_i 的显式相位表 oracle（Jordan 缩放）。"""
    dimension = len(coefficients)
    grid_points = 1 << grid_bits
    angles = []
    for value in range(1 << (dimension * grid_bits)):
        phase = 0.0
        for i, coefficient in enumerate(coefficients):
            chunk = (value >> (i * grid_bits)) & (grid_points - 1)
            phase += coefficient * chunk / grid_points
        angles.append(2 * math.pi * grid_points * phase)
    return gate_phase_oracle(
        dimension * grid_bits, angles, phase_scale=grid_points, name=name
    )


def readout_distribution(state):
    """target 为唯一根寄存器时的读出分布。"""
    return {key[0]: abs(amplitude) ** 2 for key, amplitude in state.amplitudes.items()}


def mode(distribution):
    return max(distribution, key=distribution.get)


class GradientTests(unittest.TestCase):
    def test_linear_function_exact_two_dimensions(self):
        # 梯度分量取网格可精确表示的值，读出寄存器应确定性地落在编码值上。
        grid_bits = 4
        oracle = linear_oracle((3 / 16, -2 / 16), grid_bits)
        operation = gradient_estimation(oracle, dimension=2, grid_bits=grid_bits)
        attrs = dict(operation.module.attributes)
        self.assertEqual(attrs["algorithm"], "jordan_gradient")
        self.assertEqual(attrs["oracle_queries"], 1)
        measured = readout_distribution(simulate(operation.program()))
        self.assertEqual(len(measured), 1)
        self.assertAlmostEqual(measured[mode(measured)], 1.0, places=12)
        self.assertEqual(
            gradient_from_readout(mode(measured), dimension=2, grid_bits=grid_bits),
            (3 / 16, -2 / 16),
        )

    def test_perturbed_linear_concentrates_with_grid_bits(self):
        # f(x) = a·x + x²/N²（N = 2**m）：相位扰动为 o(1/N)，正确读出的概率单调上升。
        a_numerator, a_denominator = 3, 8
        probabilities = []
        for grid_bits in (3, 4, 5):
            grid_points = 1 << grid_bits
            angles = [
                2
                * math.pi
                * grid_points
                * (
                    a_numerator / a_denominator * value / grid_points
                    + (value / grid_points) ** 2 / grid_points**2
                )
                for value in range(grid_points)
            ]
            oracle = gate_phase_oracle(grid_bits, angles, phase_scale=grid_points)
            operation = gradient_estimation(oracle, dimension=1, grid_bits=grid_bits)
            measured = readout_distribution(simulate(operation.program()))
            exact = a_numerator * grid_points // a_denominator
            self.assertEqual(mode(measured), exact)
            probabilities.append(measured[exact])
        self.assertLess(probabilities[0], probabilities[1])
        self.assertLess(probabilities[1], probabilities[2])

    def test_abstract_oracle_binds_to_gate_implementation(self):
        grid_bits = 3
        oracle = abstract_phase_oracle("JordanPhase", 2 * grid_bits, phase_scale=1 << grid_bits)
        operation = gradient_estimation(oracle, dimension=2, grid_bits=grid_bits)
        program = operation.program()
        self.assertEqual([r.name for r in unresolved(program)], ["JordanPhase"])
        bound = bind(
            program,
            {"JordanPhase": linear_oracle((1 / 8, 1 / 4), grid_bits).operation},
        )
        self.assertFalse(unresolved(bound))
        measured = readout_distribution(simulate(bound))
        self.assertEqual(
            gradient_from_readout(mode(measured), dimension=2, grid_bits=grid_bits),
            (1 / 8, 1 / 4),
        )

    def test_function_phase_oracle_from_mathfunc(self):
        # mathfunc 路径：f(x) = 0.25·x 在定点格式下精确，读出应恢复 1/4。
        # 梯度分量只在 mod 1 意义下可分辨（相位 e^{2πi·N·a·j} 对 a 与 a±1 相同），
        # 见证因此取 |a| < 1/2 的分量。
        grid_bits = 2
        oracle = function_phase_oracle(
            "def f(x0):\n    return 0.25 * x0\n",
            dimension=1,
            grid_bits=grid_bits,
        )
        self.assertEqual(oracle.phase_scale, 1 << grid_bits)
        operation = gradient_estimation(oracle, dimension=1, grid_bits=grid_bits)
        measured = readout_distribution(simulate(operation.program()))
        self.assertEqual(
            gradient_from_readout(mode(measured), dimension=1, grid_bits=grid_bits),
            (0.25,),
        )

    def test_invalid_inputs_fail_at_generation(self):
        grid_bits = 3
        oracle = linear_oracle((1 / 8,), grid_bits)
        cases = [
            lambda: gradient_estimation(oracle, dimension=2, grid_bits=grid_bits),
            lambda: gradient_estimation(oracle, dimension=1, grid_bits=grid_bits + 1),
            lambda: gradient_estimation(oracle, dimension=0, grid_bits=grid_bits),
            lambda: gradient_estimation(oracle, dimension=1, grid_bits=0),
            lambda: gate_phase_oracle(2, [0.0, 0.0], phase_scale=4),
            lambda: gate_phase_oracle(2, [0.0] * 4, phase_scale=0),
            lambda: abstract_phase_oracle("Bad", 65, phase_scale=1),
            lambda: function_phase_oracle(
                "def f(x0):\n    return x0\n", dimension=1, grid_bits=0
            ),
            lambda: gradient_from_readout(-1, dimension=1, grid_bits=2),
        ]
        for case in cases:
            with self.subTest(case=case):
                with self.assertRaises(ValidationError):
                    case()


if __name__ == "__main__":
    unittest.main()
