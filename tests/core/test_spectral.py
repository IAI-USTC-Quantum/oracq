"""谱原语的语义验证：傅里叶相位、稀疏谱块编码与谱态制备对照 DFT 参考。"""

import cmath
import math
import unittest

from pyqecclang import Builder, estimate_resources, simulate
from pyqecclang.algorithms.input_model.oracles import gate_state_prep
from pyqecclang.algorithms.input_model.spectral import (
    fourier_phase,
    fourier_phase_encoding,
    frequency_amplitudes,
    pruned_state_prep,
    spectral_diagonal,
    spectral_state_prep,
)


def applied_state(operation, initial):
    """把指定寄存器置为经典初值后调用 operation，返回参考执行器的振幅。"""
    b = Builder("column_probe", {r.name: r.type for r in operation.module.registers})
    for key, value in initial.items():
        for bit in range(b[key].width):
            if (value >> bit) & 1:
                b.x(b[key][bit])
    b.call(operation, **{r.name: b[r.name] for r in operation.module.registers})
    return simulate(b.finish().program()).amplitudes


def reference_value(spectrum, width, index):
    return sum(
        coefficient * cmath.exp(2j * math.pi * k * index / (1 << width))
        for k, coefficient in spectrum.items()
    )


class FourierPhaseTests(unittest.TestCase):
    def test_matches_dft_diagonal(self):
        for width in (2, 3):
            for k in (0, 1, 2, -1, -3, 5):
                operation = fourier_phase(width, k)
                dim = 1 << width
                for column in range(dim):
                    amplitudes = applied_state(operation, {"target": column})
                    expected = cmath.exp(2j * math.pi * k * column / dim)
                    self.assertAlmostEqual(
                        amplitudes[(column,)], expected, places=12, msg=f"w={width} k={k} J={column}"
                    )

    def test_encoding_declares_unit_alpha(self):
        be = fourier_phase_encoding(3, 5)
        self.assertEqual(be.alpha, 1.0)
        self.assertEqual(be.signal_qubits, 0)


class SpectralDiagonalTests(unittest.TestCase):
    CASES = [
        {0: 0.6 + 0.1j, 1: -0.3 + 0.2j, 2: 0.4},
        {-2: 0.5, -1: 0.25j, 0: -0.35, 1: 0.2},
        {3: 1.0},
    ]

    def _corner(self, be, width):
        dim = 1 << width
        corner = {}
        for column in range(dim):
            amplitudes = applied_state(be.operation, {"target": column})
            corner[column] = amplitudes[(column, 0)]
        return corner

    def test_sequential_matches_reference(self):
        for spectrum in self.CASES:
            width = 3
            be = spectral_diagonal(spectrum, width)
            alpha = sum(abs(c) for c in spectrum.values())
            corner = self._corner(be, width)
            for column, amplitude in corner.items():
                expected = reference_value(spectrum, width, column) / alpha
                self.assertAlmostEqual(amplitude, expected, places=11)

    def test_variants_agree(self):
        for spectrum in self.CASES:
            width = 3
            sequential = self._corner(spectral_diagonal(spectrum, width), width)
            naive = self._corner(spectral_diagonal(spectrum, width, variant="naive"), width)
            for column in sequential:
                self.assertAlmostEqual(sequential[column], naive[column], places=11)

    def test_resources_favor_sequential(self):
        spectrum = {k: 0.3 + 0.02j * k for k in range(-15, 16)}
        sequential = estimate_resources(spectral_diagonal(spectrum, 6).operation.program())
        naive = estimate_resources(
            spectral_diagonal(spectrum, 6, variant="naive").operation.program()
        )
        self.assertLess(len(sequential.rotations), len(naive.rotations))
        self.assertLess(sequential.atoms.get("toffoli", 0), naive.atoms.get("toffoli", 0))


class SpectralStatePrepTests(unittest.TestCase):
    def test_matches_normalized_field(self):
        spectrum = {-1: 0.4 + 0.2j, 0: 0.9, 1: -0.3j, 2: 0.25}
        width = 3
        prep = spectral_state_prep(spectrum, width)
        amplitudes = applied_state(prep.operation, {})
        norm = math.sqrt(sum(abs(reference_value(spectrum, width, j)) ** 2 for j in range(1 << width)))
        for index in range(1 << width):
            expected = reference_value(spectrum, width, index) / norm
            self.assertAlmostEqual(amplitudes[(index, 0)], expected, places=11)

    def test_pruned_matches_gate_prep(self):
        amplitudes = [0.5, 0, 0, 0.5j, 0, -0.5, 0, 0.5]
        pruned = applied_state(pruned_state_prep(amplitudes).operation, {})
        plain = applied_state(gate_state_prep(amplitudes).operation, {})
        self.assertEqual(set(pruned), set(plain))
        for key, value in plain.items():
            self.assertAlmostEqual(pruned[key], value, places=12)


class InterfaceTests(unittest.TestCase):
    def test_rejects_disjoint_band(self):
        from pyqecclang.infrastructure.ir import ValidationError

        with self.assertRaises(ValidationError):
            spectral_diagonal({0: 1.0, 5: 0.5}, 3)

    def test_frequency_amplitudes_aliases_modulo(self):
        amplitudes = frequency_amplitudes((0.3, 0.4), -1, 3)
        self.assertAlmostEqual(abs(amplitudes[7]) ** 2 + abs(amplitudes[0]) ** 2, 1.0)


if __name__ == "__main__":
    unittest.main()
