"""实际后端执行算法协议生成的门，比较完整复幅度。"""

import unittest

from pyqecclang import Bits, Builder, identity, run_originir, run_pysparq, simulate
from pyqecclang.algorithms.common.hamiltonian import PauliHamiltonian, hamiltonian_simulation
from pyqecclang.algorithms.input_model.block_encoding import lcu


class AlgorithmContractsNativeTests(unittest.TestCase):
    def compare_backends(self, program):
        expected = simulate(program).amplitudes
        native = run_pysparq(program).amplitudes
        vector = run_originir(program)
        for key in expected.keys() | native.keys():
            self.assertAlmostEqual(expected.get(key, 0), native.get(key, 0), places=10)
        widths = [r.type.width for r in program.main.registers]
        for index, amplitude in enumerate(vector):
            cursor, key = 0, []
            for width in widths:
                key.append((index >> cursor) & ((1 << width) - 1))
                cursor += width
            self.assertAlmostEqual(expected.get(tuple(key), 0), complex(amplitude), places=10)

    def test_unitary_input_and_trotter_protocol(self):
        b = Builder("PhaseGateInput", {"q": Bits(1)})
        b.h(b["q"])
        b.gate("phase", b["q"], 0.3)
        self.compare_backends(lcu(((1j, b.finish()), (2, identity(1)))).operation.program())
        h = PauliHamiltonian(((0.3, "I"), (0.7, "X")))
        self.compare_backends(hamiltonian_simulation(h, 0.4, steps=3).operation.program())
