"""Prepared ODE instances checked on real PySparQ and UnifiedQuantum backends."""

import math
import unittest
from functools import partial

import oracq
from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.qode.cbmd import ContourPlan
from oracq.algorithms.qode.lchs import QuadraturePlan
from oracq.algorithms.qode.schrodingerization import SchrodingerPlan


class ODEInstanceNativeTests(unittest.TestCase):
    def assert_backends(self, instance):
        instance.prepare()
        circuit = instance.circuit()
        expected = oracq.simulate(circuit.program, circuit.memory)
        sparse, originir = instance.run_pysparq(), instance.run_originir_ext()
        self.assertEqual((sparse.backend, originir.backend), ("pysparq", "originir_ext"))
        for result in (sparse, originir):
            for key in expected.amplitudes.keys() | result.state.amplitudes.keys():
                self.assertAlmostEqual(result.state.amplitudes.get(key, 0), expected.amplitudes.get(key, 0), places=10)
        self.assertAlmostEqual(sparse.success_probability, originir.success_probability, places=10)
        return sparse

    def test_coupled_generator_physical_values_on_both_backends(self):
        p = oracq.ODEProblem(oracq.MatrixInput([[-1, 2], [0, -2]]), oracq.ArrayInput([2, -1]), 0.1)
        result = self.assert_backends(oracq.qode_solve(p))
        for actual, expected in zip(result.physical_values(), (1.64, -0.82), strict=True):
            self.assertAlmostEqual(actual, expected, places=10)

    def test_complex_scalar_and_padding_on_both_backends(self):
        p = oracq.ODEProblem(oracq.MatrixInput([[1j]]), oracq.ArrayInput([2 - 1j]), 0.1)
        result = self.assert_backends(oracq.qode_solve(p))
        self.assertAlmostEqual(result.physical_values()[0], (2 - 1j) * (1 + 0.1j - 0.005), places=10)

    def test_qram_matrix_and_initial_snapshots_on_both_backends(self):
        config = oracq.QRAMMatrixConfig(2, offset=4)
        cells = config.encode_matrix([[-0.5, 0.5], [0, -0.5]], 3)
        bank = oracq.register_qram(cells, 4, 3, name="generator")
        p = oracq.ODEProblem(oracq.QRAMMatrixInput(bank, config), oracq.ArrayInput([1, 1], encoding="qram", angle_width=3), 0.1)
        result = self.assert_backends(oracq.qode_solve(p, config=oracq.ODEConfig(taylor_degree=1)))
        expected = [1 + 0.1 * sum(math.cos(math.pi * cells[4 + 2 * row + col] / 8) for col in range(2)) for row in range(2)]
        for actual, value in zip(result.physical_values(), expected, strict=True):
            self.assertAlmostEqual(actual, value, places=10)

    def test_all_existing_named_algorithms_execute_on_both_backends(self):
        p = oracq.ODEProblem(oracq.MatrixInput([[-0.5, 0], [0, -0.5]]), oracq.UniformInput(), 0.01, dissipative=True)
        for method, plan in [("lchs", QuadraturePlan((0.5,), (1,))),
                             ("cbmd", ContourPlan(cutoff=0)),
                             ("schrodingerization", SchrodingerPlan(auxiliary_width=1, selected_index=0))]:
            with self.subTest(method=method):
                config = oracq.ODEConfig(plan=plan, hamiltonian_function=partial(taylor_hamiltonian, degree=1))
                result = self.assert_backends(oracq.qode_solve(p, method=method, config=config))
                self.assertGreater(result.success_probability, 0)

    def test_zero_solution_on_both_backends(self):
        p = oracq.ODEProblem(oracq.MatrixInput([[-1, 0], [0, -1]]), oracq.UniformInput(0))
        result = self.assert_backends(oracq.qode_solve(p))
        self.assertEqual(result.physical_values(), (0j, 0j))
        self.assertEqual(result.success_probability, 0)
