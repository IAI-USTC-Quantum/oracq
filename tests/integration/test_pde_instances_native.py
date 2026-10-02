"""Prepared PDE instances cross-checked with real PySparQ and OriginIR execution."""

import unittest

import oracq


class PDEInstanceNativeTests(unittest.TestCase):
    def assert_backends(self, instance):
        instance.prepare()
        circuit = instance.circuit()
        expected = oracq.simulate(circuit.program, circuit.memory)
        sparse = instance.run_pysparq()
        originir = instance.run_originir_ext()
        self.assertEqual(sparse.backend, "pysparq")
        self.assertEqual(originir.backend, "originir_ext")
        for result in (sparse, originir):
            for key in expected.amplitudes.keys() | result.state.amplitudes.keys():
                self.assertAlmostEqual(result.state.amplitudes.get(key, 0), expected.amplitudes.get(key, 0), places=10)
        self.assertAlmostEqual(sparse.success_probability, originir.success_probability, places=10)
        return sparse

    def test_heat_evolution_and_physical_magnitude_on_both_backends(self):
        p = oracq.PDEProblem(grid=oracq.UniformGrid1D(0, 2, 2), initial=oracq.ArrayInput([0.2, 0.1]), type="heat", final_time=0.01)
        result = self.assert_backends(oracq.qpde_solve(p))
        for actual, expected in zip(result.physical_values()["u"], (0.1998, 0.1002), strict=True):
            self.assertAlmostEqual(actual, expected, places=10)

    def test_qram_initial_snapshot_is_bound_on_both_backends(self):
        p = oracq.PDEProblem(grid=oracq.UniformGrid1D(0, 2, 2), initial=oracq.ArrayInput([0.2, 0.2], encoding="qram", angle_width=4), type="heat", final_time=0.01)
        result = self.assert_backends(oracq.qpde_solve(p))
        for actual in result.physical_values()["u"]:
            self.assertAlmostEqual(actual, 0.2, places=10)

    def test_qram_mesh_queries_use_the_same_operator_on_both_backends(self):
        d = (("x", 2),)
        config = oracq.UnstructuredGridConfig(2, {d: 0}, value_scale=2)
        qram = oracq.register_qram(config.encode_matrices({d: [[-1, 1], [1, -1]]}, 3), 2, 3, name="mesh")
        grid = oracq.UnstructuredGrid(qram, config)
        p = oracq.PDEProblem(grid=grid, initial=oracq.ArrayInput([0.2, 0.1]), type="heat", final_time=0.01)
        result = self.assert_backends(oracq.qpde_solve(p))
        expected = [value + 0.001 * sum(weight * (0.2, 0.1)[col] for col, weight in grid.derivative_row(d, row)) for row, value in enumerate((0.2, 0.1))]
        for actual, value in zip(result.physical_values()["u"], expected, strict=True):
            self.assertAlmostEqual(actual, value, places=10)
