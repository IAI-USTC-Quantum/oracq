"""Input composition, physical scaling, memory lineage, and prepared PDE contracts."""

import math
import subprocess
import sys
import unittest
from dataclasses import replace

import oracq
from oracq.algorithms.input_model.oracles import abstract_state_prep, gate_state_prep
from oracq.applications.qham import Field, PolynomialPDE


class PDEInstanceTests(unittest.TestCase):
    def reference_result(self, instance):
        circuit = instance.circuit()
        return oracq.PDEResult(
            oracq.simulate(circuit.program, circuit.memory),
            oracq.InitialLayout(instance.problem.grid, instance.problem.equation.fields),
            "reference", instance.prepared.amplitude_scale,
        )

    def test_register_qram_validates_and_snapshots_words(self):
        data = {1: 3}
        bank = oracq.register_qram(data, 4, 5, name="samples")
        data[1] = 9
        self.assertEqual(bank.snapshot(), {1: 3})
        self.assertEqual(bank.type, oracq.QRAM(4, 5))
        query = bank.database()
        b = oracq.Builder("query", {"address": oracq.Bits(4), "data": oracq.Bits(5)}, {"samples": bank.type})
        b.x(b["address"][0])
        b.call(query.operation, resources={"samples": "samples"}, address=b["address"], data=b["data"])
        self.assertEqual(oracq.simulate(b.finish().program(), {"samples": bank.snapshot()}).amplitudes, {(1, 3): 1})
        for data, a, w in [({-1: 1}, 2, 3), ({4: 1}, 2, 3), ([8], 2, 3), ([True], 2, 3), ([0], 0, 3)]:
            with self.subTest(data=data, a=a, w=w), self.assertRaises(oracq.ValidationError):
                oracq.register_qram(data, a, w)

    def test_uniform_interval_boundary_and_physical_padding(self):
        grid = oracq.UniformGrid1D(0, 10, 10)
        self.assertEqual((grid.size, grid.spatial_width, grid.spacing, grid.coordinate(9)), (10, 4, 1, 9))
        self.assertEqual(dict(grid.derivative_row((("x", 1),), 0)), {1: 0.5, 9: -0.5})
        initial = oracq.UniformInput(2).prepare_initial(oracq.InitialLayout(grid, ("u",)))
        state = oracq.simulate(initial.preparation.operation.program())
        self.assertAlmostEqual(initial.norm, 2 * math.sqrt(10))
        for i in range(16):
            self.assertAlmostEqual(state.amplitudes.get((i, 0), 0), 1 / math.sqrt(10) if i < 10 else 0)
        interior = oracq.UniformGrid1D(0, 1, 3, boundary="dirichlet_zero")
        self.assertEqual([interior.coordinate(i) for i in range(3)], [0.25, 0.5, 0.75])
        with self.assertRaises(oracq.ValidationError):
            oracq.UniformGrid1D(1, 0, 2)
        with self.assertRaises(oracq.ValidationError):
            oracq.UniformGrid1D(0, 1, 33).derivative_encoding((("x", 1),))

    def test_spatial_encoding_keeps_padding_zero(self):
        grid = oracq.UniformGrid1D(0, 1, 3, boundary="dirichlet_zero")
        derivative = (("x", 2),)
        be = grid.derivative_encoding(derivative)
        for column in (0, 2, 3):
            b = oracq.Builder("column", {r.name: r.type for r in be.operation.module.registers})
            for bit in range(grid.spatial_width):
                if column >> bit & 1:
                    b.x(b["target"][bit])
            b.call(be.operation, target=b["target"], signal=b["signal"])
            state = oracq.simulate(b.finish().program())
            for row in range(4):
                self.assertAlmostEqual(state.amplitudes.get((row, 0), 0) * be.alpha, dict(grid.derivative_row(derivative, row)).get(column, 0))

    def test_burgers_recipe_prepares_a_configurable_ten_node_grid(self):
        p = oracq.PDEProblem(type="burgers", grid=oracq.UniformGrid1D(0, 10, 10), initial=oracq.UniformInput(), final_time=0.01)
        instance = oracq.qpde_solve(p, method="QHAM").prepare()
        self.assertEqual(instance.prepared.state.width, 4)
        self.assertEqual(instance.prepared.model.plan.pde.degree, 2)
        self.assertTrue(instance.prepared.report.ok)
        self.assertTrue(instance.circuit().resource_estimate().complete)

    def test_lifecycle_and_finite_evolution_physical_values(self):
        p = oracq.PDEProblem(grid=oracq.UniformGrid1D(0, 2, 2), initial=oracq.ArrayInput([0.2, 0.1]), type="heat", final_time=0.01)
        instance = oracq.qpde_solve(p)
        for action in (instance.circuit, instance.run_pysparq, instance.run_originir_ext):
            with self.assertRaisesRegex(oracq.ValidationError, "prepare"):
                action()
        self.assertIs(instance.prepare(), instance)
        circuit = instance.circuit()
        instance.prepare()
        self.assertIs(instance.circuit(), circuit)
        self.assertTrue(instance.prepared.report.ok)
        self.assertEqual(oracq.loads(circuit.dumps()), circuit.program)
        self.assertIn("DEF ", circuit.originir_ext().text)
        result = self.reference_result(instance)
        for actual, expected in zip(result.physical_values()["u"], (0.1998, 0.1002), strict=True):
            self.assertAlmostEqual(actual, expected)
        self.assertAlmostEqual(sum(abs(v) ** 2 for v in result.physical_amplitudes()["u"]), 1)

    def test_initial_inputs_are_independent_of_grid_and_algorithm(self):
        grid = oracq.UniformGrid1D(0, 2, 2)
        p = oracq.PDEProblem(grid=grid, initial=oracq.UniformInput(0.2), type="heat", final_time=0)
        inputs = [oracq.UniformInput(0.2), oracq.ArrayInput([0.2, 0.2]),
                  oracq.ArrayInput([0.2, 0.2], encoding="qram", angle_width=6),
                  oracq.OracleInput(gate_state_prep([1, 1]), math.sqrt(0.08))]
        for initial in inputs:
            with self.subTest(initial=type(initial).__name__):
                instance = oracq.qpde_solve(replace(p, initial=initial), config=oracq.QHAMConfig(order=0, taylor_degree=0)).prepare()
                for value in self.reference_result(instance).physical_values()["u"]:
                    self.assertAlmostEqual(value, 0.2)

    def test_qram_rotation_tree_offsets_and_module_resource_lineage(self):
        grid = oracq.UniformGrid1D(0, 2, 2)
        # The root angle pi/2 prepares a uniform state. Tree address 3 sits in
        # a wider external bank, demonstrating independent memory widths.
        bank = oracq.register_qram({3: 16}, 4, 6, name="initial_angles")
        initial = oracq.QRAMInput(bank, oracq.QRAMInputConfig(1, math.sqrt(0.08), offset=3))
        p = oracq.PDEProblem(grid=grid, initial=initial, type="heat", final_time=0)
        instance = oracq.qpde_solve(p, config=oracq.QHAMConfig(order=0, taylor_degree=0)).prepare()
        circuit = instance.circuit()
        self.assertTrue(circuit.memory)
        self.assertTrue(all(name.endswith("initial_angles") for name in circuit.memory))
        self.assertTrue(all(resource.type == bank.type for resource in circuit.program.main.resources))
        copy = circuit.memory
        next(iter(copy.values()))[3] = 0
        self.assertTrue(all(cells == {3: 16} for cells in circuit.memory.values()))
        self.assertGreater(circuit.resource_estimate().qram_total, 0)
        for value in self.reference_result(instance).physical_values()["u"]:
            self.assertAlmostEqual(value, 0.2)

    def test_unstructured_qram_encoding_is_a_nonsymmetric_discrete_operator(self):
        derivative = (("x", 2),)
        cfg = oracq.UnstructuredGridConfig(2, {derivative: 4}, value_scale=2)
        bank = oracq.register_qram(cfg.encode_matrices({derivative: [[-1, 1], [0.5, -0.5]]}, 8), 4, 8, name="mesh")
        grid = oracq.UnstructuredGrid(bank, cfg)
        be = grid.derivative_encoding(derivative)
        for column in range(2):
            b = oracq.Builder("column", {r.name: r.type for r in be.operation.module.registers}, {"mesh": bank.type})
            if column:
                b.x(b["target"])
            b.call(be.operation, resources={"mesh": "mesh"}, target=b["target"], signal=b["signal"])
            state = oracq.simulate(b.finish().program(), {"mesh": bank.snapshot()})
            for row in range(2):
                self.assertAlmostEqual(state.amplitudes.get((row, 0), 0) * be.alpha, dict(grid.derivative_row(derivative, row))[column])
        p = oracq.PDEProblem(grid=grid, initial=oracq.ArrayInput([0.2, 0.1], encoding="qram", angle_width=8), type="heat", final_time=0.01)
        instance = oracq.qpde_solve(p).prepare()
        self.assertTrue(any(name.endswith("mesh") for name in instance.circuit().memory))
        self.assertTrue(any(not name.endswith("mesh") for name in instance.circuit().memory))
        # Check the actual quantized spatial and initial data against a direct
        # finite-polynomial witness rather than the unquantized input arrays.
        layout = oracq.InitialLayout(grid, ("u",))
        prepared = p.initial.prepare_initial(layout)
        initial = oracq.simulate(prepared.preparation.operation.program(), {b.name: b.snapshot() for b in prepared.qrams})
        values = [initial.amplitudes.get((i, 0), 0) * prepared.norm for i in range(2)]
        expected = [values[row] + 0.001 * sum(weight * values[col] for col, weight in grid.derivative_row(derivative, row)) for row in range(2)]
        for actual, value in zip(self.reference_result(instance).physical_values()["u"], expected, strict=True):
            self.assertAlmostEqual(actual, value)

    def test_negative_angle_endpoint_does_not_wrap_to_positive(self):
        d = (("x", 1),)
        cfg = oracq.UnstructuredGridConfig(2, {d: 0})
        cells = cfg.encode_matrices({d: [[-1, 0], [0, 1]]}, 4)
        self.assertEqual((cells[0], cells[1], cells[3]), (15, 8, 0))
        with self.assertRaises(oracq.ValidationError):
            oracq.UnstructuredGrid(oracq.register_qram([], 2, 4), oracq.UnstructuredGridConfig(2, {d: 1}))
        with self.assertRaises(oracq.ValidationError):
            oracq.UnstructuredGrid(oracq.register_qram([], 4, 4), oracq.UnstructuredGridConfig(2, {d: 0, (("x", 2),): 2}))

    def test_incompatible_oracle_contracts_and_missing_memory_fail_prepare(self):
        p = oracq.PDEProblem(grid=oracq.UniformGrid1D(0, 2, 2), type="heat")
        with self.assertRaises(oracq.ContractError):
            oracq.qpde_solve(replace(p, initial=oracq.OracleInput(gate_state_prep([1, 0, 0, 0]), 1))).prepare()
        with self.assertRaisesRegex(oracq.ValidationError, "unbound"):
            oracq.qpde_solve(replace(p, initial=oracq.OracleInput(abstract_state_prep("initial", 1), 1))).prepare()
        bank = oracq.register_qram([16], 1, 6, name="angles")
        prepared = oracq.QRAMInput(bank, oracq.QRAMInputConfig(1, 1)).prepare_initial(oracq.InitialLayout(p.grid, ("u",)))
        with self.assertRaisesRegex(oracq.ValidationError, "snapshot"):
            oracq.qpde_solve(replace(p, initial=oracq.OracleInput(prepared.preparation, 1))).prepare()
        with self.assertRaises(oracq.ValidationError):
            oracq.qpde_solve(p, method="unknown")

    def test_field_major_layout_and_zero_initial_solution(self):
        u, v = Field("u"), Field("v")
        equation = PolynomialPDE.from_equations({"u": -u, "v": u - v})
        grid = oracq.UniformGrid1D(0, 1, 3)
        layout = oracq.InitialLayout(grid, equation.fields)
        self.assertEqual(layout.encode({"u": [1, 2, 3], "v": [4, 5, 6]}), (1, 2, 3, 0, 4, 5, 6, 0))
        p = oracq.PDEProblem(grid=grid, initial=oracq.UniformInput(0))
        result = self.reference_result(oracq.qpde_solve(p).prepare())
        self.assertEqual(result.physical_values(), {"u": (0j, 0j, 0j)})
        with self.assertRaises(oracq.ValidationError):
            result.physical_amplitudes()

    def test_custom_method_owns_contract_and_reuses_instance_lifecycle(self):
        class CustomMethod:
            calls = 0

            def prepare(self, problem):
                self.calls += 1
                # A custom method may restrict inputs independently of QHAM.
                if not isinstance(problem.initial, oracq.UniformInput):
                    raise oracq.ValidationError("custom method requires uniform input")
                return oracq.QHAMMethod(oracq.QHAMConfig(order=0, taylor_degree=0)).prepare(problem)

        method = CustomMethod()
        p = oracq.PDEProblem(grid=oracq.UniformGrid1D(0, 2, 2))
        instance = oracq.qpde_solve(p, method=method).prepare()
        instance.prepare()
        self.assertEqual(method.calls, 1)
        with self.assertRaises(oracq.ValidationError):
            oracq.qpde_solve(replace(p, initial=oracq.ArrayInput([1, 1])), method=method).prepare()

    def test_external_linear_solver_keeps_magnitude_recovery_explicit(self):
        from functools import partial

        from oracq.algorithms.input_model.qham import taylor_qode

        p = oracq.PDEProblem(grid=oracq.UniformGrid1D(0, 2, 2), type="heat")
        cfg = oracq.QHAMConfig(order=0, linear_solver=partial(taylor_qode, degree=0))
        result = self.reference_result(oracq.qpde_solve(p, config=cfg).prepare())
        self.assertGreater(result.success_probability, 0)
        self.assertAlmostEqual(sum(abs(v) ** 2 for v in result.physical_amplitudes()["u"]), 1)
        with self.assertRaisesRegex(oracq.ValidationError, "magnitude recovery"):
            result.physical_values()

    def test_resource_basis_conversion_preserves_symbolic_repeat(self):
        b = oracq.Builder("cost", {"target": oracq.Bits(3)})
        with b.repeat(10**12), b.control(b["target"][:2]):
            b.x(b["target"][2])
        circuit = oracq.PDECircuit(b.finish().program())
        native = circuit.resource_estimate("toffoli+clifford+t+qram")
        ct = circuit.resource_estimate("clifford+t+qram")
        self.assertEqual((native.toffoli, ct.toffoli, ct.t_exact, ct.clifford), (10**12, 0, 7 * 10**12, 8 * 10**12))
        self.assertIn("Repeat", circuit.dumps())
        self.assertLess(len(circuit.dumps()), 2000)
        with self.assertRaises(oracq.ValidationError):
            circuit.resource_estimate("anything")

    def test_prepare_and_export_do_not_import_native_backends(self):
        subprocess.run([sys.executable, "-c", """
import sys
import oracq
p = oracq.PDEProblem(grid=oracq.UniformGrid1D(0, 2, 2), type='heat')
s = oracq.qpde_solve(p, config=oracq.QHAMConfig(order=0, taylor_degree=0)).prepare()
s.circuit().originir_ext()
s.circuit().resource_estimate()
assert 'pysparq' not in sys.modules
assert 'uniqc' not in sys.modules
"""], check=True, capture_output=True)
