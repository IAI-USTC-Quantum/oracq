"""ODE physical witnesses, interchangeable inputs, method contracts, and lifecycle."""

import math
import subprocess
import sys
import unittest
from dataclasses import replace

import oracq
from oracq.algorithms.input_model.oracles import abstract_state_prep, gate_state_prep
from oracq.algorithms.qode.cbmd import ContourPlan
from oracq.algorithms.qode.lchs import QuadraturePlan
from oracq.algorithms.qode.ode import linear_qode
from oracq.algorithms.qode.schrodingerization import SchrodingerPlan


class ODEInstanceTests(unittest.TestCase):
    def reference_result(self, instance):
        circuit = instance.circuit()
        return oracq.ODEResult(oracq.simulate(circuit.program, circuit.memory),
                               oracq.ODELayout(instance.problem.size, instance.problem.components),
                               "reference", instance.prepared.amplitude_scale)

    def test_taylor_matches_independent_coupled_matrix_polynomial(self):
        matrix, values, time = [[-1, 2], [0, -2]], [2, -1], 0.1
        p = oracq.ODEProblem(oracq.MatrixInput(matrix), oracq.ArrayInput(values), time)
        instance = oracq.qode_solve(p)
        for action in (instance.circuit, instance.run_pysparq, instance.run_originir_ext):
            with self.assertRaisesRegex(oracq.ValidationError, "prepare"):
                action()
        instance.prepare()
        circuit = instance.circuit()
        instance.prepare()
        self.assertIs(circuit, instance.circuit())
        self.assertTrue(instance.prepared.report.ok)
        self.assertEqual(oracq.loads(circuit.dumps()), circuit.program)
        self.assertIn("DEF ", circuit.originir_ext().text)
        self.assertTrue(circuit.resource_estimate().complete)
        first = [sum(matrix[i][j] * values[j] for j in range(2)) for i in range(2)]
        second = [sum(matrix[i][j] * first[j] for j in range(2)) for i in range(2)]
        result = self.reference_result(instance)
        expected = [values[i] + time * first[i] + time**2 / 2 * second[i] for i in range(2)]
        for actual, value in zip(result.physical_values(), expected, strict=True):
            self.assertAlmostEqual(actual, value)
        self.assertAlmostEqual(sum(abs(v)**2 for v in result.physical_amplitudes()), 1)
        self.assertAlmostEqual(result.success_probability, sum(abs(v)**2 for v in expected) / result.amplitude_scale**2)

    def test_scalar_complex_and_three_component_padding(self):
        for matrix, values in [([[1j]], [2 - 1j]), ([[-1, 0, 0], [0, 2, 0], [0, 0, 1j]], [1, -2, 3j])]:
            with self.subTest(size=len(values)):
                p = oracq.ODEProblem(oracq.MatrixInput(matrix), oracq.ArrayInput(values), 0.1)
                instance = oracq.qode_solve(p, config=oracq.ODEConfig(taylor_degree=1)).prepare()
                result = self.reference_result(instance)
                for i, actual in enumerate(result.physical_values()):
                    self.assertAlmostEqual(actual, values[i] * (1 + 0.1 * matrix[i][i]))
                self.assertEqual(len(result.physical_values()), len(values))

    def test_generator_padding_cannot_couple_through_unphysical_components(self):
        # X couples a scalar to its unused address; PGP is the zero generator.
        p = oracq.ODEProblem(oracq.OracleGeneratorInput(oracq.pauli_x(1), 1), oracq.ArrayInput([3]), 0.2)
        result = self.reference_result(oracq.qode_solve(p).prepare())
        self.assertAlmostEqual(result.physical_values()[0], 3)

    def test_large_oracle_input_does_not_materialize_component_names_or_samples(self):
        p = oracq.ODEProblem(oracq.OracleGeneratorInput(oracq.identity(40), 1 << 40))
        instance = oracq.qode_solve(p, config=oracq.ODEConfig(taylor_degree=0)).prepare()
        self.assertEqual(p.components, ())
        self.assertEqual(instance.prepared.state.width, 40)
        self.assertLess(len(instance.circuit().dumps()), 10000)

    def test_matrix_and_initial_sources_snapshot_physical_data(self):
        matrix, values = [[-1, 0], [0, -2]], [1, 2]
        p = oracq.ODEProblem(oracq.MatrixInput(matrix), oracq.ArrayInput(values), 0)
        matrix[0][0], values[0] = 99, 99
        result = self.reference_result(oracq.qode_solve(p).prepare())
        self.assertEqual(p.generator.values[0][0], -1)
        self.assertAlmostEqual(result.physical_values()[0], 1)
        layout = oracq.ODELayout(3, ("x", "y", "z"))
        self.assertEqual(layout.encode({"x": [1], "y": [2], "z": [3]}), (1, 2, 3, 0))

    def test_uniform_qram_and_external_initial_sources_are_interchangeable(self):
        p = oracq.ODEProblem(oracq.MatrixInput([[-1, 0], [0, -1]]), final_time=0.1)
        bank = oracq.register_qram({3: 2}, 4, 3, name="initial")
        sources = [oracq.UniformInput(2), oracq.ArrayInput([2, 2], encoding="qram", angle_width=3),
                   oracq.QRAMInput(bank, oracq.QRAMInputConfig(1, math.sqrt(8), offset=3)),
                   oracq.OracleInput(gate_state_prep([1, 1]), math.sqrt(8))]
        for initial in sources:
            with self.subTest(source=type(initial).__name__):
                instance = oracq.qode_solve(replace(p, initial=initial)).prepare()
                for value in self.reference_result(instance).physical_values():
                    self.assertAlmostEqual(value, 2 * (1 - 0.1 + 0.1**2 / 2))

    def test_quantized_nonsymmetric_qram_generator_and_resource_aliases(self):
        config = oracq.QRAMMatrixConfig(2, value_scale=2, offset=4)
        cells = config.encode_matrix([[-1, 1], [0.5, -2]], 4)
        bank = oracq.register_qram(cells, 4, 4, name="operator")
        p = oracq.ODEProblem(oracq.QRAMMatrixInput(bank, config), oracq.ArrayInput([1, 2], encoding="qram", angle_width=3), 0.1)
        instance = oracq.qode_solve(p, config=oracq.ODEConfig(taylor_degree=1)).prepare()
        layout = oracq.ODELayout(2)
        initial = p.initial.prepare_initial(layout)
        state = oracq.simulate(initial.preparation.operation.program(), {b.name: b.snapshot() for b in initial.qrams})
        values = [state.amplitudes.get((i, 0), 0) * initial.norm for i in range(2)]
        matrix = [[2 * math.cos(math.pi * cells[4 + 2 * r + c] / 16) for c in range(2)] for r in range(2)]
        expected = [values[r] + 0.1 * sum(matrix[r][c] * values[c] for c in range(2)) for r in range(2)]
        result = self.reference_result(instance)
        for actual, value in zip(result.physical_values(), expected, strict=True):
            self.assertAlmostEqual(actual, value)
        self.assertTrue(any(name.endswith("operator") for name in instance.circuit().memory))
        self.assertGreater(instance.circuit().resource_estimate().qram_total, 0)
        memory = instance.circuit().memory
        next(iter(memory.values())).clear()
        self.assertTrue(all(instance.circuit().memory.values()))

    def test_qram_matrix_array_adapter_padding_and_signed_endpoint(self):
        matrix = [[-1, 1, 0], [0, 0.5, 0], [0, 0, 1]]
        p = oracq.ODEProblem(oracq.MatrixInput(matrix, encoding="qram", angle_width=3), oracq.ArrayInput([1, 2, 3]), 0.1)
        instance = oracq.qode_solve(p, config=oracq.ODEConfig(taylor_degree=0)).prepare()
        for actual, value in zip(self.reference_result(instance).physical_values(), [1, 2, 3], strict=True):
            self.assertAlmostEqual(actual, value)
        config = oracq.QRAMMatrixConfig(3)
        cells = config.encode_matrix(matrix, 3)
        self.assertEqual((cells[0], cells[1], cells[3], cells[15]), (7, 0, 4, 4))

    def test_methods_enforce_dissipativity_and_own_plan_types(self):
        p = oracq.ODEProblem(oracq.MatrixInput([[-1, 0], [0, -1]]), final_time=0)
        for name, plan in [("LCHS", QuadraturePlan((0,), (1,))), ("CBMD", ContourPlan(cutoff=0)),
                           ("schrodingerization", SchrodingerPlan(auxiliary_width=1, selected_index=0))]:
            with self.subTest(method=name):
                cfg = oracq.ODEConfig(plan=plan)
                if name != "schrodingerization":
                    with self.assertRaises(oracq.ContractError):
                        oracq.qode_solve(p, method=name, config=cfg).prepare()
                instance = oracq.qode_solve(replace(p, dissipative=True), method=name, config=cfg).prepare()
                self.assertTrue(instance.prepared.report.ok)
                self.assertGreater(self.reference_result(instance).success_probability, 0)
                with self.assertRaisesRegex(oracq.ValidationError, "magnitude recovery"):
                    self.reference_result(instance).physical_values()
        with self.assertRaises(oracq.ValidationError):
            oracq.qode_solve(p, method="lchs", config=oracq.ODEConfig(plan=ContourPlan()))
        with self.assertRaises(oracq.ValidationError):
            oracq.qode_solve(p, config=oracq.ODEConfig(plan=QuadraturePlan((0,), (1,))))

    def test_external_protocol_and_custom_method_keep_lifecycle_and_scale_explicit(self):
        p = oracq.ODEProblem(oracq.MatrixInput([[-1, 0], [0, -1]]), final_time=0, dissipative=True)
        external = linear_qode("lchs", plan=QuadraturePlan((0,), (1,)))
        instance = oracq.qode_solve(p, method=external).prepare()
        self.assertIsNone(instance.prepared.amplitude_scale)

        class CustomMethod:
            calls = 0

            def prepare(self, problem):
                self.calls += 1
                return oracq.qode_solve(problem).prepare().prepared

        method = CustomMethod()
        instance = oracq.qode_solve(p, method=method).prepare()
        instance.prepare()
        self.assertEqual(method.calls, 1)
        with self.assertRaises(oracq.ValidationError):
            oracq.qode_solve(p, method=external, config=oracq.ODEConfig())

    def test_bad_inputs_contracts_and_missing_snapshots_fail_before_execution(self):
        p = oracq.ODEProblem(oracq.MatrixInput([[-1, 0], [0, -1]]))
        for initial in [oracq.ArrayInput([1]), oracq.OracleInput(gate_state_prep([1, 0, 0, 0]), 1),
                        oracq.OracleInput(abstract_state_prep("unbound", 1), 1)]:
            with self.subTest(initial=type(initial).__name__), self.assertRaises(oracq.ValidationError):
                oracq.qode_solve(replace(p, initial=initial)).prepare()
        encoded = oracq.MatrixInput([[-1, 0], [0, -1]], encoding="qram", angle_width=3).prepare_generator(oracq.ODELayout(2))
        with self.assertRaisesRegex(oracq.ValidationError, "snapshot"):
            oracq.qode_solve(replace(p, generator=oracq.OracleGeneratorInput(encoded.encoding, 2))).prepare()
        for factory in [lambda: oracq.MatrixInput([[1, 2]]), lambda: oracq.MatrixInput([[float("nan")]]),
                        lambda: replace(p, final_time=-1), lambda: oracq.ODELayout(2, ("x", "x")),
                        lambda: oracq.QRAMMatrixInput(oracq.register_qram([], 1, 4), oracq.QRAMMatrixConfig(2)),
                        lambda: oracq.qode_solve(p, method="unknown")]:
            with self.assertRaises(oracq.ValidationError):
                factory()

    def test_zero_solution_has_no_conditional_state(self):
        p = oracq.ODEProblem(oracq.MatrixInput([[-1, 0], [0, -1]]), oracq.UniformInput(0))
        result = self.reference_result(oracq.qode_solve(p).prepare())
        self.assertEqual(result.physical_values(), (0j, 0j))
        self.assertEqual(result.success_probability, 0)
        with self.assertRaises(oracq.ValidationError):
            result.physical_amplitudes()

    def test_shared_resource_basis_adds_costs_and_preserves_large_repeats(self):
        b = oracq.Builder("mixed_cost", {"target": oracq.Bits(3)})
        b.h(b["target"][0])
        with b.repeat(10**12), b.control(b["target"][:2]):
            b.x(b["target"][2])
        circuit = oracq.ODECircuit(b.finish().program())
        cost = circuit.resource_estimate()
        self.assertEqual(cost.clifford, 8 * 10**12 + 1)
        self.assertEqual(cost.t_exact, 7 * 10**12)
        self.assertIn("Repeat", circuit.dumps())
        self.assertLess(len(circuit.dumps()), 2000)

    def test_prepare_export_and_estimate_import_no_native_backend(self):
        subprocess.run([sys.executable, "-c", """
import sys
import oracq
p = oracq.ODEProblem(oracq.MatrixInput([[-1, 0], [0, -1]]))
s = oracq.qode_solve(p).prepare()
s.circuit().originir_ext()
s.circuit().resource_estimate()
assert 'uniqc' not in sys.modules
assert 'pysparq' not in sys.modules
"""], check=True, capture_output=True)
