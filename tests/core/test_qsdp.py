"""Numerical witnesses for the QSDP framework: trace-estimation circuits, MMW driver convergence, and iterative circuit generation."""

import math
import unittest

from oracq import ValidationError, simulate, unresolved
from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding
from oracq.algorithms.input_model.density import gate_purification
from oracq.algorithms.qml.qsdp import (
    SdpInstance,
    classical_estimator,
    iteration_circuits,
    penalty_hamiltonian,
    qsdp_gibbs_solve,
    trace_estimate_circuit,
    trace_from_probe,
)

Z = ((1.0 + 0j, 0j), (0j, -1.0 + 0j))
X = ((0j, 1.0 + 0j), (1.0 + 0j, 0j))
NEG_Z = tuple(tuple(-v for v in row) for row in Z)
NEG_X = tuple(tuple(-v for v in row) for row in X)


def probe_probability(operation):
    state = simulate(operation.program())
    return sum(
        abs(a) ** 2 for k, a in state.amplitudes.items() if k[-1] == 1
    )


class TraceEstimateTests(unittest.TestCase):
    def test_pure_state_observables(self):
        # ρ = |0⟩⟨0|: Tr(Zρ) = 1, Tr(Xρ) = 0.
        purification = gate_purification(((1.0 + 0j, 0j), (0j, 0j)))
        circuit = trace_estimate_circuit(purification, matrix_pauli_encoding(Z))
        attrs = dict(circuit.module.attributes)
        self.assertEqual(attrs["algorithm"], "trace_estimate")
        alpha = attrs["be_alpha"]
        self.assertAlmostEqual(trace_from_probe(probe_probability(circuit), alpha), 1.0, places=9)
        circuit = trace_estimate_circuit(purification, matrix_pauli_encoding(X))
        self.assertAlmostEqual(
            trace_from_probe(probe_probability(circuit), attrs["be_alpha"]), 0.0, places=9
        )

    def test_mixed_state_observable(self):
        # ρ = diag(0.75, 0.25): Tr(Zρ) = 0.5.
        purification = gate_purification(((0.75 + 0j, 0j), (0j, 0.25 + 0j)))
        be = matrix_pauli_encoding(Z)
        circuit = trace_estimate_circuit(purification, be)
        self.assertAlmostEqual(
            trace_from_probe(probe_probability(circuit), be.alpha), 0.5, places=9
        )


class DriverTests(unittest.TestCase):
    def test_penalty_hamiltonian(self):
        instance = SdpInstance(((Z, 0.2), (X, 0.1)))
        hamiltonian = penalty_hamiltonian(instance, (1.0, 0.0))
        expected = ((0.8 + 0j, 0j), (0j, -1.2 + 0j))
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(hamiltonian[i][j], expected[i][j], places=12)

    def test_equality_via_doubled_constraints(self):
        # The equalities r_z = 0.2 and r_x = 0.1 split into two inequality pairs each; MMW should converge within ε.
        instance = SdpInstance(((Z, 0.2), (NEG_Z, -0.2), (X, 0.1), (NEG_X, -0.1)))
        result = qsdp_gibbs_solve(instance, epsilon=0.08)
        self.assertTrue(result["converged"])
        for violation in result["violations"]:
            self.assertLessEqual(violation, 0.08)
        rho = result["rho"]
        self.assertAlmostEqual(rho[0][0] + rho[1][1], 1.0, places=9)
        # The Bloch vector approaches the target (0.1, 0, 0.2).
        r_z = (rho[0][0] - rho[1][1]).real
        r_x = (2 * rho[0][1]).real
        self.assertAlmostEqual(r_z, 0.2, delta=0.12)
        self.assertAlmostEqual(r_x, 0.1, delta=0.12)

    def test_infeasible_instance_does_not_converge(self):
        # r_z ≤ 0.2 and r_z ≥ 0.5 are contradictory.
        instance = SdpInstance(((Z, 0.2), (NEG_Z, -0.5)))
        result = qsdp_gibbs_solve(instance, epsilon=0.05, max_iterations=4000)
        self.assertFalse(result["converged"])

    def test_classical_estimator_matches_gibbs(self):
        instance = SdpInstance(((Z, 0.0),))
        rho, estimates = classical_estimator(Z, 1.0, instance)
        z = math.exp(1.0) + math.exp(-1.0)
        self.assertAlmostEqual(estimates[0], (math.exp(-1.0) - math.exp(1.0)) / z, places=10)
        self.assertAlmostEqual(rho[0][0] + rho[1][1], 1.0, places=12)


class IterationCircuitsTests(unittest.TestCase):
    def test_generation_resolves(self):
        instance = SdpInstance(((Z, 0.2), (NEG_Z, -0.2)))
        purification, traces = iteration_circuits(
            instance, (0.6, 0.4), 0.5, hamiltonian_encoding=matrix_pauli_encoding, error=0.1
        )
        self.assertEqual(purification.attributes["algorithm"], "gibbs_purification")
        self.assertEqual(len(traces), 2)
        for circuit in traces:
            self.assertFalse(unresolved(circuit.program()))
            names = [r.name for r in circuit.module.registers]
            self.assertEqual(
                names, ["system", "environment", "signal", "probe", "purification_signal"]
            )

    def test_invalid_inputs_fail_at_generation(self):
        purification = gate_purification(((1.0 + 0j, 0j), (0j, 0j)))
        instance = SdpInstance(((Z, 0.2),))
        cases = [
            lambda: SdpInstance(()),
            lambda: SdpInstance((((1.0, 0.0), (0.0, 0.1)), 0.0)),  # not Hermitian
            lambda: penalty_hamiltonian(instance, (0.5, 0.5)),
            lambda: trace_estimate_circuit(purification, matrix_pauli_encoding(Z), component="abs"),
            lambda: trace_estimate_circuit(object(), matrix_pauli_encoding(Z)),
            lambda: trace_from_probe(1.5, 1.0),
            lambda: qsdp_gibbs_solve(instance, epsilon=0.0),
            lambda: iteration_circuits(instance, (1.0,), 0.0, hamiltonian_encoding=matrix_pauli_encoding),
        ]
        for case in cases:
            with self.subTest(case=case):
                with self.assertRaises(ValidationError):
                    case()


if __name__ == "__main__":
    unittest.main()
