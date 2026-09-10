"""SDE/Fokker–Planck 输入模型的数值见证与求解器契约测试。"""

import math
import unittest
from functools import partial

from pyqecclang import dumps, loads, simulate
from pyqecclang.algorithms.hamiltonian import taylor_hamiltonian
from pyqecclang.algorithms.lchs import QuadraturePlan
from pyqecclang.algorithms.ode import linear_qode
from pyqecclang.algorithms.ode_models import LinearODE
from pyqecclang.algorithms.oracles import StatePreparation
from pyqecclang.algorithms.sde import (
    FokkerPlanckProblem,
    boltzmann_distribution,
    distribution_moments,
    evolve_distribution,
    matrix_exponential,
    sde_state_angles,
    sde_state_preparation,
    stationary_distribution,
)
from pyqecclang.infrastructure.ir import ValidationError


def ou_problem(size=16, half_width=3.0, theta=1.0, diffusion=0.5):
    """Ornstein–Uhlenbeck 过程 a(x) = -theta*x、D 常量的小网格实例。"""
    h = 2 * half_width / size
    points = [-half_width + (i + 0.5) * h for i in range(size)]
    return FokkerPlanckProblem([-theta * x for x in points], [diffusion] * size, points)


def gaussian_distribution(points, mean, variance):
    weights = [math.exp(-((x - mean) ** 2) / (2 * variance)) for x in points]
    total = math.fsum(weights)
    return [w / total for w in weights]


class GeneratorDiscretizationTests(unittest.TestCase):
    def test_column_sums_vanish(self):
        g = ou_problem(size=8).generator_matrix()
        for column in range(8):
            self.assertAlmostEqual(sum(g[row][column] for row in range(8)), 0.0, places=12)

    def test_evolution_conserves_probability(self):
        problem = ou_problem(size=8)
        g = problem.generator_matrix()
        initial = gaussian_distribution(problem.points, 0.3, 0.2)
        exact = evolve_distribution(g, initial, 0.7)
        self.assertAlmostEqual(math.fsum(exact), 1.0, places=9)
        euler = evolve_distribution(g, initial, 0.7, steps=20000)
        self.assertAlmostEqual(math.fsum(euler), 1.0, places=9)
        for left, right in zip(exact, euler, strict=True):
            self.assertAlmostEqual(left, right, delta=1e-4)

    def test_matrix_exponential_matches_euler_limit(self):
        problem = ou_problem(size=8)
        g = problem.generator_matrix()
        initial = gaussian_distribution(problem.points, -0.2, 0.3)
        propagated = matrix_exponential(g, 0.25)
        exact = [
            math.fsum(row[j] * initial[j] for j in range(len(initial))) for row in propagated
        ]
        coarse = evolve_distribution(g, initial, 0.25, steps=2000)
        fine = evolve_distribution(g, initial, 0.25, steps=4000)
        # 显式 Euler 一阶：细步误差约为粗步一半。
        err_coarse = max(abs(a - b) for a, b in zip(coarse, exact, strict=True))
        err_fine = max(abs(a - b) for a, b in zip(fine, exact, strict=True))
        self.assertLess(err_fine, 0.6 * err_coarse)

    def test_ou_moments_match_closed_form(self):
        theta, diffusion, time = 1.0, 0.5, 0.4
        problem = ou_problem(size=16, half_width=3.0, theta=theta, diffusion=diffusion)
        mean0, var0 = 0.4, 0.16
        initial = gaussian_distribution(problem.points, mean0, var0)
        evolved = evolve_distribution(problem.generator_matrix(), initial, time)
        mean, second = distribution_moments(problem.points, evolved)
        variance = second - mean * mean
        self.assertAlmostEqual(mean, mean0 * math.exp(-theta * time), delta=2e-3)
        expected_var = var0 * math.exp(-2 * theta * time) + diffusion / theta * (
            1 - math.exp(-2 * theta * time)
        )
        self.assertAlmostEqual(variance, expected_var, delta=2e-2)

    def test_stationary_approaches_boltzmann(self):
        problem = ou_problem(size=16, half_width=2.0, theta=1.0, diffusion=1.0)
        uniform = [1.0 / problem.size] * problem.size
        relaxed = evolve_distribution(problem.generator_matrix(), uniform, 30.0)
        stationary = stationary_distribution(problem)
        for left, right in zip(relaxed, stationary, strict=True):
            self.assertAlmostEqual(left, right, delta=1e-6)
        boltzmann = boltzmann_distribution(problem)
        for left, right in zip(stationary, boltzmann, strict=True):
            self.assertAlmostEqual(left, right, delta=2e-2)


class StatePreparationTests(unittest.TestCase):
    def test_gate_preparation_amplitudes_are_root_probabilities(self):
        probabilities = [0.1, 0.2, 0.3, 0.4]
        preparation = sde_state_preparation(probabilities)
        self.assertIsInstance(preparation, StatePreparation)
        self.assertEqual(preparation.width, 2)
        state = simulate(preparation.operation.program())
        for index, p in enumerate(probabilities):
            self.assertAlmostEqual(
                abs(state.amplitudes.get((index, 0), 0)) ** 2, p, places=10
            )

    def test_qram_preparation_declares_angle_table(self):
        probabilities = [0.25, 0.25, 0.25, 0.25]
        preparation = sde_state_preparation(probabilities, implementation="qram", angle_width=6)
        self.assertEqual(preparation.width, 2)
        resources = {r.name: r.type for r in preparation.operation.module.resources}
        self.assertIn("angles", resources)
        angles = sde_state_angles(probabilities, angle_width=6)
        self.assertTrue(angles)
        self.assertTrue(all(0 <= v < 64 for v in angles.values()))

    def test_invalid_probabilities_rejected(self):
        with self.assertRaises(ValidationError):
            sde_state_preparation([0.5, -0.1, 0.3, 0.3])
        with self.assertRaises(ValidationError):
            sde_state_preparation([0.0, 0.0, 0.0, 0.0])
        with self.assertRaises(ValidationError):
            sde_state_preparation([0.5, 0.25, 0.25])
        with self.assertRaises(ValidationError):
            sde_state_preparation([0.5, 0.5], implementation="sparse")


class SolverContractTests(unittest.TestCase):
    def test_qode_problem_accepted_by_lchs(self):
        problem = ou_problem(size=4, half_width=2.0, theta=1.0, diffusion=1.0)
        initial = sde_state_preparation(boltzmann_distribution(problem))
        qode = problem.qode_problem(initial=initial)
        self.assertTrue(qode.dissipative)
        self.assertEqual(qode.generator.width, 2)
        protocol = linear_qode(
            "lchs",
            hamiltonian_function=partial(taylor_hamiltonian, degree=1),
            plan=QuadraturePlan.cauchy(cutoff=0),
        )
        protocol.check(qode, time=0.1).require()
        result = protocol.solve(qode, 0.1)
        self.assertEqual(result.width, 2)
        attributes = dict(result.operation.module.attributes)
        self.assertTrue(attributes["qode_dissipative_promise"])
        self.assertEqual(loads(dumps(result.operation.program())), result.operation.program())

    def test_linear_ode_view_matches_contract(self):
        problem = ou_problem(size=4, half_width=2.0)
        model = problem.linear_ode()
        self.assertIsInstance(model, LinearODE)
        self.assertEqual(model.parts.hermitian.width, model.initial.width)
        self.assertEqual(model.parts.hermitian.width, model.parts.h.width)

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValidationError):
            FokkerPlanckProblem([0.0] * 4, [1.0, -0.5, 1.0, 1.0], [0.0, 1.0, 2.0, 3.0])
        with self.assertRaises(ValidationError):
            FokkerPlanckProblem([0.0] * 4, [1.0] * 4, [0.0, 1.0, 1.5, 2.5])
        with self.assertRaises(ValidationError):
            FokkerPlanckProblem([0.0] * 3, [1.0] * 4, [0.0, 1.0, 2.0, 3.0])
        with self.assertRaises(ValidationError):
            FokkerPlanckProblem([0.0] * 4, [1.0] * 4, [0.0, 1.0, 0.5, 3.0])
        odd = FokkerPlanckProblem([0.0] * 3, [1.0] * 3, [0.0, 1.0, 2.0])
        with self.assertRaises(ValidationError):
            odd.generator_encoding()
        problem = ou_problem(size=4)
        with self.assertRaises(ValidationError):
            problem.qode_problem(initial=sde_state_preparation([0.5, 0.5]))


if __name__ == "__main__":
    unittest.main()
