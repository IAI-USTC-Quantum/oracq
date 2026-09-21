"""DM input model 与 Gibbs 态制备的数值见证。"""

import math
import unittest

from pyqecclang import ValidationError, bind, simulate, unresolved
from pyqecclang.algorithms.input_model.block_encoding import matrix_pauli_encoding
from pyqecclang.algorithms.input_model.density import (
    PurificationAccess,
    abstract_purification,
    gate_purification,
    gibbs_purification,
    gibbs_state,
    maximally_mixed_purification,
    partial_trace,
    trace_distance,
)
from pyqecclang.algorithms.input_model.oracles import gate_state_prep


def dense_state(state, width, *, signal_width=0):
    """把 simulate 的稀疏幅度字典展开成稠密态向量（丢弃 signal != 0 分支）。

    基态下标约定为 system | (environment << system_width) | (signal << (2*width))。
    """
    size = 1 << (2 * width)
    vector = [0j] * size
    for key, amplitude in state.amplitudes.items():
        system, environment = key[0], key[1]
        if signal_width and key[2] != 0:
            continue
        vector[system + (environment << width)] += amplitude
    return vector


def assert_matrix_close(case, actual, expected, places=7):
    case.assertEqual(len(actual), len(expected))
    for row_a, row_b in zip(actual, expected, strict=True):
        for a, b in zip(row_a, row_b, strict=True):
            case.assertAlmostEqual(a, b, places=places)


def renormalized(reduced):
    trace = sum(reduced[i][i] for i in range(len(reduced)))
    return tuple(
        tuple(value / trace for value in row) for row in reduced
    ), trace


RHO = ((0.7 + 0j, 0.1 - 0.05j), (0.1 + 0.05j, 0.3 + 0j))


class PurificationTests(unittest.TestCase):
    def test_gate_purification_recovers_rho(self):
        access = gate_purification(RHO)
        self.assertEqual(access.width, 1)
        vector = dense_state(simulate(access.operation.program()), 1)
        reduced = partial_trace(vector, 1, 1)
        assert_matrix_close(self, reduced, RHO)

    def test_maximally_mixed_purification(self):
        access = maximally_mixed_purification(2)
        vector = dense_state(simulate(access.operation.program()), 2)
        reduced = partial_trace(vector, 2, 2)
        expected = tuple(
            tuple((0.25 if i == j else 0.0) + 0j for j in range(4)) for i in range(4)
        )
        assert_matrix_close(self, reduced, expected)
        # Schmidt 结构：仅 system == environment 的基态有幅度。
        for key, amplitude in simulate(access.operation.program()).amplitudes.items():
            if abs(amplitude) > 1e-12:
                self.assertEqual(key[0], key[1])

    def test_pure_state_adapter(self):
        preparation = gate_state_prep([math.sqrt(0.3), math.sqrt(0.7)])
        access = PurificationAccess.from_state_preparation(preparation)
        vector = dense_state(simulate(access.operation.program()), 1)
        reduced = partial_trace(vector, 1, 1)
        expected = ((0.3 + 0j, math.sqrt(0.21) + 0j), (math.sqrt(0.21) + 0j, 0.7 + 0j))
        assert_matrix_close(self, reduced, expected)

    def test_abstract_purification_binds(self):
        abstract = abstract_purification("RhoAccess", 1)
        self.assertEqual(
            [r.name for r in unresolved(abstract.operation.program())], ["RhoAccess"]
        )
        bound = bind(
            abstract.operation.program(),
            {"RhoAccess": gate_purification(RHO).operation},
        )
        self.assertFalse(unresolved(bound))
        vector = dense_state(simulate(bound), 1)
        reduced = partial_trace(vector, 1, 1)
        assert_matrix_close(self, reduced, RHO)

    def test_classical_tools(self):
        # 对角 Gibbs 参考与迹距离的解析值。
        hamiltonian = ((1.0 + 0j, 0j), (0j, -1.0 + 0j))
        rho = gibbs_state(hamiltonian, 2.0)
        z = math.exp(-2.0) + math.exp(2.0)
        assert_matrix_close(
            self,
            rho,
            ((math.exp(-2.0) / z + 0j, 0j), (0j, math.exp(2.0) / z + 0j)),
        )
        pure = ((1.0 + 0j, 0j), (0j, 0j))
        mixed = ((0.5 + 0j, 0j), (0j, 0.5 + 0j))
        self.assertAlmostEqual(trace_distance(pure, mixed), 0.5, places=12)
        self.assertAlmostEqual(trace_distance(pure, pure), 0.0, places=12)


class GibbsTests(unittest.TestCase):
    def test_beta_zero_is_maximally_mixed(self):
        be = matrix_pauli_encoding(((1.0, 0.0), (0.0, -1.0)))
        access = gibbs_purification(be, 0.0)
        attrs = access.attributes
        self.assertEqual(attrs["algorithm"], "gibbs_purification")
        self.assertEqual(attrs["beta"], 0.0)
        vector = dense_state(simulate(access.operation.program()), 1)
        reduced = partial_trace(vector, 1, 1)
        assert_matrix_close(
            self, reduced, ((0.5 + 0j, 0j), (0j, 0.5 + 0j))
        )

    def test_diagonal_hamiltonian_gibbs(self):
        hamiltonian = ((1.0 + 0j, 0j), (0j, -0.5 + 0j))
        beta = 0.6
        be = matrix_pauli_encoding(hamiltonian)
        access = gibbs_purification(be, beta, error=0.05)
        program = access.operation.program()
        vector = dense_state(simulate(program), 1, signal_width=access.signal_qubits)
        reduced, _ = renormalized(partial_trace(vector, 1, 1))
        expected = gibbs_state(hamiltonian, beta)
        distance = trace_distance(reduced, expected)
        self.assertLess(distance, 0.1)

    def test_non_diagonal_hamiltonian_gibbs(self):
        hamiltonian = ((0.3 + 0j, 0.2 - 0.1j), (0.2 + 0.1j, -0.4 + 0j))
        beta = 0.8
        be = matrix_pauli_encoding(hamiltonian)
        access = gibbs_purification(be, beta, error=0.05)
        vector = dense_state(
            simulate(access.operation.program()), 1, signal_width=access.signal_qubits
        )
        reduced, _ = renormalized(partial_trace(vector, 1, 1))
        expected = gibbs_state(hamiltonian, beta)
        self.assertLess(trace_distance(reduced, expected), 0.1)

    def test_error_convergence_decreases(self):
        # 实测口径：error 控制多项式一致截断误差，是迹距离的上界但远不紧——
        # 各档实测距离 5.6e-4 / 8.7e-5 / 8.7e-5（error=0.2 与 0.1 落到同一
        # 截断度数，距离相同），因此断言单调不增而非严格下降，并要求每档 ≤ error。
        hamiltonian = ((1.0 + 0j, 0j), (0j, -0.5 + 0j))
        beta = 0.8
        errors = (0.4, 0.2, 0.1)
        distances = []
        for error in errors:
            be = matrix_pauli_encoding(hamiltonian)
            access = gibbs_purification(be, beta, error=error)
            vector = dense_state(
                simulate(access.operation.program()), 1, signal_width=access.signal_qubits
            )
            reduced, _ = renormalized(partial_trace(vector, 1, 1))
            distances.append(trace_distance(reduced, gibbs_state(hamiltonian, beta)))
        for error, distance in zip(errors, distances, strict=True):
            self.assertLessEqual(
                distance, error, msg=f"error={error} 档迹距离 {distance} 超过上界；各档 {distances}"
            )
        self.assertLessEqual(distances[1], distances[0], msg=f"各档距离 {distances}")
        self.assertLessEqual(distances[2], distances[1], msg=f"各档距离 {distances}")

    def test_error_bound_uniform_in_beta(self):
        # 固定 error=0.1 扫 β：各点迹距离 ≤ error，误差不随 β 恶化到越界。
        # 实测距离 2.7e-6 / 1.5e-5 / 1.9e-4，随 c=βα/2 增长但始终远小于 error。
        hamiltonian = ((1.0 + 0j, 0j), (0j, -0.5 + 0j))
        error = 0.1
        distances = []
        for beta in (0.2, 0.5, 1.0):
            be = matrix_pauli_encoding(hamiltonian)
            access = gibbs_purification(be, beta, error=error)
            vector = dense_state(
                simulate(access.operation.program()), 1, signal_width=access.signal_qubits
            )
            reduced, _ = renormalized(partial_trace(vector, 1, 1))
            distances.append(trace_distance(reduced, gibbs_state(hamiltonian, beta)))
        for beta, distance in zip((0.2, 0.5, 1.0), distances, strict=True):
            self.assertLessEqual(
                distance, error, msg=f"β={beta} 迹距离 {distance} 超过 {error}；各点 {distances}"
            )

    def test_invalid_inputs_fail_at_generation(self):
        be = matrix_pauli_encoding(((1.0, 0.0), (0.0, -1.0)))
        cases = [
            lambda: gate_purification(((0.5, 0.0), (0.0, 0.4))),
            lambda: gate_purification(((0.5, 0.2), (0.0, 0.5))),
            lambda: gate_purification(((1.0, 0.0), (0.0, 0.0), (0.0, 0.0))),
            lambda: gibbs_purification(be, -1.0),
            lambda: gibbs_purification(be, 1.0, error=0.0),
            lambda: gibbs_purification(be, 1.0, error=1.5),
            lambda: gibbs_purification(object(), 1.0),
            lambda: maximally_mixed_purification(0),
            lambda: partial_trace([1.0, 0.0], 1, 1),
            lambda: trace_distance(((1.0, 0.0), (0.0, 0.0)), ((1.0, 0.0),)),
        ]
        for case in cases:
            with self.subTest(case=case):
                with self.assertRaises(ValidationError):
                    case()


if __name__ == "__main__":
    unittest.main()
