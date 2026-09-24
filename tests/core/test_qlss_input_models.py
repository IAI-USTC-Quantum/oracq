"""输入模型、缩放、稀疏适配与数据结构的针对性见证。"""

import math
import unittest

from oracq import FixedFormat, ValidationError, simulate, unresolved
from oracq.algorithms.input_model.oracles import (
    SparseAccess,
    basis_state,
    gate_database,
    qram_state_angles,
    qram_state_prep,
    sparse_entry,
    sparse_location_gate,
)
from oracq.algorithms.input_model.sparse import chebyshev_block
from oracq.algorithms.qlss.qlss import (
    BlockSystem,
    CKSConfig,
    CostaConfig,
    LinearSystem,
    SparseSystem,
    SpectralPromise,
    make_cks_qlss,
    make_costa_qlss,
)
from oracq.applications.flow_data import RoeFlowData
from oracq.applications.qfvm import roe_qfvm_inputs, roe_qfvm_problem, roe_qfvm_step
from oracq.infrastructure.execution import events
from oracq.infrastructure.ir import Load


def small_problem():
    fmt = FixedFormat(4, 2)
    matrix = [[0.5, -0.25], [-0.25, 0.5]]
    location = sparse_location_gate(1, [[0, 1], [0, 1]], work_width=0)
    entry = sparse_entry(
        gate_database(
            2, 4, {r + (c << 1): fmt.encode(matrix[r][c]) for r in range(2) for c in range(2)}
        ),
        1,
    )
    sparse = SparseSystem(
        SparseAccess(location, entry, 1, 4, 2),
        fmt,
        1,
        basis_state(1),
        SpectralPromise(0.75, 0.25),
        diagonal_nonnegative=True,
        hermitian=True,
    )
    return LinearSystem(sparse=sparse, rhs_norm=2), matrix


class QLSSInputTests(unittest.TestCase):
    def test_signed_sparse_encoding_and_chebyshev(self):
        problem, matrix = small_problem()
        encoded, _ = problem.block_input()
        be = encoded.encoding
        self.assertEqual(be.alpha, 2)
        for col in range(2):
            state = simulate(be.operation.program(), initial={"target": col})
            for row in range(2):
                self.assertAlmostEqual(
                    state.amplitudes.get((row, 0), 0) * be.alpha, matrix[row][col], places=11
                )
        h = [[v / be.alpha for v in row] for row in matrix]
        cube = [
            [
                sum(h[i][k] * h[k][inner] * h[inner][j] for k in range(2) for inner in range(2))
                for j in range(2)
            ]
            for i in range(2)
        ]
        for col in range(2):
            state = simulate(chebyshev_block(be, 3).operation.program(), initial={"target": col})
            for row in range(2):
                self.assertAlmostEqual(
                    state.amplitudes.get((row, 0), 0),
                    4 * cube[row][col] - 3 * h[row][col],
                    places=11,
                )

    def test_protocols_consume_different_models_and_scale_kappa(self):
        problem, _ = small_problem()
        costa = make_costa_qlss(CostaConfig(steps=1, kappa=999))
        a = costa(problem)
        b = make_cks_qlss(CKSConfig(order=2))(problem)
        self.assertEqual(a.encoded_inverse_bound, 8)
        self.assertEqual(a.state.width, b.state.width)
        self.assertEqual(a.input_alpha, b.input_alpha)
        self.assertEqual(dict(a.operation.module.attributes)["qlss_input_model"], "block_encoding")
        self.assertEqual(dict(b.operation.module.attributes)["qlss_input_model"], "sparse")
        underlying = next(
            m
            for m in a.operation.program().modules
            if dict(m.attributes).get("algorithm") == "costa_qlss"
        )
        self.assertEqual(dict(underlying.attributes)["encoded_inverse_bound"], 8)
        self.assertFalse(unresolved(a.operation.program()))
        self.assertFalse(unresolved(b.operation.program()))

    def test_costa_rhs_reflection_is_independent_of_unitary_extension(self):
        from oracq import Bits, Builder, identity
        from oracq.algorithms.input_model.oracles import StatePreparation, annotate
        from oracq.algorithms.qlss.qlss import costa_walk

        first = basis_state(1, work_width=1)
        builder = Builder("alternate_rhs_extension", {"target": Bits(1), "work": Bits(1)})
        builder.swap(builder["target"], builder["work"])
        second = StatePreparation(
            annotate(builder.finish(), "state_prep_isometry", zero_input=True, clean_work=True)
        )
        wa, wb = costa_walk(identity(1), first, 0.3), costa_walk(identity(1), second, 0.3)
        for target in range(2):
            for work in range(2):
                a = simulate(wa.program(), initial={"target": target, "signal": work}).amplitudes
                b = simulate(wb.program(), initial={"target": target, "signal": work}).amplitudes
                for key in set(a) | set(b):
                    self.assertAlmostEqual(a.get(key, 0), b.get(key, 0), places=11)

    def test_no_implicit_be_to_sparse(self):
        problem, _ = small_problem()
        block, _ = problem.block_input()
        only = LinearSystem(block=BlockSystem(block.encoding, block.rhs, block.spectrum))
        with self.assertRaisesRegex(ValidationError, "不能从一般 BE"):
            make_cks_qlss()(only)

    def test_norm_recovery_uses_conditional_matrix_probe(self):
        problem, _ = small_problem()
        result = make_cks_qlss()(problem)
        self.assertAlmostEqual(result.recover_norm(0.5, 0.125), 2.0)
        self.assertEqual(result.norm_probe.width, result.state.width)
        self.assertIn(
            "matrix_norm_probe", dict(result.norm_probe.operation.module.attributes)["algorithm"]
        )
        with self.assertRaises(ValidationError):
            result.recover_norm(0.1, 0.2)

    def test_matrix_probe_recovers_scalar_system_norm(self):
        from dataclasses import replace

        problem, _ = small_problem()
        scalar_entry = sparse_entry(gate_database(2, 4, {0: 2, 3: 2}), 1)
        access = replace(problem.sparse.access, entry=scalar_entry)
        result = make_cks_qlss()(replace(problem, sparse=replace(problem.sparse, access=access)))
        solution = simulate(result.operation.program())
        probe = simulate(result.norm_probe.operation.program())
        p_solution = sum(abs(v) ** 2 for k, v in solution.amplitudes.items() if k[1] == 0)
        p_joint = sum(abs(v) ** 2 for k, v in probe.amplitudes.items() if k[1] == 0)
        self.assertAlmostEqual(result.recover_norm(p_solution, p_joint), 4.0, places=10)

    def test_zero_rhs_is_not_a_state_preparation_problem(self):
        problem, _ = small_problem()
        from dataclasses import replace

        with self.assertRaisesRegex(ValidationError, "零右端"):
            make_cks_qlss()(replace(problem, rhs_norm=0))

    def test_qfvm_preserves_sparse_input_for_both_solvers(self):
        inputs = roe_qfvm_inputs(fmt=FixedFormat(4, 1), angle_width=4)
        problem = roe_qfvm_problem(inputs, spectrum=SpectralPromise(10, 0.5), rhs_norm=1, amax=4)
        self.assertEqual(problem.sparse.access.sparsity, 9)
        self.assertEqual(problem.physical_high_value, 1)
        for protocol in (make_cks_qlss(), make_costa_qlss(CostaConfig(steps=1))):
            result = protocol(problem)
            self.assertEqual(result.state.width, inputs.cell_width + 2)
            self.assertEqual(result.encoded_inverse_bound, 72)
            self.assertEqual(
                {r.name for r in unresolved(result.operation.program())},
                {"RoeQfvmGeometry", "RoeQfvmRho", "RoeQfvmMomentum", "RoeQfvmEnergy", "RoeQfvmRhs"},
            )
        with self.assertRaisesRegex(ValidationError, "SpectralPromise"):
            roe_qfvm_step(inputs, make_costa_qlss())

    def test_tree_preparation_queries_each_layer_coherently(self):
        prep = qram_state_prep(4, 6)
        self.assertEqual(
            sum(isinstance(node, Load) for node, _, _ in events(prep.operation.program())), 8
        )
        amplitudes = [
            0.0,
            1.0,
            2.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
        memory = {"angles": qram_state_angles(amplitudes, 6)}
        state = simulate(prep.operation.program(), memory)
        self.assertTrue(all(key[1] == 0 for key in state.amplitudes))
        expected_norm = math.sqrt(sum(v * v for v in amplitudes))
        for index, value in enumerate(amplitudes):
            self.assertAlmostEqual(
                abs(state.amplitudes.get((index, 0), 0)), value / expected_norm, delta=0.07
            )

    def test_local_update_does_not_scan_entire_flow(self):
        flow = RoeFlowData([(1, 0.125, 2)] * 16, fmt=FixedFormat(10, 5))

        class LocalOnly(list):
            def __iter__(self):
                raise AssertionError("局部更新不应复制或扫描全部流场")

        flow.states = LocalOnly(flow.states)
        patch = flow.update({5: (1.25, 0.25, 2.5)})
        self.assertEqual(patch.recomputed_faces, (4, 5))
        self.assertEqual(patch.recomputed_cells, (4, 5, 6))
        for cell in patch.recomputed_cells:
            for component, value in enumerate(flow.residuals[cell]):
                self.assertEqual(
                    flow.store.banks["rhs_values"].get(4 * cell + component, 0),
                    flow.fmt.encode(value),
                )
