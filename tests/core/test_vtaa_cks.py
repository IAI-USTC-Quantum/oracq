"""VTAA-CKS 变时求解器的 GPE 判决、频带逆与端到端见证。"""

import math
import unittest
from dataclasses import replace

from oracq import FixedFormat, ValidationError, dumps, loads, simulate, unresolved
from oracq.algorithms.common.qsvt import qsp_response
from oracq.algorithms.input_model.oracles import (
    SparseAccess,
    gate_database,
    gate_state_prep,
    sparse_entry,
    sparse_location_gate,
)
from oracq.algorithms.qlss.qlss import LinearSystem, SparseSystem, SpectralPromise
from oracq.algorithms.qlss.vtaa_cks import (
    VTAAConfig,
    band_inverse_step,
    gapped_phase_estimation,
    gpe_fire_phases,
    make_vtaa_cks_qlss,
    tunable_rounds,
)
from oracq.infrastructure.ir import Adjoint, Control, Repeat


def chebyshev_value(coefficients, x):
    odd, even, total = x, 1.0, 0.0
    for coefficient in coefficients:
        total += coefficient * odd
        odd, even = 2 * x * odd - even, odd
        odd, even = 2 * x * odd - even, odd
    return total


def diagonal_problem(eigenvalues, entry_bound=1.5):
    """对角谱问题的稀疏输入；特征值必须落在定点网格上，entry_bound 留出松弛。

    谱 {1, 1/2}：kappa_phys = 2 对应两级时钟；x = lambda/alpha = 2/3 与 1/3，
    分别落在 band1 的 fire 验证边界与 band2 的 fire 验证边界。
    """
    fmt = FixedFormat(4, 2, signed=False)
    if any(abs(fmt.encode(v) / (1 << fmt.fraction) - v) > 1e-12 for v in eigenvalues):
        raise AssertionError("测试特征值必须在定点网格上")
    matrix = [[float(eigenvalues[r]) if r == c else 0.0 for c in range(2)] for r in range(2)]
    location = sparse_location_gate(1, [[0, 1], [1, 0]], work_width=0)
    entry = sparse_entry(
        gate_database(
            2, 4, {r + (c << 1): fmt.encode(matrix[r][c]) for r in range(2) for c in range(2)}
        ),
        1,
    )
    sparse = SparseSystem(
        SparseAccess(location, entry, 1, 4, 1),
        fmt,
        entry_bound,
        gate_state_prep([1.0, 1.0]),
        SpectralPromise(max(eigenvalues), min(eigenvalues)),
        diagonal_nonnegative=True,
        hermitian=True,
    )
    return LinearSystem(sparse=sparse, rhs_norm=math.sqrt(2))


def solution_amplitudes(result):
    state = simulate(result.state.operation.program()).amplitudes
    success = {k[0]: v for k, v in state.items() if k[1] == 0}
    mass = math.sqrt(sum(abs(v) ** 2 for v in success.values()))
    return {index: value / mass for index, value in success.items()}


class GappedPhaseEstimationTests(unittest.TestCase):
    def test_marker_circuit_matches_qsp_response_exactly(self):
        problem = diagonal_problem([1.0, 0.5])
        block, _ = problem.block_input()
        a = block.encoding
        self.assertEqual(a.alpha, 1.5)
        x_edge = block.spectrum.norm_upper / a.alpha
        # 寄存器顺序 target/signal/decision：fire 判决振幅精确等于 |P(x)|。
        for eigenvalue in (1.0, 0.5):
            x = eigenvalue / a.alpha
            for threshold in (x_edge, x_edge / 2):
                phases, _ = gpe_fire_phases(threshold, x_edge, 0.005, 40)
                gpe = gapped_phase_estimation(
                    a, threshold, x_edge, epsilon=0.005
                )
                state = simulate(
                    gpe.program(), initial={"target": 0 if eigenvalue == 1.0 else 1}
                ).amplitudes
                fire = sum(abs(v) ** 2 for k, v in state.items() if k[2] == 1)
                self.assertAlmostEqual(fire, abs(qsp_response(x, phases)) ** 2, places=11)
                self.assertTrue(
                    all(k[0] == (0 if eigenvalue == 1.0 else 1) for k in state)
                )

    def test_fire_side_hard_bound_and_junk_retained(self):
        problem = diagonal_problem([1.0, 0.5])
        block, _ = problem.block_input()
        a = block.encoding
        x_edge = block.spectrum.norm_upper / a.alpha
        gpe1 = gapped_phase_estimation(a, x_edge, x_edge, epsilon=0.005)
        state = simulate(gpe1.program(), initial={"target": 0}).amplitudes
        self.assertGreater(sum(abs(v) ** 2 for k, v in state.items() if k[2] == 1), 0.99)
        # 延后分支按论文保留 |gamma> 垃圾（signal 寄存器不复净）。
        small = simulate(gpe1.program(), initial={"target": 1}).amplitudes
        self.assertGreater(sum(abs(v) ** 2 for k, v in small.items() if k[1] != 0), 0.1)

    def test_gpe_rejects_tight_geometry(self):
        problem = diagonal_problem([1.0, 0.5])
        block, _ = problem.block_input()
        a = block.encoding
        with self.assertRaisesRegex(ValidationError, "上限"):
            gapped_phase_estimation(a, 0.01, 0.6, epsilon=0.005, degree_cap=12)


class BandInverseTests(unittest.TestCase):
    def test_band_inverse_step_matches_chebyshev_polynomial(self):
        problem = diagonal_problem([1.0, 0.5])
        a = problem.block_input()[0].encoding
        coefficients = (1.25, -0.25)
        program = band_inverse_step(a, coefficients, 1.5).program()
        for column, lam in ((0, 1.0), (1, 0.5)):
            state = simulate(program, initial={"target": column}).amplitudes
            x = lam / a.alpha
            expected = chebyshev_value(coefficients, x) / 1.5
            self.assertAlmostEqual(state.get((column, 0, 1), 0), expected, places=11)
            flagged = sum(abs(v) ** 2 for k, v in state.items() if k[1] == 0 and k[2] == 1)
            self.assertAlmostEqual(flagged, expected**2, places=11)
        with self.assertRaisesRegex(ValidationError, "alpha_max"):
            band_inverse_step(a, coefficients, 1.0)


class VTAAEndToEndTests(unittest.TestCase):
    def test_uniform_spectrum_single_band_is_exact(self):
        # kappa_phys=1：单频带，判决确定，h(x)=x 作用于谱 {x_edge}，
        # 均匀化旋转退化为 X，成功分支精确携带 |b> 方向。
        problem = diagonal_problem([1.0, 1.0])
        result = make_vtaa_cks_qlss(
            VTAAConfig(order=1, marker_epsilon=0.005)
        )(problem)
        amplitudes = solution_amplitudes(result)
        for index in (0, 1):
            self.assertAlmostEqual(amplitudes.get(index, 0), 1 / math.sqrt(2), delta=0.02)

    def test_variable_time_clock_separates_bands(self):
        # 谱 {1, 1/2}：lambda=1 在 band1 求逆（order 1，h(x)=x）；lambda=1/2 的
        # band1 判决落在过渡带（CKS 未承诺带），主要路径为延后至 band2（order 2）。
        # 逐本征值的路径幅值由判决响应精确给出，端到端断言取幅值区间。
        problem = diagonal_problem([1.0, 0.5])
        config = VTAAConfig(order=1, marker_epsilon=0.005)
        result = make_vtaa_cks_qlss(config)(problem)
        amplitudes = solution_amplitudes(result)
        block, _ = problem.block_input()
        a = block.encoding
        x_edge = block.spectrum.norm_upper / a.alpha
        phases1, _ = gpe_fire_phases(x_edge, x_edge, 0.005, 40)
        phases2, _ = gpe_fire_phases(x_edge / 2, x_edge, 0.005, 40)
        alpha_max = 1.5
        paths = []
        for eigenvalue, h in ((1.0, chebyshev_value((1.0,), 1.0 / a.alpha)),
                              (0.5, chebyshev_value((1.25, -0.25), 0.5 / a.alpha))):
            x = eigenvalue / a.alpha
            fire1 = abs(qsp_response(x, phases1))
            defer1 = math.sqrt(1 - fire1**2)
            fire2 = abs(qsp_response(x, phases2))
            defer2 = math.sqrt(1 - fire2**2)
            paths.append([
                fire1 * h / alpha_max,
                defer1 * fire2 * h / alpha_max,
                defer1 * defer2,
            ])
        dominant = [max(p) for p in paths]
        slack = [sum(p) - max(p) for p in paths]
        # 归一化把两个本征值的幅值耦合起来；区间按另一分量取极值端点计算。
        low = [max(0.0, d - s) for d, s in zip(dominant, slack, strict=True)]
        high = [d + s for d, s in zip(dominant, slack, strict=True)]
        for index in (0, 1):
            other = 1 - index
            mass = abs(amplitudes.get(index, 0))
            lower = low[index] / math.hypot(low[index], high[other])
            upper = high[index] / math.hypot(high[index], low[other])
            self.assertGreaterEqual(mass, lower - 1e-6)
            self.assertLessEqual(mass, upper + 1e-6)
        # 变时结构可见：band1 直接求逆的分量占优。
        self.assertGreater(dominant[0] / sum(dominant), 0.5)
        # 日程无关性：放大轮数只改成功率，条件解态幅值不变。
        amplified = solution_amplitudes(
            make_vtaa_cks_qlss(replace(config, rounds=(1, 0)))(problem)
        )
        for index in (0, 1):
            self.assertAlmostEqual(abs(amplitudes[index]), abs(amplified[index]), delta=0.03)


class VTAAStructureTests(unittest.TestCase):
    def test_protocol_output_structure_and_attributes(self):
        problem = diagonal_problem([1.0, 0.5])
        result = make_vtaa_cks_qlss(VTAAConfig(order=1))(problem)
        program = result.state.operation.program()
        attributes = dict(result.state.operation.module.attributes)
        self.assertEqual(attributes["qlss_protocol"], "vtaa_cks")
        self.assertEqual(attributes["qlss_input_model"], "sparse")
        self.assertFalse(unresolved(program))
        loads(dumps(program))
        module_attributes = {m.name: dict(m.attributes) for m in program.modules}
        kernel = next(
            attrs for attrs in module_attributes.values() if attrs.get("algorithm") == "vtaa_cks"
        )
        self.assertEqual(kernel["clock_steps"], 2)
        algorithms = {attrs.get("algorithm") for attrs in module_attributes.values()}
        for expected in (
            "vtaa_variable_step",
            "vtaa_uncompute_step",
            "gapped_phase_estimation",
            "vtaa_band_inverse",
            "vtaa_amplified_stage",
            "vtaa_prefix",
        ):
            self.assertIn(expected, algorithms)

    def amplified_repeat_counts(self, program):
        counts = []
        for module in program.modules:
            if dict(module.attributes).get("algorithm") != "vtaa_amplified_stage":
                continue

            def scan(nodes):
                for node in nodes:
                    if isinstance(node, Repeat) and node.count:
                        counts.append(node.count)
                    body = getattr(node, "body", ())
                    if isinstance(node, (Control, Adjoint)) and scan(body):
                        pass
                    if isinstance(node, Repeat):
                        scan(node.body)

            scan(module.body or ())
        return counts

    def test_repeat_rounds_survive_generation(self):
        problem = diagonal_problem([1.0, 0.5])
        flat = make_vtaa_cks_qlss(VTAAConfig(order=1))(problem)
        boosted = make_vtaa_cks_qlss(VTAAConfig(order=1, rounds=(0, 1)))(problem)
        self.assertEqual(self.amplified_repeat_counts(flat.state.operation.program()), [])
        self.assertEqual(
            self.amplified_repeat_counts(boosted.state.operation.program()), [1]
        )

    def test_estimate_counts_are_symbolic(self):
        problem = diagonal_problem([1.0, 0.5])
        result = make_vtaa_cks_qlss(VTAAConfig(order=1))(problem)
        estimate = result.state.operation.estimate()
        self.assertGreater(sum(estimate.atoms.values()), 0)
        self.assertGreater(estimate.qubits, problem.sparse.access.width)

    def test_rejects_invalid_configurations(self):
        with self.assertRaisesRegex(ValidationError, "0.2"):
            VTAAConfig(marker_epsilon=0.5)
        with self.assertRaisesRegex(ValidationError, "rounds"):
            VTAAConfig(rounds=())
        problem = diagonal_problem([1.0, 0.5])
        with self.assertRaisesRegex(ValidationError, "clock_steps"):
            make_vtaa_cks_qlss(VTAAConfig(clock_steps=1))(problem)
        with self.assertRaisesRegex(ValidationError, "上限"):
            make_vtaa_cks_qlss(VTAAConfig(clock_steps=3, degree_cap=14))(problem)
        with self.assertRaisesRegex(ValidationError, "rounds"):
            make_vtaa_cks_qlss(VTAAConfig(rounds=(0,)))(problem)

    def test_tunable_rounds_formula(self):
        self.assertEqual(tunable_rounds([1.0]), (0,))
        self.assertEqual(tunable_rounds([0.1]), (2,))
        self.assertEqual(tunable_rounds([0.2, 0.2]), (1, 1))
        with self.assertRaisesRegex(ValidationError, "阶段范数"):
            tunable_rounds([0.0])
        with self.assertRaisesRegex(ValidationError, "阈值"):
            tunable_rounds([0.5], thresholds=[0.1, 0.2])


if __name__ == "__main__":
    unittest.main()
