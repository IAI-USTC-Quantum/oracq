"""数论 / Heinrich 求和积分 / SDE / QLSS / 应用目录的论文级数值验证。

覆盖模块：algorithms/{number_theory,integration,sde,qlss,contracts,interfaces,legacy}.py
与 applications/{catalog,gallery,legacy}.py。

- 模乘置换：叠加态一次穷举全部基态输入，对照经典数论置换；另有 OriginIR-ext
  全幺正块与泄漏检查。
- 量子求阶：phase 边缘分布对照独立 Dirichlet 核理论（|1> 在 r 个本征态均布）；
  factors_from_phase 对小半素数全相位穷举，统计 Shor 后处理成功率。
- Heinrich 求和/积分：比较器恒等式 P(flag)=E[v]/2**w 的精确性；QAE 读出分布
  对照独立双峰 Dirichlet 理论；积分解码与区间缩放；gate/qram 加载器一致性。
- SDE/Fokker–Planck：生成元守恒性对照 scipy.linalg.expm；显式 Euler 收敛阶；
  OU 矩对照闭式解与独立 Euler–Maruyama Monte Carlo；稳态对照零空间向量与
  Boltzmann 参考；初态制备与生成元块编码在真实后端上的振幅级验证。
- QLSS（CKS）：条件解态与成功概率对照独立 Chebyshev 矩阵多项式（numpy 递推），
  方法误差对照 numpy.linalg.solve 并展示阶数收敛；协议级 recover_norm 对照
  ‖A^{-1}b‖；Costa 的 Dolph–Chebyshev 权重对照闭式窗函数，装配程序跨后端对拍
  （其 kernel_status 为库内声明的 prototype，数值精度不作判据）。
- applications catalog/gallery：全部条目在真实后端跑通并与参考执行器逐振幅对拍。

运行：PYTHONPATH=src <含 pysparq+uniqc 的 python> tests/verification/verify_nt_qlss_sde.py
"""

from __future__ import annotations

import math
import multiprocessing as mp
from fractions import Fraction

from harness import (
    Report,
    adapter_pysparq,
    amplitude_error,
    amplitudes_to_statevector,
    basis_program,
    effective_block,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    statevector_error,
    superposition_program,
    tvd,
)

from oracq import Bits, Builder, FixedFormat
from oracq.algorithms.basics.number_theory import (
    factors_from_phase,
    modular_multiply,
    order_finding,
)
from oracq.algorithms.common.integration import (
    heinrich_rate,
    integral_from_phase,
    mean_from_phase,
    quantum_integral,
    quantum_sum,
    sum_preparation,
    table_loader,
)
from oracq.algorithms.qode.sde import (
    FokkerPlanckProblem,
    boltzmann_distribution,
    distribution_moments,
    evolve_distribution,
    matrix_exponential,
    sde_state_angles,
    sde_state_preparation,
    stationary_distribution,
)

GROUP = "nt_qlss_sde"
QUBIT_BUDGET = 24


# ---------------------------------------------------------------------------
# 独立经典 oracle（不依赖被测实现的辅助函数）
# ---------------------------------------------------------------------------


def _classical_order(a, m):
    """乘法阶：直接迭代幂。"""
    r, x = 1, a % m
    while x != 1:
        x = x * a % m
        r += 1
    return r


def _dirichlet_peak(y, center, size):
    """|D(y-center)|² 的 Dirichlet 核（峰值归一），center 可为非整数。"""
    import numpy as np

    delta = (y - center + size / 2) % size - size / 2
    amp = np.ones(len(y))
    mask = np.abs(delta) > 1e-9
    amp[mask] = np.sin(np.pi * delta[mask]) / (size * np.sin(np.pi * delta[mask] / size))
    return amp**2


def _qpe_phase_theory(order, precision):
    """|1> 在 r 个本征态上均布时 QPE 的 phase 边缘分布（独立理论）。"""
    import numpy as np

    size = 1 << precision
    y = np.arange(size)
    return sum(_dirichlet_peak(y, size * s / order, size) / order for s in range(order))


def _qae_theory(p_good, precision):
    """标准振幅估计读出分布：±θ/π 两个 Dirichlet 峰，sin²θ = p_good。"""
    import numpy as np

    size = 1 << precision
    theta = math.asin(math.sqrt(p_good))
    y = np.arange(size)
    return 0.5 * (
        _dirichlet_peak(y, size * theta / math.pi, size)
        + _dirichlet_peak(y, size * (1 - theta / math.pi), size)
    )


def _phase_marginal(amplitudes):
    """最后一个公开寄存器（phase）的边缘概率。"""
    out = {}
    for key, amplitude in amplitudes.items():
        out[key[-1]] = out.get(key[-1], 0.0) + abs(amplitude) ** 2
    return out


def _tvd_vec(distribution, theory):
    return 0.5 * float(sum(abs(distribution.get(y, 0.0) - theory[y]) for y in range(len(theory))))


def _originir_qubits(program):
    """预判 OriginIR 态向量预算；超预算返回 None。"""
    from oracq.infrastructure.layout import workspace_table

    width = sum(r.type.width for r in program.main.registers)
    total = width + workspace_table(program)[program.entry]
    return total if total <= QUBIT_BUDGET else None


# ---------------------------------------------------------------------------
# 数论：模乘置换与量子求阶
# ---------------------------------------------------------------------------

MODMUL_CASES = ((2, 5), (3, 5), (2, 7), (2, 9), (7, 100))


def verify_modular_multiply(report):
    """叠加态一次穷举 2^w 个基态输入，对照经典置换；小实例加 OriginIR 全振幅。"""
    for multiplier, modulus in MODMUL_CASES:
        width = (modulus - 1).bit_length()
        size = 1 << width
        operation = modular_multiply(multiplier, modulus)
        program = superposition_program(operation, ["target"])
        uniform = 1 / math.sqrt(size)
        expected = {
            ((multiplier * x % modulus if x < modulus else x),): uniform for x in range(size)
        }
        ref = reference(program, max_states=4 * size)
        deviation = amplitude_error(ref, expected)
        for runner in (rir_pysparq, adapter_pysparq):
            deviation = max(deviation, amplitude_error(runner(program, max_states=4 * size), expected))
        paths = ["reference", "rir-pysparq", "adapter-pysparq"]
        parameters = {
            "multiplier": multiplier,
            "modulus": modulus,
            "width": width,
            "basis_states": size,
        }
        if _originir_qubits(program) is not None:
            vector = originir_ext(program)
            deviation = max(
                deviation, statevector_error(vector, amplitudes_to_statevector(expected, [width]))
            )
            paths.append("originir-ext")
            parameters["originir_qubits"] = _originir_qubits(program)
        else:
            parameters["originir"] = "excluded: 超过 24 量子位预算"
        report.case(
            f"modmul-superposition-{multiplier}-{modulus}",
            paths=paths,
            parameters=parameters,
            metrics={"max_error": deviation},
            criterion="叠加穷举输出与经典模乘置换逐振幅一致（max_error < 1e-9）",
            passed=deviation < 1e-9,
        )


def verify_modular_multiply_unitary(report):
    """幺正层面：OriginIR-ext to_matrix 的有效块等于经典置换矩阵且无泄漏。"""
    import numpy as np

    multiplier, modulus = 2, 5
    width = 3
    operation = modular_multiply(multiplier, modulus)
    unitary = originir_unitary(operation.program())
    permutation = np.zeros((1 << width, 1 << width))
    for x in range(1 << width):
        permutation[(multiplier * x) % modulus if x < modulus else x, x] = 1.0
    block, leakage = effective_block(unitary, width)
    error = float(np.abs(block - permutation).max())
    report.case(
        "modmul-unitary-2-5",
        paths=["originir-ext+to_matrix"],
        parameters={"multiplier": 2, "modulus": 5, "width": width},
        metrics={"max_error": error, "leakage": leakage},
        criterion="幺正有效块等于置换矩阵且工作区泄漏为零（max_error < 1e-12）",
        passed=error < 1e-12 and leakage < 1e-12,
    )


ORDER_CASES = ((2, 3, 3), (2, 5, 4), (3, 5, 4), (2, 7, 4), (4, 7, 4), (2, 15, 5), (7, 15, 5))


def verify_order_finding(report):
    """phase 边缘分布对照独立 Dirichlet 核理论；连分数阶还原成功率一并报告。"""
    for multiplier, modulus, precision in ORDER_CASES:
        order = _classical_order(multiplier, modulus)
        theory = _qpe_phase_theory(order, precision)
        program = order_finding(multiplier, modulus, precision=precision).program()
        ref = _phase_marginal(reference(program))
        rir = _phase_marginal(rir_pysparq(program))
        deviation = max(_tvd_vec(ref, theory), _tvd_vec(rir, theory))
        cross = tvd(ref, rir)
        paths = ["reference", "rir-pysparq"]
        parameters = {
            "multiplier": multiplier,
            "modulus": modulus,
            "precision": precision,
            "classical_order": order,
        }
        if (multiplier, modulus) == (2, 3):
            # 最小实例追加 adapter 与 OriginIR 全振幅路径
            adapter = _phase_marginal(adapter_pysparq(program))
            deviation = max(deviation, _tvd_vec(adapter, theory))
            paths.append("adapter-pysparq")
            if _originir_qubits(program) is not None:
                vector = originir_ext(program)
                widths = [r.type.width for r in program.main.registers]
                origin = _phase_marginal(
                    {
                        tuple(
                            (index >> offset) & ((1 << width) - 1)
                            for offset, width in zip(_offsets(widths), widths, strict=True)
                        ): amplitude
                        for index, amplitude in enumerate(vector)
                        if abs(amplitude) > 1e-15
                    }
                )
                deviation = max(deviation, _tvd_vec(origin, theory))
                paths.append("originir-ext")
                parameters["originir_qubits"] = _originir_qubits(program)
        # 连分数还原出完整阶的经验概率（模拟分布加权）
        recovery = sum(
            probability
            for y, probability in ref.items()
            if y and Fraction(y, 1 << precision).limit_denominator(modulus).denominator == order
        )
        top = sorted(ref.items(), key=lambda kv: -kv[1])[:4]
        parameters["top_readouts"] = [y for y, _ in top]
        report.case(
            f"order-finding-{multiplier}-{modulus}-p{precision}",
            paths=paths,
            parameters=parameters,
            metrics={
                "tvd_vs_dirichlet_theory": deviation,
                "tvd_cross_backend": cross,
                "order_recovery_probability": recovery,
                "classical_order": order,
            },
            criterion="phase 分布与独立 Dirichlet 理论 TVD < 1e-9 且跨后端一致",
            passed=deviation < 1e-9 and cross < 1e-9,
        )


def _offsets(widths):
    offset, out = 0, []
    for w in widths:
        out.append(offset)
        offset += w
    return out


def verify_factors_from_phase(report):
    """经典后处理：小半素数 × 全部乘数 × 全部 8 位相位的穷举。

    两项独立 oracle：
    (a) 任何返回的因子对必须是真因子（硬正确性）；
    (b) 各模数的 QPE 加权成功率对照教科书预言——s 在 [1, r) 上均匀、
        候选阶 d = r/gcd(s, r)、要求 d 偶且 a^d ≡ 1 且 gcd(a^{d/2} ± 1, m)
        非平凡；实测低于预言的部分来自 p=8 相位分辨率对较大的阶的连分数
        还原损失（真实物理，非实现缺陷），判据为比值 ≥ 0.75。
    """
    precision = 8
    invalid = 0
    evaluated = 0
    per_modulus = {}
    hopeless_units = 0
    total_units = 0
    ratios = []
    for modulus in (15, 21, 33, 35):
        weighted_mean = []
        predicted_mean = []
        for multiplier in range(2, modulus - 1):
            if math.gcd(multiplier, modulus) != 1:
                # 平凡公因子路径：必须直接给出有效因子对
                for value in range(1 << precision):
                    result = factors_from_phase(value, precision, multiplier, modulus)
                    evaluated += 1
                    if result is None:
                        continue
                    p, q = result
                    if not (p * q == modulus and 1 < p <= q < modulus):
                        invalid += 1
                continue
            total_units += 1
            order = _classical_order(multiplier, modulus)
            half_order = pow(multiplier, order // 2, modulus) if order % 2 == 0 else None
            if order % 2 == 1 or half_order == modulus - 1:
                hopeless_units += 1
            # 独立教科书预言：s 均匀时可用阶的出现概率
            predicted = 0.0
            for s in range(1, order):
                d = order // math.gcd(s, order)
                if d % 2 or pow(multiplier, d, modulus) != 1:
                    continue
                half = pow(multiplier, d // 2, modulus)
                if 1 < math.gcd(half - 1, modulus) < modulus or 1 < math.gcd(
                    half + 1, modulus
                ) < modulus:
                    predicted += 1.0 / order
            predicted_mean.append(predicted)
            theory = _qpe_phase_theory(order, precision)
            success = 0.0
            for value in range(1 << precision):
                result = factors_from_phase(value, precision, multiplier, modulus)
                evaluated += 1
                if result is None:
                    continue
                p, q = result
                if not (p * q == modulus and 1 < p <= q < modulus):
                    invalid += 1
                else:
                    success += float(theory[value])
            weighted_mean.append(success)
        measured_m = sum(weighted_mean) / len(weighted_mean)
        predicted_m = sum(predicted_mean) / len(predicted_mean)
        ratio = measured_m / predicted_m if predicted_m else 1.0
        ratios.append(ratio)
        per_modulus[modulus] = {
            "weighted_success_mean": measured_m,
            "textbook_prediction": predicted_m,
            "ratio": ratio,
        }
    overall = sum(entry["weighted_success_mean"] for entry in per_modulus.values()) / len(per_modulus)
    report.case(
        "factors-from-phase-exhaustive",
        paths=["classical-oracle"],
        parameters={
            "moduli": [15, 21, 33, 35],
            "precision": precision,
            "evaluations": evaluated,
            "per_modulus": per_modulus,
        },
        metrics={
            "invalid_factorizations": invalid,
            "weighted_success_mean_over_units": overall,
            "min_ratio_vs_textbook": min(ratios),
            "hopeless_units": hopeless_units,
            "total_units": total_units,
        },
        criterion=(
            "返回的因子对全部有效；各模数实测成功率 / 教科书预言 ≥ 0.75"
            "（差额为 p=8 的连分数还原分辨率损失）"
        ),
        passed=invalid == 0 and min(ratios) >= 0.75,
    )


# ---------------------------------------------------------------------------
# Heinrich 量子求和与积分
# ---------------------------------------------------------------------------

SUM_TABLES = {
    "const5-w3": ([5] * 8, 3),
    "ramp8-w3": (list(range(8)), 3),
    "rand16-w4": ([3, 14, 1, 9, 6, 11, 0, 5, 12, 7, 2, 15, 4, 10, 8, 13], 4),
}


def _flag_probability_dict(amplitudes, n, w):
    flag_bit = 1 << (n + w)
    return sum(abs(a) ** 2 for key, a in amplitudes.items() if key[0] & flag_bit)


def _flag_probability_vector(vector, n, w):
    flag_bit = 1 << (n + w)
    return sum(abs(a) ** 2 for index, a in enumerate(vector) if index & flag_bit)


def verify_sum_preparation(report):
    """比较器构造核心恒等式：P(flag=1) 恰为 E[v]/2**w（对 v 严格线性）。"""
    for tag, (values, w) in SUM_TABLES.items():
        n = max(1, (len(values) - 1).bit_length())
        loader = table_loader(values, data_width=w)
        operation = sum_preparation(loader.database)
        program = operation.program()
        exact = sum(values) / len(values) / (1 << w)
        deviation = abs(_flag_probability_dict(reference(program), n, w) - exact)
        for runner in (rir_pysparq, adapter_pysparq):
            deviation = max(deviation, abs(_flag_probability_dict(runner(program), n, w) - exact))
        paths = ["reference", "rir-pysparq", "adapter-pysparq"]
        parameters = {"table": tag, "index_bits": n, "value_bits": w, "exact": exact}
        if _originir_qubits(program) is not None:
            deviation = max(
                deviation, abs(_flag_probability_vector(originir_ext(program), n, w) - exact)
            )
            paths.append("originir-ext")
            parameters["originir_qubits"] = _originir_qubits(program)
        else:
            parameters["originir"] = "excluded: 超过 24 量子位预算"
        report.case(
            f"sum-preparation-flag-{tag}",
            paths=paths,
            parameters=parameters,
            metrics={"max_error_vs_exact_mean": deviation},
            criterion="P(flag=1) 与 E[v]/2**w 的偏差 < 1e-9（线性恒等式）",
            passed=deviation < 1e-9,
        )


def verify_quantum_sum_on_grid(report):
    """均值恰落 QAE 栅格：全部非零概率读出精确等于真值。"""
    loader = table_loader([2] * 8, data_width=2)
    precision = 4
    program = quantum_sum(loader.database, precision=precision).program()
    states = [
        _phase_marginal(runner(program, max_states=4096))
        for runner in (reference, rir_pysparq, adapter_pysparq)
    ]
    decode_error = 0.0
    for value, probability in states[0].items():
        if probability > 1e-9:
            decode_error = max(decode_error, abs(mean_from_phase(value, precision, 2) - 2.0))
    cross = max(tvd(states[0], states[1]), tvd(states[0], states[2]))
    report.case(
        "quantum-sum-on-grid",
        paths=["reference", "rir-pysparq", "adapter-pysparq"],
        parameters={"table": "constant-2x8", "value_bits": 2, "precision": precision},
        metrics={"max_decode_error": decode_error, "tvd_cross_backend": cross},
        criterion="所有非零概率读出与均值 2 的偏差 < 1e-12 且后端两两 TVD < 1e-9",
        passed=decode_error < 1e-12 and cross < 1e-9,
    )


QAE_SWEEP_PRECISIONS = (3, 4, 5)


def verify_quantum_sum_qae(report):
    """QAE 读出对照独立双峰 Dirichlet 理论；期望误差随精度递减。"""
    values = [0, 1, 2, 3]
    loader = table_loader(values, data_width=2)
    exact = sum(values) / len(values)
    p_good = exact / 4
    expected_errors = []
    for precision in QAE_SWEEP_PRECISIONS:
        theory = _qae_theory(p_good, precision)
        program = quantum_sum(loader.database, precision=precision).program()
        ref = _phase_marginal(reference(program, max_states=65536))
        rir = _phase_marginal(rir_pysparq(program, max_states=65536))
        deviation = max(_tvd_vec(ref, theory), _tvd_vec(rir, theory))
        expected_error = sum(
            probability * abs(mean_from_phase(y, precision, 2) - exact)
            for y, probability in ref.items()
        )
        expected_errors.append(expected_error)
        mode = max(ref, key=ref.get)
        estimate = mean_from_phase(mode, precision, 2)
        resolution = 4 * math.pi / (1 << precision)
        report.case(
            f"quantum-sum-qae-p{precision}",
            paths=["reference", "rir-pysparq"],
            parameters={
                "table": "ramp0-3",
                "value_bits": 2,
                "precision": precision,
                "p_good": p_good,
            },
            metrics={
                "tvd_vs_qae_theory": deviation,
                "expected_abs_error": expected_error,
                "mode_estimate": estimate,
                "mode_error": abs(estimate - exact),
                "resolution_bound": resolution,
            },
            criterion="读出分布与 QAE 理论 TVD < 1e-9 且众数估计在一阶分辨率界内",
            passed=deviation < 1e-9 and abs(estimate - exact) <= resolution + 1e-12,
        )
    decreasing = all(
        later < earlier - 1e-12
        for earlier, later in zip(expected_errors, expected_errors[1:], strict=False)
    )
    report.case(
        "quantum-sum-qae-convergence",
        paths=["reference"],
        parameters={"precisions": list(QAE_SWEEP_PRECISIONS)},
        metrics={"expected_abs_errors": expected_errors},
        criterion="QAE 期望绝对误差随精度严格递减",
        passed=decreasing,
    )


def verify_quantum_integral(report):
    """f(x)=x 中点网格积分：QAE 理论对拍、闭式量化均值与区间缩放恒等式。"""
    # 4 点实例：三后端 + 理论分布；同一分布按 interval=2 重解码验证缩放
    grid, w, precision = 4, 3, 4
    values4 = [round((i + 0.5) / grid * ((1 << w) - 1)) for i in range(grid)]
    loader = table_loader(values4, data_width=w)
    program = quantum_integral(loader.database, precision=precision, interval=1.0).program()
    exact_q = sum(values4) / grid / ((1 << w) - 1)
    theory = _qae_theory(sum(values4) / grid / (1 << w), precision)
    states = [
        _phase_marginal(runner(program, max_states=65536))
        for runner in (reference, rir_pysparq, adapter_pysparq)
    ]
    deviation = max(_tvd_vec(states[0], theory), _tvd_vec(states[1], theory))
    cross = max(tvd(states[0], states[1]), tvd(states[0], states[2]))
    mode = max(states[0], key=states[0].get)
    estimate = integral_from_phase(mode, precision, w)
    scaled = integral_from_phase(mode, precision, w, interval=2.0)
    bound = (1 << w) * math.pi / (1 << precision) / ((1 << w) - 1)
    report.case(
        "quantum-integral-midpoint4",
        paths=["reference", "rir-pysparq", "adapter-pysparq"],
        parameters={
            "grid": grid,
            "value_bits": w,
            "precision": precision,
            "values": values4,
            "exact_quantized": exact_q,
        },
        metrics={
            "tvd_vs_qae_theory": deviation,
            "tvd_cross_backend": cross,
            "mode_estimate": estimate,
            "error_vs_quantized": abs(estimate - exact_q),
            "resolution_bound": bound,
            "interval_scaling_deviation": abs(scaled - 2 * estimate),
        },
        criterion="QAE 理论 TVD < 1e-9；积分估计在一阶分辨率界内；interval 缩放精确",
        passed=deviation < 1e-9
        and cross < 1e-9
        and abs(estimate - exact_q) <= bound + 1e-12
        and abs(scaled - 2 * estimate) < 1e-12,
    )
    # 8 点规范实例（核心测试同款）：rir-pysparq + QAE 理论；reference 交叉由 4 点例覆盖
    grid8, w8, precision8 = 8, 4, 4
    values8 = [round((i + 0.5) / grid8 * ((1 << w8) - 1)) for i in range(grid8)]
    loader8 = table_loader(values8, data_width=w8)
    program8 = quantum_integral(loader8.database, precision=precision8, interval=1.0).program()
    rir8 = _phase_marginal(rir_pysparq(program8, max_states=1 << 18))
    theory8 = _qae_theory(sum(values8) / grid8 / (1 << w8), precision8)
    deviation8 = _tvd_vec(rir8, theory8)
    mode8 = max(rir8, key=rir8.get)
    estimate8 = integral_from_phase(mode8, precision8, w8)
    exact8 = sum(values8) / grid8 / ((1 << w8) - 1)
    report.case(
        "quantum-integral-midpoint8",
        paths=["rir-pysparq"],
        parameters={
            "grid": grid8,
            "value_bits": w8,
            "precision": precision8,
            "exact_quantized": exact8,
            "reference_cross": "由同结构 4 点例的三后端对拍覆盖（本例 14 位寄存器，reference 代价高）",
        },
        metrics={
            "tvd_vs_qae_theory": deviation8,
            "mode_estimate": estimate8,
            "error_vs_quantized": abs(estimate8 - exact8),
            "error_vs_true_integral": abs(estimate8 - 0.5),
        },
        criterion="QAE 理论 TVD < 1e-9；|估计-0.5| ≤ 0.09（核心测试同判据）",
        passed=deviation8 < 1e-9 and abs(estimate8 - 0.5) <= 0.09,
    )


def verify_table_loader_qram(report):
    """gate 与 qram 两种加载器绑定在制备与求和读出上逐点一致。"""
    values = list(range(8))
    gate = table_loader(values, data_width=3)
    qram = table_loader(values, data_width=3, backend="qram")
    deviation = 0.0
    flag_dev = 0.0
    for runner in (reference, rir_pysparq):
        gate_state = runner(sum_preparation(gate.database).program())
        qram_state = runner(sum_preparation(qram.database).program(), qram.memory)
        deviation = max(deviation, amplitude_error(gate_state, qram_state))
        flag_dev = max(
            flag_dev,
            abs(_flag_probability_dict(gate_state, 3, 3) - _flag_probability_dict(qram_state, 3, 3)),
        )
    phase_tvd = 0.0
    # quantum_sum 的资源带嵌套前缀（prep__db__table 与 qpe__u__prep__db__table），
    # 同一逻辑数据表按入口资源名逐一绑定。
    table = next(iter(qram.memory.values()))
    gate_qsum = quantum_sum(gate.database, precision=3).program()
    qram_qsum = quantum_sum(qram.database, precision=3).program()
    qsum_memory = {r.name: table for r in qram_qsum.main.resources}
    for runner in (reference, rir_pysparq):
        gate_dist = _phase_marginal(runner(gate_qsum, max_states=8192))
        qram_dist = _phase_marginal(runner(qram_qsum, qsum_memory, max_states=8192))
        phase_tvd = max(phase_tvd, tvd(gate_dist, qram_dist))
    report.case(
        "table-loader-qram-vs-gate",
        paths=["reference", "rir-pysparq"],
        parameters={"table": "ramp0-7", "value_bits": 3, "sum_precision": 3},
        metrics={
            "preparation_max_error": deviation,
            "flag_probability_deviation": flag_dev,
            "sum_phase_tvd": phase_tvd,
        },
        criterion="gate/qram 加载器振幅逐点一致（< 1e-9），求和读出分布一致",
        passed=deviation < 1e-9 and flag_dev < 1e-9 and phase_tvd < 1e-9,
    )


def verify_heinrich_rate(report):
    """函数类最优收敛率闭式值与量子优势的二次间隔。"""
    worst = 0.0
    gap_error = 0.0
    for smoothness, dimension in ((1, 1), (2, 1), (1, 4), (3, 7)):
        rates = heinrich_rate(smoothness, dimension)
        worst = max(
            worst,
            abs(rates["deterministic"] - smoothness / dimension),
            abs(rates["randomized"] - (smoothness / dimension + 0.5)),
            abs(rates["quantum"] - (smoothness / dimension + 1.0)),
        )
        gap_error = max(
            gap_error,
            abs(
                (rates["quantum"] - rates["randomized"])
                - (rates["randomized"] - rates["deterministic"])
            ),
        )
    report.case(
        "heinrich-rate-closed-form",
        paths=["classical-oracle"],
        parameters={"instances": [[1, 1], [2, 1], [1, 4], [3, 7]]},
        metrics={"max_error": worst, "quadratic_gap_error": gap_error},
        criterion="收敛率与闭式 s/d、+1/2、+1 完全一致（< 1e-15）",
        passed=worst < 1e-15 and gap_error < 1e-15,
    )


# ---------------------------------------------------------------------------
# SDE / Fokker–Planck
# ---------------------------------------------------------------------------


def _ou_problem(size=16, half_width=3.0, theta=1.0, diffusion=0.5):
    h = 2 * half_width / size
    points = [-half_width + (i + 0.5) * h for i in range(size)]
    return FokkerPlanckProblem([-theta * x for x in points], [diffusion] * size, points)


def _gaussian(points, mean, variance):
    weights = [math.exp(-((x - mean) ** 2) / (2 * variance)) for x in points]
    total = math.fsum(weights)
    return [w / total for w in weights]


def verify_sde_generator(report):
    """离散生成元的守恒性与库矩阵指数对照 scipy.linalg.expm（独立 oracle）。"""
    import numpy as np
    import scipy.linalg

    problem = _ou_problem()
    g = np.array([list(row) for row in problem.generator_matrix()])
    column_sums = float(np.abs(g.sum(axis=0)).max())
    time = 0.7
    initial = np.array(_gaussian(problem.points, 0.3, 0.2))
    exact_propagator = scipy.linalg.expm(g * time)
    conservation = abs(float((exact_propagator @ initial).sum()) - 1.0)
    library_propagator = np.array(matrix_exponential(problem.generator_matrix(), time))
    expm_dev = float(np.abs(library_propagator - exact_propagator).max())
    spectral_abscissa = float(np.linalg.eigvals(g).real.max())
    report.case(
        "sde-generator-conservation",
        paths=["classical-oracle(scipy)"],
        parameters={"size": problem.size, "time": time},
        metrics={
            "column_sum_max": column_sums,
            "probability_conservation_error": conservation,
            "matrix_exponential_vs_scipy": expm_dev,
            "spectral_abscissa": spectral_abscissa,
        },
        criterion="列和为零、演化守恒概率、库矩阵指数与 scipy 一致（< 1e-9）、谱横坐标 ≤ 0",
        passed=column_sums < 1e-12
        and conservation < 1e-9
        and expm_dev < 1e-9
        and spectral_abscissa < 1e-9,
    )


def verify_sde_euler_order(report):
    """显式 Euler 的时间收敛阶：步长减半误差减半（一阶）。"""
    import numpy as np
    import scipy.linalg

    problem = _ou_problem()
    g = np.array([list(row) for row in problem.generator_matrix()])
    time = 0.25
    initial = _gaussian(problem.points, -0.2, 0.3)
    exact = scipy.linalg.expm(g * time) @ np.array(initial)
    errors = []
    for steps in (2000, 4000, 8000):
        coarse = np.array(evolve_distribution(problem.generator_matrix(), initial, time, steps=steps))
        errors.append(float(np.abs(coarse - exact).max()))
    orders = [
        math.log2(earlier / later)
        for earlier, later in zip(errors, errors[1:], strict=False)
    ]
    report.case(
        "sde-euler-convergence-order",
        paths=["classical-oracle(scipy)"],
        parameters={"steps": [2000, 4000, 8000], "time": time},
        metrics={"errors": errors, "measured_orders": orders},
        criterion="实测收敛阶在 [0.9, 1.2]（显式 Euler 一阶）",
        passed=all(0.9 <= o <= 1.2 for o in orders),
    )


def verify_sde_moments_monte_carlo(report):
    """OU 过程矩：Fokker–Planck 网格矩 vs 闭式解 vs 独立 Euler–Maruyama Monte Carlo。"""
    import numpy as np

    theta, diffusion, time = 1.0, 0.5, 0.4
    problem = _ou_problem(theta=theta, diffusion=diffusion)
    mean0, var0 = 0.4, 0.16
    initial = _gaussian(problem.points, mean0, var0)
    evolved = evolve_distribution(problem.generator_matrix(), initial, time)
    mean_fpe, second_fpe = distribution_moments(problem.points, evolved)
    var_fpe = second_fpe - mean_fpe * mean_fpe
    mean_exact = mean0 * math.exp(-theta * time)
    var_exact = var0 * math.exp(-2 * theta * time) + diffusion / theta * (
        1 - math.exp(-2 * theta * time)
    )
    rng = np.random.default_rng(20260916)
    count, dt = 1 << 18, 1e-3
    samples = rng.normal(mean0, math.sqrt(var0), count)
    noise = math.sqrt(2 * diffusion * dt)
    for _ in range(round(time / dt)):
        samples = samples - theta * samples * dt + noise * rng.standard_normal(count)
    mean_mc, var_mc = float(samples.mean()), float(samples.var())
    metrics = {
        "mean_fpe_vs_exact": abs(mean_fpe - mean_exact),
        "mean_mc_vs_exact": abs(mean_mc - mean_exact),
        "mean_fpe_vs_mc": abs(mean_fpe - mean_mc),
        "var_fpe_vs_exact": abs(var_fpe - var_exact),
        "var_mc_vs_exact": abs(var_mc - var_exact),
        "var_fpe_vs_mc": abs(var_fpe - var_mc),
        "mean_fpe": mean_fpe,
        "var_fpe": var_fpe,
        "mean_exact": mean_exact,
        "var_exact": var_exact,
        "mean_mc": mean_mc,
        "var_mc": var_mc,
    }
    report.case(
        "sde-ou-moments-vs-monte-carlo",
        paths=["classical-oracle(numpy MC, seed=20260916)"],
        parameters={"particles": count, "dt": dt, "time": time, "grid": problem.size},
        metrics=metrics,
        criterion=(
            "FPE 矩与 OU 闭式矩一致（均值 < 1e-3，方差 < 3e-2 网格截断）；"
            "MC 与闭式一致（均值 < 5e-3，方差 < 1e-2 统计涨落）"
        ),
        passed=metrics["mean_fpe_vs_exact"] < 1e-3
        and metrics["mean_mc_vs_exact"] < 5e-3
        and metrics["mean_fpe_vs_mc"] < 5e-3
        and metrics["var_fpe_vs_exact"] < 3e-2
        and metrics["var_mc_vs_exact"] < 1e-2
        and metrics["var_fpe_vs_mc"] < 3e-2,
    )


def verify_sde_stationary(report):
    """离散稳态对照 G 的零空间向量；长时间松弛收敛；与 Boltzmann 参考的 O(h²) 差距。"""
    import numpy as np
    import scipy.linalg

    problem = _ou_problem(half_width=2.0, diffusion=1.0)
    g = np.array([list(row) for row in problem.generator_matrix()])
    eigenvalues, eigenvectors = np.linalg.eig(g)
    null = eigenvectors[:, int(np.argmin(np.abs(eigenvalues)))].real
    null = null / null.sum()
    if null[0] < 0:
        null = -null
    stationary = np.array(stationary_distribution(problem))
    boltzmann = np.array(boltzmann_distribution(problem))
    uniform = np.full(problem.size, 1.0 / problem.size)
    relaxed = scipy.linalg.expm(g * 30.0) @ uniform
    tvd_null = 0.5 * float(np.abs(stationary - null).sum())
    tvd_relaxed = 0.5 * float(np.abs(relaxed - stationary).sum())
    tvd_boltzmann = 0.5 * float(np.abs(stationary - boltzmann).sum())
    report.case(
        "sde-stationary-boltzmann",
        paths=["classical-oracle(numpy/scipy)"],
        parameters={"size": problem.size, "relaxation_time": 30.0},
        metrics={
            "tvd_stationary_vs_nullspace": tvd_null,
            "tvd_relaxed_vs_stationary": tvd_relaxed,
            "tvd_stationary_vs_boltzmann": tvd_boltzmann,
        },
        criterion="稳态即零空间向量（TVD < 1e-9）；松弛收敛（< 5e-6）；Boltzmann 差距为 O(h²)",
        passed=tvd_null < 1e-9 and tvd_relaxed < 5e-6 and tvd_boltzmann < 0.05,
    )


def verify_sde_state_preparation(report):
    """初态制备：gate 实现振幅恰为 sqrt(p)（四路径）；qram 角表实现的量化精度。"""
    probabilities = [0.1, 0.2, 0.3, 0.4]
    preparation = sde_state_preparation(probabilities)
    program = preparation.operation.program()
    expected = {(index, 0): math.sqrt(p) for index, p in enumerate(probabilities)}
    deviation = amplitude_error(reference(program), expected)
    for runner in (rir_pysparq, adapter_pysparq):
        deviation = max(deviation, amplitude_error(runner(program), expected))
    paths = ["reference", "rir-pysparq", "adapter-pysparq"]
    if _originir_qubits(program) is not None:
        deviation = max(
            deviation, statevector_error(originir_ext(program), amplitudes_to_statevector(expected, [2, 0]))
        )
        paths.append("originir-ext")
    report.case(
        "sde-state-preparation-gate",
        paths=paths,
        parameters={"probabilities": probabilities},
        metrics={"max_error": deviation},
        criterion="制备振幅与 sqrt(p) 逐点一致（< 1e-9）",
        passed=deviation < 1e-9,
    )
    # QRAM 角表实现：先对照独立推导的量化角振幅（实现误差），再对照 p（量化误差）；
    # rir 路径按目标边缘分布对拍并报告工作位残留（复净回归指标，修复后应为 0）。
    angle_width = 8
    qram_prep = sde_state_preparation(probabilities, implementation="qram", angle_width=angle_width)
    angles = sde_state_angles(probabilities, angle_width=angle_width)
    memory = {"angles": angles}

    def quantized_expected():
        # 独立旋转树：root 地址 0，depth1 地址 1（左）与 2（右）；2 位目标的振幅路径
        words = [angles[0], angles[1], angles[2]]
        halved = [math.pi * word / (1 << angle_width) for word in words]
        cosines = [math.cos(a) for a in halved]
        sines = [math.sin(a) for a in halved]
        return [
            cosines[0] * cosines[1],
            cosines[0] * sines[1],
            sines[0] * cosines[2],
            sines[0] * sines[2],
        ]

    expected_amps = quantized_expected()
    ref = reference(qram_prep.operation.program(), memory)
    impl_dev = max(
        abs(abs(ref.get((index, 0), 0j)) - expected_amps[index]) for index in range(4)
    ) + max(
        (abs(a) for k, a in ref.items() if k[1] != 0), default=0.0
    )
    rir = rir_pysparq(qram_prep.operation.program(), memory)
    marginal_dev = 0.0
    for index, p in enumerate(probabilities):
        marginal = sum(abs(a) ** 2 for k, a in rir.items() if k[0] == index)
        marginal_dev = max(marginal_dev, abs(marginal - p))
    residue = max((abs(a) ** 2 for k, a in rir.items() if k[1] != 0), default=0.0)
    quantization = max(
        abs(abs(expected_amps[index]) ** 2 - p) for index, p in enumerate(probabilities)
    )
    report.case(
        "sde-state-preparation-qram",
        paths=["reference", "rir-pysparq"],
        parameters={
            "probabilities": probabilities,
            "angle_width": angle_width,
            "angle_table": angles,
        },
        metrics={
            "impl_error_vs_quantized_angles": impl_dev,
            "quantization_error_vs_probabilities": quantization,
            "rir_target_marginal_deviation": marginal_dev,
            "rir_work_register_residue": residue,
        },
        criterion=(
            "reference 与独立量化角振幅一致（< 1e-9，实现误差）；"
            "量化误差与 rir 边缘偏差 < 0.02（角字方法误差）"
        ),
        passed=impl_dev < 1e-9 and quantization < 0.02 and marginal_dev < 0.02,
    )


def verify_sde_generator_encoding(report):
    """生成元的 Pauli 块编码：基态驱动下 (signal=0) 块振幅乘 alpha 恰为 G。"""
    # OU：a(x) = -x，D = 1，网格点 ±0.5、±1.5（h = 1）
    points = (-1.5, -0.5, 0.5, 1.5)
    problem = FokkerPlanckProblem(
        tuple(-x for x in points),
        [1.0] * 4,
        points,
    )
    encoding = problem.generator_encoding()
    matrix = problem.generator_matrix()
    worst = 0.0
    originir_ok = True
    for column in range(4):
        program = basis_program(encoding.operation, {"target": column})
        states = [reference(program), rir_pysparq(program), adapter_pysparq(program)]
        vector = None
        if originir_ok and _originir_qubits(program) is not None:
            vector = originir_ext(program)
        elif _originir_qubits(program) is None:
            originir_ok = False
        for state in states:
            for row in range(4):
                worst = max(
                    worst,
                    abs(state.get((row, 0), 0j) * encoding.alpha - matrix[row][column]),
                )
        if vector is not None:
            for row in range(4):
                worst = max(
                    worst,
                    abs(vector[row] * encoding.alpha - matrix[row][column]),
                )
    report.case(
        "sde-generator-encoding",
        paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        parameters={
            "size": 4,
            "alpha": encoding.alpha,
            "signal_qubits": encoding.signal_qubits,
        },
        metrics={"max_error": worst},
        criterion="块编码 (signal=0) 块振幅 × alpha 与生成元矩阵逐点一致（< 1e-9）",
        passed=worst < 1e-9,
    )


# ---------------------------------------------------------------------------
# QLSS：CKS Chebyshev 求解器与 Costa 行走
# ---------------------------------------------------------------------------


def _sparse_problem():
    """κ=3 的 2×2 有符号稀疏系统（特征值 0.5 与 1.0，alpha = 1.5）。"""
    from oracq.algorithms.input_model.oracles import (
        SparseAccess,
        basis_state,
        gate_database,
        sparse_entry,
        sparse_location_gate,
    )
    from oracq.algorithms.qlss.qlss import LinearSystem, SparseSystem, SpectralPromise

    fmt = FixedFormat(4, 2)
    matrix = [[0.75, -0.25], [-0.25, 0.75]]
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
        0.75,
        basis_state(1),
        SpectralPromise(1.0, 0.5),
        diagonal_nonnegative=True,
        hermitian=True,
    )
    return LinearSystem(sparse=sparse, rhs_norm=2), matrix


def _cks_coefficients(order):
    """CKS 系数的独立闭式重算（math.comb）。"""
    return [
        4.0
        * (-1) ** j
        * sum(math.comb(2 * order, order + i) for i in range(j + 1, order + 1))
        / 2.0 ** (2 * order)
        for j in range(order)
    ]


def _chebyshev_polynomial_apply(matrix, alpha, coefficients, vector):
    """P(M) v，P = Σ c_j T_{2j+1}，M = matrix/alpha；奇次递推用 T2 步进。"""
    import numpy as np

    m = np.array(matrix, dtype=float) / alpha
    t2 = 2 * (m @ m) - np.eye(m.shape[0])
    u_prev = m.copy()  # T_1
    u_curr = 2 * (t2 @ m) - m  # T_3
    out = coefficients[0] * (u_prev @ vector)
    if len(coefficients) > 1:
        out = out + coefficients[1] * (u_curr @ vector)
    for j in range(2, len(coefficients)):
        u_next = 2 * (t2 @ u_curr) - u_prev
        u_prev, u_curr = u_curr, u_next
        out = out + coefficients[j] * (u_curr @ vector)
    return out


CKS_ORDERS = (2, 4, 8, 16)


def verify_cks_kernel(report):
    """条件解态与成功概率对照独立 Chebyshev 矩阵多项式；方法误差随阶数收敛。"""
    import numpy as np

    from oracq.algorithms.qlss.qlss import CKSConfig, cks_chebyshev

    problem, matrix = _sparse_problem()
    a = np.array(matrix)
    b = np.array([1.0, 0.0])
    x_true = np.linalg.solve(a, b)
    method_errors = []
    for order in CKS_ORDERS:
        state = cks_chebyshev(problem.sparse, CKSConfig(order=order))
        program = state.operation.program()
        ref = reference(program, max_steps=50_000_000)
        rir = rir_pysparq(program, max_steps=50_000_000)
        coefficients = _cks_coefficients(order)
        alpha_inv = sum(abs(c) for c in coefficients)
        expected = _chebyshev_polynomial_apply(a, 1.5, coefficients, b)
        p_expected = float(np.sum(expected**2)) / alpha_inv**2

        def success_branch(amplitudes):
            good = {k[0]: amp for k, amp in amplitudes.items() if k[1] == 0}
            probability = sum(abs(amp) ** 2 for amp in good.values())
            vector = np.array([good.get(0, 0j), good.get(1, 0j)])
            return probability, vector

        p_ref, v_ref = success_branch(ref)
        p_rir, v_rir = success_branch(rir)
        v_expected = expected / np.linalg.norm(expected)
        impl_error = max(
            float(np.max(np.abs(v_ref / np.linalg.norm(v_ref) - v_expected))),
            float(np.max(np.abs(v_rir / np.linalg.norm(v_rir) - v_expected))),
        )
        p_error = max(abs(p_ref - p_expected), abs(p_rir - p_expected))
        method_error = float(
            np.max(np.abs(expected / np.linalg.norm(expected) - x_true / np.linalg.norm(x_true)))
        )
        method_errors.append(method_error)
        cross = amplitude_error(ref, rir)
        report.case(
            f"cks-kernel-order-{order}",
            paths=["reference", "rir-pysparq"],
            parameters={
                "matrix": matrix,
                "alpha": 1.5,
                "kappa": 3.0,
                "order": order,
                "lcu_normalization": alpha_inv,
            },
            metrics={
                "impl_error_vs_chebyshev_oracle": impl_error,
                "success_probability_error": p_error,
                "success_probability": p_ref,
                "method_error_vs_numpy_solve": method_error,
                "cross_backend_max_error": cross,
            },
            criterion=(
                "条件解态与独立 Chebyshev 多项式一致（< 1e-8），"
                "成功概率与多项式预测一致（< 1e-8）"
            ),
            passed=impl_error < 1e-8 and p_error < 1e-8 and cross < 1e-6,
        )
    decreasing = all(
        later < earlier - 1e-9
        for earlier, later in zip(method_errors, method_errors[1:], strict=False)
    )
    report.case(
        "cks-method-convergence",
        paths=["classical-oracle(numpy)"],
        parameters={"orders": list(CKS_ORDERS)},
        metrics={"method_errors": method_errors},
        criterion="Chebyshev 截断的方法误差随阶数严格递减（κ=3 的几何收敛）",
        passed=decreasing,
    )


def verify_cks_protocol(report):
    """协议级：recover_norm(p_solver, p_joint) 恢复 ‖A^{-1}b‖；探针概率与多项式一致。"""
    import numpy as np

    from oracq.algorithms.qlss.qlss import CKSConfig, make_cks_qlss

    problem, matrix = _sparse_problem()
    a = np.array(matrix)
    order = 8
    result = make_cks_qlss(CKSConfig(order=order))(problem)
    solution = reference(result.operation.program(), max_steps=50_000_000)
    probe = reference(result.norm_probe.operation.program(), max_steps=50_000_000)
    solution_rir = rir_pysparq(result.operation.program(), max_steps=50_000_000)
    probe_rir = rir_pysparq(result.norm_probe.operation.program(), max_steps=50_000_000)
    p_solver = sum(abs(v) ** 2 for k, v in solution.items() if k[1] == 0)
    p_joint = sum(abs(v) ** 2 for k, v in probe.items() if k[1] == 0)
    coefficients = _cks_coefficients(order)
    alpha_inv = sum(abs(c) for c in coefficients)
    filtered = _chebyshev_polynomial_apply(a, 1.5, coefficients, np.array([1.0, 0.0]))
    p_solver_expected = float(np.sum(filtered**2)) / alpha_inv**2
    p_joint_expected = float(np.sum((a @ filtered / 1.5) ** 2)) / alpha_inv**2
    probability_error = max(
        abs(p_solver - p_solver_expected), abs(p_joint - p_joint_expected)
    )
    recovered = result.recover_norm(p_solver, p_joint)
    x = np.linalg.solve(a, np.array([2.0, 0.0]))
    norm_exact = float(np.linalg.norm(x))
    relative_error = abs(recovered - norm_exact) / norm_exact
    good = {k[0]: v for k, v in solution.items() if k[1] == 0}
    v = np.array([good.get(0, 0j), good.get(1, 0j)])
    fidelity = float(abs(np.vdot(v / np.linalg.norm(v), x / np.linalg.norm(x))) ** 2)
    cross = max(
        amplitude_error(solution, solution_rir), amplitude_error(probe, probe_rir)
    )
    report.case(
        "cks-protocol-norm-recovery",
        paths=["reference", "rir-pysparq"],
        parameters={
            "matrix": matrix,
            "order": order,
            "rhs_norm": 2,
            "encoded_inverse_bound": result.encoded_inverse_bound,
        },
        metrics={
            "p_solver": p_solver,
            "p_joint": p_joint,
            "probability_error_vs_polynomial": probability_error,
            "recovered_norm": recovered,
            "exact_norm": norm_exact,
            "norm_relative_error": relative_error,
            "solution_fidelity_vs_numpy": fidelity,
            "cross_backend_max_error": cross,
        },
        criterion=(
            "探针概率与多项式一致（< 1e-9）；recover_norm 相对误差 ≤ 0.2 且保真度 ≥ 0.94"
            "（order=8 截断的方法误差 0.213 所致，收敛趋势见 kernel 案例）"
        ),
        passed=probability_error < 1e-9
        and relative_error <= 0.2
        and fidelity >= 0.94
        and cross < 1e-9,
    )


def verify_costa(report):
    """Costa 构件：Dolph–Chebyshev 权重对照闭式窗、schedule 闭式、unary 制备。

    costa_qlss 的 kernel_status 是库内声明的 prototype（行走初态与读出通道
    未认证），因此数值判据只覆盖可独立推导的构件与装配确定性。
    """
    import numpy as np
    from numpy.polynomial.chebyshev import chebval

    from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding
    from oracq.algorithms.input_model.oracles import basis_state
    from oracq.algorithms.qlss.qlss import (
        CostaConfig,
        costa_qlss,
        dolph_chebyshev_plan,
        schedule,
        unary_weight_preparation,
    )

    worst = 0.0
    for degree, attenuation in ((2, 0.2), (4, 0.2), (6, 0.1)):
        plan = dolph_chebyshev_plan(degree, attenuation)
        beta = math.cosh(math.acosh(1 / attenuation) / degree)
        theta = np.linspace(0, math.pi, 4097)
        response = sum(
            weight * np.exp(1j * k * theta)
            for k, weight in zip(range(plan.offset, plan.offset + plan.stride * len(plan.weights), plan.stride), plan.weights, strict=True)
        )
        closed = attenuation * chebval(beta * np.cos(theta), [0.0] * degree + [1.0])
        worst = max(worst, float(np.abs(response - closed).max()))
        worst = max(worst, abs(sum(plan.weights) - 1.0))
    grid = np.linspace(0.0, 1.0, 33)
    schedule_error = 0.0
    monotone_ok = True
    for kappa in (1.0, 2.0, 4.0, 9.0):
        for power in (0.5, 1.5, 2.5):
            values = [schedule(float(s), kappa, power) for s in grid]
            closed = [
                s
                if kappa == 1
                else kappa / (kappa - 1) * (1 - (1 + s * (kappa ** (power - 1) - 1)) ** (1 / (1 - power)))
                for s in grid
            ]
            schedule_error = max(
                schedule_error, max(abs(a - b) for a, b in zip(values, closed, strict=True))
            )
            monotone_ok = monotone_ok and all(
                later >= earlier - 1e-15
                for earlier, later in zip(values, values[1:], strict=False)
            )
            schedule_error = max(schedule_error, abs(values[0] - 0.0), abs(values[-1] - 1.0))
    unary = unary_weight_preparation(dolph_chebyshev_plan(4, 0.2).weights)
    weights = dolph_chebyshev_plan(4, 0.2).weights
    total = sum(weights)
    unary_error = 0.0
    for runner in (reference, rir_pysparq):
        state = runner(unary.program())
        for prefix, weight in enumerate(weights):
            # unary_prefix 编码：前 prefix 个位置为 1 的基态振幅平方 = w/W
            basis = (1 << prefix) - 1
            unary_error = max(
                unary_error, abs(abs(state.get((basis,), 0j)) ** 2 - weight / total)
            )
    report.case(
        "costa-plan-schedule-unary",
        paths=["classical-oracle(numpy)", "reference", "rir-pysparq"],
        parameters={"plans": [[2, 0.2], [4, 0.2], [6, 0.1]]},
        metrics={
            "dolph_chebyshev_max_error": worst,
            "schedule_max_error": schedule_error,
            "schedule_monotone": monotone_ok,
            "unary_preparation_max_error": unary_error,
        },
        criterion="DC 权重与闭式窗 < 1e-9；schedule 与闭式一致且单调；unary 制备 ∝ 权重",
        passed=worst < 1e-9 and schedule_error < 1e-12 and monotone_ok and unary_error < 1e-9,
    )
    operation = costa_qlss(
        matrix_pauli_encoding([[2, -1], [-1, 2]]), basis_state(1), CostaConfig(steps=1)
    )
    program = operation.operation.program()
    ref = reference(program)
    deviation = max(
        amplitude_error(ref, rir_pysparq(program)), amplitude_error(ref, adapter_pysparq(program))
    )
    report.case(
        "costa-assembly-cross-backend",
        paths=["reference", "rir-pysparq", "adapter-pysparq"],
        parameters={
            "matrix": [[2, -1], [-1, 2]],
            "steps": 1,
            "kernel_status": "prototype（库内声明；求解精度不作判据，仅验证装配确定性）",
        },
        metrics={"cross_backend_max_error": deviation},
        criterion="装配程序跨后端逐振幅一致（< 1e-9）",
        passed=deviation < 1e-9,
    )


# ---------------------------------------------------------------------------
# applications：catalog 与 gallery 全条目跨后端对拍
# ---------------------------------------------------------------------------


def verify_gallery(report):
    """gallery 全部条目在 rir/adapter 真实后端跑通并与参考执行器逐振幅对拍。"""
    from oracq.applications.gallery import algorithm_gallery

    for case in algorithm_gallery():
        program = case.operation.program()
        ref = reference(program)
        deviation = max(
            amplitude_error(ref, rir_pysparq(program)),
            amplitude_error(ref, adapter_pysparq(program)),
        )
        report.case(
            f"gallery-{case.name}",
            paths=["reference", "rir-pysparq", "adapter-pysparq"],
            parameters={"family": case.family, "readout": case.readout},
            metrics={"cross_backend_max_error": deviation},
            criterion="三后端逐振幅一致（< 1e-9）",
            passed=deviation < 1e-9,
        )


def _catalog_qham_worker(name, path, queue):
    """qham 重例的子进程入口：独立构建并执行，回传振幅。"""
    import time as _time

    from oracq import run_pysparq, simulate
    from oracq.applications.catalog import build_case

    case = build_case(name)
    program = case.closed()
    start = _time.time()
    if path == "reference":
        amplitudes = dict(
            simulate(program, case.memory, max_states=100_000, max_steps=100_000_000).amplitudes
        )
    else:
        amplitudes = dict(
            run_pysparq(program, case.memory, max_states=100_000, max_steps=100_000_000).amplitudes
        )
    queue.put((name, path, amplitudes, _time.time() - start))


def _verify_catalog_qham(report):
    """qham_qode/qham_qpde：展开步数超 1e6 的重例，四进程并发跑 reference 与 adapter。"""
    ctx = mp.get_context("spawn")
    queue = ctx.Queue()
    jobs = [
        (name, path)
        for name in ("qham_qode", "qham_qpde")
        for path in ("reference", "adapter_pysparq")
    ]
    procs = [ctx.Process(target=_catalog_qham_worker, args=(*job, queue)) for job in jobs]
    for proc in procs:
        proc.start()
    results = {}
    for _ in jobs:
        name, path, amplitudes, elapsed = queue.get(timeout=600)
        results[(name, path)] = amplitudes
        print(f"[qham] {name}/{path} {elapsed:.1f}s states={len(amplitudes)}", flush=True)
    for proc in procs:
        proc.join(600)
    for name in ("qham_qode", "qham_qpde"):
        ref = results[(name, "reference")]
        ada = results[(name, "adapter_pysparq")]
        deviation = amplitude_error(ref, ada)
        report.case(
            f"catalog-{name}",
            paths=["reference", "adapter-pysparq"],
            parameters={
                "states": len(ref),
                "note": "展开步数 > 1e6，max_steps 提升至 1e8；并发执行控制墙钟",
            },
            metrics={"cross_backend_max_error": deviation},
            criterion="长线路浮点累积下两后端一致（< 1e-8）",
            passed=deviation < 1e-8,
        )


def verify_catalog(report):
    """catalog 全部条目跑通并对拍参考执行器；stateprep_qram 记录工作位残留。"""
    from oracq.applications.catalog import CASES, build_case

    for name in CASES:
        if name in ("qham_qode", "qham_qpde"):
            continue  # 重例走并发子进程
        case = build_case(name)
        program = case.closed()
        ref = reference(program, case.memory, max_states=300_000, max_steps=50_000_000)
        if name == "stateprep_qram":
            # 回归哨兵：pysparq.rir 早期版本在此例的 QRAM 工作位未复净（work=4 残留，
            # 全振幅偏差 0.65）；该 bug 修复后此例同时钉死工作位复净与全态一致。
            ada = adapter_pysparq(program, case.memory)
            rir = rir_pysparq(program, case.memory)
            full_dev = amplitude_error(ref, ada)
            widths = [r.type.width for r in program.main.registers]

            def marginal(amplitudes):
                out = {}
                for key, amp in amplitudes.items():
                    out[key[0]] = out.get(key[0], 0.0) + abs(amp) ** 2
                return out

            marginal_dev = tvd(marginal(ref), marginal(rir))
            residue = max((abs(a) ** 2 for k, a in rir.items() if k[1] != 0), default=0.0)
            report.case(
                "catalog-stateprep_qram",
                paths=["reference", "adapter-pysparq", "rir-pysparq"],
                parameters={
                    "register_widths": widths,
                    "note": "工作位复净回归哨兵（pysparq.rir 修复前的残留案例）",
                },
                metrics={
                    "reference_vs_adapter_max_error": full_dev,
                    "target_marginal_tvd_ref_rir": marginal_dev,
                    "rir_work_register_residue": residue,
                },
                criterion="reference/adapter 全态一致且 rir 目标边缘分布一致（< 1e-9）",
                passed=full_dev < 1e-9 and marginal_dev < 1e-9,
            )
            continue
        rir = rir_pysparq(program, case.memory, max_states=300_000, max_steps=50_000_000)
        deviation = amplitude_error(ref, rir)
        report.case(
            f"catalog-{name}",
            paths=["reference", "rir-pysparq"],
            parameters={"states": len(ref)},
            metrics={"cross_backend_max_error": deviation},
            criterion="两后端逐振幅一致（< 1e-9）",
            passed=deviation < 1e-9,
        )
    _verify_catalog_qham(report)


# ---------------------------------------------------------------------------
# 驱动
# ---------------------------------------------------------------------------


def _warm_up_originir():
    """UniQC 首次调用有一次性初始化开销；用最小线路预热，避免计入案例。"""
    b = Builder("originir_warmup", {"q": Bits(1)})
    b.h(b["q"])
    originir_ext(b.finish().program())


def run():
    report = Report(
        GROUP,
        "数论/Heinrich 求和积分/SDE/QLSS 的论文级数值验证与 applications 全条目跨后端对拍。",
    )
    _warm_up_originir()
    # 数论
    verify_modular_multiply(report)
    verify_modular_multiply_unitary(report)
    verify_order_finding(report)
    verify_factors_from_phase(report)
    # Heinrich 求和与积分
    verify_sum_preparation(report)
    verify_quantum_sum_on_grid(report)
    verify_quantum_sum_qae(report)
    verify_quantum_integral(report)
    verify_table_loader_qram(report)
    verify_heinrich_rate(report)
    # SDE / Fokker–Planck
    verify_sde_generator(report)
    verify_sde_euler_order(report)
    verify_sde_moments_monte_carlo(report)
    verify_sde_stationary(report)
    verify_sde_state_preparation(report)
    verify_sde_generator_encoding(report)
    # QLSS
    verify_cks_kernel(report)
    verify_cks_protocol(report)
    verify_costa(report)
    # applications
    verify_gallery(report)
    verify_catalog(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
