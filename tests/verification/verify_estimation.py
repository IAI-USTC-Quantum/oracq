"""估计类算法（estimation.py / gradient.py）的论文级数值验证。

覆盖五个读出算法，全部以真实后端跑出的分布/期望值对照独立经典闭式：

- QPE：已知本征相位的幺正（相位门、add_const 的 Fourier 本征态），相位寄存器
  直方图对照 Dirichlet 核闭式分布，多精度位数扫描；另有 OriginIR 全幺正矩阵
  对照 numpy 独立构造的 QPE 幺正。
- QAE：均匀/非均匀制备下的幅度估计，相位分布对照 Grover 双峰闭式，峰值解码
  误差对照 Brassard 界。
- Hadamard test：probe 的 Z 期望对照由角度表/态矢量独立计算的 Re/Im <ψ|U|ψ>。
- Swap test：probe=0 概率对照 (1+|<a|c>|^2)/2，重叠由 numpy 独立计算。
- Jordan 梯度：线性函数（gate 相位表与 mathfunc 相位 oracle）读出精确等于
  梯度，并与中心有限差分对照；扰动线性的失败概率衰减率随网格位数下降。

运行：PYTHONPATH=src <含 pysparq+uniqc 的 python> tests/verification/verify_estimation.py
"""

from __future__ import annotations

import math

from harness import (
    Report,
    adapter_pysparq,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    tvd,
)

from pyqecclang import Bits, Builder, UInt
from pyqecclang.algorithms.common.estimation import (
    amplitude_estimation,
    hadamard_test,
    phase_estimation,
    swap_test,
)
from pyqecclang.algorithms.common.fourier import qft
from pyqecclang.algorithms.input_model.oracles import gate_state_prep, uniform_state
from pyqecclang.algorithms.optimization.gradient import (
    function_phase_oracle,
    gate_phase_oracle,
    gradient_estimation,
    gradient_from_readout,
)
from pyqecclang.infrastructure.layout import workspace_table

# ---------------------------------------------------------------------------
# 通用辅助
# ---------------------------------------------------------------------------


def _marginal_from_amplitudes(amplitudes, registers, name):
    """字典稀疏态（寄存器元组→振幅）对指定寄存器的边际分布。"""
    index = [r.name for r in registers].index(name)
    result = {}
    for key, amplitude in amplitudes.items():
        result[key[index]] = result.get(key[index], 0.0) + abs(amplitude) ** 2
    return result


def _marginal_from_statevector(vector, registers, name):
    """OriginIR 全振幅态向量的寄存器边际；寄存器按声明顺序占据低位量子位。"""
    widths = [r.type.width for r in registers]
    index = [r.name for r in registers].index(name)
    offset = sum(widths[:index])
    mask = (1 << widths[index]) - 1
    result = {}
    for basis, amplitude in enumerate(vector):
        value = (basis >> offset) & mask
        result[value] = result.get(value, 0.0) + abs(amplitude) ** 2
    return result


def _readout_marginals(program, readout, *, use_originir=True):
    """各后端路径的读出寄存器边际分布；返回 {路径名: 分布}。"""
    marginals = {}
    for label, runner in (
        ("reference", reference),
        ("rir-pysparq", rir_pysparq),
        ("adapter-pysparq", adapter_pysparq),
    ):
        marginals[label] = _marginal_from_amplitudes(
            runner(program), program.main.registers, readout
        )
    if use_originir:
        marginals["originir-ext"] = _marginal_from_statevector(
            originir_ext(program), program.main.registers, readout
        )
    return marginals


def _dirichlet_kernel(precision, delta):
    """QPE 单峰闭式：D(δ) = sin²(π·2^p·δ) / (4^p·sin²(π·δ))，δ≡0 (mod 1) 时为 1。"""
    points = 1 << precision
    numerator = math.sin(math.pi * points * delta)
    denominator = math.sin(math.pi * delta)
    if abs(denominator) < 1e-12:
        return 1.0
    return (numerator / (points * denominator)) ** 2


def _qpe_closed_form(phase, precision):
    """本征相位 φ 的 QPE 相位寄存器闭式分布。"""
    points = 1 << precision
    return {y: _dirichlet_kernel(precision, y / points - phase) for y in range(points)}


def _qae_closed_form(amplitude, precision):
    """好状态概率 a 的 QAE 相位闭式分布：±θ/π 两个 Dirichlet 峰的等权混合。"""
    points = 1 << precision
    theta = math.asin(math.sqrt(amplitude))
    return {
        y: 0.5
        * (
            _dirichlet_kernel(precision, y / points - theta / math.pi)
            + _dirichlet_kernel(precision, y / points + theta / math.pi)
        )
        for y in range(points)
    }


def _distribution_checks(marginals, closed_form, extra_metrics):
    """对每条路径的读出分布做 TVD 与交叉对拍，返回指标字典与闭式峰值。"""
    worst_tvd = 0.0
    for distribution in marginals.values():
        worst_tvd = max(worst_tvd, tvd(distribution, closed_form))
    paths = list(marginals)
    pairwise = 0.0
    for i, left in enumerate(paths):
        for right in paths[i + 1 :]:
            pairwise = max(pairwise, tvd(marginals[left], marginals[right]))
    peak = max(closed_form, key=closed_form.get)
    metrics = {
        "max_tvd_vs_closed_form": worst_tvd,
        "max_pairwise_tvd": pairwise,
        "peak_value": peak,
        "peak_probability": max(d[peak] for d in marginals.values()),
        **extra_metrics,
    }
    return metrics, peak


# ---------------------------------------------------------------------------
# QPE
# ---------------------------------------------------------------------------


def _phase_unitary(phase):
    """单比特对角相位门 U|1> = exp(2πi·phase)|1>。"""
    label = str(phase).replace(".", "p").replace("-", "m")
    b = Builder(f"qpe_witness_{label}", {"target": Bits(1)})
    b.gate("phase", b["target"], 2 * math.pi * phase)
    return b.finish()


def _qpe_driver(operation, precision, prepare, name):
    """制备本征态后调用 QPE 的驱动程序。"""
    qpe = phase_estimation(operation, precision=precision)
    b = Builder(name, {r.name: r.type for r in qpe.module.registers})
    prepare(b)
    b.call(qpe, **{r.name: b[r.name] for r in qpe.module.registers})
    return b.finish().program()


def verify_qpe_ongrid(report):
    """栅格上的本征相位：读出应确定性地等于 2^p·φ，多精度扫描。"""
    for precision in range(2, 7):
        phase = 0.625  # φ = 5/8，对 p ≥ 3 恰在栅格上；p = 2 时取 φ' = 0.5 分支
        if precision == 2:
            phase = 0.5
        program = _qpe_driver(
            _phase_unitary(phase), precision, lambda b: b.x(b["target"]), f"qpe_grid_{precision}"
        )
        marginals = _readout_marginals(program, "phase")
        closed = _qpe_closed_form(phase, precision)
        expected = round((1 << precision) * phase)
        metrics, peak = _distribution_checks(marginals, closed, {"expected_peak": expected})
        report.case(
            f"qpe-phase-ongrid-w{precision}",
            paths=list(marginals),
            parameters={"phase": phase, "precision": precision},
            metrics=metrics,
            criterion="峰值 == 2^p·φ 且概率为 1，分布 TVD < 1e-9",
            passed=peak == expected
            and metrics["peak_probability"] > 1 - 1e-9
            and metrics["max_tvd_vs_closed_form"] < 1e-9,
        )


def verify_qpe_offgrid(report):
    """栅格外相位 φ = 0.3：峰值在最近栅格点，概率下界 4/π²，分布对照 Dirichlet 核。"""
    phase = 0.3
    for precision in (3, 4, 5, 6):
        program = _qpe_driver(
            _phase_unitary(phase), precision, lambda b: b.x(b["target"]), f"qpe_off_{precision}"
        )
        marginals = _readout_marginals(program, "phase")
        closed = _qpe_closed_form(phase, precision)
        expected = round((1 << precision) * phase) % (1 << precision)
        metrics, peak = _distribution_checks(marginals, closed, {"expected_peak": expected})
        report.case(
            f"qpe-phase-offgrid-w{precision}",
            paths=list(marginals),
            parameters={"phase": phase, "precision": precision},
            metrics=metrics,
            criterion="峰值在最近栅格点且概率 ≥ 4/π²，分布 TVD < 1e-9",
            passed=peak == expected
            and metrics["peak_probability"] >= 4 / math.pi**2 - 1e-9
            and metrics["max_tvd_vs_closed_form"] < 1e-9,
        )


def verify_qpe_fourier_eigenstate(report):
    """非对角幺正：add_const(1) 的 Fourier 本征态 |φ̃_j>，本征相位 (-j/4) mod 1。

    经典参考用 numpy 独立构造置换矩阵与 Fourier 态求本征值，不复用库内实现。
    """
    import numpy as np

    permutation = np.zeros((4, 4), dtype=complex)
    for x in range(4):
        permutation[(x + 1) % 4, x] = 1.0
    add = Builder("qpe_add_w2", {"w": UInt(2)})
    add.add_const(add["w"], 1)
    add_operation = add.finish()
    for j, precision in ((1, 3), (3, 3), (1, 4)):
        fourier_state = np.exp(2j * math.pi * j * np.arange(4) / 4) / 2
        eigenvalue = fourier_state.conjugate() @ permutation @ fourier_state
        phase = (math.atan2(eigenvalue.imag, eigenvalue.real) / (2 * math.pi)) % 1.0

        def prepare(b, j=j):
            for bit in range(2):
                if (j >> bit) & 1:
                    b.x(b["w"][bit])
            b.call(qft(2), target=b["w"].reinterpret("bits"))

        program = _qpe_driver(add_operation, precision, prepare, f"qpe_fourier_{j}_{precision}")
        marginals = _readout_marginals(program, "phase")
        closed = _qpe_closed_form(phase, precision)
        expected = round((1 << precision) * phase) % (1 << precision)
        metrics, peak = _distribution_checks(
            marginals, closed, {"expected_peak": expected, "eigenvalue_phase": phase}
        )
        report.case(
            f"qpe-addconst-fourier-j{j}-w{precision}",
            paths=list(marginals),
            parameters={"fourier_index": j, "precision": precision},
            metrics=metrics,
            criterion="numpy 本征相位确定性读出（峰值概率 1，TVD < 1e-9）",
            passed=peak == expected
            and metrics["peak_probability"] > 1 - 1e-9
            and metrics["max_tvd_vs_closed_form"] < 1e-9,
        )


def verify_qpe_unitary_matrix(report):
    """幺正层面：OriginIR 全幺正对照 numpy 独立组装的 QPE 幺正（含制备 X）。"""
    import numpy as np

    phase, precision = 0.625, 3
    points = 1 << precision
    program = _qpe_driver(
        _phase_unitary(phase), precision, lambda b: b.x(b["target"]), "qpe_matrix"
    )
    actual = originir_unitary(program)
    # numpy 独立组装：H^⊗p（相位寄存器）→ 受控 U^z → 逆 DFT；量子位 0 为 target
    u = np.diag([1.0, np.exp(2j * math.pi * phase)])
    identity2 = np.eye(2)
    hadamard_full = np.empty((points, points))
    for y in range(points):
        for z in range(points):
            hadamard_full[y, z] = (-1) ** bin(y & z).count("1") / math.sqrt(points)
    controlled = np.zeros((2 * points, 2 * points), dtype=complex)
    power = np.eye(2, dtype=complex)
    for z in range(points):
        controlled[2 * z : 2 * z + 2, 2 * z : 2 * z + 2] = power
        power = power @ u
    inverse_dft = np.empty((points, points), dtype=complex)
    for y in range(points):
        for z in range(points):
            inverse_dft[y, z] = np.exp(-2j * math.pi * y * z / points) / math.sqrt(points)
    expected = (
        np.kron(inverse_dft, identity2)
        @ controlled
        @ np.kron(hadamard_full, identity2)
        @ np.kron(np.eye(points), np.array([[0, 1], [1, 0]]))
    )
    error = float(np.abs(actual - expected).max())
    report.case(
        "qpe-unitary-matrix-w3",
        paths=["originir-ext+to_matrix"],
        parameters={"phase": phase, "precision": precision, "matrix_dim": actual.shape[0]},
        metrics={"max_error": error},
        criterion="线路幺正与 numpy 独立组装逐元素一致（max_error < 1e-12）",
        passed=error < 1e-12,
    )


# ---------------------------------------------------------------------------
# QAE
# ---------------------------------------------------------------------------


def _brassard_bound(amplitude, precision):
    """Brassard 等人的峰值解码误差界 |â - a| ≤ 2π√(a(1-a))/2^p + π²/4^p。"""
    points = 1 << precision
    return 2 * math.pi * math.sqrt(amplitude * (1 - amplitude)) / points + math.pi**2 / points**2


def _verify_qae(report, name, preparation, marked, amplitude, precision, exact):
    operation = amplitude_estimation(preparation, marked, precision=precision)
    program = operation.program()
    marginals = _readout_marginals(program, "phase")
    closed = _qae_closed_form(amplitude, precision)
    points = 1 << precision
    peak = max(max(m.values(), default=0.0) for m in marginals.values())
    # 峰值解码：对每条路径取各自 argmax，误差取最劣（镜像峰 y 与 2^p-y 解码相同）
    worst_decode = 0.0
    for distribution in marginals.values():
        mode = max(distribution, key=distribution.get)
        estimate = math.sin(math.pi * mode / points) ** 2
        worst_decode = max(worst_decode, abs(estimate - amplitude))
    metrics, _ = _distribution_checks(
        marginals,
        closed,
        {
            "decode_error": worst_decode,
            "brassard_bound": _brassard_bound(amplitude, precision),
            "mode_peak_probability": peak,
        },
    )
    if exact:
        criterion = "a 恰在解码栅格上：decode_error == 0 且分布 TVD < 1e-9"
        passed = worst_decode < 1e-12 and metrics["max_tvd_vs_closed_form"] < 1e-9
    else:
        criterion = "分布 TVD < 1e-9 且峰值解码误差 ≤ Brassard 界"
        passed = (
            metrics["max_tvd_vs_closed_form"] < 1e-9
            and worst_decode <= metrics["brassard_bound"]
        )
    report.case(
        name,
        paths=list(marginals),
        parameters={
            "amplitude": amplitude,
            "precision": precision,
            "marked": list(marked),
        },
        metrics=metrics,
        criterion=criterion,
        passed=passed,
    )


def verify_qae(report):
    # a = 1/2 恰在栅格上（θ/π = 1/4）：确定性双峰，解码精确
    for precision in (3, 4, 5):
        _verify_qae(
            report,
            f"qae-uniform-half-w{precision}",
            uniform_state(1),
            (1,),
            0.5,
            precision,
            exact=True,
        )
    # a = 3/8 不在栅格上：对照 Grover 双峰闭式与 Brassard 界
    for precision in (4, 5, 6):
        _verify_qae(
            report,
            f"qae-uniform-3over8-w{precision}",
            uniform_state(3),
            (0, 5, 7),
            3 / 8,
            precision,
            exact=False,
        )
    # 非均匀制备 a = 0.3：gate_state_prep 幅度已知，闭式独立
    _verify_qae(
        report,
        "qae-nonuniform-0p3-w5",
        gate_state_prep([math.sqrt(0.7), math.sqrt(0.3)]),
        (1,),
        0.3,
        5,
        exact=False,
    )


# ---------------------------------------------------------------------------
# Hadamard test 与 swap test
# ---------------------------------------------------------------------------


def _diagonal_unitary(angles, name):
    """对角幺正 U|x> = exp(i·angles[x])|x>（受控全局相位实现）。"""
    width = (len(angles) - 1).bit_length()
    b = Builder(name, {"target": Bits(width)})
    for value, angle in enumerate(angles):
        if angle:
            with b.control(b["target"], value):
                b.global_phase(angle)
    return b.finish()


def verify_hadamard(report):
    # 单比特相位门作用于 |1>：期望 e^{iθ}
    for angle in (0.6, -1.1):
        unitary = _phase_unitary(angle / (2 * math.pi))
        from pyqecclang.algorithms.input_model.oracles import basis_state

        for component, expected in (
            ("real", math.cos(angle)),
            ("imag", math.sin(angle)),
        ):
            program = hadamard_test(
                unitary, basis_state(1, 1), component=component
            ).program()
            marginals = _readout_marginals(program, "probe")
            worst = max(
                abs(dist.get(0, 0.0) - dist.get(1, 0.0) - expected)
                for dist in marginals.values()
            )
            report.case(
                f"hadamard-phase-{component}-{angle}",
                paths=list(marginals),
                parameters={"angle": angle, "component": component},
                metrics={"z_expectation_error": worst, "expected": expected},
                criterion="probe Z 期望 == " + ("cos θ" if component == "real" else "sin θ") + "（误差 < 1e-9）",
                passed=worst < 1e-9,
            )
    # 双比特对角幺正 + 复幅度制备：期望 Σ_x |ψ_x|² e^{iθ_x}
    angles = (0.35, -0.9, 1.7, 0.55)
    amplitudes = (0.5, 0.5j, 0.5, -0.5)
    unitary = _diagonal_unitary(angles, "hadamard_diag_2q")
    preparation = gate_state_prep(amplitudes)
    expected_value = sum(
        abs(a) ** 2 * math.cos(t) + 1j * abs(a) ** 2 * math.sin(t)
        for a, t in zip(amplitudes, angles, strict=True)
    )
    for component, expected in (
        ("real", expected_value.real),
        ("imag", expected_value.imag),
    ):
        program = hadamard_test(unitary, preparation, component=component).program()
        marginals = _readout_marginals(program, "probe")
        worst = max(
            abs(dist.get(0, 0.0) - dist.get(1, 0.0) - expected)
            for dist in marginals.values()
        )
        report.case(
            f"hadamard-diagonal-2q-{component}",
            paths=list(marginals),
            parameters={"angles": list(angles), "component": component},
            metrics={"z_expectation_error": worst, "expected": expected},
            criterion="probe Z 期望 == Σ_x |ψ_x|² e^{iθ_x} 的分量（误差 < 1e-9）",
            passed=worst < 1e-9,
        )


def verify_swap(report):
    # (幅度向量1, 幅度向量2)，重叠 F 由 numpy 独立计算
    pairs = [
        ("same", (0.6, 0.8), (0.6, 0.8)),
        ("orthogonal", (0.6, 0.8), (0.8, -0.6)),
        ("partial", (1.0, 0.0), (0.6, 0.8)),
        ("two-qubit", (1 / math.sqrt(2), 1 / math.sqrt(2), 0.0, 0.0), (0.0, 0.5, 0.5, 1 / math.sqrt(2))),
    ]
    for name, first, second in pairs:
        import numpy as np

        overlap = abs(np.vdot(np.asarray(first, dtype=complex), np.asarray(second, dtype=complex))) ** 2
        expected = float((1 + overlap) / 2)
        program = swap_test(gate_state_prep(first), gate_state_prep(second)).program()
        marginals = _readout_marginals(program, "probe")
        worst = max(abs(dist.get(0, 0.0) - expected) for dist in marginals.values())
        report.case(
            f"swap-test-{name}",
            paths=list(marginals),
            parameters={
                "first": [str(complex(a)) for a in first],
                "second": [str(complex(a)) for a in second],
            },
            metrics={
                "p0_error": worst,
                "expected_p0": expected,
                "fidelity_overlap": float(overlap),
            },
            criterion="P(probe=0) == (1+F)/2（误差 < 1e-9）",
            passed=worst < 1e-9,
        )


# ---------------------------------------------------------------------------
# Jordan 梯度
# ---------------------------------------------------------------------------


def _linear_phase_oracle(coefficients, grid_bits):
    """线性函数 f(x) = Σ_i c_i·x_i 的显式相位表（Jordan 缩放约定）。"""
    dimension = len(coefficients)
    points = 1 << grid_bits
    angles = []
    for value in range(1 << (dimension * grid_bits)):
        phase = 0.0
        for i, coefficient in enumerate(coefficients):
            chunk = (value >> (i * grid_bits)) & (points - 1)
            phase += coefficient * chunk / points
        angles.append(2 * math.pi * points * phase)
    return gate_phase_oracle(dimension * grid_bits, angles, phase_scale=points)


def _finite_difference(function, point, step):
    """中心有限差分梯度（独立经典参考）。"""
    return tuple(
        (
            function(*[x + (step if i == j else 0) for i, x in enumerate(point)])
            - function(*[x - (step if i == j else 0) for i, x in enumerate(point)])
        )
        / (2 * step)
        for j in range(len(point))
    )


def _verify_jordan_exact(
    report, name, oracle, dimension, grid_bits, gradient, paths_note=None, use_originir=True
):
    """线性函数：读出分布应确定性落在编码值上，解码精确等于梯度。"""
    operation = gradient_estimation(oracle, dimension=dimension, grid_bits=grid_bits)
    program = operation.program()
    marginals = _readout_marginals(program, "target", use_originir=use_originir)
    worst_decode = 0.0
    worst_peak = 1.0
    for distribution in marginals.values():
        mode = max(distribution, key=distribution.get)
        decoded = gradient_from_readout(mode, dimension=dimension, grid_bits=grid_bits)
        worst_decode = max(
            worst_decode, max(abs(g - e) for g, e in zip(decoded, gradient, strict=True))
        )
        worst_peak = min(worst_peak, distribution[mode])
    pairwise = 0.0
    labels = list(marginals)
    for i, left in enumerate(labels):
        for right in labels[i + 1 :]:
            pairwise = max(pairwise, tvd(marginals[left], marginals[right]))
    parameters = {"dimension": dimension, "grid_bits": grid_bits, "gradient": list(gradient)}
    if paths_note:
        parameters["paths_note"] = paths_note
    report.case(
        name,
        paths=labels,
        parameters=parameters,
        metrics={
            "max_decode_error": worst_decode,
            "min_peak_probability": worst_peak,
            "max_pairwise_tvd": pairwise,
        },
        criterion="读出确定性（峰概率 1）且解码精确等于梯度（误差 < 1e-12）",
        passed=worst_decode < 1e-12 and worst_peak > 1 - 1e-9,
    )


def verify_jordan_gate_linear(report):
    # 一维：含负分量与不同栅格，梯度分量取网格精确值
    for grid_bits, coefficient in ((3, 3 / 8), (4, -5 / 16), (5, 9 / 32)):
        _verify_jordan_exact(
            report,
            f"jordan-linear-gate-d1-w{grid_bits}",
            _linear_phase_oracle((coefficient,), grid_bits),
            1,
            grid_bits,
            (coefficient,),
        )
    # 多维：d = 2 与 d = 3
    _verify_jordan_exact(
        report,
        "jordan-linear-gate-d2-w4",
        _linear_phase_oracle((3 / 16, -2 / 16), 4),
        2,
        4,
        (3 / 16, -2 / 16),
    )
    _verify_jordan_exact(
        report,
        "jordan-linear-gate-d3-w3",
        _linear_phase_oracle((1 / 8, -3 / 8, 1 / 4), 3),
        3,
        3,
        (1 / 8, -3 / 8, 1 / 4),
    )


def verify_jordan_function_oracle(report):
    """mathfunc 相位 oracle 路径：定点算术 + 相位踢回后的精确线性读出。

    默认定点格式的工作位使总量子位远超 OriginIR 24 位预算，只走寄存器级路径。
    """
    # d = 1：f(x) = 0.25x
    oracle1 = function_phase_oracle(
        "def f(x0):\n return 0.25 * x0", dimension=1, grid_bits=4
    )
    program1 = gradient_estimation(oracle1, dimension=1, grid_bits=4).program()
    budget1 = sum(r.type.width for r in program1.main.registers) + workspace_table(program1)[
        program1.entry
    ]
    _verify_jordan_exact(
        report,
        "jordan-function-oracle-d1-w4",
        oracle1,
        1,
        4,
        (0.25,),
        paths_note=f"originir-ext 跳过：总量子位 {budget1} 超出 24 位预算",
        use_originir=False,
    )
    # d = 2：f(x0, x1) = 0.25·x0 - 0.125·x1，对照中心有限差分
    source_gradient = (0.25, -0.125)
    oracle2 = function_phase_oracle(
        "def f(x0, x1):\n return 0.25 * x0 - 0.125 * x1", dimension=2, grid_bits=3
    )
    program2 = gradient_estimation(oracle2, dimension=2, grid_bits=3).program()
    budget2 = sum(r.type.width for r in program2.main.registers) + workspace_table(program2)[
        program2.entry
    ]
    finite_difference = _finite_difference(
        lambda x0, x1: 0.25 * x0 - 0.125 * x1, (0.5, 0.5), 1 / 8
    )
    fd_gap = max(
        abs(a - b) for a, b in zip(finite_difference, source_gradient, strict=True)
    )
    _verify_jordan_exact(
        report,
        "jordan-function-oracle-d2-w3",
        oracle2,
        2,
        3,
        source_gradient,
        paths_note=f"originir-ext 跳过：总量子位 {budget2} 超出 24 位预算；"
        f"中心有限差分与真值差距 {fd_gap:.3e}（线性函数为 0）",
        use_originir=False,
    )


def verify_jordan_perturbed(report):
    """扰动线性 f(x) = a·x + x²/N²：峰位恒为真值，失败概率近似二次衰减。

    信息性指标同时给出中心有限差分（f 非线性时与线性系数有 O(1/N²) 差距）。
    """
    a_numerator, a_denominator = 3, 8
    success = []
    for grid_bits in (3, 4, 5):
        points = 1 << grid_bits
        angles = [
            2
            * math.pi
            * points
            * (
                a_numerator / a_denominator * value / points
                + (value / points) ** 2 / points**2
            )
            for value in range(points)
        ]
        oracle = gate_phase_oracle(grid_bits, angles, phase_scale=points)
        program = gradient_estimation(oracle, dimension=1, grid_bits=grid_bits).program()
        marginals = _readout_marginals(program, "target")
        exact = a_numerator * points // a_denominator
        peak_probability = min(d[exact] for d in marginals.values())
        mode_ok = all(max(d, key=d.get) == exact for d in marginals.values())
        success.append(peak_probability)
        # 中心有限差分：f'(1/2) = a + 1/N²（信息性）
        fd = _finite_difference(
            lambda x, points=points: a_numerator / a_denominator * x + x**2 / points**2,
            (0.5,),
            1 / points,
        )[0]
        report.case(
            f"jordan-perturbed-w{grid_bits}",
            paths=list(marginals),
            parameters={"dimension": 1, "grid_bits": grid_bits, "linear_coefficient": 3 / 8},
            metrics={
                "success_probability": peak_probability,
                "failure_probability": 1 - peak_probability,
                "peak_matches_truth": mode_ok,
                "finite_difference_at_center": fd,
                "peak_gradient": exact / points,
            },
            criterion="峰位 == round(N·a)，成功概率随网格位数单调上升",
            passed=mode_ok
            and (not success[:-1] or peak_probability > success[-2]),
        )
    ratios = [
        (1 - success[i + 1]) / (1 - success[i]) for i in range(len(success) - 1)
    ]
    report.case(
        "jordan-perturbed-decay-rate",
        paths=["derived"],
        parameters={"grid_bits_sequence": [3, 4, 5]},
        metrics={
            "success_probabilities": success,
            "failure_decay_ratios": ratios,
        },
        criterion="失败概率衰减率 q_{m+1}/q_m < 0.34（近似二次收敛）",
        passed=all(ratio < 0.34 for ratio in ratios),
    )


# ---------------------------------------------------------------------------
# 驱动
# ---------------------------------------------------------------------------


def run():
    report = Report(
        "estimation",
        "QPE/QAE/Hadamard/swap/Jordan 梯度的分布级与期望值级数值验证，"
        "对照 Dirichlet 核闭式、numpy 独立参考与有限差分。",
    )
    verify_qpe_ongrid(report)
    verify_qpe_offgrid(report)
    verify_qpe_fourier_eigenstate(report)
    verify_qpe_unitary_matrix(report)
    verify_qae(report)
    verify_hadamard(report)
    verify_swap(report)
    verify_jordan_gate_linear(report)
    verify_jordan_function_oracle(report)
    verify_jordan_perturbed(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
