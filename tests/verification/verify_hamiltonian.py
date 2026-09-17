"""Hamiltonian 模拟与 QSVT 标准变换的论文级数值验证。

覆盖 src/pyqecclang/algorithms/hamiltonian.py 与 src/pyqecclang/algorithms/qsvt.py：

- trotter_hamsim / hamiltonian_simulation：线路幺正（UniQC ``Circuit.to_matrix``）
  与乘积公式的独立经典矩阵逐元素对拍（单项精确、含恒等项全局相位、1–4 量子位）；
  以 scipy.linalg.expm 为独立 oracle 拟合一阶 Lie–Trotter 的收敛阶（误差 ∝ 1/steps）。
- qsp_phases / qsvt_sequence / qubitization_walk：相位合成往返（numpy 按文档约定
  独立实现的 QSP 响应）、电路零信号块与目标多项式在后端四条路径上逐点对拍、
  qubitization 的 W^n 零信号块 = T_n(x) 递推见证。
- qsvt_hamiltonian_simulation / qsvt_matrix_inversion / eigenstate_filter /
  gibbs_purification（density.py，QSVT 消费端）：零信号块与 scipy expm /
  numpy.linalg.inv / 解析过滤多项式 / 经典 Gibbs 态对拍；报告实现误差与方法误差、
  恢复逆矩阵的条件数对比、通阻带压制、迹距离与成功概率。
- oblivious_amplification：库算子 W = U·[R U† R U]（历史缺陷已修复：原迭代体
  [R U† R U] 缺收尾 U，零信号块退化为 2B†B − I、不执行放大）。修复后零信号块
  满足切比雪夫放大恒等式 ΠWΠ = B(4B†B − 3I)（精确验证）；脚本按文献序列
  U R U† R U 独立组装，断言库算子与之逐振幅一致（回归钉），V/2 → −V 放大
  语义见 docs/manual/algorithms/oblivious-amplification.md 的数值验证节。

经典 oracle 全部独立：numpy/scipy/math 闭式与矩阵例程均在本脚本内直接构造，
不复用被测模块的内部辅助函数。

运行：PYTHONPATH=src <含 pysparq+uniqc 的 python> tests/verification/verify_hamiltonian.py
"""

from __future__ import annotations

import cmath
import math
import random

import numpy as np
from harness import (
    Report,
    adapter_pysparq,
    amplitudes_to_statevector,
    basis_program,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    statevector_error,
    superposition_program,
)
from scipy.linalg import expm as scipy_expm
from scipy.special import jv as bessel_j

from pyqecclang.algorithms.block_encoding import matrix_pauli_encoding, reflect_zero
from pyqecclang.algorithms.density import gibbs_purification
from pyqecclang.algorithms.hamiltonian import (
    EncodedOperator,
    PauliHamiltonian,
    hamiltonian_simulation,
    trotter_hamsim,
)
from pyqecclang.algorithms.operators import block_encoding, linear_combination, zero
from pyqecclang.algorithms.oracles import invoke, resources_for
from pyqecclang.algorithms.qsvt import (
    eigenstate_filter,
    qsp_phases,
    qsvt_hamiltonian_simulation,
    qsvt_matrix_inversion,
)
from pyqecclang.algorithms.transforms import (
    oblivious_amplification,
    qsvt_sequence,
    qubitization_walk,
)
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits
from pyqecclang.infrastructure.layout import workspace_table

# ---------------------------------------------------------------------------
# 独立经典 oracle：Pauli 矩阵、乘积公式、QSP 响应、目标多项式。
# ---------------------------------------------------------------------------

_PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_matrix(word):
    """Pauli 词的稠密矩阵；word[i] 作用于第 i 量子位（权重 2^i，与 RIR 寄存器一致）。"""
    result = _PAULI[word[0]]
    for letter in word[1:]:
        result = np.kron(_PAULI[letter], result)
    return result


def hamiltonian_matrix(terms):
    return sum(coefficient * pauli_matrix(word) for coefficient, word in terms)


def term_unitary(coefficient, word, time):
    """exp(−i·c·t·P) 的闭式：P² = I 故为 cos(ct)·I − i·sin(ct)·P。"""
    theta = coefficient * time
    dim = 1 << len(word)
    return math.cos(theta) * np.eye(dim) - 1j * math.sin(theta) * pauli_matrix(word)


def product_formula(terms, final_time, steps):
    """一阶 Lie–Trotter 乘积公式的经典矩阵（项序与电路作用顺序一致）。"""
    dim = 1 << len(terms[0][1])
    step = np.eye(dim, dtype=complex)
    for coefficient, word in terms:
        step = term_unitary(coefficient, word, final_time / steps) @ step
    return np.linalg.matrix_power(step, steps)


def qsp_response_ref(x, phases):
    """按模块 docstring 的反射约定用 numpy 独立实现的 QSP 响应。

    p(x) = [S(φ_0) W(x) S(φ_1) W(x) … S(φ_d)]_00，S(φ) = diag(e^{iφ}, e^{−iφ})，
    W(x) = [[x, s], [s, −x]]。与被测 qsvt.qsp_response 互为独立实现。
    """
    s = math.sqrt(max(0.0, 1.0 - x * x))
    walk = np.array([[x, s], [s, -x]])
    total = np.eye(2, dtype=complex)
    for i, phi in enumerate(phases):
        if i:
            total = walk @ total
        total = np.diag([cmath.exp(1j * phi), cmath.exp(-1j * phi)]) @ total
    return total[0, 0]


def chebyshev_coeffs(n):
    """T_n 的升幂单项式系数（三重递推，独立实现）。"""
    if n == 0:
        return (1.0,)
    if n == 1:
        return (0.0, 1.0)
    a, b = (1.0,), (0.0, 1.0)
    for _ in range(2, n + 1):
        a, b = b, tuple(
            2.0 * (b[i - 1] if i else 0.0) - (a[i] if i < len(a) else 0.0)
            for i in range(len(b) + 1)
        )
    return b


def poly_eval(coeffs, x):
    result = 0.0
    for c in reversed(coeffs):
        result = result * x + c
    return result


def matrix_polynomial(coeffs, mat):
    """升幂系数多项式的矩阵 Horner 求值（独立 oracle）。"""
    result = np.zeros_like(mat)
    for c in reversed(coeffs):
        result = result @ mat + c * np.eye(mat.shape[0])
    return result


def chebyshev_value(n, x):
    """T_n(x)，|x| > 1 时用双曲延拓。"""
    if abs(x) <= 1.0:
        return math.cos(n * math.acos(x))
    return math.cosh(n * math.acosh(abs(x))) * (1.0 if x > 0 or n % 2 == 0 else -1.0)


def filter_polynomial(gap, degree, x):
    """Lin–Tong 过滤多项式 f(x) = T_d(g(x²))/T_d(r) 的独立求值。"""
    r = (1.0 + gap**2) / (1.0 - gap**2)
    g = 2.0 * (x * x - gap**2) / (1.0 - gap**2) - 1.0
    return chebyshev_value(degree, g) / math.cosh(degree * math.acosh(r))


def jacobi_anger_tail(t, kmax):
    """e^{itx} 的 Jacobi–Anger 截断尾部上界 2·Σ_{j>K} |J_j(t)|（scipy 独立求值）。"""
    return float(2.0 * sum(abs(bessel_j(j, t)) for j in range(kmax + 1, kmax + 200)))


def gibbs_reference(matrix, beta):
    """经典 Gibbs 态 e^{−βH}/Tr（scipy expm，独立 oracle）。"""
    weights = scipy_expm(-beta * np.asarray(matrix, dtype=complex))
    return weights / np.trace(weights)


def trace_distance(rho, sigma):
    return 0.5 * float(np.abs(np.linalg.eigvalsh(rho - sigma)).sum())


# ---------------------------------------------------------------------------
# 后端执行辅助：列读出零信号块、预算预判、统一条目解码。
# ---------------------------------------------------------------------------


def budget_qubits(program):
    """寄存器 + 工作区总量（OriginIR 态向量预算 24 量子位）。"""
    return sum(r.type.width for r in program.main.registers) + workspace_table(program)[
        program.entry
    ]


def embed_unitary(circuit, dim):
    """把 UniQC to_matrix 的幺正嵌入到 dim 维。

    UniQC 按最高被引用量子位定矩阵维数，会裁掉尾部没有任何门的量子位
    （如 Pauli 词尾随的 I）；被裁部分作用为恒等，故右 Kronecker 恒等因子。
    本验证中所有程序的已用量子位均为前缀，嵌入是精确的。
    """
    current = circuit.shape[0]
    if current == dim:
        return circuit
    if current > dim or dim % current:
        raise AssertionError(f"幺正维数 {current} 无法嵌入 {dim}")
    return np.kron(np.eye(dim // current, dtype=complex), circuit)


def _register_widths(operation):
    return [r.type.width for r in operation.module.registers]


def _decode_index(index, widths):
    values = []
    for w in widths:
        values.append(index & ((1 << w) - 1))
        index >>= w
    return tuple(values)


def _iter_entries(state, widths):
    """把 dict 稀疏态与 originir 稠密态向量统一为 (寄存器值元组, 振幅) 迭代。"""
    if isinstance(state, dict):
        yield from state.items()
    else:
        for index, amplitude in enumerate(state):
            if amplitude:
                yield _decode_index(index, widths), amplitude


def run_column(operation, runner, column):
    """X 门制备 target 基态列后执行（四条路径统一走 driver，起点均为 |0>）。"""
    return runner(basis_program(operation, {"target": column}, name=f"column_{column}"))


def zero_signal_block(operation, runner, dim):
    """逐列读出 signal == 0 分支的 target 块与该列成功概率。"""
    widths = _register_widths(operation)
    block = np.zeros((dim, dim), dtype=complex)
    success = []
    for column in range(dim):
        state = run_column(operation, runner, column)
        total = 0.0
        for key, amplitude in _iter_entries(state, widths):
            if key[1] == 0:
                block[key[0], column] = amplitude
                total += abs(amplitude) ** 2
        success.append(total)
    return block, success


# ---------------------------------------------------------------------------
# A. Trotter 线路语义：幺正 vs 乘积公式独立矩阵。
# ---------------------------------------------------------------------------


def verify_trotter_single_term(report):
    """单项 Pauli 演化（steps=1 无 Trotter 误差）：幺正与闭式逐元素一致。"""
    instances = [
        (0.7, "X", 0.4),
        (-0.3, "Y", 1.1),
        (0.55, "Z", -0.8),
        (0.25, "I", 0.9),
        (0.9, "ZZ", 0.6),
        (-0.4, "XY", 1.2),
        (0.35, "YZI", -0.5),
        (-0.8, "ZXI", 0.7),
        (0.45, "XYZI", 0.3),
        (-0.6, "IYXZ", -1.0),
    ]
    worst = 0.0
    per_word = {}
    for coefficient, word, time in instances:
        circuit = embed_unitary(
            originir_unitary(trotter_hamsim([(coefficient, word)], time, steps=1).program()),
            1 << len(word),
        )
        error = float(np.abs(circuit - term_unitary(coefficient, word, time)).max())
        per_word[f"{coefficient}·{word}·t={time}"] = error
        worst = max(worst, error)
    report.case(
        "trotter-single-term-exact",
        paths=["originir-ext+to_matrix"],
        parameters={"instances": len(instances)},
        metrics={"max_error": worst, "per_word": per_word},
        criterion="单项演化幺正等于 cos(ct)I − i·sin(ct)P（max_error < 1e-12）",
        passed=worst < 1e-12,
    )


def verify_trotter_product_formula(report):
    """多项非对易分解：幺正与乘积公式经典矩阵逐步数/时间逐元素一致。"""
    instances = [
        (
            "1q",
            [(0.3, "I"), (0.7, "X"), (-0.2, "Z")],
            [1, 3],
            [0.4, -0.7],
        ),
        (
            "2q",
            [(0.9, "ZZ"), (0.6, "XI"), (0.45, "IX"), (0.25, "II")],
            [1, 2, 5],
            [0.4, 1.3],
        ),
        (
            "3q",
            [(0.9, "ZZI"), (0.9, "IZZ"), (0.6, "XII"), (0.6, "IXI"), (0.6, "IIX")],
            [2, 3],
            [0.4, -0.7],
        ),
        (
            "4q",
            [(0.4, "ZIZI"), (-0.7, "XIII"), (0.2, "IXII"), (0.35, "IIXI"), (-0.15, "IIIX")],
            [3],
            [0.9],
        ),
    ]
    for label, terms, steps_list, times in instances:
        worst = 0.0
        gap = 0.0
        exact = scipy_expm(-1j * hamiltonian_matrix(terms) * times[0])
        for steps in steps_list:
            for time in times:
                program = trotter_hamsim(terms, time, steps=steps).program()
                circuit = embed_unitary(originir_unitary(program), 1 << len(terms[0][1]))
                expected = product_formula(terms, time, steps)
                worst = max(worst, float(np.abs(circuit - expected).max()))
                if time == times[0]:
                    gap = max(gap, float(np.linalg.norm(circuit - exact, 2)))
        report.case(
            f"trotter-product-formula-{label}",
            paths=["originir-ext+to_matrix"],
            parameters={
                "terms": len(terms),
                "steps": steps_list,
                "times": times,
                "qubits": budget_qubits(program),
            },
            metrics={"max_error": worst, "trotter_gap_vs_expm": gap},
            criterion="线路幺正等于乘积公式矩阵（max_error < 1e-12）",
            passed=worst < 1e-12,
        )


def verify_trotter_superposition_cross(report):
    """3 量子位 TFIM 在均匀叠加输入下的四路径对拍（一次运行覆盖全部基态）。"""
    terms = [(0.9, "ZZI"), (0.9, "IZZ"), (0.6, "XII"), (0.6, "IXI"), (0.6, "IIX")]
    operation = trotter_hamsim(terms, 0.9, steps=3)
    program = superposition_program(operation, ["target"])
    expected = product_formula(terms, 0.9, 3) @ np.full(8, 1.0 / math.sqrt(8))
    widths = [3, 0]
    deviations = {}
    for name, runner in (
        ("reference", reference),
        ("rir-pysparq", rir_pysparq),
        ("adapter-pysparq", adapter_pysparq),
        ("originir-ext", originir_ext),
    ):
        state = runner(program)
        vector = (
            amplitudes_to_statevector(state, widths) if isinstance(state, dict) else list(state)
        )
        deviations[name] = float(statevector_error(vector, list(expected)))
    worst = max(deviations.values())
    report.case(
        "trotter-superposition-cross-3q",
        paths=list(deviations),
        parameters={"steps": 3, "time": 0.9, "qubits": budget_qubits(program)},
        metrics={"max_deviation": worst, "per_path": deviations},
        criterion="四路径叠加态与乘积公式逐振幅一致（max_deviation < 1e-9）",
        passed=worst < 1e-9,
    )


def _convergence_scan(terms, time, grid):
    """对固定 t 扫步数 r，返回谱范数误差序列（幺正路径，vs scipy expm）。"""
    exact = scipy_expm(-1j * hamiltonian_matrix(terms) * time)
    errors = []
    for steps in grid:
        circuit = embed_unitary(
            originir_unitary(trotter_hamsim(terms, time, steps=steps).program()),
            1 << len(terms[0][1]),
        )
        errors.append(float(np.linalg.norm(circuit - exact, 2)))
    return errors


def verify_trotter_convergence(report):
    """一阶 Lie–Trotter 收敛阶拟合：log–log 斜率应接近 −1（误差 ∝ t²/r）。"""
    instances = [
        ("2q", [(0.9, "ZZ"), (0.6, "XI"), (0.45, "IX")], 1.0, [1, 2, 4, 8, 16, 32, 64]),
        (
            "3q",
            [(0.9, "ZZI"), (0.9, "IZZ"), (0.6, "XII"), (0.6, "IXI"), (0.6, "IIX")],
            1.0,
            [1, 2, 4, 8, 16, 32],
        ),
    ]
    for label, terms, time, grid in instances:
        errors = _convergence_scan(terms, time, grid)
        fit_r = np.log([r for r in grid if r >= 2])
        fit_e = np.log([e for r, e in zip(grid, errors, strict=True) if r >= 2])
        slope = float(np.polyfit(fit_r, fit_e, 1)[0])
        report.case(
            f"trotter-convergence-{label}",
            paths=["originir-ext+to_matrix"],
            parameters={"time": time, "steps_grid": grid, "terms": len(terms)},
            metrics={
                "errors": {str(r): e for r, e in zip(grid, errors, strict=True)},
                "fitted_order": -slope,
                "error_at_max_steps": errors[-1],
            },
            criterion="拟合收敛阶在 [0.8, 1.3]（一阶乘积公式理论值 1）",
            passed=0.8 <= -slope <= 1.3,
        )


# ---------------------------------------------------------------------------
# B. hamiltonian_simulation 协议层（Trotter 路由与 QSP 注入）。
# ---------------------------------------------------------------------------


def verify_protocol_trotter(report):
    """PauliHamiltonian 经协议路由到 Trotter：幺正精确等于乘积公式。"""
    terms = [(0.9, "ZZ"), (0.6, "XI"), (0.45, "IX")]
    operator = PauliHamiltonian(terms)
    exact = scipy_expm(-1j * hamiltonian_matrix(terms))
    worst_formula = 0.0
    gaps = {}
    for steps in (2, 16):
        result = hamiltonian_simulation(operator, 1.0, steps=steps)
        circuit = embed_unitary(originir_unitary(result.operation.program()), 4)
        worst_formula = max(
            worst_formula, float(np.abs(circuit - product_formula(terms, 1.0, steps)).max())
        )
        gaps[str(steps)] = float(np.linalg.norm(circuit - exact, 2))
    report.case(
        "hamsim-protocol-trotter-2q",
        paths=["originir-ext+to_matrix"],
        parameters={"time": 1.0, "steps": [2, 16], "alpha": 1.0},
        metrics={
            "max_error_vs_product_formula": worst_formula,
            "trotter_gap_vs_expm": gaps,
        },
        criterion="协议组央幺正等于乘积公式（< 1e-12）且 r=16 的 expm 偏差小于 r=2",
        passed=worst_formula < 1e-12 and gaps["16"] < gaps["2"],
    )


def verify_protocol_qsp_injection(report):
    """EncodedOperator 经 auto 路由到注入的 QSVT 内核：零信号块 vs scipy expm。"""
    terms = [(0.6, "Z"), (0.4, "X")]
    ham = PauliHamiltonian(terms)
    operator = EncodedOperator(ham.block_encoding(), True)
    time = 0.5

    def injected_qsp(be, evolution_time):
        # e^{iτA/α} 中取 τ = −t·α 即得 e^{−itH}
        return qsvt_hamiltonian_simulation(be, -evolution_time * be.alpha, error=0.02)

    result = hamiltonian_simulation(operator, time, method="auto", qsp=injected_qsp)
    scale = dict(result.operation.module.attributes)["sim_scale"]
    expected = scipy_expm(-1j * hamiltonian_matrix(terms) * time) / scale
    worst = 0.0
    success = []
    for runner in (reference, rir_pysparq, originir_ext):
        block, probs = zero_signal_block(result.operation, runner, 2)
        worst = max(worst, float(np.abs(block - expected).max()))
        success = probs
    report.case(
        "hamsim-protocol-qsp-injection",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"time": time, "qsp_error": 0.02, "sim_scale": scale},
        metrics={"max_error": worst, "success_probability": success},
        criterion="零信号块 ≈ e^{−itH}/sim_scale（max_error < 3e-2，含 QSVT 逼近预算）",
        passed=worst < 3e-2,
    )


# ---------------------------------------------------------------------------
# C. QSP 相位合成、QSVT 序列约定与 qubitization。
# ---------------------------------------------------------------------------


def verify_phase_synthesis(report):
    """qsp_phases 往返：合成相位经独立响应求值回到目标多项式（401 点网格）。"""
    grid = [-1.0 + 2.0 * i / 400 for i in range(401)]
    targets = {f"chebyshev-T{n}": (chebyshev_coeffs(n), None) for n in range(1, 7)}
    targets["explicit-imag"] = ((0.0, 0.5), (0.0, math.sqrt(0.75)))
    errors = {}
    for name, (coeffs, imag) in targets.items():
        phases = qsp_phases(coeffs, imag=imag)
        errors[name] = float(
            max(
                abs(
                    qsp_response_ref(x, phases)
                    - poly_eval(coeffs, x)
                    - (1j * poly_eval(imag, x) if imag else 0.0)
                )
                for x in grid
            )
        )
    worst = max(errors.values())
    report.case(
        "qsp-phase-synthesis-roundtrip",
        paths=["numpy-independent-response"],
        parameters={"grid_points": len(grid), "targets": len(targets)},
        metrics={"max_error": worst, "per_target": errors},
        criterion="合成相位往返误差 < 1e-8（模块声称典型 1e-9 量级）",
        passed=worst < 1e-8,
    )


def verify_qsvt_sequence_convention(report):
    """随机相位下 qsvt_sequence 电路零信号块与独立响应逐点一致（四路径）。"""
    rng = random.Random(20240901)
    phases = tuple(rng.uniform(-math.pi, math.pi) for _ in range(6))
    matrix = [[0.85, 0.0], [0.0, -0.55]]
    be = matrix_pauli_encoding(matrix)
    operation = qsvt_sequence(be, phases)
    xs = [matrix[0][0] / be.alpha, matrix[1][1] / be.alpha]
    deviations = {}
    for name, runner in (
        ("reference", reference),
        ("rir-pysparq", rir_pysparq),
        ("adapter-pysparq", adapter_pysparq),
        ("originir-ext", originir_ext),
    ):
        block, _ = zero_signal_block(operation, runner, 2)
        deviations[name] = float(
            max(abs(block[col, col] - qsp_response_ref(xs[col], phases)) for col in (0, 1))
        )
    worst = max(deviations.values())
    report.case(
        "qsvt-sequence-convention-1q",
        paths=list(deviations),
        parameters={"phases": len(phases), "spectral_points": xs},
        metrics={"max_deviation": worst, "per_path": deviations},
        criterion="四路径零信号块与独立 QSP 响应一致（max_deviation < 1e-9）",
        passed=worst < 1e-9,
    )


def verify_qsvt_sequence_matrix_block(report):
    """2 量子位非对角 BE：qsvt_sequence 的完整零信号块等于矩阵多项式 T_4(A/α)。"""
    terms = [(0.45, "ZZ"), (0.25, "XI"), (-0.15, "IZ"), (0.10, "II")]
    matrix = hamiltonian_matrix(terms)
    be = matrix_pauli_encoding(matrix)
    coeffs = chebyshev_coeffs(4)
    operation = qsvt_sequence(be, qsp_phases(coeffs))
    expected = matrix_polynomial(coeffs, matrix / be.alpha)
    worst = 0.0
    used = []
    for name, runner in (
        ("reference", reference),
        ("rir-pysparq", rir_pysparq),
        ("originir-ext", originir_ext),
    ):
        block, _ = zero_signal_block(operation, runner, 4)
        worst = max(worst, float(np.abs(block - expected).max()))
        used.append(name)
    report.case(
        "qsvt-sequence-matrix-block-2q",
        paths=used,
        parameters={"degree": 4, "alpha": be.alpha, "qubits": 4},
        metrics={"max_error": worst},
        criterion="零信号块等于矩阵多项式 T_4(A/α)（max_error < 1e-9）",
        passed=worst < 1e-9,
    )


def verify_qubitization(report):
    """qubitization_walk：W^n 的零信号块等于 T_n(x)（n = 1..5，双本征态）。"""
    matrix = [[0.85, 0.0], [0.0, -0.55]]
    be = matrix_pauli_encoding(matrix)
    walk = qubitization_walk(be)
    xs = [matrix[0][0] / be.alpha, matrix[1][1] / be.alpha]
    worst = 0.0
    for n in (1, 2, 3, 5):
        for col in (0, 1):
            b = Builder(
                f"walk_pow_{n}_{col}",
                {r.name: r.type for r in walk.module.registers},
                resources_for(("w", walk)),
            )
            for bit in range(b["target"].width):
                if (col >> bit) & 1:
                    b.x(b["target"][bit])
            with b.repeat(n):
                invoke(b, walk, "w", target=b["target"], signal=b["signal"])
            program = b.finish().program()
            widths = [r.type.width for r in program.main.registers]
            expected = chebyshev_value(n, max(-1.0, min(1.0, xs[col])))
            for runner in (reference, rir_pysparq, originir_ext):
                state = runner(program)
                amplitude = next(
                    (a for k, a in _iter_entries(state, widths) if k[0] == col and k[1] == 0),
                    0j,
                )
                worst = max(worst, abs(amplitude - expected))
    report.case(
        "qubitization-chebyshev-recurrence",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"powers": [1, 2, 3, 5], "spectral_points": xs},
        metrics={"max_error": worst},
        criterion="W^n 零信号块等于 T_n(x)（max_error < 1e-12）",
        passed=worst < 1e-12,
    )


# ---------------------------------------------------------------------------
# D. QSVT 标准变换：HamSim、矩阵求逆、特征态过滤。
# ---------------------------------------------------------------------------


def verify_qsvt_hamsim(report):
    """qsvt_hamiltonian_simulation：零信号块 vs scipy expm；方法误差与实现误差分列。"""
    matrix = np.array([[0.5, 0.2], [0.2, -0.3]])
    be = matrix_pauli_encoding(matrix)
    for time in (0.7, 2.0):
        result = qsvt_hamiltonian_simulation(be, time, error=0.01)
        attrs = dict(result.operation.module.attributes)
        scale, degree = attrs["sim_scale"], attrs["qsp_degree"]
        expected = scipy_expm(1j * time * matrix / be.alpha) / scale
        worst = 0.0
        for runner in (reference, rir_pysparq, originir_ext):
            block, _ = zero_signal_block(result.operation, runner, 2)
            worst = max(worst, float(np.abs(block - expected).max()))
        tail = jacobi_anger_tail(time, degree)
        report.case(
            f"qsvt-hamsim-t{time}",
            paths=["reference", "rir-pysparq", "originir-ext"],
            parameters={"time": time, "error": 0.01, "qsp_degree": degree, "sim_scale": scale},
            metrics={
                "impl_error": worst,
                "method_tail_bound": tail,
                "total_error_bound": worst + tail,
            },
            criterion="块·sim_scale 与 e^{itA/α} 的总偏差 < 2×error（0.02）",
            passed=worst + tail < 0.02,
        )


def _inversion_case(report, name, matrix, kappa, error, runners):
    be = matrix_pauli_encoding(matrix)
    result = qsvt_matrix_inversion(be, kappa, error=error)
    attrs = dict(result.operation.module.attributes)
    scale, degree = attrs["inverse_scale"], attrs["qsp_degree"]
    matrix = np.asarray(matrix, dtype=complex)
    exact_inverse = np.linalg.inv(matrix)
    expected = scale * exact_inverse
    dim = matrix.shape[0]
    worst = 0.0
    success = []
    for runner in runners:
        block, probs = zero_signal_block(result.operation, runner, dim)
        worst = max(worst, float(np.abs(block - expected).max()))
        success = probs
    recovered = np.linalg.svd(block / scale, compute_uv=False)
    kappa_recovered = float(recovered[0] / recovered[-1])
    relative = float(
        np.linalg.norm(block / scale - exact_inverse, 2) / np.linalg.norm(exact_inverse, 2)
    )
    # 每个奇异值带 ≤ error 的相对偏差，故恢复条件数落在 κ·[(1−e)/(1+e), (1+e)/(1−e)]
    kappa_exact = float(np.linalg.cond(matrix))
    band = (kappa_exact * (1 - error) / (1 + error), kappa_exact * (1 + error) / (1 - error))
    report.case(
        name,
        paths=[r.__name__ for r in runners],
        parameters={"kappa": kappa, "error": error, "qsp_degree": degree, "inverse_scale": scale},
        metrics={
            "max_error": worst,
            "relative_error": relative,
            "kappa_recovered": kappa_recovered,
            "kappa_exact": kappa_exact,
            "kappa_band": list(band),
            "success_probability": success,
        },
        criterion="恢复逆的相对谱误差 ≤ 1.05×error 且恢复条件数落在 κ·(1±e)/(1∓e) 区间",
        passed=relative <= 1.05 * error and band[0] <= kappa_recovered <= band[1],
    )


def verify_qsvt_inversion(report):
    """qsvt_matrix_inversion：块/缩放 vs numpy 精确逆；条件数对比；error 扫描。"""
    matrix_k2 = [[0.6, -0.2], [-0.2, 0.6]]  # 本征值 0.4/0.8，κ = 2
    for error in (0.15, 0.10):
        _inversion_case(
            report,
            f"qsvt-inversion-kappa2-error{error}",
            matrix_k2,
            2.0,
            error,
            (reference, rir_pysparq),
        )
    _inversion_case(
        report,
        "qsvt-inversion-kappa3-1q",
        [[0.4, -0.2], [-0.2, 0.4]],  # 本征值 0.2/0.6，κ = 3
        3.0,
        0.4,
        (reference, rir_pysparq, originir_ext),
    )
    _inversion_case(
        report,
        "qsvt-inversion-kappa3-2q",
        [[0.3, 0, 0, 0], [0, 0.4, 0, 0], [0, 0, 0.6, 0], [0, 0, 0, 0.9]],
        3.0,
        0.4,
        (reference, rir_pysparq, originir_ext),
    )


def verify_eigenstate_filter(report):
    """eigenstate_filter：零信号块对角幅度 vs 解析过滤多项式；通阻带压制。"""
    # 1 量子位带中心平移（docs 案例）：谱变量平移后一通一阻。
    be = matrix_pauli_encoding([[0.05, 0.0], [0.0, 0.5]])
    gap, degree, center = 0.2, 8, 0.1
    flt = eigenstate_filter(be, gap, degree, center=center)
    attrs = dict(flt.operation.module.attributes)
    alpha_shifted = be.alpha + abs(center)
    xs = [(0.05 - center) / alpha_shifted, (0.5 - center) / alpha_shifted]
    expected = [filter_polynomial(gap, degree, x) for x in xs]
    worst = 0.0
    for runner in (reference, rir_pysparq, originir_ext):
        block, _ = zero_signal_block(flt.operation, runner, 2)
        worst = max(worst, float(max(abs(block[c, c] - expected[c]) for c in (0, 1))))
    report.case(
        "eigenstate-filter-centered-1q",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"gap": gap, "degree": degree, "center": center, "shifted_x": xs},
        metrics={
            "max_error": worst,
            "passband_amplitude": expected[0],
            "stopband_amplitude": expected[1],
            "suppression_attribute": attrs["suppression"],
        },
        criterion="块对角幅度等于解析过滤多项式（max_error < 1e-5，合成度数 16 的剥离噪声量级）"
        "且阻带 ≤ suppression",
        passed=worst < 1e-5 and abs(expected[1]) <= attrs["suppression"] + 1e-9,
    )
    # 2 量子位四个本征值：两个通带两个阻带。
    matrix = [[0.02, 0, 0, 0], [0, -0.03, 0, 0], [0, 0, 0.4, 0], [0, 0, 0, 0.55]]
    be2 = matrix_pauli_encoding(matrix)
    gap2, degree2 = 0.2, 6
    flt2 = eigenstate_filter(be2, gap2, degree2)
    attrs2 = dict(flt2.operation.module.attributes)
    eigenvalues = [0.02, -0.03, 0.4, 0.55]
    xs2 = [v / be2.alpha for v in eigenvalues]
    expected2 = [filter_polynomial(gap2, degree2, x) for x in xs2]
    worst2 = 0.0
    for runner in (reference, rir_pysparq, originir_ext):
        block, _ = zero_signal_block(flt2.operation, runner, 4)
        worst2 = max(worst2, float(max(abs(block[c, c] - expected2[c]) for c in range(4))))
    stopband = max(abs(v) for v in expected2[2:])
    report.case(
        "eigenstate-filter-2q",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"gap": gap2, "degree": degree2, "alpha": be2.alpha, "spectral_x": xs2},
        metrics={
            "max_error": worst2,
            "passband_min": min(abs(v) for v in expected2[:2]),
            "stopband_max": stopband,
            "suppression_attribute": attrs2["suppression"],
        },
        criterion="全部四个本征态幅度与解析多项式一致（max_error < 1e-5，剥离自检容差）"
        "且阻带 ≤ suppression",
        passed=worst2 < 1e-5 and stopband <= attrs2["suppression"] + 1e-9,
    )


# ---------------------------------------------------------------------------
# E. Oblivious 振幅放大与 Gibbs 态制备。
# ---------------------------------------------------------------------------


def verify_oaa_iterate_identity(report):
    """库算子 W = U·[R U† R U] 的代数恒等式：ΠWΠ = B(4B†B − 3I)（对任意 BE 精确）。

    历史缺陷已修复：原迭代体 [R U† R U] 的零信号块退化为 2B†B − I，不执行
    文献中的 OAA 放大（见组报告与 oaa-standard-sequence 案例）；修复后的
    算子补上收尾 U，零信号块为切比雪夫放大 B(4B†B − 3I)——对零信号块为
    V/2 的输入恰为 −V。本案例把库函数实际实现的语义钉死到机器精度。
    """
    matrix = np.array([[0.6, -0.2], [-0.2, 0.6]])
    be = matrix_pauli_encoding(matrix)
    source_block = matrix / be.alpha
    expected = source_block @ (4.0 * source_block.conj().T @ source_block - 3.0 * np.eye(2))
    oaa = oblivious_amplification(be, iterations=1)
    worst = 0.0
    for runner in (reference, rir_pysparq, originir_ext):
        block, _ = zero_signal_block(oaa, runner, 2)
        worst = max(worst, float(np.abs(block - expected).max()))
    report.case(
        "oaa-iterate-block-identity",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"alpha": be.alpha, "iterations": 1},
        metrics={"max_error": worst},
        criterion="ΠWΠ 等于 B(4B†B − 3I)（max_error < 1e-12）",
        passed=worst < 1e-12,
    )


def verify_oaa_standard_sequence(report):
    """标准 OAA 序列 U R U† R U：V/2 块编码经一次迭代恢复 −V（文献语义）。

    序列用库公开组件（invoke + reflect_zero + adjoint）在脚本内组装；修复后
    的库 oblivious_amplification 与之逐振幅一致——本案例同时断言两者相等，
    作为"库 = 标准三查询序列"的回归钉。
    """
    vb = Builder("oaa_v", {"target": Bits(1), "signal": Bits(0)})
    vb.ry(vb["target"][0], 0.9)
    vb.rz(vb["target"][0], 0.4)
    vbe = block_encoding(vb.finish())
    theta_y, theta_z = 0.9, 0.4
    unitary_v = np.array(
        [[cmath.exp(-1j * theta_z / 2), 0], [0, cmath.exp(1j * theta_z / 2)]]
    ) @ np.array(
        [
            [math.cos(theta_y / 2), -math.sin(theta_y / 2)],
            [math.sin(theta_y / 2), math.cos(theta_y / 2)],
        ]
    )
    be = linear_combination(1.0, vbe, 1.0, zero(1))  # 零信号块 = V/2，α = 2
    b = Builder(
        "oaa_standard",
        {"target": Bits(1), "signal": Bits(be.signal_qubits)},
        resources_for(("a", be.operation)),
    )
    # 文献标准序列 U R U† R U：与修复后的库算子逐振幅一致
    invoke(b, be.operation, "a", target=b["target"], signal=b["signal"])
    reflect_zero(b, b["signal"])
    with b.adjoint():
        invoke(b, be.operation, "a", target=b["target"], signal=b["signal"])
    reflect_zero(b, b["signal"])
    invoke(b, be.operation, "a", target=b["target"], signal=b["signal"])
    operation = b.finish()
    library = oblivious_amplification(be, iterations=1)
    before, _ = zero_signal_block(be.operation, reference, 2)
    worst = 0.0
    lib_vs_script = 0.0
    measured = None
    for runner in (reference, rir_pysparq, originir_ext):
        block, _ = zero_signal_block(operation, runner, 2)
        lib_block, _ = zero_signal_block(library, runner, 2)
        measured = block
        worst = max(worst, float(np.abs(block + unitary_v).max()))
        lib_vs_script = max(lib_vs_script, float(np.abs(block - lib_block).max()))
    amplification = float(abs(measured[0, 0]) / abs(before[0, 0]))
    report.case(
        "oaa-standard-sequence-amplification",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"iterations": 1, "input_block": "V/2 (V 幺正)"},
        metrics={
            "max_error_vs_minus_V": worst,
            "library_vs_script_error": lib_vs_script,
            "amplitude_before": float(abs(before[0, 0])),
            "amplitude_after": float(abs(measured[0, 0])),
            "amplification": amplification,
        },
        criterion=(
            "零信号块等于 −V（max_error < 1e-12），幅度放大 0.5 → 1.0；"
            "库算子与脚本组装的标准序列逐振幅一致（library_vs_script_error < 1e-12）"
        ),
        passed=worst < 1e-12 and lib_vs_script < 1e-12,
    )


def _gibbs_case(report, name, matrix, beta, error, runners):
    matrix = np.asarray(matrix, dtype=complex)
    be = matrix_pauli_encoding(matrix)
    gibbs = gibbs_purification(be, beta, error=error)
    operation = gibbs.operation
    widths = _register_widths(operation)
    n_sys, n_env = widths[0], widths[1]
    expected = gibbs_reference(matrix, beta)
    worst = 0.0
    success = 0.0
    for runner in runners:
        state = runner(basis_program(operation, {}, name="gibbs_run"))
        psi = np.zeros((1 << n_sys, 1 << n_env), dtype=complex)
        total = 0.0
        for key, amplitude in _iter_entries(state, widths):
            if key[2] == 0:
                psi[key[0], key[1]] = amplitude
                total += abs(amplitude) ** 2
        rho = psi @ psi.conj().T
        rho = rho / np.trace(rho)
        worst = max(worst, trace_distance(rho, expected))
        success = total
    report.case(
        name,
        paths=[r.__name__ for r in runners],
        parameters={"beta": beta, "error": error, "qubits": sum(widths)},
        metrics={"trace_distance": worst, "success_probability": success},
        criterion=f"约化密度矩阵与经典 Gibbs 态的迹距离 < 3×error（{3 * error:.2f}）",
        passed=worst < 3 * error,
    )


def verify_gibbs(report):
    """gibbs_purification：后选择偏迹与经典 Gibbs 态的迹距离（对角/非对角/β=0）。"""
    _gibbs_case(
        report,
        "gibbs-trace-distance-nondiag",
        [[0.5, 0.2], [0.2, -0.3]],
        1.2,
        0.02,
        (reference, rir_pysparq, originir_ext),
    )
    _gibbs_case(
        report,
        "gibbs-trace-distance-diag",
        [[0.8, 0.0], [0.0, -0.4]],
        0.6,
        0.02,
        (reference, rir_pysparq),
    )
    _gibbs_case(
        report,
        "gibbs-beta0-maximally-mixed",
        [[0.5, 0.2], [0.2, -0.3]],
        0.0,
        0.02,
        (reference, rir_pysparq, originir_ext),
    )


def run():
    report = Report(
        "hamiltonian",
        "Trotter/协议层/QSVT 标准变换的幺正-乘积公式对拍、收敛阶拟合、"
        "多项式逐点对拍、求逆条件数对比与 Gibbs 迹距离，覆盖 1–4 量子位。",
    )
    verify_trotter_single_term(report)
    verify_trotter_product_formula(report)
    verify_trotter_superposition_cross(report)
    verify_trotter_convergence(report)
    verify_protocol_trotter(report)
    verify_protocol_qsp_injection(report)
    verify_phase_synthesis(report)
    verify_qsvt_sequence_convention(report)
    verify_qsvt_sequence_matrix_block(report)
    verify_qubitization(report)
    verify_qsvt_hamsim(report)
    verify_qsvt_inversion(report)
    verify_eigenstate_filter(report)
    verify_oaa_iterate_identity(report)
    verify_oaa_standard_sequence(report)
    verify_gibbs(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
