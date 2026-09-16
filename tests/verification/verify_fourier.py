"""Fourier 与 transforms 模块的论文级数值验证。

覆盖两个源文件：

- ``algorithms/fourier.py``（qft / qft_with_work / inverse_qft / fourier_add）：
  * n ≤ 4 全幺正与独立 DFT 矩阵逐元素对比（OriginIR-ext → UniQC ``to_matrix``
    加 ``effective_block`` 提取）；
  * n = 6/8/12 叠加态相位逐点对比与 Fourier 模式聚焦（rir-pysparq）；
  * fourier_add 与整数加法真值表（n ≤ 4 幺正与态向量穷举、n = 8 叠加穷举、
    n = 12 确定性采样基态对）；
  * 逆 QFT 往返恒等（幺正级 n = 3/4，寄存器级 n = 8/12）。
- ``algorithms/transforms.py``（qubitization_walk / qsvt_sequence /
  oblivious_amplification）：输入为开放块编码接口，按任务约定以 gate 级绑定的
  小实例（matrix_pauli_encoding 的 1 比特目标厄米矩阵、门级 X/2 夹具、
  signal 为 0 宽的 Pauli 字）做端到端验证：
  * 行走算子全幺正 = (2Π−I)U，谱转角 = ±arccos(λ/α)；
  * QSVT 序列全幺正 = 交替 S(φ)·U/U† 乘积；零信号块 = p(A/α)，经典参考为
    独立实现的 2×2 乘积公式与 Chebyshev 矩阵递推（不复用被测组装逻辑）；
  * OAA 迭代全幺正 = [R U† R U]^it；并量测零块为 X/2 的最小实例的迭代行为，
    与标准三查询 OAA 叙事（sin θ → sin 3θ）的差距作为信息性指标记录。

经典参考全部独立构造：numpy/cmath 的 DFT 矩阵与置换矩阵、逐比特相位门制备
（不经过被测 qft 组装）、显式 Pauli 展开求 α、2×2 QSP 乘积与 Chebyshev 递推。

运行：PYTHONPATH=src <含 pysparq+uniqc 的 python> tests/verification/verify_fourier.py
"""

from __future__ import annotations

import cmath
import math
import random

import numpy as np
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
    sampled_inputs,
    statevector_error,
)

from pyqecclang import Bits, Builder
from pyqecclang.algorithms.block_encoding import matrix_pauli_encoding, pauli_word
from pyqecclang.algorithms.fourier import fourier_add, inverse_qft, qft, qft_with_work
from pyqecclang.algorithms.operators import BlockEncoding
from pyqecclang.algorithms.oracles import annotate
from pyqecclang.algorithms.qsvt import qsp_phases
from pyqecclang.algorithms.transforms import (
    oblivious_amplification,
    qsvt_sequence,
    qubitization_walk,
)

# ---------------------------------------------------------------------------
# 独立经典参考
# ---------------------------------------------------------------------------


def _dft_matrix(n):
    """正号 DFT 矩阵：F[y, x] = exp(2πi x y / 2^n) / sqrt(2^n)。"""
    size = 1 << n
    scale = 1.0 / math.sqrt(size)
    return np.array(
        [
            [cmath.exp(2j * math.pi * x * y / size) * scale for x in range(size)]
            for y in range(size)
        ]
    )


def _add_permutation(n):
    """fourier_add 的经典真值表置换：|a, b> → |a, (a+b) mod 2^n>（a 在低位）。"""
    size = 1 << n
    perm = np.zeros((size * size, size * size), dtype=complex)
    for a in range(size):
        for b in range(size):
            perm[a | (((a + b) % size) << n), a | (b << n)] = 1.0
    return perm


def _pauli_alpha(matrix):
    """显式 Pauli 展开求块编码归一化 α = Σ_P |Tr(P M)| / d（独立于库实现）。"""
    m = np.asarray(matrix, dtype=complex)
    d = m.shape[0]
    n = (d - 1).bit_length()
    paulis = [
        np.eye(2),
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]),
        np.array([[1, 0], [0, -1]]),
    ]
    alpha = 0.0
    for index in range(4**n):
        word = np.array([[1.0 + 0j]])
        rest = index
        for _ in range(n):
            word = np.kron(paulis[rest % 4], word)
            rest //= 4
        alpha += abs(np.trace(word @ m)) / d
    return alpha


def _signal_projector(dim, data_dim):
    """Π = |0><0|_signal ⊗ I（入口布局中 signal 为高位，零块即低 data_dim 维）。"""
    proj = np.zeros((dim, dim))
    proj[:data_dim, :data_dim] = np.eye(data_dim)
    return proj


def _signal_phase(phi, dim, data_dim):
    """QSVT 相位算子 S(φ)：信号零分支 e^{+iφ}，其余分支 e^{-iφ}。"""
    diag = np.full(dim, cmath.exp(-1j * phi), dtype=complex)
    diag[:data_dim] = cmath.exp(1j * phi)
    return np.diag(diag)


def _qsp_response_independent(x, phases):
    """文档约定的 2×2 乘积参考：p(x) = [S(φ_d) W … W S(φ_0)]_00，右端最先作用。"""
    s = math.sqrt(max(0.0, 1.0 - x * x))
    w = np.array([[x, s], [s, -x]], dtype=complex)
    total = np.diag([cmath.exp(1j * phases[0]), cmath.exp(-1j * phases[0])])
    for phi in phases[1:]:
        total = np.diag([cmath.exp(1j * phi), cmath.exp(-1j * phi)]) @ w @ total
    return total[0, 0]


def _chebyshev_matrix(degree, m):
    """T_d(M) 的独立递推：T_0 = I，T_1 = M，T_{k+1} = 2 M T_k − T_{k−1}。"""
    if degree == 0:
        return np.eye(m.shape[0], dtype=complex)
    t0, t1 = np.eye(m.shape[0], dtype=complex), np.asarray(m, dtype=complex).copy()
    for _ in range(2, degree + 1):
        t0, t1 = t1, 2 * m @ t1 - t0
    return t1


def _angle_distance(a, b):
    """两个转角在圆周上的距离（[0, π]）。"""
    d = abs(a - b) % (2 * math.pi)
    return min(d, 2 * math.pi - d)


# ---------------------------------------------------------------------------
# 驱动程序与夹具
# ---------------------------------------------------------------------------


def _fourier_mode_program(n, k):
    """逐比特相位门制备负号 Fourier 模式后接 qft。

    制备态 |φ_k⟩ = Σ_x exp(-2πi k x / 2^n) |x⟩ / sqrt(2^n) 是可分解乘积态，
    只用 H 与单比特相位门（不经过被测 qft 组装）；它是正号 DFT 的第 -k 行，
    故 QFT 后全部振幅应聚焦到基态索引 k。
    """
    size = 1 << n
    b = Builder(f"fourier_mode_{n}_{k}", {"target": Bits(n)})
    for j in range(n):
        b.h(b["target"][j])
        b.gate("phase", b["target"][j], -2.0 * math.pi * k * (1 << j) / size)
    b.call(qft(n), target=b["target"])
    return b.finish().program()


def _roundtrip_unitary_program(n):
    """qft 后接 inverse_qft 的裸复合（幺正级往返检查，不含态制备）。"""
    b = Builder(f"qft_roundtrip_unitary_{n}", {"target": Bits(n)})
    b.call(qft(n), target=b["target"])
    b.call(inverse_qft(n), target=b["target"])
    return b.finish().program()


def _roundtrip_program(n, *, value=None):
    """qft 与 inverse_qft 的往返程序；value 为 None 时以全叠加出发。"""
    b = Builder(f"qft_roundtrip_{n}_{value}", {"target": Bits(n)})
    if value is None:
        b.h(b["target"])
    else:
        for bit in range(n):
            if (value >> bit) & 1:
                b.x(b["target"][bit])
    b.call(qft(n), target=b["target"])
    b.call(inverse_qft(n), target=b["target"])
    return b.finish().program()


def _fourier_add_fixed_b(n, b0):
    """a 全叠加、b 固定为已知常量后调用 fourier_add：分支 (a, (a+b0) mod 2^n)。"""
    b = Builder(f"fadd_superposed_a_{n}_{b0}", {"a": Bits(n), "b": Bits(n)})
    b.h(b["a"])
    for bit in range(n):
        if (b0 >> bit) & 1:
            b.x(b["b"][bit])
    b.call(fourier_add(n), a=b["a"], b=b["b"])
    return b.finish().program()


# transforms 的 gate 级绑定夹具：1 比特目标厄米矩阵（谱值归一后落在 (-1, 1)）
MATRIX = np.array([[0.5, 0.3], [0.3, -0.1]])


def _fixture_be():
    """matrix_pauli_encoding 属 block_encoding 模块，此处仅作绑定实例来源。"""
    return matrix_pauli_encoding(MATRIX.tolist())


def _half_x_be():
    """门级最小 OAA 夹具：U = M ⊗ X，零信号块恰为 X/2（sin θ = 1/2 的 V/2 实例）。"""
    b = Builder("be_half_x", {"target": Bits(1), "signal": Bits(1)})
    b.x(b["target"])
    b.z(b["signal"])
    b.ry(b["signal"], 2.0 * math.pi / 3.0)
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))


def _unitary_from_runs(operation, target_width, signal_width):
    """经 reference 路径逐基态列组装小实例幺正（与 OriginIR 路径交叉对拍）。"""
    data = 1 << target_width
    dim = 1 << (target_width + signal_width)
    matrix = np.zeros((dim, dim), dtype=complex)
    for s in range(1 << signal_width):
        for t in range(data):
            state = reference(basis_program(operation, {"target": t, "signal": s}))
            for (tv, sv), amp in state.items():
                matrix[tv | (sv << target_width), t | (s << target_width)] = amp
    return matrix


# ---------------------------------------------------------------------------
# fourier.py：幺正级验证（n ≤ 4，OriginIR-ext + UniQC to_matrix）
# ---------------------------------------------------------------------------


def verify_qft_unitary(report):
    for n in (1, 2, 3, 4):
        unitary = originir_unitary(qft(n).program())
        block, leakage = effective_block(unitary, n)
        error = float(np.abs(block - _dft_matrix(n)).max())
        report.case(
            f"qft-unitary-dft-n{n}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"width": n, "dimension": 1 << n},
            metrics={"max_error": error, "workspace_leakage": leakage},
            criterion="全幺正与正号 DFT 矩阵逐元素一致（max_error < 1e-12，泄漏为 0）",
            passed=error < 1e-12 and leakage < 1e-12,
        )


def verify_inverse_qft_unitary(report):
    for n in (1, 2, 3, 4):
        unitary = originir_unitary(inverse_qft(n).program())
        block, leakage = effective_block(unitary, n)
        error = float(np.abs(block - _dft_matrix(n).conj().T).max())
        report.case(
            f"inverse-qft-unitary-n{n}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"width": n, "dimension": 1 << n},
            metrics={"max_error": error, "workspace_leakage": leakage},
            criterion="全幺正与 DFT 的伴随矩阵逐元素一致（max_error < 1e-12）",
            passed=error < 1e-12 and leakage < 1e-12,
        )
    for n in (3, 4):
        unitary = originir_unitary(_roundtrip_unitary_program(n))
        error = float(np.abs(unitary - np.eye(1 << n)).max())
        report.case(
            f"inverse-qft-roundtrip-unitary-n{n}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"width": n},
            metrics={"max_error": error},
            criterion="qft 后接 inverse_qft 的复合幺正为恒等（max_error < 1e-12）",
            passed=error < 1e-12,
        )


# ---------------------------------------------------------------------------
# fourier.py：寄存器级验证（rir-pysparq，n = 6/8/12）
# ---------------------------------------------------------------------------


def verify_qft_basis_rows(report):
    """基态 |x> 经 QFT 的全部 2^n 个振幅与 DFT 行逐点对比。"""
    for n in (6, 8, 12):
        size = 1 << n
        uniform = 1.0 / math.sqrt(size)
        inputs = sorted({0, 1, size // 2, size - 1, 0x5A5A % size})
        worst = 0.0
        for x in inputs:
            state = rir_pysparq(basis_program(qft(n), {"target": x}), max_states=1 << (n + 1))
            for y in range(size):
                expected = cmath.exp(2j * math.pi * x * y / size) * uniform
                worst = max(worst, abs(state.get((y,), 0j) - expected))
        report.case(
            f"qft-basis-row-pointwise-n{n}",
            paths=["rir-pysparq"],
            parameters={"width": n, "inputs": inputs, "amplitudes_per_input": size},
            metrics={"max_error": worst},
            criterion="每个输入的全部 2^n 个振幅与 DFT 行逐点一致（max_error < 1e-9）",
            passed=worst < 1e-9,
        )


def verify_qft_fourier_mode_focus(report):
    """独立门级制备的 Fourier 模式经 QFT 后聚焦到单个基态索引。"""
    for n in (6, 8, 12):
        size = 1 << n
        modes = sorted({1, size // 3, size - 1})
        min_probability = 1.0
        worst_leak = 0.0
        for k in modes:
            state = rir_pysparq(_fourier_mode_program(n, k), max_states=1 << (n + 1))
            min_probability = min(min_probability, abs(state.get((k,), 0j)) ** 2)
            leak = max(
                (abs(v) for key, v in state.items() if key != (k,)),
                default=0.0,
            )
            worst_leak = max(worst_leak, leak)
        report.case(
            f"qft-fourier-mode-focus-n{n}",
            paths=["rir-pysparq"],
            parameters={"width": n, "modes": modes},
            metrics={
                "min_success_probability": min_probability,
                "max_leaked_amplitude": worst_leak,
            },
            criterion="聚焦成功概率 ≥ 1 − 1e-12 且模式外泄漏振幅 < 1e-6",
            passed=min_probability >= 1 - 1e-12 and worst_leak < 1e-6,
        )


def verify_qft_cross_path(report):
    """n = 4 富相位态（QFT 作用于基态）四路径对拍，附 qft_with_work 适配等价。"""
    program = basis_program(qft(4), {"target": 11})
    ref = reference(program)
    deviation = max(
        amplitude_error(ref, rir_pysparq(program)),
        amplitude_error(ref, adapter_pysparq(program)),
        statevector_error(originir_ext(program), amplitudes_to_statevector(ref, [4])),
    )
    legacy = reference(basis_program(qft_with_work(4), {"target": 11}))
    legacy_deviation = max(
        abs(legacy.get((y, 0), 0j) - ref.get((y,), 0j)) for y in range(16)
    )
    deviation = max(deviation, legacy_deviation)
    report.case(
        "qft-cross-path-n4",
        paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        parameters={"width": 4, "input": 11},
        metrics={"max_pairwise_deviation": deviation},
        criterion="四路径及 qft_with_work 适配两两偏差 < 1e-9",
        passed=deviation < 1e-9,
    )


# ---------------------------------------------------------------------------
# fourier.py：Fourier 加法真值表
# ---------------------------------------------------------------------------


def verify_fourier_add_unitary(report):
    for n in (1, 2, 3, 4):
        unitary = originir_unitary(fourier_add(n).program())
        block, leakage = effective_block(unitary, 2 * n)
        error = float(np.abs(block - _add_permutation(n)).max())
        report.case(
            f"fourier-add-unitary-n{n}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"width": n, "dimension": 1 << (2 * n)},
            metrics={"max_error": error, "workspace_leakage": leakage},
            criterion="全幺正与整数模加法置换矩阵一致（max_error < 1e-12）",
            passed=error < 1e-12 and leakage < 1e-12,
        )


def verify_fourier_add_truth_table(report):
    """n ≤ 4：固定 b0、叠加 a 逐次覆盖全部 4^n 个 (a, b) 输入（双路径对拍）。"""
    for n in (1, 2, 3, 4):
        size = 1 << n
        uniform = 1.0 / math.sqrt(size)
        worst = 0.0
        for b0 in range(size):
            program = _fourier_add_fixed_b(n, b0)
            expected = {(a, (a + b0) % size): uniform for a in range(size)}
            worst = max(worst, amplitude_error(reference(program), expected))
            vector = originir_ext(program)
            worst = max(
                worst,
                statevector_error(vector, amplitudes_to_statevector(expected, [n, n])),
            )
        report.case(
            f"fourier-add-truthtable-n{n}",
            paths=["reference", "originir-ext"],
            parameters={"width": n, "input_pairs": size * size},
            metrics={"max_error": worst},
            criterion="全部 (a, b) 输入对映射到 |a, (a+b) mod 2^n> 且振幅均匀（max_error < 1e-9）",
            passed=worst < 1e-9,
        )


def verify_fourier_add_superposed_a_n8(report):
    """n = 8：rir-pysparq 一次叠加穷举 a 的全部 256 个分支（中间态峰值 2^16）。"""
    n, size = 8, 256
    uniform = 1.0 / math.sqrt(size)
    worst = 0.0
    for b0 in (0, 1, 0x55, 0x80, 0xFF):
        state = rir_pysparq(_fourier_add_fixed_b(n, b0), max_states=1 << 17)
        expected = {(a, (a + b0) % size): uniform for a in range(size)}
        worst = max(worst, amplitude_error(state, expected))
    report.case(
        "fourier-add-superposed-a-n8",
        paths=["rir-pysparq"],
        parameters={"width": n, "constants": [0, 1, 0x55, 0x80, 0xFF], "branches": size},
        metrics={"max_error": worst},
        criterion="每个 b0 的 256 个分支逐振幅一致（max_error < 1e-9）",
        passed=worst < 1e-9,
    )


def verify_fourier_add_sampled_n12(report):
    """n = 12：采样基态对的确定性输出检查（中间态峰值 2^12，叠加 a 会超预算）。"""
    n = 12
    values, _ = sampled_inputs(n)
    pairs = list(zip(values[:12], values[1:13], strict=True))
    operation = fourier_add(n)
    failures = 0
    for a, b0 in pairs:
        state = rir_pysparq(
            basis_program(operation, {"a": a, "b": b0}), max_states=1 << 13
        )
        if set(state) != {(a, (a + b0) % (1 << n))} or abs(
            next(iter(state.values())) - 1
        ) > 1e-12:
            failures += 1
    report.case(
        "fourier-add-sampled-basis-n12",
        paths=["rir-pysparq"],
        parameters={"width": n, "sampled_pairs": len(pairs)},
        metrics={"failures": failures},
        criterion="采样输入的输出恰为 |a, (a+b) mod 2^12>（failures == 0）",
        passed=failures == 0,
    )


# ---------------------------------------------------------------------------
# fourier.py：逆 QFT 往返恒等（寄存器级）
# ---------------------------------------------------------------------------


def verify_inverse_qft_roundtrip_wide(report):
    for n in (8, 12):
        size = 1 << n
        uniform = 1.0 / math.sqrt(size)
        state = rir_pysparq(_roundtrip_program(n), max_states=1 << (n + 1))
        expected = {(x,): uniform for x in range(size)}
        uniform_error = amplitude_error(state, expected)
        failures = 0
        for x in (0, 1, size // 2, size - 1, 0x5A5A % size):
            single = rir_pysparq(_roundtrip_program(n, value=x), max_states=1 << (n + 1))
            if set(single) != {(x,)} or abs(next(iter(single.values())) - 1) > 1e-12:
                failures += 1
        report.case(
            f"inverse-qft-roundtrip-wide-n{n}",
            paths=["rir-pysparq"],
            parameters={"width": n, "branches": size, "basis_samples": 5},
            metrics={"uniform_max_error": uniform_error, "basis_failures": failures},
            criterion="全叠加往返恢复均匀态（< 1e-9）且采样基态往返恒等（failures == 0）",
            passed=uniform_error < 1e-9 and failures == 0,
        )


# ---------------------------------------------------------------------------
# transforms.py：门级绑定小实例的端到端验证
# ---------------------------------------------------------------------------


def verify_qubitization_walk(report):
    be = _fixture_be()
    tw, sw = be.width, be.signal_qubits
    data, dim = 1 << tw, 1 << (tw + sw)
    u_runs = _unitary_from_runs(be.operation, tw, sw)
    u_origin = originir_unitary(be.operation.program())
    crosscheck = float(np.abs(u_runs - u_origin).max())
    proj = _signal_projector(dim, data)
    walk_expected = (2 * proj - np.eye(dim)) @ u_runs
    walk_actual = originir_unitary(qubitization_walk(be).program())
    unitary_error = float(np.abs(walk_actual - walk_expected).max())
    alpha = _pauli_alpha(MATRIX)
    zero_block_error = float(np.abs(walk_actual[:data, :data] - MATRIX / alpha).max())
    # 谱性质：每个谱值 x = λ/α 对应一对本征角 ±arccos(x)，其余空间转角 ∈ {0, π}
    eigenangles = [abs(float(a)) for a in np.angle(np.linalg.eigvals(walk_actual))]
    allowed = {0.0, math.pi}
    for x in np.linalg.eigvalsh(MATRIX / alpha):
        allowed.add(math.acos(float(np.clip(x, -1.0, 1.0))))
    spectrum_dev = max(
        min(_angle_distance(angle, expect) for expect in allowed) for angle in eigenangles
    )
    report.case(
        "qubitization-walk-unitary-spectrum",
        paths=["originir-ext+uniqc-to_matrix", "reference"],
        parameters={
            "matrix": MATRIX.tolist(),
            "alpha": alpha,
            "signal_qubits": sw,
            "dimension": dim,
        },
        metrics={
            "unitary_max_error": unitary_error,
            "be_path_crosscheck": crosscheck,
            "zero_block_max_error": zero_block_error,
            "spectrum_max_deviation": spectrum_dev,
        },
        criterion=(
            "行走幺正 = (2Π−I)U（< 1e-12），零信号块 = A/α，"
            "本征角落入 {±arccos(λ/α)} ∪ {0, π}（角度偏差 < 1e-9）"
        ),
        passed=unitary_error < 1e-12
        and zero_block_error < 1e-12
        and spectrum_dev < 1e-9,
    )


def verify_qsvt_sequence(report):
    be = _fixture_be()
    tw, sw = be.width, be.signal_qubits
    data, dim = 1 << tw, 1 << (tw + sw)
    u = originir_unitary(be.operation.program())
    rng = random.Random(20260916)
    phases = tuple(rng.uniform(-math.pi, math.pi) for _ in range(7))
    actual = originir_unitary(qsvt_sequence(be, phases).program())
    # 全幺正参考：时间正序交替 S(φ)·U/U†（U 取 BE 的实测幺正）
    expected = _signal_phase(phases[0], dim, data)
    for k in range(len(phases) - 1):
        expected = (u if k % 2 == 0 else u.conj().T) @ expected
        expected = _signal_phase(phases[k + 1], dim, data) @ expected
    unitary_error = float(np.abs(actual - expected).max())
    # 零信号块谱语义：p(A/α)，p 来自独立 2×2 乘积公式
    alpha = _pauli_alpha(MATRIX)
    values, vectors = np.linalg.eigh(MATRIX / alpha)
    response = [complex(_qsp_response_independent(float(x), phases)) for x in values]
    expected_block = (vectors * response) @ vectors.conj().T
    block_error = float(np.abs(actual[:data, :data] - expected_block).max())
    report.case(
        "qsvt-sequence-random-phases",
        paths=["originir-ext+uniqc-to_matrix"],
        parameters={"degree": len(phases) - 1, "seed": 20260916},
        metrics={
            "unitary_max_error": unitary_error,
            "zero_block_max_error": block_error,
        },
        criterion="序列幺正 = 交替 S(φ)·U/U† 乘积且零块 = p(A/α)（均 < 1e-12）",
        passed=unitary_error < 1e-12 and block_error < 1e-12,
    )
    # Chebyshev 多项式端到端：qsp_phases 仅作相位输入生成，参考为独立 T_d 递推
    for degree in (4, 5):
        coeffs = tuple(float(c) for c in np.polynomial.chebyshev.cheb2poly([0] * degree + [1]))
        t_phases = qsp_phases(coeffs)
        actual = originir_unitary(qsvt_sequence(be, t_phases).program())
        expected_block = _chebyshev_matrix(degree, MATRIX / alpha)
        block_error = float(np.abs(actual[:data, :data] - expected_block).max())
        report.case(
            f"qsvt-sequence-chebyshev-t{degree}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"degree": degree, "alpha": alpha},
            metrics={"zero_block_max_error": block_error},
            criterion=f"零信号块 = T_{degree}(A/α)（独立 Chebyshev 递推，< 1e-12）",
            passed=block_error < 1e-12,
        )


def verify_qsvt_degenerate(report):
    """signal_qubits == 0 的退化分支：相位退化为无条件全局相位。"""
    be = pauli_word("X")
    x_matrix = np.array([[0, 1], [1, 0]], dtype=complex)
    for phases in ((0.3, -0.7, 1.1), (0.3, -0.7, 1.1, 0.5)):
        degree = len(phases) - 1
        base = np.eye(2, dtype=complex) if degree % 2 == 0 else x_matrix
        expected = cmath.exp(1j * sum(phases)) * base
        actual = originir_unitary(qsvt_sequence(be, phases).program())
        error = float(np.abs(actual - expected).max())
        report.case(
            f"qsvt-sequence-degenerate-signal0-d{degree}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"degree": degree, "be": "pauli_word(X)"},
            metrics={"max_error": error},
            criterion="退化序列 = e^{iΣφ}·X^d（max_error < 1e-12）",
            passed=error < 1e-12,
        )


def verify_oaa_unitary(report):
    be = _fixture_be()
    tw, sw = be.width, be.signal_qubits
    data, dim = 1 << tw, 1 << (tw + sw)
    u = originir_unitary(be.operation.program())
    proj = _signal_projector(dim, data)
    reflect = np.eye(dim) - 2 * proj
    iterate = reflect @ u.conj().T @ reflect @ u
    expected = np.eye(dim, dtype=complex)
    for iterations in (1, 2, 3):
        expected = expected @ iterate
        actual = originir_unitary(oblivious_amplification(be, iterations).program())
        error = float(np.abs(actual - expected).max())
        report.case(
            f"oaa-unitary-it{iterations}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"iterations": iterations, "dimension": dim},
            metrics={"max_error": error},
            criterion="迭代幺正 = [R U† R U]^it（R = I − 2Π，max_error < 1e-12）",
            passed=error < 1e-12,
        )


def verify_oaa_half_block(report):
    """零块为 X/2 的最小实例：量测实现算子的零块行为并记录与标准 OAA 的差距。

    注：实现的迭代 R U† R U 的零块为 2B†B − I（B = X/2 时为 −I/2），
    与标准三查询 OAA（sin θ → sin 3θ，θ = π/6 时应恢复为 X）不同；
    本案例按文档给出的算子公式判定，并把与标准 OAA 的偏差列为信息性指标。
    """
    x_matrix = np.array([[0, 1], [1, 0]], dtype=complex)
    b_block = x_matrix / 2
    actual = originir_unitary(oblivious_amplification(_half_x_be(), 1).program())
    zero_block = actual[:2, :2]
    implemented_error = float(
        np.abs(zero_block - (2 * b_block.conj().T @ b_block - np.eye(2))).max()
    )
    standard_deviation = float(np.abs(zero_block - x_matrix).max())
    report.case(
        "oaa-half-block-behavior",
        paths=["originir-ext+uniqc-to_matrix"],
        parameters={"fixture": "zero block = X/2 (sin theta = 1/2)", "iterations": 1},
        metrics={
            "implemented_operator_error": implemented_error,
            "standard_oaa_deviation": standard_deviation,
        },
        criterion=(
            "零块与文档算子 R U† R U 的代数结果 2B†B−I 一致（< 1e-12）；"
            "standard_oaa_deviation 为信息性指标，记录与标准 OAA 叙事的差距"
        ),
        passed=implemented_error < 1e-12,
    )


def run():
    report = Report(
        "fourier",
        "QFT/逆 QFT/Fourier 加法的幺正-叠加-真值表三级验证（n ≤ 12），"
        "以及 transforms 块编码变换（行走/QSVT/OAA）的门级绑定端到端验证。",
    )
    verify_qft_unitary(report)
    verify_inverse_qft_unitary(report)
    verify_qft_basis_rows(report)
    verify_qft_fourier_mode_focus(report)
    verify_qft_cross_path(report)
    verify_fourier_add_unitary(report)
    verify_fourier_add_truth_table(report)
    verify_fourier_add_superposed_a_n8(report)
    verify_fourier_add_sampled_n12(report)
    verify_inverse_qft_roundtrip_wide(report)
    verify_qubitization_walk(report)
    verify_qsvt_sequence(report)
    verify_qsvt_degenerate(report)
    verify_oaa_unitary(report)
    verify_oaa_half_block(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
