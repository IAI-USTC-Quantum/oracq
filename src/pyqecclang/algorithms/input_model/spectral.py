"""谱输入/输出原语：傅里叶基对角块编码、稀疏谱块编码与谱态制备。

对应 arXiv:2509.08807 的 Lemma C.8（SS-BE）与层级谱编码的输入侧。
傅里叶基函数 F_k(J)=exp(2*pi*1j*k*J/2**w) 的对角矩阵按 J 的二进制位分解
为单比特相位门的张量积（论文 Eq. C58），受控乘积形式给出零 Toffoli 的
谱块编码（Eq. C60）；谱态制备在同一个格点寄存器上稀疏制备频域幅度，
再经 QFT 完成 S 维谱空间到 N 维格点空间的 Hilbert 空间放大。

约定：f(J) = sum_k c_k * exp(2*pi*1j*k*J/2**width)，频率 k 为整数
（负频率按 mod 2**width 理解）；谱块编码的 alpha = sum_k abs(c_k)。
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Mapping, Sequence

from pyqecclang.algorithms.common.fourier import qft
from pyqecclang.algorithms.input_model.block_encoding import lcu
from pyqecclang.algorithms.input_model.operators import BlockEncoding, _name
from pyqecclang.algorithms.input_model.oracles import (
    StatePreparation,
    _state_angles,
    annotate,
    gate_state_prep,
)
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, Ref, ValidationError


def normalized_spectrum(spectrum: Mapping[int, complex]) -> tuple[int, tuple[complex, ...]]:
    """校验谱系数并规范成连续频带。

    输入 ``{k: c_k}`` 表示 f(J) = sum_k c_k exp(2*pi*1j*k*J/2**width)。
    受控乘积分解要求频带连续；返回 ``(k_min, coefficients)``，coefficients
    覆盖 ``[k_min, k_max]`` 每个整数频率（带内全零项允许，计零系数）。

    Args:
        spectrum: 频率到复系数的非空映射，须构成连续频带且系数不全为零。

    Returns:
        tuple[int, tuple[complex, ...]]: ``(k_min, coefficients)``，即频带
        下限与按升频排列、含带内零系数的元组。
    """
    if not isinstance(spectrum, dict) or not spectrum:
        raise ValidationError("谱系数必须是非空字典 {频率: 复系数}")
    ks = sorted(spectrum)
    if any(type(k) is not int for k in ks):
        raise ValidationError("谱频率必须是整数")
    if ks[-1] - ks[0] + 1 != len(spectrum):
        raise ValidationError("谱频率必须构成连续频带；间断频带请拆分后用 LCU 组合")
    k_min = ks[0]
    coefficients = tuple(spectrum[k] for k in range(ks[0], ks[-1] + 1))
    if all(c == 0 for c in coefficients):
        raise ValidationError("谱系数不能全为零")
    for c in coefficients:
        if not (math.isfinite(c.real) and math.isfinite(c.imag)):
            raise ValidationError("谱系数必须有限")
    return k_min, coefficients


def _spectrum_register_width(coefficients: Sequence[complex]) -> int:
    """频带项数 ``S`` 所需的谱寄存器位宽 ``ceil(log2 S)``，至少一位。"""
    return max(1, (len(coefficients) - 1).bit_length())


def frequency_amplitudes(
    coefficients: Sequence[complex], k_min: int, width: int, *, normalize: bool = True
) -> list[complex]:
    """谱系数到 2**width 维幅度向量（索引为 k mod 2**width）。

    Args:
        coefficients: 频带内按升频排列的复系数序列。
        k_min: 频带下限频率。
        width: 格点寄存器位宽，决定幅度向量的维数 2**width。
        normalize: 为 True 时除以系数的欧几里得范数。

    Returns:
        list[complex]: 长度 2**width 的频域幅度向量。
    """
    scale = math.sqrt(sum(abs(c) ** 2 for c in coefficients)) if normalize else 1.0
    if scale == 0:
        raise ValidationError("谱系数范数为零")
    amplitudes = [0j] * (1 << width)
    for j, c in enumerate(coefficients):
        if c:
            amplitudes[(k_min + j) % (1 << width)] = c / scale
    return amplitudes


def pruned_state_prep(
    amplitudes: Sequence[complex], *, name: str | None = None
) -> StatePreparation:
    """零角剪枝的多路旋转态制备；与 gate_state_prep 酉等价但省去零旋转。

    稀疏谱输入（2**width 维向量只有 S 个非零幅度）下，全零子树的所有
    旋转角为零，剪枝后资源随 S 而不是 2**width 增长。

    Args:
        amplitudes: 稀疏复幅度向量，长度须为二的幂。
        name: 生成操作的名称；缺省由幅度内容派生。

    Returns:
        StatePreparation: 剪枝多路旋转实现的制备视图。
    """
    values, n, nodes = _state_angles(amplitudes)
    if n == 0:
        raise ValidationError("态制备需要至少一个目标位")
    b = Builder(
        name or _name("pruned_state", values),
        {"target": Bits(n), "work": Bits(0)},
    )
    for depth, prefix, bit, angle in nodes:
        if angle == 0:
            continue
        if not depth:
            b.ry(b["target"][bit], angle)
        else:
            with b.control(b["target"][bit + 1 :], prefix):
                b.ry(b["target"][bit], angle)
    for index, value in enumerate(values):
        if value and cmath.phase(value):
            with b.control(b["target"], index):
                b.global_phase(cmath.phase(value))
    return StatePreparation(
        annotate(
            b.finish(),
            "state_prep_isometry",
            zero_input=True,
            clean_work=True,
            implementation="multiplexed_rotations_pruned",
        )
    )


def fourier_phase(width: int, k: int) -> Operation:
    """diag_J exp(2*pi*1j*k*J/2**width) 的精确电路（论文 Eq. C58）。

    e^{2*pi*1j*k*J/2**width} 按 J 的二进制位分解为张量积，每位一个单比特
    相位门；角度是 pi 的二进制幂倍数，落在精确角度网格上。只含 target
    寄存器，作为子模块被调用。

    Args:
        width: 目标寄存器位宽，范围 1..64。
        k: 整数频率，按 ``mod 2**width`` 理解。

    Returns:
        Operation: 单比特相位门张量积组成的对角相位操作。
    """
    if type(width) is not int or not 1 <= width <= 64:
        raise ValidationError("fourier_phase.width 必须处于 1..64")
    if type(k) is not int:
        raise ValidationError("fourier_phase.k 必须是整数")
    k %= 1 << width
    b = Builder(f"fourier_phase_{width}_{k}", {"target": Bits(width)})
    for q in range(width):
        angle = 2 * math.pi * k * (1 << q) / (1 << width)
        if angle % (2 * math.pi):
            b.gate("phase", b["target"][q], angle)
    return b.finish()


def fourier_phase_encoding(width: int, k: int) -> BlockEncoding:
    """fourier_phase 的 (1, 0, 0) 块编码：零宽 signal 接口。

    Args:
        width: 目标寄存器位宽，范围 1..64。
        k: 整数频率，按 ``mod 2**width`` 理解。

    Returns:
        BlockEncoding: 归一化常数为 1、signal 零宽的块编码。
    """
    if type(width) is not int or not 1 <= width <= 64:
        raise ValidationError("fourier_phase_encoding.width 必须处于 1..64")
    if type(k) is not int:
        raise ValidationError("fourier_phase_encoding.k 必须是整数")
    k %= 1 << width
    b = Builder(
        _name("fourier_phase_be", width, k),
        {"target": Bits(width), "signal": Bits(0)},
        attributes={"matrix_interpretation": f"diag(exp(2*pi*1j*{k}*J/2^{width}))"},
    )
    for q in range(width):
        angle = 2 * math.pi * k * (1 << q) / (1 << width)
        if angle % (2 * math.pi):
            b.gate("phase", b["target"][q], angle)
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))


def _phase_grid(
    builder: Builder,
    control_builder: Ref,
    k_min: int,
    coefficients: Sequence[complex],
    width: int,
    spectrum_width: int,
) -> None:
    """谱寄存器控制的相位网格（论文 Eq. C60 的受控乘积部分）。

    对谱寄存器的每一位 bit 与目标寄存器的每一位 qubit 施加受控相位
    2*pi*2**bit*2**q/2**width；频带下限 k_min 的贡献是无条件相位。
    角度为 2*pi 整数倍的项被跳过（恒等）。
    """
    target, spectrum = builder["target"], control_builder
    for q in range(width):
        if k_min:
            angle = 2 * math.pi * k_min * (1 << q) / (1 << width)
            if angle % (2 * math.pi):
                builder.gate("phase", target[q], angle)
    for bit in range(spectrum_width):
        for q in range(width):
            exponent = (1 << bit) * (1 << q)
            if exponent % (1 << width):
                angle = 2 * math.pi * exponent / (1 << width)
                with builder.control(spectrum[bit], 1):
                    builder.gate("phase", target[q], angle)


def spectral_diagonal(
    spectrum: Mapping[int, complex], width: int, *, variant: str = "sequential"
) -> BlockEncoding:
    """稀疏谱块编码：BE(diag f)，f 由连续频带傅里叶和给出。

    variant="sequential"（论文 Lemma C.8 / Eq. (7)）：P_L 在 s=log S 位谱
    寄存器上制备幅度 sqrt(abs(c_j)/alpha)，谱寄存器各位控制目标寄存器的相位
    网格，每模相位由受控 gphase 补回，P_R=P_L^dagger 复净。无 Toffoli。
    variant="naive"：每个基函数作为独立酉算子交给通用 LCU，作为结构化
    构造（顺序形式）收益的对照基线。

    Args:
        spectrum: 频率到复系数的连续频带映射。
        width: 格点寄存器位宽。
        variant: ``sequential`` 或 ``naive``，分别走受控乘积分解或通用 LCU。

    Returns:
        BlockEncoding: ``BE(diag f)``，alpha = sum_k abs(c_k)。
    """
    if variant not in {"sequential", "naive"}:
        raise ValidationError("variant 必须是 sequential 或 naive")
    k_min, coefficients = normalized_spectrum(spectrum)
    if variant == "naive":
        return lcu(
            [
                (c, fourier_phase_encoding(width, (k_min + j) % (1 << width)))
                for j, c in enumerate(coefficients)
                if c
            ]
        )
    s = _spectrum_register_width(coefficients)
    alpha = sum(abs(c) for c in coefficients)
    amplitudes: list[complex] = []
    for c in coefficients:
        amplitudes.append(math.sqrt(abs(c) / alpha) * (c / abs(c) if c else 0))
    amplitudes += [0.0] * ((1 << s) - len(amplitudes))
    prep = gate_state_prep(amplitudes)
    b = Builder(
        _name("spectral_diagonal", coefficients, k_min, width),
        {"target": Bits(width), "signal": Bits(s)},
        attributes={
            "spectral_sparsity": sum(1 for c in coefficients if c),
            "band_low": k_min,
            "implementation": "sequential_fourier_spectral",
        },
    )
    b.call(prep.operation, target=b["signal"], work=b["signal"][:0])
    _phase_grid(b, b["signal"], k_min, coefficients, width, s)
    for j, c in enumerate(coefficients):
        if c and cmath.phase(c):
            with b.control(b["signal"], j):
                b.global_phase(cmath.phase(c))
    with b.adjoint():
        b.call(prep.operation, target=b["signal"], work=b["signal"][:0])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=alpha))


def spectral_state_prep(spectrum: Mapping[int, complex], width: int) -> StatePreparation:
    """谱输入（层级谱编码的输入侧）：制备归一化态 ``|f⟩``。

    在同一个 width 位格点寄存器上稀疏制备频域幅度 sum_k c_k / norm(c) 于基态 ``|k⟩``，
    再做正号 QFT 放大到格点空间：QFT 把基态映射为逐位相位均匀态，线性组合
    得到 sum_J f(J)*ket(J)/sqrt(N*sum(abs(c)^2))，即归一化态（Parseval）。

    Args:
        spectrum: 频率到复系数的连续频带映射，频带项数不超过 2**width。
        width: 格点寄存器位宽。

    Returns:
        StatePreparation: 经频域稀疏制备加 QFT 放大实现的制备视图。
    """
    k_min, coefficients = normalized_spectrum(spectrum)
    if len(coefficients) > (1 << width):
        raise ValidationError("谱频带宽度超过格点寄存器可表示的范围")
    amplitudes = frequency_amplitudes(coefficients, k_min, width)
    prep = pruned_state_prep(amplitudes)
    b = Builder(
        _name("spectral_state_prep", coefficients, k_min, width),
        {"target": Bits(width), "work": Bits(0)},
        attributes={
            "implementation": "hierarchy_spectral_encoding",
            "spectral_sparsity": sum(1 for c in coefficients if c),
            "zero_input": True,
            "clean_work": True,
        },
    )
    b.call(prep.operation, target=b["target"], work=b["target"][:0])
    b.call(qft(width), target=b["target"])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry"))
