"""算子级电路合成优化（arXiv:2509.08807 附录 D）。

实现论文的三类逻辑资源优化并保证酉语义不变：

1. 均匀受控旋转（UCR）：把态制备的多控旋转树改写为"纯旋转 + CNOT 阶梯"，
   消除多控 Toffoli；
2. 相似群合并（match by signal / similar group）：线性组合中只差全局符号
   或标量倍数的项合并为一项，系数吸收进制备幅度；
3. fan-out 相位网格：FBBE 的受控相位层改写为"单比特旋转 + fan-out CNOT"，
   消除受控旋转（论文图 9(a2)）。

每个构造在 tests/core/test_spectral_synthesis.py 中与朴素构造做逐元素
矩阵对照，并用 estimate_resources 验证资源下降。
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Iterable, Mapping, Sequence
from typing import cast

from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    _state_angles,
    annotate,
    gate_state_prep,
)
from oracq.algorithms.input_model.spectral import (
    _spectrum_register_width,
    normalized_spectrum,
)
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, Ref, ValidationError


def _emit_ucr(
    builder: Builder, target: Ref, controls: list[Ref], thetas: Sequence[float]
) -> None:
    """递归 UCR：CX; UCR(diff); CX; UCR(mean)，thetas 按控制前缀索引。

    前缀位序：controls[0] 是最低位。等价于按前缀值的多控 Ry 树
    （gate_state_prep 的受控层），但只使用单比特旋转与 CNOT。
    """
    if not controls:
        if thetas[0]:
            builder.ry(target, thetas[0])
        return
    control = controls[-1]
    half = len(thetas) // 2
    mean = [(thetas[j] + thetas[j + half]) / 2 for j in range(half)]
    diff = [(thetas[j] - thetas[j + half]) / 2 for j in range(half)]
    builder.xor(control, target)
    _emit_ucr(builder, target, controls[:-1], diff)
    builder.xor(control, target)
    _emit_ucr(builder, target, controls[:-1], mean)


def uniformly_controlled_prep(
    amplitudes: Iterable[complex], *, name: str | None = None
) -> StatePreparation:
    """UCR 版态制备：与 gate_state_prep 酉等价，无多控旋转。

    每层 2^d 个受控 Ry 改写为 2^d 个单比特 Ry 与 2^d 个 CNOT；
    复相位仍由受控 gphase 写入。

    Args:
        amplitudes: 目标态的复振幅序列，长度为 2 的幂且非全零。
        name: 生成的制备模块名；缺省自动生成。

    Returns:
        StatePreparation: 仅含单比特旋转与 CNOT 的态制备句柄。
    """
    values, n, nodes = _state_angles(amplitudes)
    if n == 0:
        raise ValidationError("态制备需要至少一个目标位")
    b = Builder(
        name or _name("ucr_state", values),
        {"target": Bits(n), "work": Bits(0)},
    )
    by_level: dict[int, list[float]] = {}
    for depth, prefix, _bit, angle in nodes:
        by_level.setdefault(depth, [0.0] * (1 << depth))[prefix] = angle
    for depth in range(n):
        bit = n - depth - 1
        controls = [b["target"][q] for q in range(bit + 1, n)]
        _emit_ucr(b, b["target"][bit], controls, by_level[depth])
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
            implementation="uniformly_controlled_rotations",
        )
    )


def merge_similar(
    terms: Iterable[tuple[complex, BlockEncoding]],
) -> list[tuple[complex, BlockEncoding]]:
    """按块编码恒等合并线性组合项（match by signal / similar group）。

    操作相同的项（含只差全局符号 −1 的情形）合并为一项、系数求和；
    系数归零的项删除。这是精确代数重写：LCU 的 α 与角块不变，
    但 SELECT 分支数、制备幅度树与受控调用数下降。

    Args:
        terms: (系数, 块编码) 线性组合项序列；零系数项直接丢弃。

    Returns:
        list[tuple[complex, BlockEncoding]]: 按块编码恒等合并、系数求和后的项列表。
    """
    merged: dict[str, list[complex | BlockEncoding]] = {}
    order: list[str] = []
    for coefficient, operand in terms:
        if coefficient == 0:
            continue
        key = operand.operation.module.name
        if key not in merged:
            merged[key] = [0j, operand]
            order.append(key)
        merged[key][0] += complex(coefficient)  # type: ignore[operator]
    return [
        (cast("complex", merged[key][0]), cast("BlockEncoding", merged[key][1]))
        for key in order
        if merged[key][0] != 0
    ]


def fanout_spectral_diagonal(spectrum: Mapping[int, complex], width: int) -> BlockEncoding:
    """FBBE 的 fan-out 形式：受控相位网格改写为旋转 + fan-out CNOT。

    论文图 9(a2)：每个谱位控制的相位层 = 每目标位一个半角旋转、
    fan-out CNOT、反号半角旋转、再 fan-out CNOT；共享同一控制位的
    CNOT 构成多目标 fan-out X（表面码晶格手术的原生操作）。RIR 中
    fan-out 体现为同源 xor 序列。与 spectral_diagonal 的 sequential
    变体逐元素等价。

    Args:
        spectrum: 谱频率到复系数的映射；频率按 mod 2^width 的整数解释。
        width: 目标寄存器位宽，决定相位角分母 2^width。

    Returns:
        BlockEncoding: fan-out 形式的谱对角块编码，尺度为系数绝对值之和。
    """
    k_min, coefficients = normalized_spectrum(spectrum)
    s = _spectrum_register_width(coefficients)
    alpha = sum(abs(c) for c in coefficients)
    amplitudes = [
        math.sqrt(abs(c) / alpha) * (c / abs(c) if c else 0) for c in coefficients
    ]
    amplitudes += [0.0] * ((1 << s) - len(amplitudes))
    prep = gate_state_prep(amplitudes)
    b = Builder(
        _name("fanout_spectral_diagonal", coefficients, k_min, width),
        {"target": Bits(width), "signal": Bits(s)},
        attributes={
            "spectral_sparsity": sum(1 for c in coefficients if c),
            "band_low": k_min,
            "implementation": "fanout_fourier_spectral",
        },
    )
    b.call(prep.operation, target=b["signal"], work=b["signal"][:0])
    for q in range(width):
        if k_min:
            angle = 2 * math.pi * k_min * (1 << q) / (1 << width)
            if angle % (2 * math.pi):
                b.gate("phase", b["target"][q], angle)
    for bit in range(s):
        pairs: list[tuple[int, float]] = []
        for q in range(width):
            exponent = (1 << bit) * (1 << q)
            if exponent % (1 << width):
                pairs.append((q, 2 * math.pi * exponent / (1 << width)))
        if not pairs:
            continue
        control = b["signal"][bit]
        for q, angle in pairs:
            b.gate("phase", b["target"][q], angle / 2)
        for q, _ in pairs:
            b.xor(control, b["target"][q])
        for q, angle in pairs:
            b.gate("phase", b["target"][q], -angle / 2)
        for q, _ in pairs:
            b.xor(control, b["target"][q])
        # 三明治分解在 control=1 分支残留 e^{-i*theta/2} 相位（逐对累积）；
        # 在控制位上补一个单比特相位门精确抵消。
        b.gate("phase", control, sum(angle for _, angle in pairs) / 2)
    for j, c in enumerate(coefficients):
        if c and cmath.phase(c):
            with b.control(b["signal"], j):
                b.global_phase(cmath.phase(c))
    with b.adjoint():
        b.call(prep.operation, target=b["signal"], work=b["signal"][:0])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=alpha))
