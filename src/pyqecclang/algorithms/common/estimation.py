"""相位估计、振幅估计与重叠测量电路；读出和统计在宿主侧完成。"""

from __future__ import annotations

import math
from collections.abc import Iterable

from pyqecclang.algorithms.common.fourier import qft
from pyqecclang.algorithms.input_model.contracts import positive_integer
from pyqecclang.algorithms.input_model.interfaces import (
    BlockEncodingProtocol,
    StatePreparationProtocol,
    as_block_encoding,
    checked_state_preparation,
)
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.algorithms.input_model.oracles import basis_state, invoke, resources_for
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, ValidationError


def phase_estimation(operation: Operation, *, precision: int = 2) -> Operation:
    """标准量子相位估计（QPE）。

    Args:
        operation: 支持受控调用的完整 Operation。应由调用方准备其输入态。
        precision: 相位寄存器位数，范围为 1..63。

    Returns:
        Operation: 保留原公开寄存器，并追加 phase。对本征相位 φ，读出值近似 2**precision * φ。

    Raises:
        ValidationError: 精度无效、输入占用 phase 名字或不支持所需受控调用。

    各次幂以 Repeat 保存；生成器不会按幂次展开门序列。"""
    positive_integer(precision, "phase_estimation.precision", maximum=63)
    registers = {r.name: r.type for r in operation.module.registers}
    if "phase" in registers:
        raise ValidationError("被调接口占用了 phase 参数名")
    b = Builder(
        _name("qpe", operation, precision),
        {**registers, "phase": Bits(precision)},
        resources_for(("u", operation)),
        attributes={"algorithm": "qpe"},
    )
    b.h(b["phase"])
    for bit in range(precision):
        with b.control(b["phase"][bit]):
            with b.repeat(1 << bit):
                invoke(b, operation, "u", **{name: b[name] for name in registers})
    with b.adjoint():
        invoke(b, qft(precision), target=b["phase"])
    return b.finish()


def hadamard_test(
    unitary: BlockEncodingProtocol,
    preparation: StatePreparationProtocol | None = None,
    *,
    component: str = "real",
) -> Operation:
    """生成复期望值的 Hadamard test 电路。

    Args:
        unitary: 完整酉操作，或无信号且 alpha=1 的 BE。
        preparation: 初态制备；省略时使用目标空间的零态。work 必须复净。
        component: ``real`` 或 ``imag``，指定读取期望值的实部或虚部。

    Returns:
        Operation: 寄存器为 target、work、probe。probe 的 Z 期望为相应的复期望分量。

    该操作不执行测量。实际应用需要采样 probe，并在经典侧计算概率差。"""
    encoded = as_block_encoding(unitary)
    if encoded.signal_qubits or encoded.alpha != 1:
        raise ValidationError("Hadamard test 当前要求无信号、alpha=1 的完整 unitary 输入")
    prep = (
        basis_state(encoded.width)
        if preparation is None
        else checked_state_preparation(preparation)
    )
    if prep.width != encoded.width or component not in {"real", "imag"}:
        raise ValidationError("Hadamard test 的初态宽度或分量选项无效")
    b = Builder(
        _name("hadamard_test", encoded.operation, prep.operation, component),
        {"target": Bits(prep.width), "work": Bits(prep.work_width), "probe": Bits(1)},
        resources_for(("u", encoded.operation), ("prep", prep.operation)),
        attributes={
            "algorithm": "hadamard_test",
            "readout_register": "probe",
            "component": component,
        },
    )
    invoke(b, prep.operation, "prep", target=b["target"], work=b["work"])
    b.h(b["probe"])
    with b.control(b["probe"]):
        invoke(b, encoded.operation, "u", target=b["target"], signal=b["work"][:0])
    if component == "imag":
        b.gate("phase", b["probe"], -math.pi / 2)
    b.h(b["probe"])
    return b.finish()


def swap_test(first: StatePreparationProtocol, second: StatePreparationProtocol) -> Operation:
    """生成两个纯态的重叠测量电路。

    Args:
        first: 第一个态的制备，要求零输入和干净工作区。
        second: 第二个同宽态的制备，要求零输入和干净工作区。

    Returns:
        Operation: 保留 left、right、两组 work 和 probe。probe 为零的概率是两态重叠模方加一后除以二。

    两个制备只需正向调用；条件交换由算法生成。"""
    a, c = checked_state_preparation(first), checked_state_preparation(second)
    if a.width != c.width:
        raise ValidationError("Swap test 的两个态必须同宽")
    b = Builder(
        _name("swap_test", a.operation, c.operation),
        {
            "left": Bits(a.width),
            "right": Bits(c.width),
            "left_work": Bits(a.work_width),
            "right_work": Bits(c.work_width),
            "probe": Bits(1),
        },
        resources_for(("a", a.operation), ("b", c.operation)),
        attributes={"algorithm": "swap_test", "readout_register": "probe"},
    )
    invoke(b, a.operation, "a", target=b["left"], work=b["left_work"])
    invoke(b, c.operation, "b", target=b["right"], work=b["right_work"])
    b.h(b["probe"])
    with b.control(b["probe"]):
        b.swap(b["left"], b["right"])
    b.h(b["probe"])
    return b.finish()


def amplitude_estimation(
    preparation: StatePreparationProtocol, marked: Iterable[int], *, precision: int = 3
) -> Operation:
    """生成对好状态概率进行估计的 QPE 电路。

    Args:
        preparation: 支持逆和受控调用的零输入态制备。
        marked: 目标空间中好状态的整数编号集合。
        precision: 相位寄存器位数，范围为 1..63。

    Returns:
        Operation: 寄存器为 target、work、phase。读取 phase 后调用 amplitude_from_phase 解码。

    有限精度读出可能对应多个近似概率。采样与统计处理由宿主负责。"""
    from pyqecclang.algorithms.common.search import grover_iterate

    prep = checked_state_preparation(preparation, adjoint=True, controlled=True)
    iterate = grover_iterate(prep, marked)
    qpe = phase_estimation(iterate, precision=precision)
    b = Builder(
        _name("amplitude_estimation", prep.operation, iterate, precision),
        {"target": Bits(prep.width), "work": Bits(prep.work_width), "phase": Bits(precision)},
        resources_for(("prep", prep.operation), ("qpe", qpe)),
        attributes={
            "algorithm": "amplitude_estimation",
            "readout_register": "phase",
            "decoder": "sin(pi*phase/2**precision)**2",
        },
    )
    invoke(b, prep.operation, "prep", target=b["target"], work=b["work"])
    invoke(b, qpe, "qpe", target=b["target"], work=b["work"], phase=b["phase"])
    return b.finish()


def amplitude_from_phase(value: int, precision: int) -> float:
    """将一个相位样本转换为好状态概率的估计。

    Args:
        value: phase 的无符号整数读出值。
        precision: phase 的位宽。

    Returns:
        float: ``sin(pi*value/2**precision)**2``，位于零和一之间。"""
    positive_integer(precision, "amplitude.precision", maximum=63)
    positive_integer(value, "amplitude.phase", minimum=0, maximum=(1 << precision) - 1)
    return math.sin(math.pi * value / (1 << precision)) ** 2
