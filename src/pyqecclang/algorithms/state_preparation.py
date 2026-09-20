"""初态扩展、编码作用与物理子空间选择的组合工具。"""

from __future__ import annotations

from pyqecclang.algorithms.operators import _name
from pyqecclang.algorithms.oracles import (
    StateOracle,
    StatePreparation,
    annotate,
    invoke,
    resources_for,
)
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse


def extend_initial(prep, extra_width):
    """把态制备 ``prep`` 的目标空间扩展 ``extra_width`` 个高位。

    原制备作用于扩展后 target 的低位，新增高位保持零；零输入承诺保留。

    Args:
        prep: 原始 ``StatePreparation``。
        extra_width: 追加的高位量子位数。

    Returns:
        StatePreparation: 目标宽度为 ``prep.width + extra_width`` 的制备。
    """
    b = Builder(
        _name("extend_initial", prep.operation, extra_width),
        {"target": Bits(prep.width + extra_width), "work": Bits(prep.work_width)},
        resources_for(("initial", prep.operation)),
    )
    invoke(b, prep.operation, "initial", target=b["target"][: prep.width], work=b["work"])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def select_subspace(state: StateOracle, output_width, high_value=0, *, label="selection"):
    """从态 oracle 的目标空间中选出物理子空间。

    低 ``output_width`` 位作为输出，其余高位并入 signal 并要求等于
    ``high_value``；成功条件为整个 signal 寄存器复零。

    Args:
        state: 宽度不小于 ``output_width`` 的输入 ``StateOracle``。
        output_width: 选出的物理输出宽度。
        high_value: 被丢弃高位必须等于的值，默认为零。
        label: 生成模块名与算法属性使用的标签。

    Returns:
        StateOracle: 输出宽度为 ``output_width``、signal 含原信号位与
        附加高位的态 oracle。

    Raises:
        ValidationError: ``output_width`` 超出 ``state.width`` 或
            ``high_value`` 越界。
    """
    extra = state.width - output_width
    if extra < 0 or not 0 <= high_value < 1 << extra:
        raise ValidationError("输出子空间布局无效")
    b = Builder(
        _name(label, state.operation, output_width, high_value),
        {"target": Bits(output_width), "signal": Bits(state.signal_qubits + extra)},
        resources_for(("state", state.operation)),
        attributes={
            "algorithm": label,
            "selected_high_value": high_value,
            "validation_stage": "paradigm",
            "success_condition": "signal == 0",
        },
    )
    old_signal, high = b["signal"][: state.signal_qubits], b["signal"][state.signal_qubits :]
    invoke(b, state.operation, "state", target=fuse(b["target"], high), signal=old_signal)
    for bit in range(extra):
        if (high_value >> bit) & 1:
            b.x(high[bit])
    return StateOracle(b.finish())


def apply_be_to_state(a, prep):
    """把块编码 ``a`` 作用到已制备的初态上，得到态 oracle。

    ``a`` 的信号位占 signal 低位，``prep`` 的工作区并入其高位；成功条件
    为整个 signal 寄存器复零。

    Args:
        a: 作用算子的 ``BlockEncoding``。
        prep: 与 ``a`` 同目标宽度的 ``StatePreparation``。

    Returns:
        StateOracle: 目标宽度不变、附带信号寄存器的态 oracle。

    Raises:
        ValidationError: 两者目标宽度不同。
    """
    if a.width != prep.width:
        raise ValidationError("算子与制备的目标宽度不符")
    b = Builder(
        _name("apply_be_state", a.operation, prep.operation),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits + prep.work_width)},
        resources_for(("a", a.operation), ("prep", prep.operation)),
    )
    invoke(b, prep.operation, "prep", target=b["target"], work=b["signal"][a.signal_qubits :])
    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"][: a.signal_qubits])
    return StateOracle(b.finish())
