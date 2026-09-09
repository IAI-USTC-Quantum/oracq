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
    b = Builder(
        _name("extend_initial", prep.operation, extra_width),
        {"target": Bits(prep.width + extra_width), "work": Bits(prep.work_width)},
        resources_for(("initial", prep.operation)),
    )
    invoke(b, prep.operation, "initial", target=b["target"][: prep.width], work=b["work"])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def select_subspace(state: StateOracle, output_width, high_value=0, *, label="selection"):
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
