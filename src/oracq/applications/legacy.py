"QFVM 与 QHAM 的可开放、可绑定组装案例。"

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass

from oracq.algorithms.common.state_preparation import select_subspace
from oracq.algorithms.input_model.block_encoding import lcu, tensor
from oracq.algorithms.input_model.operators import BlockEncoding, _name, identity
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    XorDatabase,
    abstract_database,
    abstract_state_prep,
    annotate,
    declare,
    invoke,
    resources_for,
)
from oracq.algorithms.qpde.pde import DiscretePDE
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Ref, ValidationError, fuse


@dataclass(frozen=True)
class QfvmInputs:
    """QFVM 组装案例的全部开放输入声明与布局常量。

    Attributes:
        position: 行坐标查询的 XOR 数据库开放声明。
        reverse_slot: 行内反向列槽查询的 XOR 数据库开放声明。
        column_position: 列坐标查询的 XOR 数据库开放声明。
        column_reverse_slot: 列内反向槽查询的 XOR 数据库开放声明。
        flow: 通量字查询的 XOR 数据库开放声明。
        boundary: 边界字查询的 XOR 数据库开放声明。
        physical_entry: 物理条目可逆函数的开放声明。
        amplitude: 振幅转导的开放声明。
        residual: 残差态制备声明。
        cell_width: 空间单元坐标位宽。
        slot_width: 列槽位宽。
        field_width: 通量与边界字宽。
        value_width: 物理条目值字宽。
    """

    position: XorDatabase
    reverse_slot: XorDatabase
    column_position: XorDatabase
    column_reverse_slot: XorDatabase
    flow: XorDatabase
    boundary: XorDatabase
    physical_entry: Operation
    amplitude: Operation
    residual: StatePreparation
    cell_width: int
    slot_width: int
    field_width: int
    value_width: int


def qfvm_inputs(
    name: str = "Qfvm",
    *,
    cell_width: int = 1,
    slot_width: int = 1,
    field_width: int = 1,
    value_width: int = 2,
) -> QfvmInputs:
    """构造 QFVM 组装案例的全部开放输入声明。

    Args:
        name: 声明名前缀，拼接到各槽位模块名前。
        cell_width: 空间单元坐标位宽。
        slot_width: 列槽位宽。
        field_width: 通量与边界字宽。
        value_width: 物理条目值字宽。

    Returns:
        QfvmInputs: 含六个 XOR 数据库、物理条目与振幅开放声明及残差制备的输入束。
    """
    n, s, w, v = cell_width, slot_width, field_width, value_width
    entry_registers = {
        "row": Bits(n),
        "slot": Bits(s),
        "column": Bits(n),
        "flow": Bits(w),
        "boundary": Bits(w),
        "value": Bits(v),
        "status": Bits(1),
    }
    return QfvmInputs(
        abstract_database(name + "Position", n + s, n),
        abstract_database(name + "ReverseSlot", n + s, s),
        abstract_database(name + "ColumnPosition", n + s, n),
        abstract_database(name + "ColumnReverseSlot", n + s, s),
        abstract_database(name + "Flow", n, w),
        abstract_database(name + "Boundary", n, w),
        declare(
            name + "PhysicalEntry",
            entry_registers,
            paradigm="reversible_function",
            attributes={"domain_contract": "frozen Roe or tabulated-entry implementation"},
        ),
        declare(
            name + "Amplitude",
            {"value": Bits(v), "amplitude": Bits(1)},
            paradigm="reversible_function",
        ),
        abstract_state_prep(name + "Residual", n),
        n,
        s,
        w,
        v,
    )


def _qfvm_prepare(inputs: QfvmInputs, *, column: bool = False) -> Operation:
    """构造行向或列向的 QFVM 访问等距制备模块。"""
    n, s, w, v = inputs.cell_width, inputs.slot_width, inputs.field_width, inputs.value_width
    position = inputs.column_position if column else inputs.position
    reverse = inputs.column_reverse_slot if column else inputs.reverse_slot
    operations = [
        ("position", position.operation),
        ("reverse", reverse.operation),
        ("flow", inputs.flow.operation),
        ("boundary", inputs.boundary.operation),
        ("entry", inputs.physical_entry),
        ("amplitude", inputs.amplitude),
    ]
    width = s + n + s + 2 * w + v + 2
    b = Builder(
        _name(
            "qfvm_prepare_column" if column else "qfvm_prepare_row", *(op for _, op in operations)
        ),
        {"target": Bits(n), "signal": Bits(width)},
        resources_for(*operations),
        attributes={"algorithm": "qfvm_access_isometry", "validation_stage": "paradigm"},
    )
    offset = 0
    views: dict[str, Ref] = {}
    for name, size in [
        ("slot", s),
        ("column", n),
        ("reverse", s),
        ("flow", w),
        ("boundary", w),
        ("value", v),
        ("status", 1),
        ("amplitude", 1),
    ]:
        views[name] = b["signal"][offset : offset + size]
        offset += size
    b.h(views["slot"])
    address = fuse(b["target"], views["slot"])
    invoke(b, position.operation, "position", address=address, data=views["column"])
    invoke(b, reverse.operation, "reverse", address=address, data=views["reverse"])
    invoke(b, inputs.flow.operation, "flow", address=b["target"], data=views["flow"])
    invoke(
        b, inputs.boundary.operation, "boundary", address=views["column"], data=views["boundary"]
    )
    arguments = {
        "row": b["target"],
        "slot": views["slot"],
        "column": views["column"],
        "flow": views["flow"],
        "boundary": views["boundary"],
        "value": views["value"],
        "status": views["status"],
    }
    invoke(b, inputs.physical_entry, "entry", **arguments)
    invoke(b, inputs.amplitude, "amplitude", value=views["value"], amplitude=views["amplitude"])
    with b.adjoint():
        invoke(b, inputs.flow.operation, "flow", address=b["target"], data=views["flow"])
        invoke(
            b,
            inputs.boundary.operation,
            "boundary",
            address=views["column"],
            data=views["boundary"],
        )
        invoke(b, inputs.physical_entry, "entry", **arguments)
    return b.finish()


def qfvm_block_encoding(inputs: QfvmInputs, *, alpha: float = 2.0) -> BlockEncoding:
    """以 ``T_L†SWAP T_R`` 结构把 QFVM 访问等距组合成块编码。

    Args:
        inputs: ``qfvm_inputs`` 构造的开放输入束。
        alpha: 写入 ``be_alpha`` 属性的归一化常数。

    Returns:
        BlockEncoding: 目标为空间单元、信号为访问图工作区的块编码视图。
    """
    right, left = _qfvm_prepare(inputs), _qfvm_prepare(inputs, column=True)
    width = right.module.registers[1].type.width
    b = Builder(
        _name("qfvm_be", right, left, alpha),
        {"target": Bits(inputs.cell_width), "signal": Bits(width)},
        resources_for(("right", right), ("left", left)),
    )
    invoke(b, right, "right", target=b["target"], signal=b["signal"])
    s, n = inputs.slot_width, inputs.cell_width
    b.swap(b["target"], b["signal"][s : s + n])
    b.swap(b["signal"][:s], b["signal"][s + n : s + n + s])
    with b.adjoint():
        invoke(b, left, "left", target=b["target"], signal=b["signal"])
    return BlockEncoding(
        annotate(
            b.finish(),
            "block_encoding",
            be_alpha=alpha,
            application="qfvm",
            construction="T_L_dagger_SWAP_T_R",
            validation_stage="paradigm",
            contract_status="access graph and layout implemented; Roe and normalization validation pending",
        )
    )


def qfvm_step(
    inputs: QfvmInputs, qlss: Callable[..., StateOracle], *, alpha: float = 2.0
) -> StateOracle:
    """把 QFVM 块编码与残差制备交给线性系统求解协议。

    Args:
        inputs: ``qfvm_inputs`` 构造的开放输入束。
        qlss: 接受块编码与态制备并返回 ``StateOracle`` 的求解协议。
        alpha: 传入块编码的归一化常数。

    Returns:
        StateOracle: 求解协议产出的 QFVM 步进解态 oracle。
    """
    return qlss(qfvm_block_encoding(inputs, alpha=alpha), inputs.residual)


def embed_block(
    a: BlockEncoding, global_width: int, row_offset: int, column_offset: int
) -> BlockEncoding:
    """把块编码嵌入到更大全局矩阵的指定对角块位置。

    先把目标回绕减去列偏移、在高位标记补齐区，再调用 ``a`` 并整体加上行偏移。

    Args:
        a: 被嵌入的块编码。
        global_width: 全局目标寄存器位宽。
        row_offset: 目标块在全局矩阵中的行偏移。
        column_offset: 目标块在全局矩阵中的列偏移。

    Returns:
        BlockEncoding: 目标宽为 ``global_width``、alpha 与 ``a`` 相同的块编码。

    Raises:
        ValidationError: 嵌入块越过全局矩阵范围。
    """
    d, size = 1 << global_width, 1 << a.width
    if not 0 <= row_offset <= d - size or not 0 <= column_offset <= d - size:
        raise ValidationError("嵌入块越过全局矩阵范围")
    b = Builder(
        _name("embed", a.operation, global_width, row_offset, column_offset),
        {"target": Bits(global_width), "signal": Bits(a.signal_qubits + 1)},
        resources_for(("a", a.operation)),
    )
    b.add_const(b["target"].reinterpret("uint"), (-column_offset) % d)
    high, flag = b["target"][a.width :], b["signal"][a.signal_qubits :]
    if high.width:
        b.x(flag)
        with b.control(high, 0):
            b.x(flag)
    invoke(
        b, a.operation, "a", target=b["target"][: a.width], signal=b["signal"][: a.signal_qubits]
    )
    b.add_const(b["target"].reinterpret("uint"), row_offset)
    return BlockEncoding(
        annotate(
            b.finish(),
            "block_encoding",
            be_alpha=a.alpha,
            row_offset=row_offset,
            column_offset=column_offset,
        )
    )


def qham_lift_m1(
    linear: BlockEncoding, fold: BlockEncoding, *, h: float = -0.25, hfun: float = 1.0
) -> BlockEncoding:
    """把线性项与折叠项提升为 QHAM m=1 封闭系统的块编码。

    线性项沿三个坐标块与两个双坐标张量积项放在对角块上，折叠项以
    ``-h*hfun`` 权重出现在两个耦合块上，整体经 LCU 组合。

    Args:
        linear: 单坐标线性算子的块编码。
        fold: 双坐标折叠算子的块编码，宽度须为 ``linear`` 的两倍。
        h: 折叠耦合常数。
        hfun: 折叠系数的缩放因子，乘积 ``h*hfun`` 作为耦合权重。

    Returns:
        BlockEncoding: LCU 组合后的提升系统块编码。

    Raises:
        ValidationError: 折叠矩阵未补齐到双坐标空间。
    """
    n = linear.width
    if fold.width != 2 * n:
        raise ValidationError("折叠矩阵必须先补齐到双坐标空间")
    d = 1 << n
    width = (3 * d + d * d - 1).bit_length()
    eta = h * hfun
    terms: list[tuple[float, BlockEncoding]] = [
        (1, embed_block(linear, width, offset, offset)) for offset in (0, d, 2 * d)
    ]
    terms += [
        (1, embed_block(tensor(linear, identity(n)), width, 3 * d, 3 * d)),
        (1, embed_block(tensor(identity(n), linear), width, 3 * d, 3 * d)),
        (-eta, embed_block(fold, width, 0, 3 * d)),
        (-eta, embed_block(fold, width, 2 * d, 3 * d)),
    ]
    return lcu(terms)


def qham_initial_vector(values: Iterable[complex]) -> tuple[complex, ...]:
    """构造 QHAM m=1 提升系统的初态振幅元组。

    按 ``[v, v, 0, v⊗v]`` 排列后补零到二的幂长度，供态制备绑定使用。

    Args:
        values: 物理通道初态振幅序列。

    Returns:
        tuple[complex, ...]: 补零到二的幂长度的复振幅元组。
    """
    values = tuple(complex(v) for v in values)
    d = len(values)
    raw = [*values, *values, *([0j] * d), *(a * b for a in values for b in values)]
    raw += [0j] * ((1 << (len(raw) - 1).bit_length()) - len(raw))
    return tuple(raw)


def qham_m1(
    linear: BlockEncoding,
    fold: BlockEncoding,
    initial: StatePreparation,
    solver: Callable[..., StateOracle],
    *,
    final_time: float = 0.1,
    h: float = -0.25,
    hfun: float = 1.0,
    via_pde: bool = False,
) -> StateOracle:
    """组装并求解 QHAM m=1 提升系统，选出物理通道子空间。

    Args:
        linear: 单坐标线性算子的块编码。
        fold: 双坐标折叠算子的块编码。
        initial: 物理通道初态制备。
        solver: 提升系统求解协议；``via_pde`` 为真时按
            ``(DiscretePDE, 终时刻)`` 调用，否则按 ``(块编码, 制备, 终时刻)``
            调用。
        final_time: 演化终时刻。
        h: 折叠耦合常数。
        hfun: 折叠系数的缩放因子。
        via_pde: 是否经 ``DiscretePDE`` 封装后交给求解协议。

    Returns:
        StateOracle: 选取物理子空间后的解态 oracle。
    """
    lifted = qham_lift_m1(linear, fold, h=h, hfun=hfun)
    state = (
        solver(DiscretePDE(lifted, initial, "qham_lifted_pde"), final_time)
        if via_pde
        else solver(lifted, initial, final_time)
    )
    return select_subspace(state, linear.width, 0, label="qham_physical_channel")
