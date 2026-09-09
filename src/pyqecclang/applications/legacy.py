"QFVM 与 QHAM 的可开放、可绑定组装案例。"

from __future__ import annotations

from dataclasses import dataclass

from pyqecclang.algorithms.block_encoding import lcu, tensor
from pyqecclang.algorithms.operators import BlockEncoding, _name, identity
from pyqecclang.algorithms.oracles import (
    StatePreparation,
    XorDatabase,
    abstract_database,
    abstract_state_prep,
    annotate,
    declare,
    invoke,
    resources_for,
)
from pyqecclang.algorithms.pde import DiscretePDE
from pyqecclang.algorithms.state_preparation import select_subspace
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse


@dataclass(frozen=True)
class QfvmInputs:
    position: XorDatabase
    reverse_slot: XorDatabase
    column_position: XorDatabase
    column_reverse_slot: XorDatabase
    flow: XorDatabase
    boundary: XorDatabase
    physical_entry: object
    amplitude: object
    residual: StatePreparation
    cell_width: int
    slot_width: int
    field_width: int
    value_width: int


def qfvm_inputs(name="Qfvm", *, cell_width=1, slot_width=1, field_width=1, value_width=2):
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


def _qfvm_prepare(inputs, *, column=False):
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
    views = {}
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


def qfvm_block_encoding(inputs: QfvmInputs, *, alpha=2.0):
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


def qfvm_step(inputs, qlss, *, alpha=2.0):
    return qlss(qfvm_block_encoding(inputs, alpha=alpha), inputs.residual)


def embed_block(a: BlockEncoding, global_width, row_offset, column_offset):
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


def qham_lift_m1(linear: BlockEncoding, fold: BlockEncoding, *, h=-0.25, hfun=1.0):
    n = linear.width
    if fold.width != 2 * n:
        raise ValidationError("折叠矩阵必须先补齐到双坐标空间")
    d = 1 << n
    width = (3 * d + d * d - 1).bit_length()
    eta = h * hfun
    terms = [(1, embed_block(linear, width, offset, offset)) for offset in (0, d, 2 * d)]
    terms += [
        (1, embed_block(tensor(linear, identity(n)), width, 3 * d, 3 * d)),
        (1, embed_block(tensor(identity(n), linear), width, 3 * d, 3 * d)),
        (-eta, embed_block(fold, width, 0, 3 * d)),
        (-eta, embed_block(fold, width, 2 * d, 3 * d)),
    ]
    return lcu(terms)


def qham_initial_vector(values):
    values = tuple(complex(v) for v in values)
    d = len(values)
    raw = [*values, *values, *([0j] * d), *(a * b for a in values for b in values)]
    raw += [0j] * ((1 << (len(raw) - 1).bit_length()) - len(raw))
    return tuple(raw)


def qham_m1(linear, fold, initial, solver, *, final_time=0.1, h=-0.25, hfun=1.0, via_pde=False):
    lifted = qham_lift_m1(linear, fold, h=h, hfun=hfun)
    state = (
        solver(DiscretePDE(lifted, initial, "qham_lifted_pde"), final_time)
        if via_pde
        else solver(lifted, initial, final_time)
    )
    return select_subspace(state, linear.width, 0, label="qham_physical_channel")
