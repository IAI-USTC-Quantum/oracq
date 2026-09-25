"Open, bindable assembly examples for QFVM and QHAM."

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
    """All open input declarations and layout constants of the QFVM assembly example.

    Attributes:
        position: Open XOR-database declaration for row-coordinate queries.
        reverse_slot: Open XOR-database declaration for reverse column-slot
            queries within a row.
        column_position: Open XOR-database declaration for column-coordinate
            queries.
        column_reverse_slot: Open XOR-database declaration for reverse slot
            queries within a column.
        flow: Open XOR-database declaration for flux-word queries.
        boundary: Open XOR-database declaration for boundary-word queries.
        physical_entry: Open declaration of the reversible physical-entry
            function.
        amplitude: Open declaration of the amplitude transduction.
        residual: Residual state preparation declaration.
        cell_width: Bit width of spatial cell coordinates.
        slot_width: Bit width of column slots.
        field_width: Width of flux and boundary words.
        value_width: Width of physical entry values.
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
    """Build all open input declarations of the QFVM assembly example.

    Args:
        name: Declaration name prefix, prepended to each slot module name.
        cell_width: Bit width of spatial cell coordinates.
        slot_width: Bit width of column slots.
        field_width: Width of flux and boundary words.
        value_width: Width of physical entry values.

    Returns:
        QfvmInputs: Input bundle with the six XOR databases, the
        physical-entry and amplitude open declarations, and the residual
        preparation.
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
    """Build the row-wise or column-wise QFVM access-isometry preparation module."""
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
    """Compose the QFVM access isometries into a block encoding with the ``T_L†SWAP T_R`` structure.

    Args:
        inputs: Open input bundle built by ``qfvm_inputs``.
        alpha: Normalization constant written into the ``be_alpha`` attribute.

    Returns:
        BlockEncoding: Block encoding view with spatial cells as the target
        and the access-graph workspace as the signal.
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
    """Hand the QFVM block encoding and residual preparation to the linear-system solver protocol.

    Args:
        inputs: Open input bundle built by ``qfvm_inputs``.
        qlss: Solver protocol accepting a block encoding and a state
            preparation and returning a ``StateOracle``.
        alpha: Normalization constant passed to the block encoding.

    Returns:
        StateOracle: QFVM stepping solution-state oracle produced by the
        solver protocol.
    """
    return qlss(qfvm_block_encoding(inputs, alpha=alpha), inputs.residual)


def embed_block(
    a: BlockEncoding, global_width: int, row_offset: int, column_offset: int
) -> BlockEncoding:
    """Embed a block encoding at a diagonal-block position of a larger global matrix.

    The target is first wrapped by subtracting the column offset and the
    padding region flagged on the high bits, then ``a`` is invoked and the row
    offset added to the whole register.

    Args:
        a: The embedded block encoding.
        global_width: Bit width of the global target register.
        row_offset: Row offset of the target block in the global matrix.
        column_offset: Column offset of the target block in the global
            matrix.

    Returns:
        BlockEncoding: Block encoding with target width ``global_width`` and
        the same alpha as ``a``.

    Raises:
        ValidationError: The embedded block crosses the global matrix
            boundary.
    """
    d, size = 1 << global_width, 1 << a.width
    if not 0 <= row_offset <= d - size or not 0 <= column_offset <= d - size:
        raise ValidationError("the embedded block crosses the global matrix boundary")
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
    """Lift the linear and folding terms into block encodings of the QHAM m=1 closed system.

    The linear terms sit on diagonal blocks along the three coordinate blocks
    and two two-coordinate tensor-product terms, the folding terms appear on
    two coupling blocks with weight ``-h*hfun``, and the whole is combined by
    LCU.

    Args:
        linear: Block encoding of the single-coordinate linear operator.
        fold: Block encoding of the two-coordinate folding operator; its
            width must be twice that of ``linear``.
        h: Folding coupling constant.
        hfun: Scaling factor of the folding coefficient; the product
            ``h*hfun`` serves as the coupling weight.

    Returns:
        BlockEncoding: LCU-combined block encoding of the lifted system.

    Raises:
        ValidationError: The folding matrix is not padded to the
            two-coordinate space.
    """
    n = linear.width
    if fold.width != 2 * n:
        raise ValidationError("the folding matrix must first be padded to the two-coordinate space")
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
    """Build the initial-state amplitude tuple of the QHAM m=1 lifted system.

    Arranged as ``[v, v, 0, v⊗v]`` and zero-padded to a power-of-two length
    for state preparation binding.

    Args:
        values: Sequence of physical-channel initial-state amplitudes.

    Returns:
        tuple[complex, ...]: Complex amplitude tuple zero-padded to a
        power-of-two length.
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
    """Assemble and solve the QHAM m=1 lifted system, selecting the physical-channel subspace.

    Args:
        linear: Block encoding of the single-coordinate linear operator.
        fold: Block encoding of the two-coordinate folding operator.
        initial: Physical-channel initial-state preparation.
        solver: Solver protocol for the lifted system; when ``via_pde`` is
            true it is called as ``(DiscretePDE, final time)``, otherwise as
            ``(block encoding, preparation, final time)``.
        final_time: Final evolution time.
        h: Folding coupling constant.
        hfun: Scaling factor of the folding coefficient.
        via_pde: Whether to wrap in ``DiscretePDE`` before handing over to
            the solver protocol.

    Returns:
        StateOracle: Solution-state oracle after selecting the physical
        subspace.
    """
    lifted = qham_lift_m1(linear, fold, h=h, hfun=hfun)
    state = (
        solver(DiscretePDE(lifted, initial, "qham_lifted_pde"), final_time)
        if via_pde
        else solver(lifted, initial, final_time)
    )
    return select_subspace(state, linear.width, 0, label="qham_physical_channel")
