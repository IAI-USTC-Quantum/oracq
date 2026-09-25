"""Direct QMem rewrite of the QFVM data path: state and geometry tables use pointer-style two-dimensional addressing, and the residual tree uses QVector.

The circuit structure is isomorphic to applications/qfvm.py (reversible Roe
arithmetic, nine-slot geometry, padding diagonal and in-place position
permutation are all preserved), but it no longer goes through abstract
database slots and bind: modules declare QRAM-form resources directly, the
three conserved variables are merged into a single (field, cell) state table,
neighbors cell±1 realize the periodic boundary through pointer arithmetic
modulo 2^cell_width, the geometry table is addressed two-dimensionally by
(column, slot), and the residual state is prepared by QVector (squared-norm
tree plus sign kickback). Semantic equivalence with the existing path is
cross-checked amplitude by amplitude on simulate by tests/core/test_qfvm_qmem.py.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import cast

from oracq.algorithms.common.arithmetic import BooleanNetwork, FixedFormat
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import StatePreparation, annotate
from oracq.algorithms.input_model.qdata import QVector
from oracq.algorithms.input_model.sparse import compare_words, value_transposition
from oracq.applications.flow_data import RoeFlowData
from oracq.applications.qfvm import RoeQfvmInputs, geometry_cells
from oracq.applications.roe import ArithmeticBuilder, Word, roe_face
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import QRAM, Adjoint, Bits, Ref, ValidationError, fuse
from oracq.infrastructure.qmem import QMem, QPtr


def _geometry_refs(ref: Ref, inputs: RoeQfvmInputs) -> list[Ref]:
    """Slice a geometry table word into seven contiguous views by their declared widths."""
    sizes = (inputs.width, 4, inputs.cell_width, 2, 2, 2, 1)
    refs: list[Ref] = []
    offset = 0
    for size in sizes:
        refs.append(ref[offset : offset + size])
        offset += size
    return refs


def qfvm_qmem_physical(
    inputs: RoeQfvmInputs,
    *,
    gamma: float = 1.4,
    entropy_delta: float = 0.125,
    mass: float = 1.0,
    dx: float = 1.0,
) -> Operation:
    """Physical computation of matrix elements: the state table is addressed two-dimensionally by (field, neighbor cell), periodic neighbors use modular-addition pointers.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``, providing the cell index width and fixed-point format.
        gamma: Ratio of specific heats.
        entropy_delta: Roe entropy fix coefficient.
        mass: The mass term, applied only to the component diagonal of the center block.
        dx: Grid step size.

    Returns:
        Operation: The matrix element computation circuit with registers
        source, row, col, band, value, status; status aggregates the failure
        flags of the arithmetic nodes.
    """
    cw, fmt = inputs.cell_width, inputs.fmt
    b = Builder(
        _name("qfvm_qmem_physical", cw, fmt, gamma, entropy_delta, mass, dx),
        {
            "source": Bits(cw),
            "row": Bits(2),
            "col": Bits(2),
            "band": Bits(2),
            "value": Bits(fmt.width),
            "status": Bits(2),
        },
        {"state": QRAM(cw + 2, fmt.width)},
        attributes={
            "algorithm": "qfvm_arithmetic_entry",
            "correctness": "pending",
            "data_access": "qmem_state_table",
            "boundary": "periodic",
            "mass": mass,
            "dx": dx,
        },
    )
    state = QMem(b, "state", shape=(3, 1 << cw))
    g = ArithmeticBuilder(b, fmt)
    words: list[list[Ref]] = []
    for offset in (-1, 0, 1):
        addr = g.local(cw)
        b.xor(b["source"], addr)
        b.add_const(addr.reinterpret("uint"), offset % (1 << cw))
        row_words: list[Ref] = []
        for field in range(3):
            word = g.local()
            cast(QPtr, state[field, addr]).load(word)
            row_words.append(word)
        words.append(row_words)
    face = roe_face(fmt=fmt, gamma=gamma, entropy_delta=entropy_delta)
    faces: list[tuple[Word, Word]] = []
    for i in range(2):
        left, right, flag = g.local(), g.local(), g.local(2)
        b.call(
            face,
            **dict(  # type: ignore[arg-type]  # dynamic keyword dispatch: mypy cannot rule out the resources parameter
                zip(
                    ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r"),
                    words[i] + words[i + 1],
                    strict=True,
                )
            ),
            row=b["row"],
            col=b["col"],
            left=left,
            right=right,
            status=flag,
        )
        g.flags.append(flag)
        faces.append((g.word(left), g.word(right)))
    equal = b.local("component_equal", Bits(1))
    net = BooleanNetwork()
    r, c = net.input("r", 2), net.input("c", 2)
    net.outputs = {"equal": [net.inv(net.any([net.xor(x, y) for x, y in zip(r, c, strict=True)]))]}
    b.call(net.operation(), r=b["row"], c=b["col"], equal=equal)
    west = -faces[0][0] / dx
    center = (faces[1][0] - faces[0][1]) / dx + g.choose(equal, mass, 0)
    east = faces[1][1] / dx
    value = g.choose3(b["band"], [west, center, east])
    return g.finish([(b["value"], value)], b["status"])


def qfvm_qmem_location(inputs: RoeQfvmInputs) -> Operation:
    """CKS in-place position permutation: nine structural slots addressed through the two-dimensional (column, slot) geometry.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``, determining the coordinate and geometry word widths.

    Returns:
        Operation: The in-place sparse position permutation circuit with registers column, index, work.
    """
    n = inputs.width
    b = Builder(
        _name("qfvm_qmem_location", n),
        {"column": Bits(n), "index": Bits(n), "work": Bits(n)},
        {"geometry": QRAM(n + 4, inputs.geometry_width)},
        attributes={
            "sparsity": 9,
            "orientation": "column",
            "padding": "identity on variable 3",
            "full_permutation_extension": True,
            "construction": "nine coherent transpositions",
            "data_access": "qmem_geometry_2d",
        },
    )
    geometry = QMem(b, "geometry", shape=(16, 1 << n))
    neighbors: list[Ref] = []
    lefts: list[Ref] = []
    for rank in range(9):
        slot = b.local("slot_" + str(rank), Bits(4))
        data = b.local("geometry_" + str(rank), Bits(inputs.geometry_width))
        neighbor = b.local("neighbor_" + str(rank), Bits(n))
        left = b.local("left_" + str(rank), Bits(n))
        for bit in range(4):
            if (rank >> bit) & 1:
                b.x(slot[bit])
        cast(QPtr, geometry[slot, b["column"]]).load(data)
        b.xor(data[:n], neighbor)
        padded = b.local("padded_neighbor_" + str(rank), Bits(n))
        b.xor(b["column"], padded)
        b.add_const(padded.reinterpret("uint"), rank)
        with b.control(b["column"][:2], 3):
            b.xor(data[:n], neighbor)
            b.xor(padded, neighbor)
        for bit in range(n):
            if (rank >> bit) & 1:
                b.x(left[bit])
        for previous in range(rank):
            b.call(value_transposition(n), index=left, a=lefts[previous], b=neighbors[previous])
        neighbors.append(neighbor)
        lefts.append(left)
    setup = tuple(b._frames[0])
    for left, right in zip(lefts, neighbors, strict=True):
        b.call(value_transposition(n), index=b["index"], a=left, b=right)
    b.emit(Adjoint(setup))
    return annotate(b.finish(), "sparse_location_inplace")


def qfvm_qmem_entry(
    inputs: RoeQfvmInputs, *, padding_value: float = 1.0, **entry_options: float
) -> Operation:
    """Arbitrary-coordinate matrix element oracle: same semantics as the entry of qfvm_sparse_access, with the data plane wired directly to QMem.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``.
        padding_value: The padding diagonal value; must be a positive number exactly representable in the fixed-point format.
        **entry_options: Arithmetic options forwarded to ``qfvm_qmem_physical`` (gamma, mass and so on).

    Returns:
        Operation: The matrix element XOR circuit with registers row, column,
        data; entries outside the structural domain are zero, entries with
        failed arithmetic are zeroed, and padded coordinates write the padding
        diagonal.
    """
    if padding_value <= 0 or inputs.fmt.decode(inputs.fmt.encode(padding_value)) != padding_value:
        raise ValidationError("The padding diagonal value must be a positive number exactly representable in the fixed-point format")
    n = inputs.width
    physical = qfvm_qmem_physical(inputs, **entry_options)
    b = Builder(
        _name("qfvm_qmem_entry", n, inputs.fmt, padding_value),
        {"row": Bits(n), "column": Bits(n), "data": Bits(inputs.fmt.width)},
        {
            "geometry": QRAM(n + 4, inputs.geometry_width),
            "state": QRAM(inputs.cell_width + 2, inputs.fmt.width),
        },
        attributes={
            "value_encoding": "signed_fixed_point",
            "value_fraction": inputs.fmt.fraction,
            "matrix": "Hermitian dilation of Roe matrix, identity on padded coordinates",
            "invalid_arithmetic": "entry totalized to zero",
            "correctness": "pending",
            "data_access": "qmem_state_geometry",
        },
    )
    geometry = QMem(b, "geometry", shape=(16, 1 << n))
    selected = b.local("selected_geometry", Bits(inputs.geometry_width))
    for rank in range(9):
        slot = b.local("slot_" + str(rank), Bits(4))
        geom = b.local("geometry_" + str(rank), Bits(inputs.geometry_width))
        match = b.local("match_" + str(rank), Bits(1))
        for bit in range(4):
            if (rank >> bit) & 1:
                b.x(slot[bit])
        cast(QPtr, geometry[slot, b["column"]]).load(geom)
        b.call(compare_words(n), a=b["row"], b=geom[:n], flag=match)
        with b.control(fuse(match, geom[inputs.geometry_width - 1])):
            b.xor(geom[: inputs.geometry_width - 1], selected[: inputs.geometry_width - 1])
            b.x(selected[inputs.geometry_width - 1])
        b.call(compare_words(n), a=b["row"], b=geom[:n], flag=match)
        cast(QPtr, geometry[slot, b["column"]]).load(geom)
    _, _, source, row, col, band, valid = _geometry_refs(selected, inputs)
    value = b.local("computed_value", Bits(inputs.fmt.width))
    status = b.local("arithmetic_status", Bits(2))
    b.call(
        physical,
        source=source,
        row=row,
        col=col,
        band=band,
        value=value,
        status=status,
        resources={"state": "state"},
    )
    same = b.local("diagonal", Bits(1))
    b.call(compare_words(n), a=b["row"], b=b["column"], flag=same)
    forward = tuple(b._frames[0])
    with b.control(valid):
        with b.control(status, 0):
            b.xor(value, b["data"])
    with b.control(fuse(same, b["column"][:2]), 7):
        raw = inputs.fmt.encode(padding_value)
        for bit in range(inputs.fmt.width):
            if (raw >> bit) & 1:
                b.x(b["data"][bit])
    b.emit(Adjoint(forward))
    return annotate(b.finish(), "sparse_entry_xor")


class RoeQmemData:
    """Classical-side flow field data: reuses the physical quantities and residuals of RoeFlowData and reorganizes the quantum plane into QMem banks."""

    def __init__(
        self,
        states: Sequence[Sequence[float]],
        *,
        fmt: FixedFormat | None = None,
        angle_width: int = 10,
        gamma: float = 1.4,
        entropy_delta: float = 0.125,
        dx: float = 1.0,
        name: str = "residual",
    ) -> None:
        """Initialize the flow field and assemble the residual QVector.

        Args:
            states: The conserved variable triple (density, momentum, energy) of each cell; the number of cells must be a power of two of at least 4.
            fmt: Fixed-point format of the conserved variables; defaults to ``FixedFormat(10, 5)``.
            angle_width: Bit width of the rotation angle words of the residual QVector.
            gamma: Ratio of specific heats.
            entropy_delta: Harten entropy fix threshold.
            dx: Grid step size; must be positive.
            name: Naming prefix of the QRAM banks of the residual QVector.

        Raises:
            ValidationError: Invalid flow field cell count, component count or dx.
        """
        fmt = fmt or FixedFormat(10, 5)
        self.flow: RoeFlowData = RoeFlowData(
            states, fmt=fmt, angle_width=angle_width, gamma=gamma,
            entropy_delta=entropy_delta, dx=dx,
        )
        self.fmt: FixedFormat = fmt
        self.angle_width: int = angle_width
        self.cell_width: int = self.flow.n.bit_length() - 1
        values: list[float] = []
        for cell in range(self.flow.n):
            residual = self.flow.residuals[cell]
            values.extend(residual[j] if j < 3 else 0.0 for j in range(4))
        self.vector: QVector = QVector(values, fmt=fmt, angle_width=angle_width, name=name)

    @property
    def state_bank(self) -> dict[int, int]:
        """(field, cell) state table: address = field*n + cell, and the word is the fixed-point encoding of the conserved variable."""
        banks = self.flow.store.snapshot()
        table: dict[int, int] = {}
        for field, key in enumerate(("rho", "momentum", "energy")):
            for address, word in banks[key].items():
                table[(field << self.cell_width) | address] = word
        return table

    def memories(self, inputs: RoeQfvmInputs) -> Mapping[str, Sequence[int] | Mapping[int, int]]:
        """Export all runtime memory tables of the QMem direct data path.

        Args:
            inputs: The slot set declared by ``roe_qfvm_inputs``, determining the geometry table layout.

        Returns:
            dict: The ``state`` state table, the ``geometry`` geometry table,
            plus the angle word and sign banks of the residual ``QVector``
            snapshot (defaulting to ``residual_angles`` and ``residual_sign``;
            the names follow the ``name`` parameter given at construction); the
            execution entry binds them by name to the QRAM resources declared
            by modules.
        """
        return {
            "state": self.state_bank,
            "geometry": geometry_cells(inputs),
            **self.vector.snapshot(signed=True),
        }


def qfvm_qmem_rhs(data: RoeQmemData) -> StatePreparation:
    """Residual state preparation: QVector squared-norm tree plus sign phase kickback (replacing the qram_state_prep plus sign combination).

    Args:
        data: Classical-side flow field data whose residual QVector provides amplitude and sign information.

    Returns:
        StatePreparation: The preparation circuit of the normalized residual
        state, with negative components carrying the sign via phase kickback.
    """
    return data.vector.preparation(signed=True)
