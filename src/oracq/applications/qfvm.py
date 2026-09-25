"""Build QFVM sparse inputs and a replaceable QLSS problem from a QRAM flow field and Roe arithmetic."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import (
    SparseAccess,
    StatePreparation,
    XorDatabase,
    abstract_database,
    abstract_state_prep,
    annotate,
    invoke,
    qram_database,
    qram_state_prep,
    resources_for,
)
from oracq.algorithms.input_model.sparse import compare_words, value_transposition
from oracq.algorithms.qlss.qlss import LinearSystem, QLSSProtocol, SolveResult, SpectralPromise
from oracq.applications.flow_data import RoeFlowData
from oracq.applications.roe import ArithmeticBuilder, Word, roe_face
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Adjoint, Bits, Program, Ref, ValidationError, fuse
from oracq.infrastructure.linking import Binding, bind


@dataclass(frozen=True)
class RoeQfvmInputs:
    """The set of abstract slots declared by the QFVM input model: raw field values, geometry, residual values and residual state preparation.

    The slots themselves carry no data; implementations are bound by
    ``bind_qfvm`` and runtime memory tables are materialized by
    ``qfvm_memories``. Matrix elements are not pre-stored in any slot but are
    computed on the fly at query time by reversible Roe arithmetic.

    Attributes:
        cell_width: Bit width of cell indices; the grid has ``2**cell_width`` cells.
        fmt: Fixed-point format used by conserved variables and residual values.
        angle_width: Bit width of rotation angle words.
        rho: Abstract declared slot of the density (first conserved variable) database.
        momentum: Abstract declared slot of the momentum (second conserved variable) database.
        energy: Abstract declared slot of the energy (third conserved variable) database.
        geometry: Slot of the geometry table database; encodes only sparse positions and indices, no matrix values.
        theta: Slot of the conversion table from residual fixed-point words to rotation angle words (reserved data plane).
        rhs: Abstract slot of the normalized residual state preparation.
        residual: Slot of the database of quantized residual values per (cell, component).
    """

    cell_width: int
    fmt: FixedFormat
    angle_width: int
    rho: XorDatabase
    momentum: XorDatabase
    energy: XorDatabase
    geometry: XorDatabase
    theta: XorDatabase
    rhs: StatePreparation
    residual: XorDatabase

    @property
    def width(self) -> int:
        """Bit width of full matrix coordinates: 2 component bits + cell_width cell bits + 1 dilation half-block flag bit."""
        return self.cell_width + 3

    @property
    def geometry_width(self) -> int:
        """Bit width of geometry table words: seven packed segments for partner coordinate, symmetric slot, source cell, row and column components, source band and the valid bit."""
        return self.width + 4 + self.cell_width + 2 + 2 + 2 + 1


def roe_qfvm_inputs(
    *,
    cell_width: int = 2,
    fmt: FixedFormat | None = None,
    angle_width: int = 8,
    prefix: str = "RoeQfvm",
) -> RoeQfvmInputs:
    """Declare the full set of abstract QFVM input slots: six databases plus one RHS state preparation.

    Everything returned is an open slot without data; data only appears at
    binding time (``bind_qfvm``) and execution time (``qfvm_memories``).

    Args:
        cell_width: Bit width of cell indices, in the range 2..25.
        fmt: Fixed-point format of conserved variables; defaults to ``FixedFormat(6, 2)``.
        angle_width: Bit width of rotation angle words, in the range 1..64.
        prefix: Prefix of the module names of the abstract slots.

    Returns:
        RoeQfvmInputs: The declared set, containing the three conserved-variable
        databases, the geometry table, the theta table, the residual value
        database and the abstract RHS preparation.

    Raises:
        ValidationError: cell_width or angle_width out of range.
    """
    fmt = fmt or FixedFormat(6, 2)
    if not 2 <= cell_width <= 25 or not 1 <= angle_width <= 64:
        raise ValidationError("QFVM requires at least four cells and geometry or angle word lengths of at most 64 bits")
    width = cell_width + 3
    geometry_width = width + 4 + cell_width + 7
    return RoeQfvmInputs(
        cell_width,
        fmt,
        angle_width,
        abstract_database(prefix + "Rho", cell_width, fmt.width),
        abstract_database(prefix + "Momentum", cell_width, fmt.width),
        abstract_database(prefix + "Energy", cell_width, fmt.width),
        abstract_database(prefix + "Geometry", width + 4, geometry_width),
        abstract_database(prefix + "Theta", fmt.width, angle_width),
        abstract_state_prep(prefix + "Rhs", width, cell_width + 2 + angle_width + 1),
        abstract_database(prefix + "RhsValues", cell_width + 2, fmt.width),
    )


def geometry_cells(inputs: RoeQfvmInputs) -> dict[int, int]:
    """Encode only sparse positions and raw data indices; contains no flow field matrix values.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``, determining coordinate widths and the geometry table layout.

    Returns:
        dict[int, int]: The geometry table contents; the address is a matrix
        coordinate concatenated with a 4-bit slot number, and the value is the
        seven-segment packed word of partner coordinate, reverse slot, source
        cell, row and column components, source band and valid bit.
    """
    cw, w = inputs.cell_width, inputs.width
    n = 1 << cw
    table: dict[int, int] = {}
    for row in range(1 << w):
        cell, var, half = (row >> 2) % n, row % 4, row >> (cw + 2)
        for slot in range(16):
            band, other_var = divmod(slot, 3)
            valid = int(slot < 9 and var < 3)
            neighbor_cell = (cell + band - 1) % n
            neighbor = other_var + 4 * neighbor_cell + ((1 - half) << (cw + 2))
            reverse = (2 - band) * 3 + var if valid else 0
            source = cell if half == 0 else neighbor_cell
            rowvar, colvar, source_band = (
                (var, other_var, band) if half == 0 else (other_var, var, 2 - band)
            )
            values = (
                (neighbor, w),
                (reverse, 4),
                (source, cw),
                (rowvar, 2),
                (colvar, 2),
                (source_band, 2),
                (valid, 1),
            )
            packed, offset = 0, 0
            for value, width in values:
                packed |= (value & ((1 << width) - 1)) << offset
                offset += width
            table[row + (slot << w)] = packed
    return table


def ptheta_cells(fmt: FixedFormat, angle_width: int, amax: float) -> dict[int, int]:
    """Generate the theta table: conversion from residual fixed-point words to rotation angle words ``2*acos(min(1, |v|/amax))``.

    Args:
        fmt: The fixed-point format; addresses cover all ``2**fmt.width`` raw words of the format.
        angle_width: Bit width of output angle words.
        amax: Magnitude upper bound of elements, used to clamp the acos argument.

    Returns:
        dict: Mapping from raw fixed-point words to quantized angle words; angle words wrap modulo ``2**angle_width``.

    Raises:
        ValueError: amax non-positive.
    """
    if amax <= 0:
        raise ValueError("amax must be positive")
    return {
        raw: round(
            2 * math.acos(min(1, abs(fmt.decode(raw)) / amax)) * (1 << angle_width) / (2 * math.pi)
        )
        % (1 << angle_width)
        for raw in range(1 << fmt.width)
    }


def roe_entry(
    inputs: RoeQfvmInputs,
    *,
    gamma: float = 1.4,
    entropy_delta: float = 0.125,
    mass: float = 1.0,
    dx: float = 1.0,
) -> Operation:
    """Build the reversible arithmetic circuit for a single Roe matrix element (without sparse addressing).

    The three conserved-variable databases are queried for ``source`` and its
    left and right neighbor cells (with periodic boundary wraparound),
    ``roe_face`` computes the flux characteristic components of the two faces,
    and the west, center and east candidate values are assembled by ``band``;
    the mass term is added to the center block diagonal only when the row and
    column components coincide.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``, providing the conserved-variable databases and fixed-point format.
        gamma: Ratio of specific heats.
        entropy_delta: Roe entropy fix coefficient.
        mass: The mass term, applied only to the component diagonal of the center block.
        dx: Grid step size.

    Returns:
        Operation: Registers are source, row, col, band, value, status; value
        writes out the selected matrix element, and status aggregates the
        failure flags of the arithmetic nodes (non-zero means overflow or
        division by zero).
    """
    cw, fmt = inputs.cell_width, inputs.fmt
    fields = (
        ("rho", inputs.rho.operation),
        ("momentum", inputs.momentum.operation),
        ("energy", inputs.energy.operation),
    )
    b = Builder(
        _name("roe_entry", cw, fmt, gamma, entropy_delta, mass, dx, *[o for _, o in fields]),
        {
            "source": Bits(cw),
            "row": Bits(2),
            "col": Bits(2),
            "band": Bits(2),
            "value": Bits(fmt.width),
            "status": Bits(2),
        },
        resources_for(*fields),
        attributes={
            "algorithm": "qfvm_arithmetic_entry",
            "correctness": "pending",
            "matrix_table": False,
            "boundary": "periodic",
            "mass": mass,
            "dx": dx,
        },
    )
    g = ArithmeticBuilder(b, fmt)
    addresses: list[Ref] = []
    words: list[list[Ref]] = []
    for offset in (-1, 0, 1):
        addr = g.local(cw)
        b.xor(b["source"], addr)
        b.add_const(addr.reinterpret("uint"), offset % (1 << cw))
        addresses.append(addr)
        state: list[Ref] = []
        for name, op in fields:
            word = g.local()
            invoke(b, op, name, address=addr, data=word)
            state.append(word)
        words.append(state)
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
    # The mass term applies only to the component diagonal of the center block.
    equal = b.local("component_equal", Bits(1))
    from oracq.algorithms.common.arithmetic import BooleanNetwork

    net = BooleanNetwork()
    r, c = net.input("r", 2), net.input("c", 2)
    net.outputs = {"equal": [net.inv(net.any([net.xor(x, y) for x, y in zip(r, c, strict=True)]))]}
    b.call(net.operation(), r=b["row"], c=b["col"], equal=equal)
    west = -faces[0][0] / dx
    center = (faces[1][0] - faces[0][1]) / dx + g.choose(equal, mass, 0)
    east = faces[1][1] / dx
    value = g.choose3(b["band"], [west, center, east])
    return g.finish([(b["value"], value)], b["status"])


def _geometry_refs(ref: Ref, inputs: RoeQfvmInputs) -> list[Ref]:
    """Slice a geometry table word into seven contiguous views by their declared widths."""
    sizes = (inputs.width, 4, inputs.cell_width, 2, 2, 2, 1)
    refs: list[Ref] = []
    offset = 0
    for size in sizes:
        refs.append(ref[offset : offset + size])
        offset += size
    return refs


def roe_qfvm_block_encoding(
    inputs: RoeQfvmInputs,
    *,
    amax: float = 8.0,
    padding_value: float = 1.0,
    **entry_options: float,
) -> BlockEncoding:
    """Explicit conversion from sparse inputs to a block encoding; QFVM itself no longer forces exposing only a block encoding.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``.
        amax: Magnitude upper bound of sparse elements, used as the scaling factor of the block encoding.
        padding_value: The padding diagonal value; must be positive and no greater than ``amax``.
        **entry_options: Arithmetic options forwarded to ``roe_entry`` (gamma, mass and so on).

    Returns:
        BlockEncoding: The real symmetric sparse block encoding of the Hermitian dilation of the Roe matrix.
    """
    from oracq.algorithms.input_model.sparse import real_symmetric_sparse_encoding

    if not 0 < padding_value <= amax:
        raise ValidationError("The padding diagonal value must be positive and no greater than the element bound")
    access = qfvm_sparse_access(inputs, padding_value=padding_value, **entry_options)
    return real_symmetric_sparse_encoding(access, inputs.fmt, amax, diagonal_nonnegative=True)


def rhs_qram_preparation(inputs: RoeQfvmInputs) -> StatePreparation:
    """Build the QRAM preparation of the normalized residual state: amplitudes come from the squared-norm tree, signs are written by phase kickback.

    Amplitudes are obtained by ``qram_state_prep`` querying the ``rhs_angles``
    angle word bank layer by layer; signs use the one-bit ``rhs_sign`` bank,
    turning negative residual components into a pi phase through Load, Z and
    un-Load. The input is the zero state and work returns clean.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``, from which coordinate and angle word widths are taken.

    Returns:
        StatePreparation: target holds the full matrix coordinate (cell_width+3
        bits); amplitudes are written only on the first cell_width+2 bits, and
        the dilation half-block flag stays zero.
    """
    n, aw = inputs.cell_width + 2, inputs.angle_width
    prep = qram_state_prep(n, aw)
    sign = qram_database(n, 1)
    b = Builder(
        _name("roe_rhs", n, aw),
        {"target": Bits(n + 1), "work": Bits(n + aw + 1)},
        resources_for(("prep", prep.operation), ("sign", sign.operation)),
        attributes={"correctness": "pending", "algorithm": "qram_signed_residual_tree"},
    )
    invoke(b, prep.operation, "prep", target=b["target"][:n], work=b["work"][: n + aw])
    flag = b["work"][n + aw :]
    invoke(b, sign.operation, "sign", address=b["target"][:n], data=flag)
    b.z(flag)
    invoke(b, sign.operation, "sign", address=b["target"][:n], data=flag)
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def qram_bindings(inputs: RoeQfvmInputs) -> dict[str, Binding]:
    """Provide the complete binding map from abstract QFVM slots to QRAM implementations.

    The three conserved-variable databases, the geometry table and the theta
    table bind to actual ``qram_database`` libraries (the table resource maps to
    the like-named classical bank); the residual value slot binds to the
    ``RoeResidualValues`` database (rhs_values); the RHS preparation slot binds
    to the ``rhs_qram_preparation`` circuit (angle word and sign resources map
    to rhs_angles and rhs_sign respectively).

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``, providing all abstract slots to bind.

    Returns:
        dict: Mapping from abstract slot module names to ``Binding`` objects, ready to hand to ``bind``.
    """
    mapping: dict[str, Binding] = {}
    for key in ("rho", "momentum", "energy", "geometry", "theta"):
        db = getattr(inputs, key)
        actual = qram_database(db.address_width, db.data_width, name="RoeMemory_" + key)
        mapping[db.operation.module.name] = Binding(actual.operation, {"table": key})
    residual = qram_database(
        inputs.residual.address_width, inputs.residual.data_width, name="RoeResidualValues"
    )
    mapping[inputs.residual.operation.module.name] = Binding(
        residual.operation, {"table": "rhs_values"}
    )
    prep = rhs_qram_preparation(inputs)
    mapping[inputs.rhs.operation.module.name] = Binding(
        prep.operation, {"prep__angles": "rhs_angles", "sign__table": "rhs_sign"}
    )
    return mapping


def bind_qfvm(program: Program, inputs: RoeQfvmInputs) -> Program:
    """Bind the still-unresolved abstract QFVM slots in a program to QRAM implementations and the signed residual tree preparation.

    Args:
        program: The program containing calls to abstract QFVM slots.
        inputs: The slot set declared by ``roe_qfvm_inputs``.

    Returns:
        Program: The program with only unresolved slots replaced; open slots unrelated to QFVM are left untouched.
    """
    from oracq.infrastructure.linking import unresolved

    missing = {r.name for r in unresolved(program)}
    return bind(program, {k: v for k, v in qram_bindings(inputs).items() if k in missing})


def qfvm_memories(
    inputs: RoeQfvmInputs, flow: RoeFlowData, *, amax: float = 8.0
) -> Mapping[str, Sequence[int] | Mapping[int, int]]:
    """Materialize the runtime memory tables of the QFVM path: a flow field snapshot plus the static geometry table and theta table.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``.
        flow: Classical-side flow field data providing the snapshot of the six banks.
        amax: Magnitude upper bound of elements used by the theta table.

    Returns:
        dict: Memory tables handed to the execution backend by name: rho,
        momentum, energy, rhs_values, rhs_sign, rhs_angles, geometry and theta.

    Raises:
        ValidationError: The flow field mismatches the declared cell count, fixed-point format or angle format.
    """
    if (
        flow.n != (1 << inputs.cell_width)
        or flow.fmt != inputs.fmt
        or flow.angle_width != inputs.angle_width
    ):
        raise ValidationError("The QFVM flow field and the oracle disagree on cell count, fixed-point format or angle format")
    result = flow.store.snapshot()
    result["geometry"] = geometry_cells(inputs)
    result["theta"] = ptheta_cells(inputs.fmt, inputs.angle_width, amax)
    return result


def roe_qfvm_problem(
    inputs: RoeQfvmInputs,
    *,
    spectrum: SpectralPromise,
    rhs_norm: float | None = None,
    amax: float = 8.0,
    padding_value: float = 1.0,
    **entry_options: float,
) -> LinearSystem:
    """D=([[0,M],[M.T,0]] on physical coordinates) + padding_value*I_pad.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``.
        spectrum: A ``SpectralPromise`` covering the spectral bounds of the quantized Roe matrix; required.
        rhs_norm: Classical right-hand-side norm; the physical magnitude of the solution can be recovered only when provided.
        amax: Magnitude upper bound of sparse elements.
        padding_value: The padding diagonal value; must be positive and no greater than ``amax``.
        **entry_options: Arithmetic options forwarded to ``roe_entry`` (gamma, mass and so on).

    Returns:
        LinearSystem: The QLSS linear system composed of the sparse access, the
        residual preparation and the extended spectral bounds with the padding
        value merged in.
    """
    from oracq.algorithms.qlss.qlss import LinearSystem, SparseSystem, SpectralPromise

    if not 0 < padding_value <= amax:
        raise ValidationError("The padding diagonal value must be positive and no greater than the element bound")
    if not isinstance(spectrum, SpectralPromise):
        raise ValidationError("The QFVM QLSS input requires an explicit SpectralPromise")
    extended = SpectralPromise(
        max(spectrum.norm_upper, padding_value),
        min(spectrum.sigma_min_lower, padding_value),
        spectrum.evidence,
    )
    sparse = SparseSystem(
        qfvm_sparse_access(inputs, padding_value=padding_value, **entry_options),
        inputs.fmt,
        amax,
        inputs.rhs,
        extended,
        diagonal_nonnegative=True,
        hermitian=True,
    )
    return LinearSystem(
        sparse=sparse,
        physical_width=inputs.cell_width + 2,
        physical_high_value=1,
        rhs_norm=rhs_norm,
        data_assumptions=(
            "coherent XOR QRAM on arbitrary address superpositions",
            "memory snapshot fixed across coherent calls, adjoints and norm probes",
            "raw conserved fields and O(Ns) geometry; no matrix-entry table",
            "signed residual value oracle and norm tree with cached rotation angles",
            "clean reversible RHS preparation; norm known separately",
            "logical local patches; native QRAM banks currently rematerialized",
        ),
    )


def roe_qfvm_step(
    inputs: RoeQfvmInputs,
    qlss: QLSSProtocol,
    *,
    spectrum: SpectralPromise | None = None,
    rhs_norm: float | None = None,
    **options: float,
) -> SolveResult:
    """Hand the QFVM problem to a replaceable QLSS protocol for solving.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``.
        qlss: A ``QLSSProtocol`` that declares ``input_model``.
        spectrum: A ``SpectralPromise`` covering the actually quantized matrix; required.
        rhs_norm: Classical right-hand-side norm; only when provided can the solution magnitude be recovered via ``SolveResult.recover_norm``.
        **options: Remaining keyword arguments forwarded to ``roe_qfvm_problem``, including amax, padding_value and matrix element arithmetic options.

    Returns:
        SolveResult: The solution state, norm probe and adaptation records produced by the protocol.

    Raises:
        ValidationError: qlss is not a ``QLSSProtocol`` declaring ``input_model``,
            or spectrum is not provided.
    """
    from oracq.algorithms.qlss.qlss import QLSSProtocol

    if not isinstance(qlss, QLSSProtocol):
        raise ValidationError(
            "Use a QLSSProtocol that declares input_model; a bare callable cannot express an input model"
        )
    if spectrum is None:
        raise ValidationError(
            "Replacing the QFVM QLSS requires minimum singular value and norm declarations: spectrum=SpectralPromise"
        )
    return qlss.solve(roe_qfvm_problem(inputs, spectrum=spectrum, rhs_norm=rhs_norm, **options))


def qfvm_sparse_access(
    inputs: RoeQfvmInputs, *, padding_value: float = 1.0, **entry_options: float
) -> SparseAccess:
    """Build the CKS sparse access for QFVM: two oracles, an in-place position permutation and an arbitrary-coordinate matrix element XOR.

    The location oracle queries the nine structural slots of each column and
    completes them into a full permutation via nine coherent value
    transpositions; the padded component has no neighbor structure, and its
    partner coordinate is rewritten as ``column + rank`` to yield the padding
    diagonal. The entry oracle selects the source cell and component with the
    geometry table and calls ``roe_entry`` to compute the element on the fly:
    entries outside the structural domain are zero, entries with failed
    arithmetic are zeroed, and the diagonal of padded coordinates writes
    padding_value.

    Args:
        inputs: The slot set declared by ``roe_qfvm_inputs``.
        padding_value: The padding diagonal value; must be a positive number exactly representable in the fixed-point format.
        **entry_options: Arithmetic options forwarded to ``roe_entry`` (gamma, mass and so on).

    Returns:
        SparseAccess: Coordinate width ``inputs.width``, value width
        ``inputs.fmt.width``, sparsity 9.

    Raises:
        ValidationError: padding_value is non-positive or not exactly representable in the fixed-point format.
    """

    if padding_value <= 0 or inputs.fmt.decode(inputs.fmt.encode(padding_value)) != padding_value:
        raise ValidationError("The padding diagonal value must be a positive number exactly representable in the fixed-point format")
    n = inputs.width
    geometry = inputs.geometry.operation
    locator = Builder(
        _name("qfvm_sparse_location", geometry, padding_value),
        {"column": Bits(n), "index": Bits(n), "work": Bits(0)},
        resources_for(("geometry", geometry)),
        attributes={
            "sparsity": 9,
            "orientation": "column",
            "padding": "identity on variable 3",
            "full_permutation_extension": True,
            "construction": "nine coherent transpositions",
        },
    )
    neighbors: list[Ref] = []
    lefts: list[Ref] = []
    for rank in range(9):
        slot = locator.local("slot_" + str(rank), Bits(4))
        data = locator.local("geometry_" + str(rank), Bits(inputs.geometry_width))
        neighbor = locator.local("neighbor_" + str(rank), Bits(n))
        left = locator.local("left_" + str(rank), Bits(n))
        for bit in range(4):
            if (rank >> bit) & 1:
                locator.x(slot[bit])
        invoke(locator, geometry, "geometry", address=fuse(locator["column"], slot), data=data)
        locator.xor(data[:n], neighbor)
        # The padded variable has diagonal padding_value; the other eight distinct positions return zero elements.
        padded = locator.local("padded_neighbor_" + str(rank), Bits(n))
        locator.xor(locator["column"], padded)
        locator.add_const(padded.reinterpret("uint"), rank)
        with locator.control(locator["column"][:2], 3):
            locator.xor(data[:n], neighbor)
            locator.xor(padded, neighbor)
        for bit in range(n):
            if (rank >> bit) & 1:
                locator.x(left[bit])
        for previous in range(rank):
            locator.call(
                value_transposition(n), index=left, a=lefts[previous], b=neighbors[previous]
            )
        neighbors.append(neighbor)
        lefts.append(left)
    setup = tuple(locator._frames[0])
    for left, right in zip(lefts, neighbors, strict=True):
        locator.call(value_transposition(n), index=locator["index"], a=left, b=right)
    locator.emit(Adjoint(setup))
    location = annotate(locator.finish(), "sparse_location_inplace")

    physical = roe_entry(inputs, **entry_options)
    b = Builder(
        _name("qfvm_sparse_entry", geometry, physical, padding_value),
        {"row": Bits(n), "column": Bits(n), "data": Bits(inputs.fmt.width)},
        resources_for(("geometry", geometry), ("physical", physical)),
        attributes={
            "value_encoding": "signed_fixed_point",
            "value_fraction": inputs.fmt.fraction,
            "matrix": "Hermitian dilation of Roe matrix, identity on padded coordinates",
            "invalid_arithmetic": "entry totalized to zero",
            "correctness": "pending",
        },
    )
    selected = b.local("selected_geometry", Bits(inputs.geometry_width))
    for rank in range(9):
        slot = b.local("slot_" + str(rank), Bits(4))
        geom = b.local("geometry_" + str(rank), Bits(inputs.geometry_width))
        match = b.local("match_" + str(rank), Bits(1))
        for bit in range(4):
            if (rank >> bit) & 1:
                b.x(slot[bit])
        invoke(b, geometry, "geometry", address=fuse(b["column"], slot), data=geom)
        b.call(compare_words(n), a=b["row"], b=geom[:n], flag=match)
        # The valid bit of the raw geometry excludes the padded variable and invalid slots.
        with b.control(fuse(match, geom[inputs.geometry_width - 1])):
            b.xor(geom[: inputs.geometry_width - 1], selected[: inputs.geometry_width - 1])
            b.x(selected[inputs.geometry_width - 1])
        b.call(compare_words(n), a=b["row"], b=geom[:n], flag=match)
        invoke(b, geometry, "geometry", address=fuse(b["column"], slot), data=geom)
    _, _, source, row, col, band, valid = _geometry_refs(selected, inputs)
    value = b.local("computed_value", Bits(inputs.fmt.width))
    status = b.local("arithmetic_status", Bits(2))
    invoke(
        b,
        physical,
        "physical",
        source=source,
        row=row,
        col=col,
        band=band,
        value=value,
        status=status,
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
    entry = annotate(b.finish(), "sparse_entry_xor")
    return SparseAccess(location, entry, n, inputs.fmt.width, 9)
