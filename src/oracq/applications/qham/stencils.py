"Structured PDE ports on periodic grids: shifted LCU, component selection, and same-point contraction."

from __future__ import annotations

import cmath
import math
from collections.abc import Callable, Sequence

from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.operators import (
    BlockEncoding,
    _name,
    identity,
    product,
    scale,
    zero,
)
from oracq.algorithms.input_model.oracles import (
    abstract_database,
    annotate,
    diagonal_block_encoding,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.algorithms.input_model.qham import PortBinding, QHAMBindings
from oracq.applications.qham.pde import EquationTerm, Monomial
from oracq.applications.qham.reference import Discretization, Grid, centered_weights
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, Ref, ValidationError, fuse


def derivative_encoding(grid: Grid, derivative: tuple[tuple[str, int], ...]) -> BlockEncoding:
    """Build a shifted-LCU block encoding of a spatial derivative on a periodic grid.

    The derivative on each axis is decomposed into cyclic shift terms by the
    centered difference stencil, with term coefficients equal to the stencil
    weights divided by ``spacing**order``; the per-axis operators compose onto
    the whole spatial address register, each acting only on its own axis's
    address bits. An axis without a derivative contributes the identity, and
    an axis whose stencil coefficients cancel entirely contributes the zero
    operator.

    Args:
        grid: A ``Grid`` with periodic boundary and axis lengths that are
            powers of 2.
        derivative: Derivative specification of the form ``((axis, order),
            ...)``.

    Returns:
        BlockEncoding: Block encoding of the difference operator whose target
        is the spatial register of width ``grid.spatial_width``.

    Raises:
        ValidationError: The grid boundary is not periodic, or some axis
            length is not a power of 2.
    """
    if grid.boundary != "periodic" or any(n & (n - 1) for n in grid.shape):
        raise ValidationError("structured shift ports require a periodic grid with power-of-two axis lengths; other boundaries may provide their own BE")
    width = grid.spatial_width
    result = identity(width)
    cursor = 0
    for axis, length, spacing in zip(grid.axes, grid.shape, grid.spacing, strict=True):
        axis_width = (length - 1).bit_length()
        order = dict(derivative).get(axis, 0)
        if order:
            grouped: dict[int, float] = {}
            for offset, coefficient in centered_weights(order):
                step = (-offset) % length
                grouped[step] = grouped.get(step, 0) + coefficient / spacing**order
            terms: list[tuple[float, BlockEncoding]] = []
            for step, coefficient in grouped.items():
                if not coefficient:
                    continue
                b = Builder(
                    _name("fd_shift", grid, axis, step), {"target": Bits(width), "signal": Bits(0)}
                )
                if axis_width:
                    b.add_const(b["target"][cursor : cursor + axis_width].reinterpret("uint"), step)
                encoded = BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))
                terms.append((coefficient, encoded))
            axis_op = lcu(terms) if terms else zero(width)
            result = product(axis_op, result)
        cursor += axis_width
    return result


def coefficient_encoding(
    discretization: Discretization, monomial: Monomial, *, max_words: int = 4096
) -> BlockEncoding:
    """Build a gate-implemented block encoding of the known-coefficient diagonal multiplier.

    The product of the known fields, including their spatial derivatives, and
    the monomial coefficient is evaluated address by address; a single signal
    qubit is rotated under address control with the complex phase compensated
    so that projecting the signal back to zero applies the diagonal value
    divided by ``alpha``; values range over the whole spatial address
    register, with padded addresses outside the grid equal to 0.

    Args:
        discretization: The ``Discretization`` providing the grid and known
            field data.
        monomial: The monomial to encode; degenerates to a constant scaling
            when there are no known fields and the grid has no padding bits.
        max_words: Maximum number of addresses allowed for the gate
            coefficient table.

    Returns:
        BlockEncoding: Block encoding whose target is the spatial register and
        whose ``alpha`` is the largest modulus of the diagonal values; returns
        the zero operator when all diagonal values are zero.

    Raises:
        ValidationError: The coefficient table exceeds the ``max_words``
            budget.
    """
    grid = discretization.grid
    width = grid.spatial_width
    if not monomial.known and grid.size == 1 << width:
        return scale(monomial.coefficient, identity(width))
    if 1 << width > max_words:
        raise ValidationError("the gate implementation of known coefficients exceeds the budget; bind a QRAM or custom coefficient BE instead")
    values = [
        discretization.known_product(monomial, row) if row < grid.size else 0j
        for row in range(1 << width)
    ]
    alpha = max((abs(v) for v in values), default=0)
    if not alpha:
        return zero(width)
    b = Builder(_name("known_multiplier", values), {"target": Bits(width), "signal": Bits(1)})
    for address, value in enumerate(values):
        with b.control(b["target"], address):
            b.ry(b["signal"], 2 * math.acos(min(1, abs(value) / alpha)))
            if value:
                b.global_phase(cmath.phase(value))
    return BlockEncoding(
        annotate(
            b.finish(), "block_encoding", be_alpha=alpha, implementation="known_diagonal_multiplier"
        )
    )


def qram_coefficient_encoding(
    discretization: Discretization,
    monomial: Monomial,
    *,
    angle_width: int = 8,
    max_words: int = 4096,
) -> BlockEncoding:
    """Open angle-database encoding of the known-coefficient diagonal; the same program can bind either a gate or a QRAM database.

    The contract matches coefficient_encoding (target on the spatial bits and
    the same alpha), but the coefficient data is not burned into gates: the
    angle words for diagonal values alpha*cos(theta_a/2) stay in open XOR
    database slots, and the runtime table is computed separately by
    qram_coefficient_memory. Only real coefficient data is accepted; angle
    quantization introduces an amplitude error of at most
    alpha*pi/2**angle_width.

    Args:
        discretization: The ``Discretization`` providing the grid and known
            field data.
        monomial: The monomial to encode; degenerates to a constant scaling
            when there are no known fields and the grid has no padding bits.
        angle_width: Bit width of an angle word, determining the quantization
            precision of diagonal values.
        max_words: Maximum number of addresses allowed for the angle
            database.

    Returns:
        BlockEncoding: Block encoding with diagonal values left in open XOR
        database slots, ``alpha`` equal to the largest modulus of the diagonal
        values; returns the zero operator when all diagonal values are zero.
    """
    grid = discretization.grid
    width = grid.spatial_width
    if not monomial.known and grid.size == 1 << width:
        return scale(monomial.coefficient, identity(width))
    if 1 << width > max_words:
        raise ValidationError("the QRAM angle table of known coefficients exceeds the budget; provide a custom coefficient BE instead")
    values = [
        discretization.known_product(monomial, row) if row < grid.size else 0j
        for row in range(1 << width)
    ]
    if any(v.imag for v in values):
        raise ValidationError("coefficient data for QRAM angle encoding must be real")
    alpha = max((abs(v.real) for v in values), default=0)
    if not alpha:
        return zero(width)
    db = abstract_database(_name("coefficient_angles", values), width, angle_width)
    return diagonal_block_encoding(db, alpha=alpha)


def qram_coefficient_memory(
    discretization: Discretization, monomial: Monomial, *, angle_width: int = 8
) -> dict[int, int]:
    """Runtime angle table matching qram_coefficient_encoding (address -> angle word).

    Args:
        discretization: The ``Discretization`` providing the grid and known
            field data.
        monomial: The same monomial used at encoding time; must agree with the
            encoding call.
        angle_width: Bit width of an angle word; must agree with the encoding
            call.

    Returns:
        dict[int, int]: Mapping from spatial address to angle word; an empty
        table when the encoding degenerates to a constant scaling or all
        diagonal values are zero.
    """
    grid = discretization.grid
    width = grid.spatial_width
    if not monomial.known and grid.size == 1 << width:
        return {}
    values = [
        discretization.known_product(monomial, row) if row < grid.size else 0j
        for row in range(1 << width)
    ]
    alpha = max((abs(v.real) for v in values), default=0)
    if not alpha:
        return {}
    step = 2 * math.pi / (1 << angle_width)
    return {
        address: round(2 * math.acos(min(1, max(-1, v.real / alpha))) / step)
        % (1 << angle_width)
        for address, v in enumerate(values)
    }


def term_encoding(
    discretization: Discretization,
    term: EquationTerm,
    *,
    max_coefficient_words: int = 4096,
    coefficient_encoder: Callable[..., BlockEncoding] = coefficient_encoding,
) -> BlockEncoding:
    """Build the multilinear block encoding of the structured difference port for a single PDE equation term.

    The port is assembled from the centered-difference derivatives of each
    factor field, the known-coefficient diagonal multiplier, component
    selection, and the outer derivative: each factor is differentiated on its
    own spatial bits, and in the multi-factor case the factors are contracted
    to the same point before the diagonal merge; the two reject bits in the
    signal flag input coordinates whose components mismatch and, for multiple
    factors, factors not at the same point, while order-zero terms flag
    nonzero addresses and prepare a uniform superposition on the spatial
    bits. The port matrix is never materialized.

    Args:
        discretization: The ``Discretization`` providing the grid, component
            layout, and known data.
        term: An ``EquationTerm`` carrying the output component and one
            monomial.
        max_coefficient_words: Address budget passed to the coefficient
            encoder.
        coefficient_encoder: Encoding function for the coefficient diagonal
            multiplier, with the same contract as ``coefficient_encoding``
            (e.g. ``qram_coefficient_encoding``).

    Returns:
        BlockEncoding: Rectangular block encoding whose target width is
        ``max(1, arity)*discretization.width``, with the arity recorded in the
        ``rectangular_arity`` attribute; ``alpha`` is the product of the
        alphas of the coefficient, the outer derivative, and each factor
        derivative, with an extra factor ``sqrt(2**spatial_width)`` for
        order-zero terms.

    Raises:
        ValidationError: The port target width exceeds 64, or the coefficient
            encoder reports an exceeded budget.
    """
    n = discretization.width
    ns = discretization.grid.spatial_width
    nc = discretization.component_width
    monomial = term.monomial
    arity = len(monomial.fields)
    width = max(1, arity) * n
    if width > 64:
        raise ValidationError("a single multilinear port exceeds the current BE target packing width")
    derivatives = [derivative_encoding(discretization.grid, a.derivative) for a in monomial.fields]
    outer = derivative_encoding(discretization.grid, monomial.outer_derivative)
    coefficient = coefficient_encoder(discretization, monomial, max_words=max_coefficient_words)
    operands = [(f"d{i}", op.operation) for i, op in enumerate(derivatives)]
    operands += [("coefficient", coefficient.operation), ("outer", outer.operation)]
    signal_width = (
        sum(op.signal_qubits for op in derivatives)
        + coefficient.signal_qubits
        + outer.signal_qubits
        + 2
    )
    b = Builder(
        _name(
            "pde_multilinear_term",
            term,
            discretization.grid,
            discretization.pde.fields,
            coefficient.operation,
            outer.operation,
        ),
        {"target": Bits(width), "signal": Bits(signal_width)},
        resources_for(*operands),
        attributes={
            "pde_arity": arity,
            "input_coordinate_order": "factor 0 low",
            "structured_stencil": True,
        },
    )
    cursor = 0
    signal_views: list[Ref] = []
    for op in derivatives:
        signal_views.append(b["signal"][cursor : cursor + op.signal_qubits])
        cursor += op.signal_qubits
    coeff_signal = b["signal"][cursor : cursor + coefficient.signal_qubits]
    cursor += coefficient.signal_qubits
    outer_signal = b["signal"][cursor : cursor + outer.signal_qubits]
    cursor += outer.signal_qubits
    reject_input, reject_diagonal = b["signal"][cursor], b["signal"][cursor + 1]
    if arity == 0:
        b.x(reject_input)
        with b.control(b["target"], 0):
            b.x(reject_input)
        b.h(b["target"][:ns])
    else:
        component_refs: list[Ref] = []
        for i, atom in enumerate(monomial.fields):
            field = discretization.pde.fields.index(atom.name)
            component = b["target"][i * n + ns : (i + 1) * n]
            for bit in range(nc):
                if (field >> bit) & 1:
                    b.x(component[bit])
            component_refs.append(component)
        if nc:
            b.x(reject_input)
            with b.control(fuse(*component_refs), 0):
                b.x(reject_input)
        for i, op in enumerate(derivatives):
            invoke(
                b,
                op.operation,
                f"d{i}",
                target=b["target"][i * n : i * n + ns],
                signal=signal_views[i],
            )
    invoke(b, coefficient.operation, "coefficient", target=b["target"][:ns], signal=coeff_signal)
    if arity > 1:
        for i in range(1, arity):
            b.xor(b["target"][:ns], b["target"][i * n : i * n + ns])
        b.x(reject_diagonal)
        with b.control(b["target"][n:], 0):
            b.x(reject_diagonal)
    output = discretization.pde.fields.index(term.output)
    for bit in range(nc):
        if (output >> bit) & 1:
            b.x(b["target"][ns + bit])
    invoke(b, outer.operation, "outer", target=b["target"][:ns], signal=outer_signal)
    alpha = coefficient.alpha * outer.alpha * math.prod(op.alpha for op in derivatives)
    if arity == 0:
        alpha *= math.sqrt(1 << ns)
    return BlockEncoding(
        annotate(
            b.finish(),
            "block_encoding",
            be_alpha=alpha,
            rectangular_arity=arity,
            implementation="shifts_components_diagonal_contraction",
        )
    )


def structured_fd_bindings(
    discretization: Discretization,
    initial: Sequence[complex],
    *,
    max_coefficient_words: int = 4096,
    coefficient_encoder: Callable[..., BlockEncoding] = coefficient_encoding,
) -> QHAMBindings:
    """Base matrices are generated from shifts and contractions, without materializing N^r x N^r port matrices.

    Args:
        discretization: The ``Discretization`` providing the grid, component
            layout, and known data.
        initial: Initial vector of length ``discretization.dimension``; when
            its norm is zero the first basis vector is prepared instead.
        max_coefficient_words: Address budget passed to the coefficient
            encoder.
        coefficient_encoder: Encoding function for the coefficient diagonal
            multiplier, with the same contract as ``coefficient_encoding``
            (e.g. ``qram_coefficient_encoding``).

    Returns:
        QHAMBindings: QHAM binding set with each port bound to a shifted-LCU
        block encoding, including the initial-state preparation and norm.
    """
    if len(initial) != discretization.dimension:
        raise ValidationError("initial values require the full register layout")
    ports: list[tuple[str, PortBinding]] = []
    for port in discretization.pde.ports:
        encoded = lcu(
            [
                (
                    1,
                    term_encoding(
                        discretization,
                        term,
                        max_coefficient_words=max_coefficient_words,
                        coefficient_encoder=coefficient_encoder,
                    ),
                )
                for term in port.terms
            ]
        )
        ports.append((port.name, PortBinding(encoded, port.arity)))
    norm = math.sqrt(sum(abs(v) ** 2 for v in initial))
    prep = (
        gate_state_prep(initial)
        if norm
        else gate_state_prep([1.0] + [0.0] * (discretization.dimension - 1))
    )
    return QHAMBindings(discretization.width, tuple(ports), prep, norm)
