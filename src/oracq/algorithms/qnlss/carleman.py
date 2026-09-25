"""Finite-order tensor lifting of polynomial ODEs, initial states, and physical channel selection."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from oracq.algorithms.common.state_preparation import select_subspace
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.interfaces import (
    as_state_preparation,
)
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.algorithms.qode._dynamics import tagged
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, ValidationError


@dataclass(frozen=True)
class PolynomialODE:
    """Polynomial ODE input model and the input declarations needed for Carleman assembly.

    Attributes:
        width: Width d of a single-level state vector (register bits).
        coefficients: Sequence of (degree p, block encoding F_p) pairs; degrees are distinct
            and the sequence is nonempty; F_p acts on the d^p-dimensional tensor power, with
            width ``max(1, p) * width``.
        initial: Preparation of the initial state u(0).
        initial_norm: Initial state norm used to weight each tensor level; defaults to 1.0.

    Raises:
        ValidationError: Raised at construction when width/degree/norm values are invalid,
            the coefficient table is empty or has duplicate degrees, an F_p width mismatches,
            or the initial-state layout does not match.
    """

    width: int
    coefficients: tuple[tuple[int, BlockEncoding], ...]
    initial: StatePreparation
    initial_norm: float = 1.0

    def __post_init__(self) -> None:
        """Validate construction-time constraints on width, the coefficient table, and the initial-state layout."""
        positive_integer(self.width, "PolynomialODE.width", maximum=64)
        object.__setattr__(self, "initial", as_state_preparation(self.initial))
        finite_real(self.initial_norm, "PolynomialODE.initial_norm", minimum=0)
        object.__setattr__(self, "coefficients", tuple(self.coefficients))
        if not self.coefficients or len({p for p, _ in self.coefficients}) != len(
            self.coefficients
        ):
            raise ValidationError("PolynomialODE coefficients must be nonempty with distinct degrees")
        if self.initial.width != self.width or self.initial_norm < 0:
            raise ValidationError("Invalid Carleman initial data layout")
        for order, coefficient in self.coefficients:
            positive_integer(order, "PolynomialODE.order", minimum=0)
            require_instance(coefficient, BlockEncoding, "PolynomialODE.coefficient")
            if order < 0 or coefficient.width != max(1, order) * self.width:
                raise ValidationError("F_p must be a d by d^p block encoding padded to max(d,d^p)")


def _carleman_term(
    coefficient: BlockEncoding,
    n: int,
    cutoff: int,
    output_level: int,
    order: int,
    position: int,
) -> BlockEncoding:
    """Assemble a single Carleman placement term: apply F_p at the selected tensor position and maintain the level count."""
    level_bits = cutoff.bit_length()
    data_width, source_level = cutoff * n, output_level + order - 1
    b = Builder(
        _name("carleman_placement", coefficient.operation, cutoff, output_level, order, position),
        {"target": Bits(data_width + level_bits), "signal": Bits(coefficient.signal_qubits + 1)},
        resources_for(("f", coefficient.operation)),
        attributes={
            "carleman_row_level": output_level,
            "carleman_column_level": source_level,
            "tensor_position": position,
            "polynomial_order": order,
            "correctness": "pending",
        },
    )
    data, level, flag = (
        b["target"][:data_width],
        b["target"][data_width:],
        b["signal"][coefficient.signal_qubits],
    )
    b.x(flag)
    with b.control(level, source_level):
        if source_level * n < data.width:
            with b.control(data[source_level * n :], 0):
                b.x(flag)
        else:
            b.x(flag)
    if order == 0:
        for j in reversed(range(position, source_level)):
            b.swap(data[j * n : (j + 1) * n], data[(j + 1) * n : (j + 2) * n])
        group = data[position * n : (position + 1) * n]
    else:
        group = data[position * n : (position + order) * n]
    invoke(
        b, coefficient.operation, "f", target=group, signal=b["signal"][: coefficient.signal_qubits]
    )
    if order > 1:
        # The F_p output lives in the lowest n bits of the group; the remaining zero rows move to the high end of the data region.
        for j in range(position + 1, output_level):
            b.swap(data[j * n : (j + 1) * n], data[(j + order - 1) * n : (j + order) * n])
    delta = source_level ^ output_level
    for bit in range(level_bits):
        if (delta >> bit) & 1:
            b.x(level[bit])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=coefficient.alpha))


def carleman_lift(problem: PolynomialODE, *, cutoff: int = 2) -> BlockEncoding:
    """Assemble the block encoding of the truncated Carleman linear embedding.

    Placement terms are enumerated for each target tensor level and polynomial degree and
    combined by an equal-weight LCU; the output attributes record ``cutoff`` and the
    padded-row zero assumption (padding rows outside the first d rows and columns outside
    d^p are zero).

    Args:
        problem: A ``PolynomialODE`` input model.
        cutoff: Carleman truncation order; must be positive.

    Returns:
        BlockEncoding: The embedding operator acting on ``cutoff * width`` data bits
        plus the level register.

    Raises:
        ValidationError: problem is not a ``PolynomialODE``, or cutoff is invalid.
    """
    require_instance(problem, PolynomialODE, "carleman.problem")
    positive_integer(cutoff, "carleman.cutoff")
    if cutoff < 1:
        raise ValidationError("The Carleman truncation order must be positive")
    terms: list[tuple[complex, BlockEncoding]] = []
    for k in range(1, cutoff + 1):
        for order, coefficient in problem.coefficients:
            source = k + order - 1
            if 0 <= source <= cutoff:
                for position in range(k):
                    terms.append(
                        (1, _carleman_term(coefficient, problem.width, cutoff, k, order, position))
                    )
    out = lcu(terms)
    return BlockEncoding(
        annotate(
            out.operation,
            "block_encoding",
            be_alpha=out.alpha,
            algorithm="carleman_lift",
            cutoff=cutoff,
            correctness="pending",
            coefficient_assumption="F_p padded rows outside first d and columns outside d^p are zero",
        )
    )


def carleman_initial(problem: PolynomialODE, *, cutoff: int = 2) -> StatePreparation:
    """Construct the zero-input preparation of the tensor initial state.

    The level register is prepared weighted by powers of ``initial_norm``, then the initial
    state is copied into the first k levels under control, forming the tensor power vector
    of u(0).

    Args:
        problem: A ``PolynomialODE`` input model.
        cutoff: Carleman truncation order; must be positive.

    Returns:
        StatePreparation: Target is ``cutoff * width`` data bits plus the level register;
        the work space is reused per level at the width of ``initial``.

    Raises:
        ValidationError: problem is not a ``PolynomialODE``, or cutoff is invalid.
    """
    require_instance(problem, PolynomialODE, "carleman.problem")
    positive_integer(cutoff, "carleman.cutoff")
    n, lb = problem.width, cutoff.bit_length()
    weights = [problem.initial_norm**k for k in range(cutoff + 1)]
    weights += [0] * ((1 << lb) - len(weights))
    levels = gate_state_prep(weights)
    b = Builder(
        _name("carleman_initial", problem.initial.operation, problem.initial_norm, cutoff),
        {"target": Bits(n * cutoff + lb), "work": Bits(problem.initial.work_width * cutoff)},
        resources_for(("initial", problem.initial.operation)),
    )
    level = b["target"][n * cutoff :]
    invoke(b, levels.operation, target=level, work=level[:0])
    for k in range(1, cutoff + 1):
        with b.control(level, k):
            for j in range(k):
                invoke(
                    b,
                    problem.initial.operation,
                    "initial",
                    target=b["target"][j * n : (j + 1) * n],
                    work=b["work"][
                        j * problem.initial.work_width : (j + 1) * problem.initial.work_width
                    ],
                )
    return StatePreparation(
        annotate(
            b.finish(),
            "state_prep_isometry",
            zero_input=True,
            algorithm="carleman_tensor_initial",
            correctness="pending",
        )
    )


def carleman_qode(
    problem: PolynomialODE,
    time: float,
    linear_solver: Callable[[BlockEncoding, StatePreparation, float], StateOracle],
    *,
    cutoff: int = 2,
) -> StateOracle:
    """Hand the Carleman embedding to a linear solver and project back to the first level.

    It first assembles ``carleman_lift`` and ``carleman_initial``, then evolves to time
    ``time`` via ``linear_solver``, and finally selects the first tensor level as the
    solution state; the effect of the truncated tail terms is recorded as pending under
    ``truncation_assumption``.

    Args:
        problem: A ``PolynomialODE`` input model.
        time: Target evolution time; non-negative.
        linear_solver: A callable linear solver of the form
            ``(generator, initial, time) -> StateOracle``.
        cutoff: Carleman truncation order; must be positive.

    Returns:
        StateOracle: The first-level solution state of width ``problem.width``.

    Raises:
        ValidationError: Input types or values are invalid, ``linear_solver`` is not
            callable, or its output is not a ``StateOracle``.
    """
    require_instance(problem, PolynomialODE, "carleman.problem")
    finite_real(time, "carleman.time", minimum=0)
    if not callable(linear_solver):
        raise ValidationError("Carleman requires a callable linear solver")
    generator, initial = (
        carleman_lift(problem, cutoff=cutoff),
        carleman_initial(problem, cutoff=cutoff),
    )
    state = linear_solver(generator, initial, time)
    require_instance(state, StateOracle, "carleman.linear_solver.output")
    selected = select_subspace(
        state, problem.width, 1 << ((cutoff - 1) * problem.width), label="carleman_level_one"
    )
    return StateOracle(
        tagged(
            selected.operation,
            "carleman_qode",
            cutoff=cutoff,
            truncation_assumption="Carleman tail pending",
        )
    )
