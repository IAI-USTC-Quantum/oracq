"""Linear ODE input model for Fokker–Planck/SDE and a pure-Python classical witness.

The probability density of the one-dimensional Ito SDE ``dx = a(x) dt + sqrt(2 D(x)) dW``
satisfies the Fokker–Planck equation ``p' = -∂x(a p) + ∂xx(D p)``. This module performs a
zero-flux finite-volume discretization on a uniform grid, yielding a discrete generator G
with zero column sums; once G becomes a BlockEncoding via an explicit Pauli expansion, a
QODEProblem/LinearODE can be assembled and handed to the existing linear ODE solvers
(LCHS etc.). The dissipativity of the generator is an input-model declaration, not proved
by the language.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding
from oracq.algorithms.input_model.contracts import finite_real, positive_integer
from oracq.algorithms.input_model.operators import BlockEncoding, scale
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    gate_state_prep,
    qram_state_angles,
    qram_state_prep,
    uniform_state,
)
from oracq.algorithms.qode.ode import QODEProblem
from oracq.algorithms.qode.ode_models import HermitianParts, LinearODE
from oracq.infrastructure.ir import ValidationError


def _coefficients(
    values: int | float | Iterable[float],
    size: int,
    path: str,
    *,
    nonnegative: bool = False,
) -> tuple[float, ...]:
    """Unify a constant or a pointwise vector into a tuple of finite reals."""
    if type(values) in (int, float):
        values = cast("Iterable[float]", (values,) * size)
    result = tuple(cast("Iterable[float]", values))
    if len(result) != size:
        raise ValidationError(path + " length must equal the number of grid points")
    for value in result:
        finite_real(value, path, minimum=0 if nonnegative else None)
    return tuple(float(v) for v in result)


def _uniform_points(grid: Iterable[float]) -> tuple[tuple[float, ...], float]:
    """Validate strictly increasing uniform grid coordinates and return (points, spacing)."""
    points = tuple(grid)
    if len(points) < 2:
        raise ValidationError("FokkerPlanckProblem.grid requires at least two grid points")
    for value in points:
        finite_real(value, "FokkerPlanckProblem.grid")
    spacing = (points[-1] - points[0]) / (len(points) - 1)
    if spacing <= 0:
        raise ValidationError("FokkerPlanckProblem.grid must be strictly increasing")
    tolerance = 1e-9 * max(1.0, abs(spacing))
    for left, right in zip(points, points[1:], strict=False):
        if abs((right - left) - spacing) > tolerance:
            raise ValidationError("Only uniform grids are currently supported; non-uniform grids are out of scope")
    return tuple(float(p) for p in points), float(spacing)


def _power_of_two_width(size: int, path: str) -> int:
    """Validate that the number of grid points is a power of two of at least 2 and return the corresponding register bit count."""
    if size < 2 or size & (size - 1):
        raise ValidationError(path + " requires a power-of-two number of grid points to match the quantum register width")
    return (size - 1).bit_length()


@dataclass(frozen=True)
class FokkerPlanckProblem:
    """Zero-flux finite-volume discretization input model of the conservative-form Fokker–Planck operator.

    Args:
        drift: Drift coefficient a(x); a constant or a pointwise vector.
        diffusion: Diffusion coefficient D(x) >= 0; a constant or a pointwise vector.
        grid: Uniformly increasing grid coordinate sequence.
    """

    drift: int | float | Iterable[float]
    diffusion: int | float | Iterable[float]
    grid: Iterable[float]

    def __post_init__(self) -> None:
        """Validate the uniform grid and normalize the drift and diffusion coefficients into pointwise tuples."""
        points, spacing = _uniform_points(self.grid)
        object.__setattr__(self, "points", points)
        object.__setattr__(self, "spacing", spacing)
        size = len(points)
        object.__setattr__(
            self, "drift", _coefficients(self.drift, size, "FokkerPlanckProblem.drift")
        )
        object.__setattr__(
            self,
            "diffusion",
            _coefficients(self.diffusion, size, "FokkerPlanckProblem.diffusion", nonnegative=True),
        )

    @property
    def size(self) -> int:
        """Number of grid points; also the dimension of the discrete generator matrix."""
        return len(self.points)  # type: ignore[attr-defined]

    @property
    def width(self) -> int:
        """Quantum register bit count for the grid, equal to ``log2(size)``.

        Raises ``ValidationError`` when the number of grid points is not a power of two."""
        return _power_of_two_width(self.size, "FokkerPlanckProblem")

    def generator_matrix(self) -> tuple[tuple[float, ...], ...]:
        """Return the discrete generator G (column-vector convention ``p' = G p``, zero column sums).

        Returns:
            tuple[tuple[float, ...], ...]: A square matrix of grid dimension, as row-nested tuples.
        """
        size, h = self.size, self.spacing  # type: ignore[attr-defined]
        g = [[0.0] * size for _ in range(size)]
        for face in range(size - 1):
            a_face = 0.5 * (
                cast("tuple[float, ...]", self.drift)[face]
                + cast("tuple[float, ...]", self.drift)[face + 1]
            )
            d_face = 0.5 * (
                cast("tuple[float, ...]", self.diffusion)[face]
                + cast("tuple[float, ...]", self.diffusion)[face + 1]
            )
            forward = (0.5 * a_face + d_face / h) / h
            backward = (0.5 * a_face - d_face / h) / h
            g[face][face] -= forward
            g[face + 1][face] += forward
            g[face][face + 1] -= backward
            g[face + 1][face + 1] += backward
        return tuple(tuple(row) for row in g)

    def generator_encoding(self) -> BlockEncoding:
        """BlockEncoding via a small-scale explicit Pauli expansion; no quantum speedup is claimed.

        Returns:
            BlockEncoding: The explicit Pauli-expansion block encoding of the discrete generator G.
        """
        if self.width > 5:
            raise ValidationError("The explicit Pauli expansion is only for small grids of at most 32 points; larger instances require an access oracle")
        return matrix_pauli_encoding(self.generator_matrix())

    def qode_problem(self, initial: StatePreparation | None = None) -> QODEProblem:
        """Assemble the QODEProblem for ``p' = G p``; dissipative is a caller declaration.

        Args:
            initial: Initial state preparation handle; defaults to the uniform distribution state over the grid width, and its width must match the grid.

        Returns:
            QODEProblem: The quantum ODE problem driven by the generator G, tagged with zero-flux finite-volume evidence.
        """
        initial = initial if initial is not None else uniform_state(self.width)
        if not isinstance(initial, StatePreparation):
            raise ValidationError("FokkerPlanckProblem.qode_problem requires a StatePreparation initial state")
        if initial.width != self.width:
            raise ValidationError("Initial state width does not match the Fokker–Planck grid width")
        return QODEProblem(
            self.generator_encoding(),
            initial,
            dissipative=True,
            evidence="fokker_planck_zero_flux_finite_volume; dissipative caller-declared",
        )

    def linear_ode(self, initial: StatePreparation | None = None) -> LinearODE:
        """A ``u' = -A u`` LinearODE view of the same problem (A = -G).

        Args:
            initial: Initial state preparation handle; defaults to the uniform distribution state over the grid width.

        Returns:
            LinearODE: The linear ODE view with A = -G, labeled ``fokker_planck_minus_A``.
        """
        initial = initial if initial is not None else uniform_state(self.width)
        return LinearODE(
            HermitianParts.from_operator(scale(-1, self.generator_encoding())),
            initial,
            label="fokker_planck_minus_A",
        )


def _probabilities(probabilities: Iterable[float], path: str) -> tuple[float, ...]:
    """Validate a non-negative probability sequence of power-of-two length and normalize it."""
    values = tuple(probabilities)
    if len(values) < 2 or len(values) & (len(values) - 1):
        raise ValidationError(path + " requires a power-of-two length with at least two points")
    for value in values:
        finite_real(value, path, minimum=0)
    total = math.fsum(values)
    if total <= 0:
        raise ValidationError(path + " must have a positive sum")
    return tuple(float(v) / total for v in values)


def sde_state_preparation(
    probabilities: Iterable[float],
    *,
    implementation: str = "gate",
    angle_width: int = 8,
    work_width: int = 0,
) -> StatePreparation:
    """Encode a discrete initial distribution as a StatePreparation with amplitudes ``sqrt(p_i)``.

    Args:
        probabilities: Non-negative discrete probabilities of power-of-two length.
        implementation: ``"gate"`` for the ordinary reused-rotation implementation, ``"qram"`` for the QRAM angle-table implementation.
        angle_width: Angle-table word width of the QRAM implementation.
        work_width: Extra work bit width of the gate implementation.

    Returns:
        StatePreparation: Preparation handle of the discrete distribution with amplitudes ``sqrt(p_i)``.
    """
    values = _probabilities(probabilities, "sde_state_preparation.probabilities")
    amplitudes = [math.sqrt(v) for v in values]
    if implementation == "gate":
        return gate_state_prep(amplitudes, work_width=work_width)
    if implementation == "qram":
        positive_integer(angle_width, "sde_state_preparation.angle_width")
        return qram_state_prep((len(values) - 1).bit_length(), angle_width)
    raise ValidationError("Unknown state preparation implementation: " + repr(implementation))


def sde_state_angles(probabilities: Iterable[float], *, angle_width: int = 8) -> dict[int, int]:
    """Angle-table binding data of the QRAM implementation; keys are rotation-tree node addresses.

    Args:
        probabilities: Non-negative discrete probabilities of power-of-two length; normalized internally.
        angle_width: Quantization bit width of the angle-table words; a positive integer.

    Returns:
        dict[int, int]: Binding data keyed by rotation-tree node address with quantized angle integers as values.
    """
    values = _probabilities(probabilities, "sde_state_angles.probabilities")
    positive_integer(angle_width, "sde_state_angles.angle_width")
    return qram_state_angles([math.sqrt(v) for v in values], angle_width)


def _square_matrix(matrix: Iterable[Iterable[float]], path: str) -> tuple[tuple[float, ...], ...]:
    """Validate a nonempty square matrix with finite entries and convert it to a float-tuple representation."""
    result = tuple(tuple(row) for row in matrix)
    size = len(result)
    if size < 1 or any(len(row) != size for row in result):
        raise ValidationError(path + " must be a nonempty square matrix")
    for row in result:
        for value in row:
            finite_real(value, path)
    return tuple(tuple(float(v) for v in row) for row in result)


def _matvec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    """Compute the matrix-vector product, accumulating row by row with ``math.fsum``."""
    return [math.fsum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def _matmul(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> list[list[float]]:
    """Compute the product of two square matrices of equal order, accumulating entry by entry with ``math.fsum``."""
    size = len(a)
    columns = [[b[i][j] for i in range(size)] for j in range(size)]
    return [[math.fsum(x * y for x, y in zip(row, col, strict=True)) for col in columns] for row in a]


def matrix_exponential(matrix: Iterable[Iterable[float]], time: float = 1.0) -> list[list[float]]:
    """Pure-Python matrix exponential by scaling-and-squaring plus Taylor; only for small-scale classical witnesses.

    Args:
        matrix: A nonempty square matrix with finite real entries.
        time: Duration scaling factor; a finite real number.

    Returns:
        list[list[float]]: The matrix exp(matrix*time), as row-nested lists.
    """
    finite_real(time, "matrix_exponential.time")
    a = [list(row) for row in _square_matrix(matrix, "matrix_exponential.matrix")]
    size = len(a)
    for i in range(size):
        for j in range(size):
            a[i][j] *= time
    norm = max((math.fsum(abs(v) for v in row) for row in a), default=0.0)
    halvings = max(0, math.ceil(math.log2(norm / 0.5))) if norm > 0.5 else 0
    if halvings:
        shrink = 2.0**halvings
        a = [[v / shrink for v in row] for row in a]
    result = [[float(i == j) for j in range(size)] for i in range(size)]
    term = [row[:] for row in result]
    for order in range(1, 200):
        term = [[v / order for v in row] for row in _matmul(term, a)]
        result = [[x + y for x, y in zip(rx, tx, strict=True)] for rx, tx in zip(result, term, strict=True)]
        if max(abs(v) for row in term for v in row) < 1e-17:
            break
    for _ in range(halvings):
        result = _matmul(result, result)
    return result


def evolve_distribution(
    matrix: Iterable[Iterable[float]],
    initial: Iterable[float],
    time: float,
    *,
    steps: int | None = None,
) -> list[float]:
    """Classical reference evolution; uses the matrix exponential when steps is None, explicit Euler otherwise.

    Args:
        matrix: Generator square matrix under the convention ``p' = matrix * p``.
        initial: Initial vector of the same dimension as the matrix.
        time: Evolution duration; a non-negative finite real number.
        steps: Explicit Euler step count; when omitted, the solve is a single matrix-exponential call.

    Returns:
        list[float]: The distribution vector evolved to the given time.
    """
    g = _square_matrix(matrix, "evolve_distribution.matrix")
    initial = tuple(initial)
    if len(initial) != len(g):
        raise ValidationError("evolve_distribution.initial length does not match the matrix")
    for value in initial:
        finite_real(value, "evolve_distribution.initial")
    finite_real(time, "evolve_distribution.time", minimum=0)
    vector = [float(v) for v in initial]
    if steps is None:
        return _matvec(matrix_exponential(g, time), vector)
    positive_integer(steps, "evolve_distribution.steps")
    dt = time / steps
    for _ in range(steps):
        delta = _matvec(g, vector)
        vector = [v + dt * dv for v, dv in zip(vector, delta, strict=True)]
    return vector


def distribution_moments(
    points: Iterable[float],
    probabilities: Iterable[float],
    orders: tuple[int, ...] = (1, 2),
) -> tuple[float, ...]:
    """Compute the moments ``<x^k>`` from grid coordinates and a probability vector; defaults to ``(<x>, <x^2>)``.

    Args:
        points: Grid coordinate sequence, of the same length as the probability vector.
        probabilities: Per-point probability weights; the sum must be positive; normalized internally.
        orders: Moment orders to compute; a tuple of positive integers.

    Returns:
        tuple[float, ...]: Normalized moments ``<x^k>`` corresponding one-to-one with the orders.
    """
    points = tuple(points)
    values = tuple(probabilities)
    if len(points) != len(values) or not points:
        raise ValidationError("distribution_moments requires nonempty coordinates and probabilities of equal length")
    for value in (*points, *values):
        finite_real(value, "distribution_moments")
    if any(type(k) is not int or k < 1 for k in orders):
        raise ValidationError("distribution_moments.orders requires positive integer orders")
    total = math.fsum(values)
    if total <= 0:
        raise ValidationError("The distribution_moments probability sum must be positive")
    return tuple(
        math.fsum(p * x**k for x, p in zip(points, values, strict=True)) / total for k in orders
    )


def _interface_ratio(problem: FokkerPlanckProblem, face: int) -> float:
    """Effective grid Péclet number at a zero-flux face, u = a h / (2 D)."""
    if not isinstance(problem, FokkerPlanckProblem):
        raise ValidationError("A FokkerPlanckProblem is required")
    a_face = 0.5 * (
        cast("tuple[float, ...]", problem.drift)[face]
        + cast("tuple[float, ...]", problem.drift)[face + 1]
    )
    d_face = 0.5 * (
        cast("tuple[float, ...]", problem.diffusion)[face]
        + cast("tuple[float, ...]", problem.diffusion)[face + 1]
    )
    if d_face <= 0:
        raise ValidationError("The stationary reference requires strictly positive diffusion coefficients at the interfaces")
    return a_face * problem.spacing / (2 * d_face)  # type: ignore[attr-defined]


def stationary_distribution(problem: FokkerPlanckProblem) -> tuple[float, ...]:
    """Exact stationary state of the zero-flux discretization; adjacent-point ratios are ``(1+u)/(1-u)``.

    Args:
        problem: A normalized Fokker–Planck discrete problem; the diffusion coefficients at all faces must be positive.

    Returns:
        tuple[float, ...]: Discrete stationary probabilities, one per grid point, summing to one.
    """
    weights = [1.0]
    for face in range(problem.size - 1):
        u = _interface_ratio(problem, face)
        if abs(u) >= 1:
            raise ValidationError("The grid Péclet number is too large; the central-difference stationary state is no longer positive")
        weights.append(weights[-1] * (1 + u) / (1 - u))
    total = math.fsum(weights)
    return tuple(w / total for w in weights)


def boltzmann_distribution(problem: FokkerPlanckProblem) -> tuple[float, ...]:
    """Continuous stationary reference ``p ∝ exp(∫ a/D dx)``, differing from the discrete stationary state by O(h^2).

    Args:
        problem: A normalized Fokker–Planck discrete problem; the diffusion coefficients at all faces must be positive.

    Returns:
        tuple[float, ...]: Continuous stationary reference probabilities, one per grid point, summing to one.
    """
    weights = [1.0]
    for face in range(problem.size - 1):
        u = _interface_ratio(problem, face)
        weights.append(weights[-1] * math.exp(2 * u))
    total = math.fsum(weights)
    return tuple(w / total for w in weights)
