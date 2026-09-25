"Classical Roe Riemann local updates and a queryable QRAM data structure; no matrix element table is stored."

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class MemoryPatch:
    """Immutable write record returned by an update.

    Attributes:
        version: Storage version after applying the writes; unchanged when nothing is actually written.
        changes: Storage items that actually changed, mapping addresses to integer values by bank name.
        recomputed_faces: Face indices whose Riemann fluxes were recomputed by this update; empty for direct storage-layer writes.
        recomputed_cells: Cell indices whose residuals were recomputed by this update; empty for direct storage-layer writes.
    """

    version: int
    changes: dict[str, dict[int, int]]
    recomputed_faces: tuple[int, ...] = ()
    recomputed_cells: tuple[int, ...] = ()


class QRAMStore:
    """Versioned multi-bank classical integer storage serving as the data source of a queryable QRAM.

    Attributes:
        banks: Storage contents organized by bank name; each bank maps addresses to integer values.
        version: The version number; incremented once per apply call that produces actual writes.
    """

    def __init__(self) -> None:
        """Initialize an empty multi-bank store with version zero."""
        self.banks: dict[str, dict[int, int]] = {}
        self.version: int = 0

    def apply(self, changes: Mapping[str, Mapping[int, int]]) -> MemoryPatch:
        """Apply a batch of updates and return the record of writes that actually happened.

        Items equal to the current stored value (unwritten addresses count as
        zero) are not recorded in the result; the version increments when actual
        writes exist.

        Args:
            changes: Writes to apply, mapping addresses to integer values by bank name.

        Returns:
            MemoryPatch: The post-write version and actual writes, excluding recomputation ranges.
        """
        actual: dict[str, dict[int, int]] = {}
        for bank, cells in changes.items():
            target = self.banks.setdefault(bank, {})
            for address, value in cells.items():
                if target.get(address, 0) != value:
                    target[address] = value
                    actual.setdefault(bank, {})[address] = value
        if actual:
            self.version += 1
        return MemoryPatch(self.version, actual)

    def snapshot(self) -> dict[str, dict[int, int]]:
        """Return an independent copy of all bank contents.

        Returns:
            dict: A copy of each bank's address-to-value mapping; modifying the snapshot does not affect the store.
        """
        return {k: dict(v) for k, v in self.banks.items()}

    def materialize_changed(
        self, patch: MemoryPatch, factory: Callable[[str, Mapping[int, int]], object]
    ) -> dict[str, object]:
        """PySparQ currently exposes no write interface: only the bank objects that changed are rebuilt, with no claim of physical in-place writes.

        Args:
            patch: The incremental write record; only banks touched by its ``changes`` are rebuilt.
            factory: Factory that builds a host-side bank object from a bank name and its address-to-value mapping.

        Returns:
            dict[str, object]: Mapping of rebuilt bank names to new bank objects; unchanged banks are absent.
        """
        return {bank: factory(bank, self.banks[bank]) for bank in patch.changes}


def _matmul(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> list[list[float]]:
    """Compute the product of two matrices."""
    return [
        [sum(x * y for x, y in zip(row, col, strict=True)) for col in zip(*b, strict=True)]
        for row in a
    ]


def _inverse3(a: Sequence[Sequence[float]]) -> list[list[float]]:
    """Invert a 3x3 matrix by Gaussian elimination with partial pivoting."""
    aug = [list(row) + [float(i == j) for j in range(3)] for i, row in enumerate(a)]
    for j in range(3):
        pivot = max(range(j, 3), key=lambda i: abs(aug[i][j]))
        aug[j], aug[pivot] = aug[pivot], aug[j]
        if abs(aug[j][j]) < 1e-15:
            raise ValidationError("The Roe eigenvector matrix is singular")
        scale = aug[j][j]
        aug[j] = [x / scale for x in aug[j]]
        for i in range(3):
            if i != j:
                scale = aug[i][j]
                aug[i] = [x - scale * y for x, y in zip(aug[i], aug[j], strict=True)]
    return [row[3:] for row in aug]


def riemann_flux(
    left: Sequence[float],
    right: Sequence[float],
    *,
    gamma: float = 1.4,
    entropy_delta: float = 0.125,
) -> tuple[float, float, float]:
    """Compute the classical Roe approximate Riemann flux on a single face.

    The eigenvector matrix ``R`` is built from the Roe average of the left and
    right conserved states; the eigenvalues ``u-c``, ``u`` and ``u+c`` take
    entropy-fixed absolute values, and ``|A| = R*diag(|lambda|)*R**-1`` is
    assembled into the flux formula.

    Args:
        left: The conserved variable triple of the left cell: density, momentum then energy.
        right: The conserved variable triple of the right cell, ordered as in ``left``.
        gamma: Ratio of specific heats.
        entropy_delta: Harten entropy fix threshold; eigenvalues with smaller absolute value use the smoothed value instead.

    Returns:
        tuple: The numerical flux of the three conserved components ``0.5*(f_l+f_r) - 0.5*|A|*(u_r-u_l)``.

    Raises:
        ValidationError: Density non-positive on either side, or the Roe-averaged squared sound speed non-positive.
    """
    def primitive(state: Sequence[float]) -> tuple[float, float, float, list[float]]:
        """Convert a one-sided conserved state to velocity, pressure, enthalpy and the flux vector."""
        rho, m, e = state
        if rho <= 0:
            raise ValidationError("The classical flow field density must be positive")
        u = m / rho
        p = (gamma - 1) * (e - 0.5 * m * u)
        return u, p, (e + p) / rho, [m, m * u + p, u * (e + p)]

    ul, _, hl, fl = primitive(left)
    ur, _, hr, fr = primitive(right)
    wl, wr = math.sqrt(left[0]), math.sqrt(right[0])
    u, h = (wl * ul + wr * ur) / (wl + wr), (wl * hl + wr * hr) / (wl + wr)
    c2 = (gamma - 1) * (h - 0.5 * u * u)
    if c2 <= 0:
        raise ValidationError("The classical Roe squared sound speed must be positive")
    c = math.sqrt(c2)
    r: list[list[float]] = [[1, 1, 1], [u - c, u, u + c], [h - u * c, 0.5 * u * u, h + u * c]]
    vals = [
        abs(x) if abs(x) >= entropy_delta else (x * x + entropy_delta**2) / (2 * entropy_delta)
        for x in (u - c, u, u + c)
    ]
    abs_a = _matmul([[x * vals[j] for j, x in enumerate(row)] for row in r], _inverse3(r))
    jump = [y - x for x, y in zip(left, right, strict=True)]
    return cast(
        "tuple[float, float, float]",
        tuple(
            0.5 * (fl[i] + fr[i]) - 0.5 * sum(abs_a[i][j] * jump[j] for j in range(3))
            for i in range(3)
        ),
    )


class RoeFlowData:
    """Periodic one-dimensional grid; three components padded to four; sign leaves plus a squared-norm and angle binary tree."""

    def __init__(
        self,
        states: Sequence[Sequence[float]],
        *,
        fmt: FixedFormat | None = None,
        gamma: float = 1.4,
        entropy_delta: float = 0.125,
        angle_width: int = 10,
        dx: float = 1.0,
    ) -> None:
        """Initialize the periodic grid and perform the first full update.

        Args:
            states: The conserved variable triple (density, momentum, energy) of each cell; the number of cells must be a power of two of at least 4.
            fmt: Fixed-point format of the conserved variables; defaults to ``FixedFormat(10, 5)``.
            gamma: Ratio of specific heats.
            entropy_delta: Harten entropy fix threshold.
            angle_width: Bit width of the rotation angle words in the residual angle tree.
            dx: Grid step size; must be positive.

        Raises:
            ValidationError: Invalid cell count, component count or dx.
        """
        fmt = fmt or FixedFormat(10, 5)
        self.states: list[tuple[float, float, float]] = [
            cast("tuple[float, float, float]", tuple(map(float, x))) for x in states
        ]
        n = len(states)
        if n < 4 or n & (n - 1) or any(len(x) != 3 for x in states) or dx <= 0:
            raise ValidationError("The flow field requires at least four power-of-two three-component cells and a positive dx")
        self.fmt: FixedFormat = fmt
        self.gamma: float = gamma
        self.delta: float = entropy_delta
        self.angle_width: int = angle_width
        self.dx: float = dx
        self.n: int = n
        self.store: QRAMStore = QRAMStore()
        self.fluxes: dict[int, tuple[float, float, float]] = {}
        self.residuals: dict[int, tuple[float, float, float]] = {}
        self.tree: list[float] = [0.0] * (8 * n)
        self.last_patch: MemoryPatch = self.update(dict(enumerate(self.states)), initialize=True)

    def update(
        self, changed: Mapping[int, Sequence[float]], *, initialize: bool = False
    ) -> MemoryPatch:
        """Perform one local update with the given new cell values and write it into the QRAM store.

        Only the Riemann fluxes on affected faces (including periodic neighbors)
        and the residuals of their adjacent cells are recomputed; the
        conserved-variable, sign, right-hand-side squared-norm tree and angle
        banks are maintained in sync.

        Args:
            changed: Mapping of cell indices to new conserved variable triples.
            initialize: When True, residuals are recomputed for all cells; used only for the initial fill at construction.

        Returns:
            MemoryPatch: The actual writes and recomputation ranges of this call, also saved as ``last_patch``.

        Raises:
            ValidationError: Cell index out of range or component count not equal to three.
        """
        for cell, values in changed.items():
            if not 0 <= cell < self.n or len(values) != 3:
                raise ValidationError("The flow field local update received an invalid address or component")
        changed = {cell: tuple(map(float, values)) for cell, values in changed.items()}

        def state_at(cell: int) -> tuple[float, float, float]:
            """Return the new value of a cell for this call, falling back to the current state when untouched."""
            return cast("tuple[float, float, float]", changed.get(cell, self.states[cell]))

        faces = sorted({face for cell in changed for face in ((cell - 1) % self.n, cell)})
        updates: dict[int, tuple[float, float, float]] = {}
        for face in faces:
            updates[face] = riemann_flux(
                state_at(face),
                state_at((face + 1) % self.n),
                gamma=self.gamma,
                entropy_delta=self.delta,
            )

        def flux_at(face: int) -> tuple[float, float, float]:
            """Return the face flux, preferring the result recomputed in this call."""
            return updates[face] if face in updates else self.fluxes[face]

        cells = sorted({cell for face in faces for cell in (face, (face + 1) % self.n)})
        if initialize:
            cells = list(range(self.n))
        changes: dict[str, dict[int, int]] = {
            bank: {}
            for bank in ("rho", "momentum", "energy", "rhs_values", "rhs_sign", "rhs_angles")
        }
        for cell in changed:
            for j, bank in enumerate(("rho", "momentum", "energy")):
                changes[bank][cell] = self.fmt.encode(changed[cell][j])
        ancestors: set[int] = set()
        for cell in cells:
            residual = cast(
                "tuple[float, float, float]",
                tuple(
                    self.fmt.decode(
                        self.fmt.encode(
                            (flux_at((cell - 1) % self.n)[j] - flux_at(cell)[j]) / self.dx
                        )
                    )
                    for j in range(3)
                ),
            )
            self.residuals[cell] = residual
            for j, value in enumerate((*residual, 0.0)):
                address = 4 * cell + j
                changes["rhs_values"][address] = self.fmt.encode(value)
                changes["rhs_sign"][address] = int(value < 0)
                node = 4 * self.n + address
                self.tree[node] = value * value
                while node > 1:
                    node //= 2
                    ancestors.add(node)
        for node in sorted(ancestors, reverse=True):
            self.tree[node] = self.tree[2 * node] + self.tree[2 * node + 1]
            angle = (
                0
                if self.tree[node] == 0
                else 2 * math.acos(math.sqrt(self.tree[2 * node] / self.tree[node]))
            )
            changes["rhs_angles"][node - 1] = round(
                angle * (1 << self.angle_width) / (2 * math.pi)
            ) % (1 << self.angle_width)
        for cell, values in changed.items():
            self.states[cell] = cast("tuple[float, float, float]", values)
        self.fluxes.update(updates)
        patch = self.store.apply(changes)
        self.last_patch = MemoryPatch(patch.version, patch.changes, tuple(faces), tuple(cells))
        return self.last_patch

    @property
    def rhs_norm(self) -> float:
        """L2 norm of the residual field, i.e. the square root of the squared-sum tree root."""
        return math.sqrt(self.tree[1])
