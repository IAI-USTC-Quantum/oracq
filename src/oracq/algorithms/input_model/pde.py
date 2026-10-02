"""Composable spatial and initial-state inputs for quantum PDE methods."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Protocol

from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding, projector
from oracq.algorithms.input_model.initial import (
    ArrayInput as ArrayInput,
)
from oracq.algorithms.input_model.initial import (
    InitialInput as InitialInput,
)
from oracq.algorithms.input_model.initial import (
    OracleInput as OracleInput,
)
from oracq.algorithms.input_model.initial import (
    PreparedInitial as PreparedInitial,
)
from oracq.algorithms.input_model.initial import (
    QRAMInput as QRAMInput,
)
from oracq.algorithms.input_model.initial import (
    QRAMInputConfig as QRAMInputConfig,
)
from oracq.algorithms.input_model.initial import (
    UniformInput as UniformInput,
)
from oracq.algorithms.input_model.operators import BlockEncoding, _name, identity, product
from oracq.algorithms.input_model.oracles import (
    annotate,
)
from oracq.applications.qham.pde import Field, PolynomialPDE
from oracq.applications.qham.reference import Grid, SpatialGrid
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, ValidationError, fuse
from oracq.infrastructure.qram_schema import RegisteredQRAM

Derivative = tuple[tuple[str, int], ...]


class PDEGrid(SpatialGrid, Protocol):
    """Spatial rows, quantum derivative encodings, and their memory snapshots."""

    @property
    def qrams(self) -> tuple[RegisteredQRAM, ...]:
        """QRAM snapshots backing the grid's derivative encodings."""

    def derivative_encoding(self, derivative: Derivative) -> BlockEncoding:
        """Return a block encoding of the sparse derivative operator."""


@dataclass(frozen=True)
class UniformGrid1D:
    """Uniform spatial samples with an explicit interval and boundary convention.

    Periodic samples cover ``[start, stop)`` with spacing ``(stop-start)/points``.
    Zero Dirichlet samples are interior points with spacing
    ``(stop-start)/(points+1)``. Power-of-two periodic grids use shift circuits;
    other grids use a bounded spatial Pauli expansion of at most five qubits.
    Only spatial derivatives are materialized, never multilinear PDE ports.
    """

    start: float
    stop: float
    points: int
    boundary: str = "periodic"
    axis: str = "x"
    _grid: Grid = field(init=False, repr=False, compare=False)
    _encodings: dict[Derivative, BlockEncoding] = field(default_factory=dict, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if (
            type(self.points) is not int or self.points < 1
            or not math.isfinite(self.start) or not math.isfinite(self.stop)
            or self.stop <= self.start or not self.axis
        ):
            raise ValidationError("a uniform grid requires a finite increasing interval and a positive point count")
        divisor = self.points if self.boundary == "periodic" else self.points + 1
        object.__setattr__(self, "_grid", Grid(
            (self.axis,), (self.points,), ((self.stop - self.start) / divisor,), self.boundary
        ))

    @property
    def axes(self) -> tuple[str, ...]:
        """Spatial axis names."""
        return self._grid.axes

    @property
    def size(self) -> int:
        """Number of physical samples, excluding register padding."""
        return self.points

    @property
    def spatial_width(self) -> int:
        """Address register width, retaining padding for arbitrary point counts."""
        return self._grid.spatial_width

    @property
    def spacing(self) -> float:
        """Physical sample spacing."""
        return self._grid.spacing[0]

    @property
    def qrams(self) -> tuple[RegisteredQRAM, ...]:
        """Uniform grids do not require a memory bank."""
        return ()

    def coordinate(self, index: int) -> float:
        """Physical coordinate of a sample."""
        if type(index) is not int or not 0 <= index < self.points:
            raise ValidationError("grid sample index is out of range")
        return self.start + (index + (self.boundary != "periodic")) * self.spacing

    def derivative_row(self, derivative: Derivative, row: int) -> tuple[tuple[int, float], ...]:
        """Finite-difference row, with padding outside the physical grid zero."""
        return self._grid.derivative_row(derivative, row)

    def derivative_encoding(self, derivative: Derivative) -> BlockEncoding:
        """Encode a derivative while preserving the physical grid size."""
        if derivative not in self._encodings:
            self._encodings[derivative] = self._derivative_encoding(derivative)
        return self._encodings[derivative]

    def _derivative_encoding(self, derivative: Derivative) -> BlockEncoding:
        from oracq.applications.qham.stencils import derivative_encoding

        if any(axis != self.axis or type(order) is not int or order < 0 for axis, order in derivative):
            raise ValidationError("the derivative must use this grid's axis and nonnegative integer orders")
        if not derivative:
            return identity(self.spatial_width) if self.points == 1 << self.spatial_width else projector(self.spatial_width, range(self.points))
        if self.boundary == "periodic" and not self.points & (self.points - 1):
            return derivative_encoding(self._grid, derivative)
        if self.spatial_width > 5:
            raise ValidationError("this uniform grid exceeds the five-qubit spatial expansion budget; supply a custom grid encoding")
        storage = 1 << self.spatial_width
        rows = [dict(self.derivative_row(derivative, row)) for row in range(storage)]
        return matrix_pauli_encoding([[rows[row].get(col, 0.0) for col in range(storage)] for row in range(storage)])


@dataclass(frozen=True)
class UnstructuredGridConfig:
    """Layout of real discrete derivative matrices in a QRAM angle bank.

    Each derivative key maps to a base address. A matrix uses
    ``base + row * 2**spatial_width + column`` and has entries
    ``value_scale*cos(pi*word/2**word_length)``. Thus a zero entry is an
    angle word ``2**(word_length-1)``, not a zero word. Boundary conditions
    and geometric discretization are already incorporated in these matrices.
    """

    size: int
    derivative_offsets: Mapping[Derivative, int]
    value_scale: float = 1.0
    axes: tuple[str, ...] = ("x",)

    def __post_init__(self) -> None:
        if type(self.size) is not int or self.size < 1 or not math.isfinite(self.value_scale) or self.value_scale <= 0:
            raise ValidationError("unstructured grid size and value scale must be positive")
        if not self.axes or len(set(self.axes)) != len(self.axes) or any(not axis for axis in self.axes):
            raise ValidationError("unstructured grid axes must be nonempty and unique")
        offsets = dict(self.derivative_offsets)
        for derivative, offset in offsets.items():
            if not derivative or any(axis not in self.axes or type(order) is not int or order < 1 for axis, order in derivative):
                raise ValidationError("stored derivative keys require known axes and positive integer orders")
            if type(offset) is not int or offset < 0:
                raise ValidationError("derivative table offsets must be nonnegative integers")
        object.__setattr__(self, "derivative_offsets", MappingProxyType(offsets))

    def encode_matrices(
        self, matrices: Mapping[Derivative, Sequence[Sequence[float]]], word_length: int,
        *, max_words: int = 1_048_576,
    ) -> dict[int, int]:
        """Encode physical derivative matrices in this layout, including zero padding.

        The amplitude error per entry is at most
        ``value_scale*pi/2**word_length``. The negative endpoint is saturated
        to the largest word instead of wrapping a full turn to zero.
        """
        if type(word_length) is not int or not 1 <= word_length <= 64:
            raise ValidationError("derivative angle word length must be an integer in 1..64")
        if set(matrices) != set(self.derivative_offsets):
            raise ValidationError("derivative matrices must match the configured table keys")
        storage = 1 << max(1, (self.size - 1).bit_length())
        if storage * storage * len(matrices) > max_words:
            raise ValidationError("derivative matrix encoding exceeds max_words; supply preencoded QRAM data")
        limit = 1 << word_length
        cells = {}
        for derivative, matrix in matrices.items():
            if len(matrix) != self.size or any(len(row) != self.size for row in matrix):
                raise ValidationError("derivative matrices must be square over physical grid nodes")
            offset = self.derivative_offsets[derivative]
            for row in range(storage):
                for col in range(storage):
                    value = matrix[row][col] if row < self.size and col < self.size else 0.0
                    if not math.isfinite(value) or abs(value) > self.value_scale:
                        raise ValidationError("derivative matrix entries must be finite and within value_scale")
                    cells[offset + row * storage + col] = min(limit - 1, round(math.acos(value / self.value_scale) * limit / math.pi))
        return cells


@dataclass(frozen=True)
class UnstructuredGrid:
    """A spatial discretization supplied as QRAM derivative matrices.

    The quantum adapter uses uniform row/column projection and QRAM angle
    transduction. Its normalization is ``2**spatial_width * value_scale``;
    it supports nonsymmetric real matrices and makes no sparse-access speedup
    claim. A different storage model can implement ``PDEGrid`` directly.
    """

    qram: RegisteredQRAM
    config: UnstructuredGridConfig
    _encodings: dict[Derivative, BlockEncoding] = field(default_factory=dict, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        length = 1 << (2 * self.spatial_width)
        intervals = sorted((offset, offset + length) for offset in self.config.derivative_offsets.values())
        if any(stop > 1 << self.qram.address_length for _, stop in intervals):
            raise ValidationError("a derivative table does not fit the registered QRAM address width")
        if any(left[1] > right[0] for left, right in zip(intervals, intervals[1:], strict=False)):
            raise ValidationError("stored derivative matrix tables must not overlap")

    @property
    def axes(self) -> tuple[str, ...]:
        """Spatial axis names."""
        return self.config.axes

    @property
    def size(self) -> int:
        """Number of physical nodes."""
        return self.config.size

    @property
    def spatial_width(self) -> int:
        """Node address width."""
        return max(1, (self.size - 1).bit_length())

    @property
    def qrams(self) -> tuple[RegisteredQRAM, ...]:
        """The registered operator memory."""
        return (self.qram,)

    def _offset(self, derivative: Derivative) -> int:
        if derivative not in self.config.derivative_offsets:
            raise ValidationError("the unstructured grid has no stored matrix for derivative " + repr(derivative))
        return self.config.derivative_offsets[derivative]

    def derivative_row(self, derivative: Derivative, row: int) -> tuple[tuple[int, float], ...]:
        """Decode a physical row for reference evaluation and known coefficients."""
        if not 0 <= row < self.size:
            return ()
        if not derivative:
            return ((row, 1.0),)
        offset, storage, cells = self._offset(derivative), 1 << self.spatial_width, self.qram.snapshot()
        return tuple((column, self.config.value_scale * math.cos(
            math.pi * cells.get(offset + row * storage + column, 0) / (1 << self.qram.word_length)
        )) for column in range(self.size))

    def derivative_encoding(self, derivative: Derivative) -> BlockEncoding:
        """Query the configured table coherently and project both padding spaces."""
        if derivative not in self._encodings:
            self._encodings[derivative] = self._derivative_encoding(derivative)
        return self._encodings[derivative]

    def _derivative_encoding(self, derivative: Derivative) -> BlockEncoding:
        n = self.spatial_width
        mask = identity(n) if self.size == 1 << n else projector(n, range(self.size))
        if not derivative:
            return mask
        offset = self._offset(derivative)
        address_width, word_width = self.qram.address_length, self.qram.word_length
        b = Builder(
            _name("qram_grid_derivative", self.qram.name, self.config, derivative),
            {"target": Bits(n), "signal": Bits(n + address_width + word_width + 1)},
            {self.qram.name: self.qram.type},
        )
        row = b["signal"][:n]
        address = b["signal"][n:n + address_width]
        angle = b["signal"][n + address_width:n + address_width + word_width]
        flag = b["signal"][n + address_width + word_width]
        b.h(row)
        b.xor(fuse(b["target"], row), address[:2 * n])
        b.add_const(address.reinterpret("uint"), offset)
        b.qram(self.qram.name, address, angle)
        for bit in range(word_width):
            with b.control(angle[bit]):
                b.ry(flag, 2 * math.pi * (1 << bit) / (1 << word_width))
        b.qram(self.qram.name, address, angle)
        b.add_const(address.reinterpret("uint"), (-offset) % (1 << address_width))
        b.xor(fuse(b["target"], row), address[:2 * n])
        b.swap(b["target"], row)
        b.h(row)
        encoded = BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=(1 << n) * self.config.value_scale))
        return encoded if self.size == 1 << n else product(mask, product(encoded, mask))


@dataclass(frozen=True)
class InitialLayout:
    """Physical fields and grid nodes, with field-major register padding."""

    grid: SpatialGrid
    fields: tuple[str, ...]
    max_words: int = 4096

    @property
    def width(self) -> int:
        """Combined component and spatial register width."""
        return self.grid.spatial_width + (len(self.fields) - 1).bit_length()

    @property
    def size(self) -> int:
        """Number of physical samples across all fields."""
        return self.grid.size * len(self.fields)

    def encode(self, values: Sequence[complex] | Mapping[str, Sequence[complex]]) -> tuple[complex, ...]:
        """Pad physical samples without treating padded addresses as grid nodes."""
        if 1 << self.width > self.max_words:
            raise ValidationError("initial sample materialization exceeds max_initial_words; provide an oracle input")
        if isinstance(values, Mapping):
            if set(values) != set(self.fields):
                raise ValidationError("initial sample fields must match the PDE fields")
            components = [tuple(map(complex, values[name])) for name in self.fields]
        else:
            raw = tuple(map(complex, values))
            if len(raw) != self.grid.size * len(self.fields):
                raise ValidationError("initial samples must cover every physical node and field")
            components = [raw[i * self.grid.size:(i + 1) * self.grid.size] for i in range(len(self.fields))]
        if any(len(values) != self.grid.size for values in components):
            raise ValidationError("each initial field must cover the physical grid")
        if any(not math.isfinite(v.real) or not math.isfinite(v.imag) for values in components for v in values):
            raise ValidationError("initial samples must be finite")
        encoded = [0j] * (1 << self.width)
        for i, values in enumerate(components):
            start = i * (1 << self.grid.spatial_width)
            encoded[start:start + self.grid.size] = values
        return tuple(encoded)


@dataclass(frozen=True)
class PDEProblem:
    """A PDE recipe whose grid and physical initial input can be composed independently.

    ``type`` accepts ``burgers``, ``heat``, or an existing ``PolynomialPDE``.
    Presets use one field ``u`` and the grid's first axis. ``known`` holds
    physical coefficient/forcing samples. Solver approximations belong to the
    method configuration, not this problem object.
    """

    grid: PDEGrid
    initial: InitialInput = field(default_factory=UniformInput)
    type: str | PolynomialPDE = "burgers"
    final_time: float = 1.0
    viscosity: float = 0.1
    known: Mapping[str, Sequence[complex]] = field(default_factory=dict)
    equation: PolynomialPDE = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.final_time) or self.final_time < 0:
            raise ValidationError("PDE final_time must be finite and nonnegative")
        if not math.isfinite(self.viscosity) or self.viscosity < 0:
            raise ValidationError("PDE viscosity must be finite and nonnegative")
        if isinstance(self.type, PolynomialPDE):
            pde = self.type.validate()
        else:
            if self.type not in {"burgers", "heat"} or len(self.grid.axes) != 1:
                raise ValidationError("PDE presets support one-dimensional burgers or heat; supply PolynomialPDE for other equations")
            u = Field("u")
            rhs = self.viscosity * u.d(self.grid.axes[0], 2)
            if self.type == "burgers":
                rhs = rhs - u * u.d(self.grid.axes[0])
            pde = PolynomialPDE.from_equations({"u": rhs}, axes=self.grid.axes, label=self.type)
        if pde.axes != self.grid.axes:
            raise ValidationError("PDE axes must match the spatial grid")
        object.__setattr__(self, "equation", pde)
        object.__setattr__(self, "known", MappingProxyType({k: tuple(v) for k, v in self.known.items()}))
