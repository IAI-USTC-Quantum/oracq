"""Composable generator and physical vector inputs for autonomous linear ODEs."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Protocol

from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding, projector
from oracq.algorithms.input_model.contracts import finite_real, positive_integer, require_instance
from oracq.algorithms.input_model.initial import InitialInput, UniformInput
from oracq.algorithms.input_model.interfaces import BlockEncodingProtocol, as_block_encoding
from oracq.algorithms.input_model.operators import BlockEncoding, _name, product
from oracq.algorithms.input_model.oracles import annotate
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, ValidationError, fuse
from oracq.infrastructure.qram_schema import RegisteredQRAM, register_qram


@dataclass(frozen=True)
class ODELayout:
    """Physical components at consecutive addresses, with optional names and zero padding."""

    size: int
    components: tuple[str, ...] = ()
    max_words: int = 4096

    def __post_init__(self) -> None:
        positive_integer(self.size, "ODE.size")
        positive_integer(self.max_words, "ODE.max_initial_words")
        components = tuple(self.components)
        if components and (len(components) != self.size or any(not isinstance(c, str) or not c for c in components) or len(set(components)) != self.size):
            raise ValidationError("ODE component names must be nonempty, unique, and match the vector size")
        object.__setattr__(self, "components", components)

    @property
    def width(self) -> int:
        """The target register width, including the scalar case."""
        return max(1, (self.size - 1).bit_length())

    def encode(self, values: Sequence[complex] | Mapping[str, Sequence[complex]]) -> tuple[complex, ...]:
        """Encode physical components using the same initial inputs as PDEs."""
        if 1 << self.width > self.max_words:
            raise ValidationError("initial vector materialization exceeds max_initial_words; provide an oracle input")
        if isinstance(values, Mapping):
            if not self.components or set(values) != set(self.components) or any(len(v) != 1 for v in values.values()):
                raise ValidationError("ODE initial mappings require one sample per named component")
            raw = tuple(complex(values[c][0]) for c in self.components)
        else:
            raw = tuple(map(complex, values))
        if len(raw) != self.size or any(not math.isfinite(v.real) or not math.isfinite(v.imag) for v in raw):
            raise ValidationError("ODE initial values must be finite and cover every physical component")
        return raw + (0j,) * ((1 << self.width) - self.size)


@dataclass(frozen=True)
class PreparedGenerator:
    """A generator block encoding and its immutable execution memory banks."""

    encoding: BlockEncoding
    qrams: tuple[RegisteredQRAM, ...] = ()


class GeneratorInput(Protocol):
    """An operator source resolved against the physical ODE vector layout."""

    @property
    def size(self) -> int:
        """Declared physical dimension of the generator's target vector."""

    def prepare_generator(self, layout: ODELayout) -> PreparedGenerator:
        """Resolve this source against the ODE layout into an encoding and memories."""


@dataclass(frozen=True)
class OracleGeneratorInput:
    """An existing block-encoding provider with a declared physical dimension."""

    encoding: BlockEncodingProtocol
    size: int
    qrams: tuple[RegisteredQRAM, ...] = ()

    def __post_init__(self) -> None:
        positive_integer(self.size, "ODE.generator.size")
        object.__setattr__(self, "qrams", tuple(self.qrams))

    def prepare_generator(self, layout: ODELayout) -> PreparedGenerator:
        """Resolve the provider; method contracts check layout and capabilities."""
        return PreparedGenerator(as_block_encoding(self.encoding), self.qrams)


def _matrix(values: Sequence[Sequence[complex]], size: int) -> tuple[tuple[complex, ...], ...]:
    result = tuple(tuple(map(complex, row)) for row in values)
    if len(result) != size or any(len(row) != size for row in result):
        raise ValidationError("ODE generator matrix must be square and match its physical dimension")
    if any(not math.isfinite(v.real) or not math.isfinite(v.imag) for row in result for v in row):
        raise ValidationError("ODE generator entries must be finite")
    return result


@dataclass(frozen=True)
class QRAMMatrixConfig:
    """Real matrix angles at ``offset + row*2**width + column``.

    A word represents ``value_scale*cos(pi*word/2**word_length)``. Absent
    words are zero words, representing positive ``value_scale`` rather than
    zero matrix entries. ``encode_matrix`` explicitly fills all padding.
    """

    size: int
    value_scale: float = 1.0
    offset: int = 0

    def __post_init__(self) -> None:
        positive_integer(self.size, "ODE.qram_matrix.size")
        finite_real(self.value_scale, "ODE.qram_matrix.value_scale", minimum=0, strict=True)
        positive_integer(self.offset, "ODE.qram_matrix.offset", minimum=0)

    def encode_matrix(self, values: Sequence[Sequence[complex]], word_length: int) -> dict[int, int]:
        """Quantize a real matrix, saturating the negative endpoint instead of wrapping."""
        positive_integer(word_length, "ODE.qram_matrix.word_length", maximum=64)
        matrix = _matrix(values, self.size)
        if any(v.imag or abs(v.real) > self.value_scale for row in matrix for v in row):
            raise ValidationError("QRAM matrix entries must be real and bounded by value_scale")
        storage, turns = 1 << max(1, (self.size - 1).bit_length()), 1 << word_length
        return {self.offset + row * storage + column: min(turns - 1, round(
            math.acos(matrix[row][column].real / self.value_scale if row < self.size and column < self.size else 0) * turns / math.pi
        )) for row in range(storage) for column in range(storage)}


@dataclass(frozen=True)
class QRAMMatrixInput:
    """Coherent real generator matrix access through a registered angle bank."""

    qram: RegisteredQRAM
    config: QRAMMatrixConfig

    def __post_init__(self) -> None:
        require_instance(self.qram, RegisteredQRAM, "ODE.qram_matrix.bank")
        require_instance(self.config, QRAMMatrixConfig, "ODE.qram_matrix.config")
        width = max(1, (self.size - 1).bit_length())
        if self.config.offset + (1 << (2 * width)) > 1 << self.qram.address_length:
            raise ValidationError("the generator matrix does not fit its QRAM address width")

    @property
    def size(self) -> int:
        """Physical generator dimension."""
        return self.config.size

    def prepare_generator(self, layout: ODELayout) -> PreparedGenerator:
        """Query and uncompute matrix angles, without materializing a gate expansion."""
        if layout.size != self.size:
            raise ValidationError("QRAM generator dimension must match the ODE layout")
        n, a, w = layout.width, self.qram.address_length, self.qram.word_length
        b = Builder(_name("qram_ode_generator", self.qram.name, self.config),
                    {"target": Bits(n), "signal": Bits(n + a + w + 1)},
                    {self.qram.name: self.qram.type})
        row = b["signal"][:n]
        address, angle = b["signal"][n:n + a], b["signal"][n + a:n + a + w]
        flag = b["signal"][n + a + w]
        b.h(row)
        b.xor(fuse(b["target"], row), address[:2 * n])
        b.add_const(address.reinterpret("uint"), self.config.offset)
        b.qram(self.qram.name, address, angle)
        for bit in range(w):
            with b.control(angle[bit]):
                b.ry(flag, 2 * math.pi * (1 << bit) / (1 << w))
        b.qram(self.qram.name, address, angle)
        b.add_const(address.reinterpret("uint"), (-self.config.offset) % (1 << a))
        b.xor(fuse(b["target"], row), address[:2 * n])
        b.swap(b["target"], row)
        b.h(row)
        encoding = BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=(1 << n) * self.config.value_scale))
        return PreparedGenerator(encoding, (self.qram,))


@dataclass(frozen=True)
class MatrixInput:
    """Finite generator matrix, encoded with small Pauli gates or real QRAM angles.

    Gate expansion supports at most five target bits. QRAM preparation still
    stores the supplied dense matrix; it does not imply a sparse-access speedup.
    """

    values: Sequence[Sequence[complex]]
    encoding: str = "gates"
    angle_width: int = 16
    value_scale: float = 1.0

    def __post_init__(self) -> None:
        positive_integer(len(self.values), "ODE.matrix.size")
        object.__setattr__(self, "values", _matrix(self.values, len(self.values)))
        if self.encoding not in {"gates", "qram"}:
            raise ValidationError("ODE matrix encoding must be gates or qram")
        positive_integer(self.angle_width, "ODE.matrix.angle_width", maximum=64)
        finite_real(self.value_scale, "ODE.matrix.value_scale", minimum=0, strict=True)

    @property
    def size(self) -> int:
        """Physical matrix dimension before register padding."""
        return len(self.values)

    def prepare_generator(self, layout: ODELayout) -> PreparedGenerator:
        """Encode G directly, retaining its sign in the convention u'=Gu."""
        if layout.size != self.size:
            raise ValidationError("generator matrix dimension must match the ODE layout")
        if self.encoding == "qram":
            config = QRAMMatrixConfig(self.size, self.value_scale)
            bank = register_qram(config.encode_matrix(self.values, self.angle_width), 2 * layout.width, self.angle_width)
            return QRAMMatrixInput(bank, config).prepare_generator(layout)
        if layout.width > 5:
            raise ValidationError("gate matrix expansion supports at most five target bits; provide a generator oracle")
        storage = 1 << layout.width
        matrix = [[self.values[r][c] if r < self.size and c < self.size else 0j for c in range(storage)] for r in range(storage)]
        return PreparedGenerator(matrix_pauli_encoding(matrix, drop_tolerance=0))


def physical_generator(encoding: BlockEncoding, layout: ODELayout) -> BlockEncoding:
    """Project padding on both sides, preventing paths through unused components."""
    if encoding.width != layout.width:
        raise ValidationError("ODE generator width must match the physical vector layout")
    if layout.size == 1 << layout.width:
        return encoding
    mask = projector(layout.width, range(layout.size))
    return product(mask, product(encoding, mask))


@dataclass(frozen=True)
class ODEProblem:
    """Autonomous homogeneous system u'=Gu with replaceable generator and initial inputs.

    Dissipativity is an explicit caller declaration, checked by methods that
    require it. Source terms, time-dependent generators, and nonlinear systems
    require separate adaptations and are not implicitly frozen or discarded.
    """

    generator: GeneratorInput
    initial: InitialInput = field(default_factory=UniformInput)
    final_time: float = 1.0
    components: tuple[str, ...] = ()
    dissipative: bool | None = None
    evidence: str = "caller_declared_unverified"

    def __post_init__(self) -> None:
        if not callable(getattr(self.generator, "prepare_generator", None)):
            raise ValidationError("ODE generator must implement prepare_generator(layout)")
        if not callable(getattr(self.initial, "prepare_initial", None)):
            raise ValidationError("ODE initial input must implement prepare_initial(layout)")
        finite_real(self.final_time, "ODE.final_time", minimum=0)
        if self.dissipative is not None and type(self.dissipative) is not bool:
            raise ValidationError("ODE dissipative must be bool or None")
        if not isinstance(self.evidence, str) or not self.evidence:
            raise ValidationError("ODE evidence requires a nonempty description")
        object.__setattr__(self, "components", ODELayout(self.generator.size, self.components).components)

    @property
    def size(self) -> int:
        """Physical vector dimension shared by generator, initial state, and readout."""
        return self.generator.size
