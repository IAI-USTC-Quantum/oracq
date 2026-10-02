"""Shared physical initial-state inputs for ODE and PDE methods."""

from __future__ import annotations

import cmath
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import Protocol

from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    annotate,
    basis_state,
    gate_state_prep,
    qram_state_angles,
)
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, ValidationError
from oracq.infrastructure.qram_schema import RegisteredQRAM, register_qram


class InitialStateLayout(Protocol):
    """Physical vector size, register width, and bounded sample encoding."""

    @property
    def size(self) -> int:
        """Number of physical samples in the layout's vector."""

    @property
    def width(self) -> int:
        """Address width in qubits indexing the physical samples."""

    @property
    def max_words(self) -> int:
        """Largest number of QRAM words the layout may use for padded samples."""

    def encode(self, values: Sequence[complex] | Mapping[str, Sequence[complex]]) -> tuple[complex, ...]:
        """Encode flat or per-component sample values into the physical vector."""


@dataclass(frozen=True)
class PreparedInitial:
    """Reversible amplitude preparation, physical norm, and memory snapshots."""

    preparation: StatePreparation
    norm: float
    qrams: tuple[RegisteredQRAM, ...] = ()

    def __post_init__(self) -> None:
        if not math.isfinite(self.norm) or self.norm < 0:
            raise ValidationError("initial physical norm must be finite and nonnegative")


class InitialInput(Protocol):
    """An initial-state source resolved against a method's physical layout."""

    def prepare_initial(self, layout: InitialStateLayout) -> PreparedInitial:
        """Resolve this source against the physical layout into a preparation."""


@dataclass(frozen=True)
class UniformInput:
    """A constant physical vector, independently of its layout implementation."""

    value: complex = 1.0

    def prepare_initial(self, layout: InitialStateLayout) -> PreparedInitial:
        """Prepare uniform amplitudes only over the layout's physical samples."""
        value = complex(self.value)
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            raise ValidationError("the uniform initial value must be finite")
        norm = abs(value) * math.sqrt(layout.size)
        if norm == 0:
            return PreparedInitial(basis_state(layout.width), norm)
        if layout.size == 1 << layout.width:
            b = Builder(_name("uniform_initial", layout.width, value), {"target": Bits(layout.width), "work": Bits(0)})
            b.h(b["target"])
            b.global_phase(cmath.phase(value))
            prep = StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True, clean_work=True))
        else:
            if 1 << layout.width > layout.max_words:
                raise ValidationError("uniform initial padding exceeds max_initial_words; provide an oracle input")
            prep = gate_state_prep(layout.encode([value] * layout.size))
        return PreparedInitial(prep, norm)


@dataclass(frozen=True)
class ArrayInput:
    """Physical initial samples, using gates or nonnegative QRAM angle preparation."""

    values: Sequence[complex] | Mapping[str, Sequence[complex]]
    encoding: str = "gates"
    angle_width: int = 16

    def __post_init__(self) -> None:
        if self.encoding not in {"gates", "qram"}:
            raise ValidationError("initial array encoding must be gates or qram")
        if type(self.angle_width) is not int or not 1 <= self.angle_width <= 64:
            raise ValidationError("initial angle width must be an integer in 1..64")
        values = MappingProxyType({k: tuple(v) for k, v in self.values.items()}) if isinstance(self.values, Mapping) else tuple(self.values)
        object.__setattr__(self, "values", values)

    def prepare_initial(self, layout: InitialStateLayout) -> PreparedInitial:
        """Convert physical samples into a normalized preparation and retain their norm."""
        values = layout.encode(self.values)
        norm = math.sqrt(sum(abs(v) ** 2 for v in values))
        amplitudes = values if norm else (1.0,) + (0.0,) * (len(values) - 1)
        if self.encoding == "gates":
            return PreparedInitial(gate_state_prep(amplitudes), norm)
        bank = register_qram(qram_state_angles(amplitudes, self.angle_width), layout.width, self.angle_width)
        return QRAMInput(bank, QRAMInputConfig(layout.width, norm)).prepare_initial(layout)


@dataclass(frozen=True)
class OracleInput:
    """An existing preparation oracle with an explicit physical vector norm.

    Signed/complex preparations and application-specific loaders can use this
    interface. Memory snapshots must use the oracle's leaf resource names.
    """

    preparation: StatePreparation
    norm: float
    qrams: tuple[RegisteredQRAM, ...] = ()

    def prepare_initial(self, layout: InitialStateLayout) -> PreparedInitial:
        """Return the supplied oracle; the consuming method checks capabilities and width."""
        return PreparedInitial(self.preparation, self.norm, self.qrams)


@dataclass(frozen=True)
class QRAMInputConfig:
    """Rotation-tree layout: node ``2**depth-1+prefix`` plus a base offset.

    Words encode Ry angles in units of ``2*pi/2**word_length``. The physical
    norm is supplied separately; this is amplitude access, not raw value access.
    """

    state_width: int
    norm: float
    offset: int = 0


@dataclass(frozen=True)
class QRAMInput:
    """Initial amplitude preparation from a registered rotation-tree bank."""

    qram: RegisteredQRAM
    config: QRAMInputConfig

    def prepare_initial(self, layout: InitialStateLayout) -> PreparedInitial:
        """Load and uncompute each angle while retaining resource names and widths."""
        width, offset = self.config.state_width, self.config.offset
        if type(width) is not int or width != layout.width:
            raise ValidationError("QRAM initial state width must match the physical state layout")
        if type(offset) is not int or offset < 0 or offset + (1 << width) - 1 > 1 << self.qram.address_length:
            raise ValidationError("the initial rotation tree does not fit its QRAM address width")
        a, w = self.qram.address_length, self.qram.word_length
        b = Builder(_name("qram_initial", self.qram.name, self.config), {"target": Bits(width), "work": Bits(a + w)}, {self.qram.name: self.qram.type})
        address, angle = b["work"][:a], b["work"][a:]
        for depth in range(width):
            bit = width - depth - 1
            if depth:
                b.xor(b["target"][bit + 1:], address[:depth])
            base = offset + (1 << depth) - 1
            b.add_const(address.reinterpret("uint"), base)
            b.qram(self.qram.name, address, angle)
            for k in range(w):
                with b.control(angle[k]):
                    b.ry(b["target"][bit], 2 * math.pi * (1 << k) / (1 << w))
            b.qram(self.qram.name, address, angle)
            b.add_const(address.reinterpret("uint"), (-base) % (1 << a))
            if depth:
                b.xor(b["target"][bit + 1:], address[:depth])
        prep = StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True, clean_work=True, implementation="qram_rotation_tree"))
        return PreparedInitial(prep, self.config.norm, (self.qram,))
