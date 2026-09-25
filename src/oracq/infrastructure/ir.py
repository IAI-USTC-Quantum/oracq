"RIR v0.1: immutable, serializable register-level operation graph."

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

VERSION = "0.1"
"""Current RIR specification version; the default ``version`` value for newly created ``Program`` objects."""


class ValidationError(ValueError):
    """Input violates a structural or semantic rule of RIR."""


@dataclass(frozen=True)
class RegType:
    """Storage type of a register, consisting of a bit pattern interpretation and a bit width.

    Attributes:
        kind: bit pattern interpretation; one of ``bits``, ``uint``, ``sint`` or ``rational``.
        width: bit width in 0..64; zero width denotes an empty interface.
    """

    kind: str
    width: int


def Bits(width: int) -> RegType:
    """Return a ``RegType`` with the ``bits`` interpretation, denoting a bit string without numeric meaning.

    Args:
        width: bit width in 0..64; zero width denotes an empty interface.

    Returns:
        RegType: storage type with kind ``bits`` and bit width equal to ``width``.
    """
    return RegType("bits", width)


def UInt(width: int) -> RegType:
    """Return a ``RegType`` with the ``uint`` interpretation, denoting an unsigned integer from 0 to 2^width-1.

    Args:
        width: bit width in 0..64; zero width denotes an empty interface.

    Returns:
        RegType: storage type with kind ``uint`` and bit width equal to ``width``.
    """
    return RegType("uint", width)


def SInt(width: int) -> RegType:
    """Return a ``RegType`` with the ``sint`` interpretation, denoting a two's-complement integer; bit operations still act on the raw word.

    Args:
        width: bit width in 0..64; zero width denotes an empty interface.

    Returns:
        RegType: storage type with kind ``sint`` and bit width equal to ``width``.
    """
    return RegType("sint", width)


def Rational(width: int) -> RegType:
    """Return a ``RegType`` with the ``rational`` interpretation, denoting an unsigned word divided by 2^width.

    Args:
        width: bit width in 0..64; also the base-2 exponent of the denominator.

    Returns:
        RegType: storage type with kind ``rational`` and bit width equal to ``width``.
    """
    return RegType("rational", width)


@dataclass(frozen=True)
class Register:
    """Formal quantum register declaration of a module.

    Attributes:
        name: register name; unique within the module and free of conflicts with resource names.
        type: storage type of the register.
    """

    name: str
    type: RegType


@dataclass(frozen=True)
class Span:
    """A contiguous qubit range within a root register.

    Attributes:
        register: root register name.
        start: index of the range start; index 0 is the least significant bit.
        width: bit width of the range.
    """

    register: str
    start: int
    width: int


@dataclass(frozen=True)
class Ref:
    """Register view concatenated from least to most significant bit, carrying no physical qubit numbers."""

    parts: tuple[Span, ...]
    type: RegType

    @property
    def width(self) -> int:
        """Total bit width of the view, i.e. the bit width of ``type``."""
        return self.type.width

    def __getitem__(self, key: int | slice) -> Ref:
        """Take a sub-view by integer index or contiguous slice, from least to most significant bit.

        Args:
            key: integer index or slice with step one; omitted start means zero, omitted stop means the bit width.

        Returns:
            Ref: extracted sub-view, interpreted as ``bits``.

        Raises:
            ValidationError: the index is not an integer, the slice step is not one, or the range is out of bounds.
        """
        if isinstance(key, int):
            start, stop = key, key + 1
        elif isinstance(key, slice) and key.step in (None, 1):
            start = 0 if key.start is None else key.start
            stop = self.width if key.stop is None else key.stop
        else:
            raise ValidationError("views only allow integer indices and contiguous slices")
        if type(start) is not int or type(stop) is not int or not 0 <= start <= stop <= self.width:
            raise ValidationError("register slice out of bounds")
        result = []
        offset = 0
        for part in self.parts:
            lo, hi = max(start, offset), min(stop, offset + part.width)
            if lo < hi:
                result.append(Span(part.register, part.start + lo - offset, hi - lo))
            offset += part.width
        return Ref(tuple(result), Bits(stop - start))

    def reinterpret(self, kind: str) -> Ref:
        """Change the interpretation of the view without altering its bit width or quantum state.

        Args:
            kind: target interpretation; one of ``bits``, ``uint``, ``sint`` or ``rational``.

        Returns:
            Ref: view reusing the original ``parts``, seen under the new interpretation.
        """
        return Ref(self.parts, RegType(kind, self.width))


def fuse(*refs: Ref) -> Ref:
    """Concatenate multiple views from least to most significant bit into a composite view that may cross root registers.

    Args:
        *refs: sequence of views occupying low to high bits in order.

    Returns:
        Ref: new view whose ``parts`` are concatenated in order, interpreted as ``bits``.
    """
    return Ref(tuple(part for ref in refs for part in ref.parts), Bits(sum(r.width for r in refs)))


@dataclass(frozen=True)
class QRAM:
    """Type declaration of a QRAM resource.

    Memory is a mapping from addresses to unsigned data words with unspecified cells
    equal to zero; contents are not written into the program JSON and are bound
    separately as execution input.

    Attributes:
        address_width: bit width of an address word in 1..64.
        data_width: bit width of a data word in 1..64.
    """

    address_width: int
    data_width: int


@dataclass(frozen=True)
class Resource:
    """Formal QRAM parameter declaration of a module.

    Attributes:
        name: resource name; unique within the module and free of conflicts with register names.
        type: ``QRAM`` type of the resource.
    """

    name: str
    type: QRAM


@dataclass(frozen=True)
class Primitive:
    """Register-level primitive instruction.

    A broadcast gate stays a single register-level operation in the IR and is not
    split for containing multiple physical gates.

    Attributes:
        op: gate name, such as ``h``, ``rx``, ``phase``, ``xor``, ``swap``, ``add_const``.
        operands: tuple of operand views.
        angle: rotation or phase angle in radians; ``None`` for gates that do not use it.
        value: modular addend of ``add_const``; ``None`` for gates that do not use it.
    """

    op: str
    operands: tuple[Ref, ...]
    angle: float | None = None
    value: int | None = None


@dataclass(frozen=True)
class Load:
    """QRAM query instruction: ``data`` is XORed with the content of the addressed cell.

    The instruction is self-inverse, does not require the target to start at zero,
    and does not modify the address or the classical memory.

    Attributes:
        resource: formal QRAM resource name being queried.
        address: address view; its width must match the resource declaration.
        data: data view; its width must match the resource declaration and it must not overlap the address.
    """

    resource: str
    address: Ref
    data: Ref


@dataclass(frozen=True)
class Store:
    """QRAM random-write instruction: at execution time assigns the classical cell to the data view.

    The address and data views must be in a definite basis state; the instruction
    changes no qubits, may only appear in a module body or a ``Repeat`` body, and
    is forbidden inside ``Control`` and ``Adjoint`` bodies.

    Attributes:
        resource: formal QRAM resource name being written.
        address: address view; its width must match the resource declaration.
        data: data view; its width must match the resource declaration and it must not overlap the address.
    """

    resource: str
    address: Ref
    data: Ref


@dataclass(frozen=True)
class Call:
    """Module call instruction; the callee definition is shared in the ``Program`` and not inlined at call sites.

    Attributes:
        module: callee module name.
        arguments: quantum argument views bound in callee signature order.
        resources: resource argument names bound in callee signature order.
    """

    module: str
    arguments: tuple[Ref, ...]
    resources: tuple[str, ...] = ()


@dataclass(frozen=True)
class Repeat:
    """Static repetition instruction: applies the body instructions in order ``count`` times.

    The generation and serialization stages do not copy the instruction body by ``count``.

    Attributes:
        count: repeat count in 0..2^63-1; zero repetitions is the identity operation.
        body: tuple of instructions being repeated; it must be structurally valid even when ``count`` is zero.
    """

    count: int
    body: tuple[Instruction, ...]


@dataclass(frozen=True)
class Control:
    """Coherent control instruction: applies the body instructions when the control view equals the given value, otherwise is the identity.

    The control register participates in coherent control as a quantum condition
    and is not measured; the control bits are protected throughout the ``body``.

    Attributes:
        register: nonempty control view.
        value: unsigned comparison value for the control view.
        body: tuple of instructions executed under control.
    """

    register: Ref
    value: int
    body: tuple[Instruction, ...]


@dataclass(frozen=True)
class Adjoint:
    """Adjoint instruction: reverses the order of the body instructions and takes the adjoint of each operation.

    Creation does not rewrite or expand its body; a double adjoint recovers the
    original operation.

    Attributes:
        body: tuple of instructions being adjointed.
    """

    body: tuple[Instruction, ...]


Instruction: TypeAlias = Primitive | Load | Store | Call | Repeat | Control | Adjoint
"""Union type alias of RIR instruction nodes."""


@dataclass(frozen=True)
class Module:
    """Module definition consisting of a call signature, QRAM resources and an ordered instruction body.

    Attributes:
        name: module name; unique within the ``Program``.
        registers: ordered list of public quantum parameters.
        resources: ordered list of formal QRAM parameters.
        body: ordered instruction body; ``None`` denotes an open declaration, i.e. an oracle slot whose implementation is not yet bound.
        attributes: ordered attribute pairs with values limited to strings, integers, finite floats or booleans.
        locals: module-private work registers with zero input and zero output, not part of the public call signature.
    """

    name: str
    registers: tuple[Register, ...]
    resources: tuple[Resource, ...]
    body: tuple[Instruction, ...] | None
    attributes: tuple[tuple[str, str | int | float | bool], ...] = ()
    locals: tuple[Register, ...] = ()


@dataclass(frozen=True)
class Program:
    """Complete program: the module table plus the entry module name.

    Attributes:
        entry: entry module name; must refer to a module defined in ``modules``.
        modules: all module definitions; definitions are shared across call sites without copying bodies.
        version: RIR version string; the reader additionally accepts historic versions ``0.1`` and ``0.2``.
    """

    entry: str
    modules: tuple[Module, ...]
    version: str = VERSION

    @property
    def module_map(self) -> dict[str, Module]:
        """Index of module definitions keyed by module name."""
        return {module.name: module for module in self.modules}

    @property
    def main(self) -> Module:
        """Definition of the entry module."""
        return self.module_map[self.entry]
