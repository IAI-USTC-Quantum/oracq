"Pointer-style reads and writes of QRAM data structures: C-style base addresses, offsets, multidimensional views and random writes."

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import QRAM, Ref, UInt, ValidationError


@dataclass(frozen=True)
class _Addr:
    """Address expression: constant + Σ(register view × positive coefficient), taken modulo 2^address width."""

    constant: int = 0
    terms: tuple[tuple[Ref, int], ...] = ()

    def plus_const(self, value: int) -> _Addr:
        """Return a new address expression with the constant term shifted by ``value``.

        Args:
            value: integer offset added to the constant term; may be positive or negative.

        Returns:
            _Addr: new address expression with the shifted constant term; the original is unchanged.
        """
        return _Addr(self.constant + value, self.terms)

    def plus_ref(self, ref: Ref, coefficient: int = 1) -> _Addr:
        """Return a new address expression with one term ``ref × coefficient`` appended.

        Args:
            ref: register view participating in the address sum.
            coefficient: integer coefficient multiplying the view; defaults to 1.

        Returns:
            _Addr: new address expression with the term appended; the original is unchanged.
        """
        return _Addr(self.constant, self.terms + ((ref, coefficient),))


def _fresh_local(builder: Builder, prefix: str, width: int) -> Ref:
    """Allocate the first free ``{prefix}_{i}`` ``UInt(width)`` local register and return its view."""
    if getattr(builder, "_closed", False):
        raise ValidationError("the builder is already closed")
    existing = {r.name for r in builder.locals} | {
        r.name for r in getattr(builder, "registers", ())
    }
    index = 1
    while f"{prefix}_{index}" in existing:
        index += 1
    return builder.local(f"{prefix}_{index}", UInt(width))


def _ripple_add(builder: Builder, target: Ref, source: Ref, carry: Ref | None) -> None:
    """target += source modulo 2^w, source unchanged; carry is zero-initialized workspace.

    Bit-by-bit carry chain: c_{i+1} = maj(target_i, source_i, c_i), sum_i = target_i ⊕ source_i ⊕ c_i.
    Carry bits stay in carry and are restored to zero by replaying under an overall Adjoint.
    """
    width = target.width
    for i in range(width):
        bit_t, bit_s = target[i], source[i]
        if carry is not None and i + 1 < width:
            with builder.control(bit_t, 1), builder.control(bit_s, 1):
                builder.x(carry[i])
            if i:
                with builder.control(bit_t, 1), builder.control(carry[i - 1], 1):
                    builder.x(carry[i])
                with builder.control(bit_s, 1), builder.control(carry[i - 1], 1):
                    builder.x(carry[i])
        builder.xor(bit_s, bit_t)
        if i:
            # i>0 implies width>1, in which case the caller guarantees carry is allocated (span>1).
            builder.xor(cast(Ref, carry)[i - 1], bit_t)


def _run_steps(
    builder: Builder, steps: list[tuple[str, Ref, Ref | None, Ref | int | None]]
) -> None:
    """Synthesize an address from a step sequence: copy is an XOR copy, add is bit-by-bit addition, const is a constant addition."""
    for kind, target, source, extra in steps:
        if kind == "copy":
            builder.xor(cast(Ref, source), target)
        elif kind == "add":
            _ripple_add(builder, target, cast(Ref, source), cast("Ref | None", extra))
        else:
            builder.add_const(target, cast(int, extra))


class QMem:
    """Multidimensional array view of a QRAM resource: row-major flattening, with indices being integers or register views."""

    def __init__(
        self,
        builder: Builder,
        resource: str,
        *,
        shape: tuple[int, ...] | None = None,
        _base: _Addr | None = None,
    ) -> None:
        """Initialize the QRAM array view.

        Args:
            builder: builder that has declared this QRAM resource.
            resource: QRAM resource name.
            shape: row-major shape tuple; defaults to a one-dimensional array covering the entire address space.
            _base: base address expression shared by the view; defaults to starting at address zero.

        Raises:
            ValidationError: the resource is not declared in the builder, shape is invalid, or the total cell count exceeds the address space.
        """
        specs = {r.name: r.type for r in builder.resources}
        if resource not in specs:
            raise ValidationError(f"builder has not declared the QRAM resource: {resource}")
        self.builder: Builder = builder
        self.resource: str = resource
        self.spec: QRAM = specs[resource]
        if shape is None:
            shape = (1 << self.spec.address_width,)
        if not isinstance(shape, tuple) or not shape:
            raise ValidationError("shape must be a nonempty tuple of positive integers")
        cells = 1
        for dim in shape:
            if type(dim) is not int or dim <= 0:
                raise ValidationError("shape must be a nonempty tuple of positive integers")
            cells *= dim
        if cells > 1 << self.spec.address_width:
            raise ValidationError("total cell count of shape exceeds the QRAM address space")
        self.shape: tuple[int, ...] = shape
        self._base: _Addr = _base if _base is not None else _Addr()

    @property
    def strides(self) -> tuple[int, ...]:
        """Row-major strides per dimension, i.e. the address increment of advancing one index along that dimension."""
        result = [1] * len(self.shape)
        for axis in range(len(self.shape) - 2, -1, -1):
            result[axis] = result[axis + 1] * self.shape[axis + 1]
        return tuple(result)

    @property
    def address_width(self) -> int:
        """Address width of the underlying QRAM resource."""
        return self.spec.address_width

    @property
    def data_width(self) -> int:
        """Data word width of the underlying QRAM resource."""
        return self.spec.data_width

    def ptr(self, base: int | Ref | None = None) -> QPtr:
        """Return a pointer; base is an integer offset or a register view holding an address, a quantum pointer.

        Args:
            base: integer offset, or a register view whose bit width does not exceed
                the address width; omitted keeps the base address of the current view.

        Returns:
            QPtr: quantum memory pointer to the shifted base address.

        Raises:
            ValidationError: the base is neither an integer nor a register view, or its bit width exceeds the limit.
        """
        if base is None:
            return QPtr(self, self._base)
        if isinstance(base, Ref):
            if base.width > self.spec.address_width:
                raise ValidationError("pointer register width exceeds the address width")
            return QPtr(self, self._base.plus_ref(base))
        if type(base) is int:
            return QPtr(self, self._base.plus_const(base))
        raise ValidationError("pointer base address must be an integer or a register view")

    def __getitem__(
        self, key: int | slice | Ref | tuple[int | slice | Ref, ...]
    ) -> QMem | QPtr:
        """Access the view by multidimensional indices.

        Args:
            key: per-dimension index as an integer, register view or slice; missing dimensions are completed as full-dimension slices.

        Returns:
            QMem | QPtr: sub-array view when a sliced dimension is present, or a pointer when every dimension resolves exactly.

        Raises:
            ValidationError: the dimension count does not match shape, a slice step is not 1, an index is out of bounds, or a type is invalid.
        """
        if not isinstance(key, tuple):
            key = (key,)
        if len(key) < len(self.shape):
            key = key + (slice(None),) * (len(self.shape) - len(key))
        if len(key) != len(self.shape):
            raise ValidationError("index dimension count does not match shape")
        addr, sliced, new_shape = self._base, False, list(self.shape)
        for axis, (index, size, stride) in enumerate(
            zip(key, self.shape, self.strides, strict=True)
        ):
            if isinstance(index, slice):
                if index.step not in (None, 1):
                    raise ValidationError("only slice steps of 1 are supported")
                start = 0 if index.start is None else index.start
                stop = size if index.stop is None else index.stop
                if type(start) is not int or type(stop) is not int or not 0 <= start < stop <= size:
                    raise ValidationError("invalid slice bounds")
                addr = addr.plus_const(start * stride)
                new_shape[axis] = stop - start
                sliced = True
            elif isinstance(index, Ref):
                if (1 << index.width) > size:
                    raise ValidationError(f"quantum index width exceeds the dimension length for axis {axis}")
                addr = addr.plus_ref(index, stride)
            elif type(index) is int:
                if not 0 <= index < size:
                    raise ValidationError("index out of bounds")
                addr = addr.plus_const(index * stride)
            else:
                raise ValidationError("indices must be integers, register views or slices")
        if sliced:
            return QMem(self.builder, self.resource, shape=tuple(new_shape), _base=addr)
        return QPtr(self, addr)

    def _emit(self, addr: _Addr, data: Ref, kind: str) -> None:
        """Dereference: materialize the address expression, then emit a Load (XOR read) or Store (random write)."""
        builder, spec = self.builder, self.spec
        if data.width != spec.data_width:
            raise ValidationError("data register width does not match the QRAM data width")
        address_width = spec.address_width
        constant = addr.constant % (1 << address_width)
        direct = (
            len(addr.terms) == 1
            and addr.terms[0][1] == 1
            and addr.terms[0][0].width == address_width
            and constant == 0
        )
        if direct:
            ref = addr.terms[0][0]
            if kind == "load":
                builder.qram(self.resource, ref, data)
            else:
                builder.store(self.resource, ref, data)
            return
        temp = _fresh_local(builder, "qmem_addr", address_width)
        steps: list[tuple[str, Ref, Ref | None, Ref | int | None]] = []
        occupied: list[tuple[int, int]] = []
        for ref, coefficient in addr.terms:
            if coefficient < 1:
                raise ValidationError("address coefficients must be positive integers")
            for shift in range(coefficient.bit_length()):
                if not (coefficient >> shift) & 1:
                    continue
                span = min(ref.width, address_width - shift)
                if span <= 0:
                    break
                target, source = temp[shift : shift + span], ref[:span]
                # When a shifted fragment falls in a not-yet-written range, that segment of temp is still zero, so XOR assigns it and no adder is needed;
                # only overlaps with an occupied range (such as a pointer + quantum offset sum) synthesize a bit-by-bit carry addition.
                if all(shift >= end or shift + span <= start for start, end in occupied):
                    steps.append(("copy", target, source, None))
                else:
                    carry = _fresh_local(builder, "qmem_carry", span - 1) if span > 1 else None
                    steps.append(("add", target, source, carry))
                occupied.append((shift, shift + span))
        if constant:
            steps.append(("const", temp, None, constant))
        _run_steps(builder, steps)
        if kind == "load":
            builder.qram(self.resource, temp, data)
        else:
            builder.store(self.resource, temp, data)
        if steps:
            with builder.adjoint():
                _run_steps(builder, steps)


class QPtr:
    """QRAM pointer: base address + offset, dereferenceable for reads (XOR-Load) or random writes (Store)."""

    def __init__(self, mem: QMem, addr: _Addr) -> None:
        """Bind the memory view the pointer belongs to and its address expression.

        Args:
            mem: array view the pointer targets.
            addr: address expression held by the pointer.
        """
        self.mem: QMem = mem
        self.addr: _Addr = addr

    def __add__(self, other: int | Ref) -> QPtr:
        """Return a new pointer advanced along the address by ``other``.

        Args:
            other: integer offset or a register view holding the offset.

        Returns:
            QPtr: new pointer after the offset; the original pointer is unchanged.

        Raises:
            ValidationError: the offset is neither an integer nor a register view.
        """
        if type(other) is int:
            return QPtr(self.mem, self.addr.plus_const(other))
        if isinstance(other, Ref):
            return QPtr(self.mem, self.addr.plus_ref(other))
        raise ValidationError("pointer offsets must be integers or register views")

    def __radd__(self, other: int | Ref) -> QPtr:
        """Support the ``integer + pointer`` form by forwarding to ``__add__``.

        Args:
            other: integer offset or a register view holding the offset.

        Returns:
            QPtr: new pointer after the offset.
        """
        return self.__add__(other)

    def __sub__(self, other: int) -> QPtr:
        """Return a new pointer moved back along the address by the integer ``other``.

        Args:
            other: integer offset.

        Returns:
            QPtr: new pointer with the offset subtracted from the base address.

        Raises:
            ValidationError: the offset is not an integer.
        """
        if type(other) is not int:
            raise ValidationError("pointer subtraction only supports integer offsets")
        return QPtr(self.mem, self.addr.plus_const(-other))

    def load(self, data: Ref) -> None:
        """XOR-Load: data ^= M[addr]; address computation is automatically restored to zero, and superposed addresses are naturally supported.

        Args:
            data: view of the same width as the data word, into which the hit content is XORed.
        """
        self.mem._emit(self.addr, data, "load")

    def store(self, data: Ref) -> None:
        """Random write: M[addr] := d. At execution time the address and data must be in a definite basis state; no gate cost is counted.

        Args:
            data: view of the same width as the data word, providing the content to write.
        """
        self.mem._emit(self.addr, data, "store")
