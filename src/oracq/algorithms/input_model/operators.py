"Block-encoding composition library; alpha is stored as an RIR module attribute and precision is not part of the language core."

from __future__ import annotations

import cmath
import hashlib
import math
from dataclasses import dataclass, replace
from typing import Protocol, cast

from oracq.algorithms.input_model.contracts import OracleView, validate_signature
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import QRAM, Bits, Ref, ValidationError
from oracq.infrastructure.serialization import dumps


class Generator(Protocol):
    """Structural protocol for quantum-operation generators.

    Any callable that is invoked with arbitrary positional and keyword arguments
    and returns an ``Operation`` structurally satisfies this protocol; the
    algorithm side thereby accepts operation factories without binding a
    concrete generation signature.
    """

    def __call__(self, *args: object, **kwargs: object) -> Operation:
        """Invoke the generator and return the operation it produces."""
        ...


def _name(kind: str, *values: object) -> str:
    """Generate a deterministic module name with the ``kind_`` prefix from the serialized content of the values."""
    serialized = [dumps(v.program()) if isinstance(v, Operation) else repr(v) for v in values]
    return kind + "_" + hashlib.sha256("\n".join(serialized).encode()).hexdigest()[:20]


@dataclass(frozen=True)
class BlockEncoding(OracleView):
    """View representing a linear operator by its zero-signal projected corner block.

    The wrapped ``Operation`` contains exactly the two bits registers ``target``
    and ``signal``, and declares a finite positive normalization constant
    through the module attribute ``be_alpha``: if the zero-signal corner block
    of the unitary ``U`` satisfies ``<0|U|0> = A/alpha``, this view treats ``U``
    as a block encoding of the operator ``A``. alpha participates in the
    composition algebra at generation time only; the executor does not scale
    quantum states by it.

    Attributes:
        operation: The wrapped ``Operation``.
    """

    oracle_kind = "block_encoding"
    operation: Operation

    def block_encoding(self) -> BlockEncoding:
        """Return self; the view-adaptation method implementing ``BlockEncodingProtocol``.

        Returns:
            BlockEncoding: This view itself.
        """
        return self

    def __post_init__(self) -> None:
        """Validate the target/signal signature and the finite-positive ``be_alpha`` attribute constraint."""
        validate_signature(self.operation, ("target", "signal"), "BlockEncoding")
        registers = {r.name: r.type for r in self.operation.module.registers}
        if set(registers) != {"target", "signal"} or any(
            t.kind != "bits" for t in registers.values()
        ):
            raise ValidationError("a BE must have exactly the two bits interfaces target and signal")
        if not registers["target"].width:
            raise ValidationError("a BE target cannot be empty")
        alpha = dict(self.operation.module.attributes).get("be_alpha")
        if (
            type(alpha) not in (int, float)
            or not math.isfinite(cast("int | float", alpha))
            or cast("int | float", alpha) <= 0
        ):
            raise ValidationError("a BE must declare a finite positive be_alpha")

    @property
    def alpha(self) -> float:
        """The normalization constant declared in the ``be_alpha`` module attribute."""
        return cast("float", dict(self.operation.module.attributes)["be_alpha"])

    @property
    def width(self) -> int:
        """Bit width of the ``target`` register."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def signal_qubits(self) -> int:
        """Bit width of the ``signal`` register."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "signal")


def block_encoding(operation: Operation, alpha: float = 1.0) -> BlockEncoding:
    """Wrap an operation with the ``target``/``signal`` signature as a block encoding.

    Writes ``be_alpha`` and ``oracle_paradigm`` into the module attributes and
    deterministically renames the module by content; the signature and alpha
    constraints are checked by the ``BlockEncoding`` constructor.

    Args:
        operation: An ``Operation`` containing exactly the ``target`` and ``signal`` bits registers.
        alpha: Finite positive normalization constant, defaulting to one.

    Returns:
        BlockEncoding: Block-encoding view with the attributes filled in.

    Raises:
        ValidationError: The operation violates the block-encoding signature or alpha constraint.
    """
    attributes = dict(operation.module.attributes)
    attributes["be_alpha"] = alpha
    attributes["oracle_paradigm"] = "block_encoding"
    module = replace(
        operation.module,
        name=_name("be", operation, alpha),
        attributes=tuple(sorted(attributes.items())),
    )
    return BlockEncoding(Operation(module, operation.dependencies))


def identity(width: int) -> BlockEncoding:
    """Build the block encoding of the identity operator; alpha is one and no signal bit is needed.

    Args:
        width: Target register bit width.

    Returns:
        BlockEncoding: Block encoding of the identity operator.
    """
    b = Builder(f"identity_{width}", {"target": Bits(width), "signal": Bits(0)})
    return block_encoding(b.finish())


def pauli_x(width: int) -> BlockEncoding:
    """Build the block encoding of the tensor power of ``width`` X gates; alpha is one with no signal bit.

    Args:
        width: Target register bit width.

    Returns:
        BlockEncoding: Block encoding of ``X^⊗width``.
    """
    b = Builder(f"pauli_x_{width}", {"target": Bits(width), "signal": Bits(0)})
    b.x(b["target"])
    return block_encoding(b.finish())


def zero(width: int) -> BlockEncoding:
    """Build the block encoding of the zero operator; a flipped signal bit makes the zero-signal corner block identically zero, with alpha one.

    Args:
        width: Target register bit width.

    Returns:
        BlockEncoding: Block encoding of the zero operator.
    """
    b = Builder(f"zero_{width}", {"target": Bits(width), "signal": Bits(1)})
    b.x(b["signal"])
    return block_encoding(b.finish())


def _resources(a: BlockEncoding, b: BlockEncoding | None = None) -> dict[str, QRAM]:
    """Collect the resources declared by each operation and rename them with the ``a__``/``b__`` prefixes."""
    result: dict[str, QRAM] = {}
    for prefix, operand in (("a", a), ("b", b)):
        if operand is not None:
            for resource in operand.operation.module.resources:
                result[prefix + "__" + resource.name] = resource.type
    return result


def _call(
    builder: Builder, operand: BlockEncoding, target: Ref, signal: Ref, prefix: str
) -> None:
    """Invoke the block-encoding operation in ``builder`` after mapping resources by prefix."""
    resources = {r.name: prefix + "__" + r.name for r in operand.operation.module.resources}
    builder.call(operand.operation, target=target, signal=signal, resources=resources)


def product(a: BlockEncoding, b: BlockEncoding) -> BlockEncoding:
    """Compose the matrix product ``A·B`` of two block encodings.

    The two operations act in turn on the shared ``target`` (``b`` first, then
    ``a``), with the signal bits concatenated at the high end for ``a``; the
    returned block encoding has alpha ``a.alpha * b.alpha`` and a signal width
    equal to the sum of both.

    Args:
        a: Left-factor block encoding.
        b: Right-factor block encoding.

    Returns:
        BlockEncoding: Block encoding of ``A·B``, with both operations kept as module calls.

    Raises:
        ValidationError: The two target widths differ.
    """
    if a.width != b.width:
        raise ValidationError("BE product has mismatched target widths")
    builder = Builder(
        _name("product", a.operation, b.operation),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits + b.signal_qubits)},
        _resources(a, b),
    )
    signal = builder["signal"]
    _call(builder, b, builder["target"], signal[a.signal_qubits :], "b")
    _call(builder, a, builder["target"], signal[: a.signal_qubits], "a")
    return block_encoding(builder.finish(), a.alpha * b.alpha)


def scale(coefficient: complex, a: BlockEncoding) -> BlockEncoding:
    """Scale the operator represented by a block encoding by a complex coefficient.

    The coefficient phase is booked with ``global_phase``, and the returned
    block encoding has alpha ``abs(coefficient) * a.alpha``; when the
    coefficient is zero, ``zero(a.width)`` is returned directly.

    Args:
        coefficient: Finite complex coefficient.
        a: The block encoding being scaled.

    Returns:
        BlockEncoding: Block encoding of ``coefficient * A``.

    Raises:
        ValidationError: The real or imaginary part of the coefficient is not finite.
    """
    if not (math.isfinite(coefficient.real) and math.isfinite(coefficient.imag)):
        raise ValidationError("BE coefficient must be finite")
    if coefficient == 0:
        return zero(a.width)
    builder = Builder(
        _name("scale", a.operation, coefficient),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        _resources(a),
    )
    builder.global_phase(cmath.phase(coefficient))
    _call(builder, a, builder["target"], builder["signal"], "a")
    return block_encoding(builder.finish(), abs(coefficient) * a.alpha)


def linear_combination(
    ca: complex, a: BlockEncoding, cb: complex, b: BlockEncoding
) -> BlockEncoding:
    """Compose the linear combination ``ca*A + cb*B`` of two block encodings.

    A single select bit branches in the weight ratio of ``abs(ca)*a.alpha`` and
    ``abs(cb)*b.alpha``; each branch compensates its coefficient phase via
    ``global_phase`` before acting on the shared ``target``. The returned block
    encoding has alpha equal to the sum of the two weights, and its signal bits
    carry one extra select bit beyond the two operations' signals. When one
    coefficient is zero it degenerates to ``scale``.

    Args:
        ca: Finite complex coefficient of ``a``.
        a: First block encoding.
        cb: Finite complex coefficient of ``b``.
        b: Second block encoding.

    Returns:
        BlockEncoding: Block encoding of ``ca*A + cb*B``.

    Raises:
        ValidationError: The target widths differ, or any coefficient contains a non-finite component.
    """
    if a.width != b.width:
        raise ValidationError("BE sum has mismatched target widths")
    if not all(math.isfinite(x) for c in (ca, cb) for x in (c.real, c.imag)):
        raise ValidationError("BE coefficients must be finite")
    if ca == 0:
        return scale(cb, b)
    if cb == 0:
        return scale(ca, a)
    alpha = abs(ca) * a.alpha + abs(cb) * b.alpha
    theta = 2 * math.acos(math.sqrt(abs(ca) * a.alpha / alpha))
    builder = Builder(
        _name("lcu", a.operation, b.operation, ca, cb),
        {"target": Bits(a.width), "signal": Bits(1 + a.signal_qubits + b.signal_qubits)},
        _resources(a, b),
    )
    signal = builder["signal"]
    select, sa, sb = signal[:1], signal[1 : 1 + a.signal_qubits], signal[1 + a.signal_qubits :]
    builder.ry(select, theta)
    with builder.control(select, 0):
        builder.global_phase(cmath.phase(ca))
        _call(builder, a, builder["target"], sa, "a")
    with builder.control(select, 1):
        builder.global_phase(cmath.phase(cb))
        _call(builder, b, builder["target"], sb, "b")
    builder.ry(select, -theta)
    return block_encoding(builder.finish(), alpha)
