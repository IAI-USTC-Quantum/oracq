"""Sparse position and entry access, plus explicit adaptation of numeric words to block encodings."""

from __future__ import annotations

import math
from collections.abc import Mapping
from functools import lru_cache

from oracq.algorithms.common.arithmetic import BooleanNetwork, FixedFormat
from oracq.algorithms.input_model.block_encoding import reflect_zero
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import (
    SparseAccess,
    XorDatabase,
    annotate,
    declare,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Ref, RegType, ValidationError, fuse


def reversible_lookup(
    inputs: Mapping[str, int],
    outputs: Mapping[str, int],
    database: XorDatabase,
    *,
    name: str | None = None,
) -> Operation:
    """Adapt the address/data interface of an ``XorDatabase`` to named input/output register groups.

    The input registers are fused, in declaration order, into an address view
    and the output registers into a data view; the query semantics remain the
    XOR database's ``data ^= memory[address]``.

    Args:
        inputs: Mapping from names to bit widths; the total width must equal ``database.address_width``.
        outputs: Mapping from names to bit widths; the total width must equal ``database.data_width``.
        database: The XOR database being adapted.
        name: Optional module name; generated from content when omitted.

    Returns:
        Operation: Registers are the union of inputs and outputs; the paradigm is ``reversible_function``.

    Raises:
        ValidationError: The input/output total widths do not match the database, or names conflict between the two sides.
    """
    if (
        sum(inputs.values()) != database.address_width
        or sum(outputs.values()) != database.data_width
    ):
        raise ValidationError("reversible lookup input and output total widths do not match")
    if set(inputs) & set(outputs):
        raise ValidationError("reversible lookup input and output names conflict")
    b = Builder(
        name
        or _name(
            "reversible_lookup", database.operation, tuple(inputs.items()), tuple(outputs.items())
        ),
        {k: Bits(v) for k, v in {**inputs, **outputs}.items()},
        resources_for(("db", database.operation)),
    )
    invoke(
        b,
        database.operation,
        "db",
        address=fuse(*(b[k] for k in inputs)),
        data=fuse(*(b[k] for k in outputs)),
    )
    return annotate(b.finish(), "reversible_function")


def word_rotation(
    value_width: int, *, scale: float | None = None, name: str | None = None
) -> Operation:
    """Reversible transduction that converts the integer value of a numeric word linearly into an Ry angle on the ``amplitude`` bit.

    With per-bit controlled superposition the total rotation angle is
    ``scale * value`` (value taken as an unsigned integer); the default
    ``scale = 2*pi / 2**value_width`` makes the whole value range sweep exactly
    one full turn.

    Args:
        value_width: Bit width of the value register.
        scale: Rotation angle per unit of integer value; defaults to ``2*pi / 2**value_width``.
        name: Optional module name; generated from content when omitted.

    Returns:
        Operation: Registers are value and the single-bit amplitude; the paradigm
        is ``reversible_function``, with module attribute ``angle_scale``
        recording the scale used.
    """
    scale = 2 * math.pi / (1 << value_width) if scale is None else scale
    b = Builder(
        name or _name("word_rotation", value_width, scale),
        {"value": Bits(value_width), "amplitude": Bits(1)},
    )
    for bit in range(value_width):
        with b.control(b["value"][bit]):
            b.ry(b["amplitude"], scale * (1 << bit))
    return annotate(b.finish(), "reversible_function", angle_scale=scale)


def sparse_block_encoding(
    access: SparseAccess, transducer: Operation | None = None, *, alpha: float | None = None
) -> BlockEncoding:
    """Historical candidate kept only for the old gallery description; the official sparse adapter is real_symmetric_sparse_encoding.

    Args:
        access: CKS sparse access bundle.
        transducer: Value-to-amplitude transduction operation; defaults to an open declaration.
        alpha: Explicitly overridden normalization constant; defaults to the sparsity.

    Returns:
        BlockEncoding: Legacy transduction construction with an unspecified contract.
    """
    n, v = access.width, access.value_width
    lw = next(r.type.width for r in access.location.module.registers if r.name == "work")
    transducer = transducer or declare(
        _name("sparse_amplitude", n, v),
        {"value": Bits(v), "amplitude": Bits(1)},
        paradigm="reversible_function",
    )
    uniform = gate_state_prep([1 if i < access.sparsity else 0 for i in range(1 << n)])
    width = n + lw + v + 1
    b = Builder(
        _name("sparse_prepare", access.location, access.entry, transducer, access.sparsity),
        {"target": Bits(n), "signal": Bits(width)},
        resources_for(
            ("position", access.location), ("entry", access.entry), ("amplitude", transducer)
        ),
        attributes={"algorithm": "sparse_isometry_extension", "validation_stage": "paradigm"},
    )
    neighbor, work = b["signal"][:n], b["signal"][n : n + lw]
    value, amplitude = b["signal"][n + lw : n + lw + v], b["signal"][width - 1 :]
    invoke(b, uniform.operation, target=neighbor, work=neighbor[:0])
    invoke(b, access.location, "position", column=b["target"], index=neighbor, work=work)
    invoke(b, access.entry, "entry", row=neighbor, column=b["target"], data=value)
    invoke(b, transducer, "amplitude", value=value, amplitude=amplitude)
    with b.adjoint():
        invoke(b, access.entry, "entry", row=neighbor, column=b["target"], data=value)
    prep = b.finish()
    out = Builder(
        _name("sparse_be", prep),
        {"target": Bits(n), "signal": Bits(width)},
        resources_for(("prep", prep)),
    )
    invoke(out, prep, "prep", target=out["target"], signal=out["signal"])
    out.swap(out["target"], out["signal"][:n])
    with out.adjoint():
        invoke(out, prep, "prep", target=out["target"], signal=out["signal"])
    return BlockEncoding(
        annotate(
            out.finish(),
            "block_encoding",
            be_alpha=float(access.sparsity if alpha is None else alpha),
            construction="Tdag_SWAP_T",
            validation_stage="paradigm",
            matrix_contract="legacy transduction unspecified; not a general sparse input adapter",
            legacy_input_model=True,
        )
    )


def batch_lookup(database: XorDatabase, count: int) -> Operation:
    """Batch circuit issuing multiple concurrent queries to the same XOR database.

    Args:
        database: The XOR database.
        count: Number of query lanes.

    Returns:
        Operation: Registers are ``address{i}`` and ``data{i}`` for i from 0 to
        count-1, with bit widths equal to the database's address/data widths
        respectively, invoking the same database operation once per lane.
    """
    registers: dict[str, RegType] = {}
    for i in range(count):
        registers[f"address{i}"] = Bits(database.address_width)
        registers[f"data{i}"] = Bits(database.data_width)
    b = Builder(
        _name("batch_lookup", database.operation, count),
        registers,
        resources_for(("db", database.operation)),
    )
    for i in range(count):
        invoke(b, database.operation, "db", address=b[f"address{i}"], data=b[f"data{i}"])
    return b.finish()


@lru_cache(maxsize=64)
def compare_words(width: int, kind: str = "eq") -> Operation:
    """Equality or less-than comparison network for two word-width unsigned integers.

    Args:
        width: Bit width of each input word.
        kind: ``"eq"`` builds an equality predicate, ``"lt"`` an unsigned less-than predicate.

    Returns:
        Operation: Input registers a and b, output single-bit flag set to 1 when
        the comparison holds. The Boolean network is compiled into a reversible
        quantum operation via compute/copy/uncompute, with the private bank zero
        in and zero out. Results are cached and reused by ``(width, kind)``.
    """
    net = BooleanNetwork()
    a, b = net.input("a", width), net.input("b", width)
    bit = (
        net.inv(net.any([net.xor(x, y) for x, y in zip(a, b, strict=True)]))
        if kind == "eq"
        else net.lt(a, b)
    )
    net.outputs = {"flag": [bit]}
    return net.operation()


@lru_cache(maxsize=64)
def value_transposition(width: int) -> Operation:
    """Swap the two bit patterns a and b within index; a and b are preserved, suitable for quantum addresses.

    Args:
        width: Bit width of each of the index and a, b registers.

    Returns:
        Operation: Operation that exchanges the basis states of index equal to a
        or b while leaving all other basis states unchanged.
    """
    net = BooleanNetwork()
    x, a, c = net.input("index", width), net.input("a", width), net.input("b", width)

    def eq(y: list[int]) -> int:
        """Output 1 if and only if ``y`` equals the input ``x`` bit by bit."""
        return net.inv(net.any([net.xor(v, w) for v, w in zip(x, y, strict=True)]))

    net.outputs = {"flag": [net.or_(eq(a), eq(c))]}
    predicate = net.operation()
    b = Builder(
        "value_transposition_" + str(width),
        {"index": Bits(width), "a": Bits(width), "b": Bits(width)},
    )
    flag = b.local("membership", Bits(1))
    b.call(predicate, index=b["index"], a=b["a"], b=b["b"], flag=flag)
    with b.control(flag):
        b.xor(b["a"], b["index"])
        b.xor(b["b"], b["index"])
    b.call(predicate, index=b["index"], a=b["a"], b=b["b"], flag=flag)
    return b.finish()


def prefix_state(width: int, count: int) -> Operation:
    """Prepare a uniform superposition over the first count basis states.

    Args:
        width: Target register bit width.
        count: Number of basis states covered by the superposition, in ``1..2**width``.

    Returns:
        Operation: Zero-input preparation on a single target register whose
        support is ``|0>`` through ``|count-1>`` with equal amplitudes.

    Raises:
        ValidationError: count is outside the ``1..2**width`` range.
    """
    if not 1 <= count <= 1 << width:
        raise ValidationError("uniform prefix range is invalid")
    b = Builder("uniform_prefix_" + str(width) + "_" + str(count), {"target": Bits(width)})

    def prepare(ref: Ref, size: int) -> None:
        """Recursively prepare a uniform superposition over the first ``size`` basis states of ``ref``."""
        if not ref.width:
            return
        if size == 1 << ref.width:
            b.h(ref)
            return
        half = 1 << (ref.width - 1)
        left, right = min(size, half), max(0, size - half)
        b.ry(ref[ref.width - 1], 2 * math.asin(math.sqrt(right / size)))
        with b.control(ref[ref.width - 1], 0):
            prepare(ref[: ref.width - 1], left)
        if right:
            with b.control(ref[ref.width - 1]):
                prepare(ref[: ref.width - 1], right)

    prepare(b["target"], count)
    return b.finish()


def magnitude_rotation(fmt: FixedFormat, amax: float) -> Operation:
    """Small plain implementation for numeric words; above 12 bits an explicit to-be-bound transducer is kept.

    Args:
        fmt: Fixed-point format of entry values.
        amax: Upper bound on entry magnitudes, a positive finite real number.

    Returns:
        Operation: Controlled Ry transduction with amplitude ``sqrt(|value|/amax)``;
        for formats wider than 12 bits an open declaration awaiting binding is
        returned.
    """
    name = _name("sparse_sqrt_rotation", fmt, amax)
    registers = {"value": Bits(fmt.width), "amplitude": Bits(1)}
    attrs: dict[str, str | int | float | bool] = {
        "entry_bound": float(amax),
        "value_fraction": fmt.fraction,
        "amplitude_contract": "good amplitude sqrt(abs(value)/entry_bound)",
    }
    if fmt.width > 12:
        return declare(name, registers, paradigm="reversible_function", attributes=attrs)
    b = Builder(name, registers, attributes=attrs)
    for raw in range(1 << fmt.width):
        magnitude = abs(fmt.decode(raw)) / amax
        angle = 2 * math.acos(math.sqrt(min(1, magnitude)))
        if angle:
            with b.control(b["value"], raw):
                b.ry(b["amplitude"], angle)
    return b.finish()


def real_symmetric_sparse_encoding(
    access: SparseAccess,
    fmt: FixedFormat,
    amax: float,
    *,
    diagonal_nonnegative: bool = False,
    rotation: Operation | None = None,
) -> BlockEncoding:
    """CKS-style T†ST: real Hermitian with non-negative diagonal; swaps the coordinates and the failure flag between the two sides.

    Args:
        access: CKS sparse access bundle.
        fmt: Fixed-point format of entry values; the bit width must match the access bundle's value width.
        amax: Upper bound on entry magnitudes, a positive finite real number.
        diagonal_nonnegative: Must be True; only matrices with a non-negative diagonal are currently supported.
        rotation: Optional amplitude transduction operation; defaults to one generated by magnitude_rotation.

    Returns:
        BlockEncoding: Block encoding of a T†ST-type real symmetric sparse matrix.
    """
    if not diagonal_nonnegative:
        raise ValidationError("the symmetric sparse adapter currently requires a non-negative diagonal; use an explicit Hermitian dilation for general matrices")
    if not math.isfinite(amax) or amax <= 0 or fmt.width != access.value_width:
        raise ValidationError("invalid sparse value format or entry bound")
    n = access.width
    rotation = rotation or magnitude_rotation(fmt, amax)
    location_work = next(r.type.width for r in access.location.module.registers if r.name == "work")
    operands = (("location", access.location), ("entry", access.entry), ("rotation", rotation))
    b = Builder(
        _name("cks_isometry", access.location, access.entry, rotation, amax),
        {"target": Bits(n), "signal": Bits(n + 2)},
        resources_for(*operands),
        attributes={
            "algorithm": "cks_real_symmetric_isometry",
            "entry_bound": float(amax),
            "sparsity": access.sparsity,
            "correctness": "pending",
            "sign_convention": "negative offdiagonal phase pi only for target < neighbor",
        },
    )
    neighbor, _first_flag, second_flag = b["signal"][:n], b["signal"][n], b["signal"][n + 1]
    value = b.local("entry_value", Bits(fmt.width))
    work = b.local("location_work", Bits(location_work))
    order = b.local("ordered", Bits(1))
    b.call(prefix_state(n, access.sparsity), target=neighbor)
    invoke(b, access.location, "location", column=b["target"], index=neighbor, work=work)
    invoke(b, access.entry, "entry", row=neighbor, column=b["target"], data=value)
    invoke(b, rotation, "rotation", value=value, amplitude=second_flag)
    if fmt.signed:
        b.call(compare_words(n, "lt"), a=b["target"], b=neighbor, flag=order)
        with b.control(fuse(order, value[fmt.width - 1])):
            b.global_phase(math.pi)
        b.call(compare_words(n, "lt"), a=b["target"], b=neighbor, flag=order)
    invoke(b, access.entry, "entry", row=neighbor, column=b["target"], data=value)
    preparation = b.finish()
    out = Builder(
        _name("cks_sparse_be", preparation),
        {"target": Bits(n), "signal": Bits(n + 2)},
        resources_for(("prep", preparation)),
    )
    invoke(out, preparation, "prep", target=out["target"], signal=out["signal"])
    out.swap(out["target"], out["signal"][:n])
    out.swap(out["signal"][n], out["signal"][n + 1])
    with out.adjoint():
        invoke(out, preparation, "prep", target=out["target"], signal=out["signal"])
    return BlockEncoding(
        annotate(
            out.finish(),
            "block_encoding",
            be_alpha=access.sparsity * amax,
            construction="CKS_Tdag_S_T",
            self_adjoint_extension=True,
            correctness="pending",
        )
    )


def chebyshev_block(a: BlockEncoding, degree: int) -> BlockEncoding:
    """Chebyshev walk power: the zero-signal block realizes the degree-th Chebyshev polynomial of the scaled matrix.

    The degree is kept as a Repeat, each step alternating between a positive
    reflection about the signal zero state and an invocation of ``a``; the
    encoded matrix enters the polynomial scaled by ``a.alpha``.

    Args:
        a: A ``BlockEncoding`` whose module attributes explicitly declare ``self_adjoint_extension``.
        degree: Chebyshev degree, a non-negative integer.

    Returns:
        BlockEncoding: ``be_alpha`` is 1.0, ``argument_scale`` records ``a.alpha``,
        and the zero-signal block is ``T_degree(A / a.alpha)``.

    Raises:
        ValidationError: The input block encoding lacks the self-adjoint unitary extension declaration, or the degree is negative.
    """
    if not dict(a.operation.module.attributes).get("self_adjoint_extension"):
        raise ValidationError("Chebyshev walk requires an explicit self-adjoint unitary extension; a plain BE does not suffice")
    if degree < 0:
        raise ValidationError("Chebyshev degree cannot be negative")
    b = Builder(
        _name("chebyshev_walk_power", a.operation, degree),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        resources_for(("a", a.operation)),
    )
    with b.repeat(degree):
        reflect_zero(b, b["signal"], positive=True)
        invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
    return BlockEncoding(
        annotate(
            b.finish(),
            "block_encoding",
            be_alpha=1.0,
            polynomial="Chebyshev",
            degree=degree,
            argument_scale=a.alpha,
            correctness="pending",
        )
    )
