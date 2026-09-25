"""Deutsch–Jozsa, Bernstein–Vazirani, and Simon query algorithms."""

from __future__ import annotations

from collections.abc import Iterable

from oracq.algorithms.input_model.contracts import positive_integer
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import XorDatabase, annotate, invoke, resources_for
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError


def deutsch_jozsa(function: XorDatabase) -> Operation:
    """Generate the D-J query circuit; the caller asserts the function is constant or balanced.

    Args:
        function: An XOR database with one-bit results, satisfying the
            constant or balanced promise.

    Returns:
        Operation: Query circuit with registers input and answer; read out
        input to decide constant versus balanced.
    """
    if function.data_width != 1:
        raise ValidationError("the D-J oracle must have exactly one result bit")
    b = Builder(
        _name("deutsch_jozsa", function.operation),
        {"input": Bits(function.address_width), "answer": Bits(1)},
        resources_for(("function", function.operation)),
        attributes={
            "algorithm": "deutsch_jozsa",
            "readout_register": "input",
            "validation_stage": "paradigm",
        },
    )
    b.x(b["answer"])
    b.h(b["answer"])
    b.h(b["input"])
    invoke(b, function.operation, "function", address=b["input"], data=b["answer"])
    b.h(b["input"])
    return b.finish()


def bernstein_vazirani(function: XorDatabase) -> Operation:
    """Generate the Bernstein–Vazirani secret string readout circuit.

    Args:
        function: An XOR database with one-bit results; the caller guarantees
            f(x)=s·x XOR c.

    Returns:
        Operation: Exposes input and answer. When the input promise holds,
        reading input yields s.

    The oracle may stay as an open declaration, with a gate or QRAM
    implementation bound later."""
    from dataclasses import replace

    from oracq.infrastructure.builder import Operation

    operation = deutsch_jozsa(function)
    attrs = dict(operation.module.attributes)
    attrs.update(algorithm="bernstein_vazirani", input_promise="f(x)=dot(s,x) XOR c over GF(2)")
    return Operation(
        replace(
            operation.module,
            name=_name("bernstein_vazirani", function.operation),
            attributes=tuple(sorted(attrs.items())),
        ),
        operation.dependencies,
    )


def affine_boolean_oracle(width: int, secret: int, *, bias: int = 0) -> XorDatabase:
    """Build a plain gate oracle for BV; bit i of secret corresponds to bit i of the address.

    Args:
        width: Address register bit width, in 1..64.
        secret: Integer encoding of the secret string, in 0..2^width−1.
        bias: Constant bias c, either 0 or 1.

    Returns:
        XorDatabase: Gate-level XOR database handle implementing f(x)=s·x
        XOR c.
    """
    positive_integer(width, "affine_boolean_oracle.width", maximum=64)
    positive_integer(secret, "affine_boolean_oracle.secret", minimum=0, maximum=(1 << width) - 1)
    positive_integer(bias, "affine_boolean_oracle.bias", minimum=0, maximum=1)
    b = Builder(
        _name("affine_boolean", width, secret, bias), {"address": Bits(width), "data": Bits(1)}
    )
    if bias:
        b.x(b["data"])
    for bit in range(width):
        if (secret >> bit) & 1:
            b.xor(b["address"][bit], b["data"])
    return XorDatabase(annotate(b.finish(), "database_xor"))


def simon_sample(function: XorDatabase) -> Operation:
    """Generate one Simon sampling circuit.

    Args:
        function: An XOR database; the caller guarantees a two-to-one map
            with a nonzero XOR period.

    Returns:
        Operation: Exposes input and output. Reading input alone yields an
        orthogonality constraint on the period.

    The output register need not be measured inside the circuit; elimination
    over repeated samples happens on the classical side."""
    if not isinstance(function, XorDatabase):
        raise ValidationError("Simon requires an XOR database")
    b = Builder(
        _name("simon_sample", function.operation),
        {"input": Bits(function.address_width), "output": Bits(function.data_width)},
        resources_for(("function", function.operation)),
        attributes={
            "algorithm": "simon_sample",
            "input_promise": "two-to-one XOR period",
            "readout_register": "input",
        },
    )
    b.h(b["input"])
    invoke(b, function.operation, "function", address=b["input"], data=b["output"])
    b.h(b["input"])
    return b.finish()


def simon_nullspace(samples: Iterable[int], width: int) -> tuple[int, ...]:
    """Compute the GF(2) null space of the Simon sample constraints.

    Args:
        samples: Integer samples encoded in little endian.
        width: Bit width of the secret string.

    Returns:
        tuple: Null space basis vectors. Only a one-dimensional null space
        determines the nonzero period directly; with insufficient samples
        several basis vectors remain."""
    positive_integer(width, "simon.width", maximum=64)
    rows = []
    for sample in samples:
        positive_integer(sample, "simon.sample", minimum=0, maximum=(1 << width) - 1)
        rows.append(sample)
    pivot_rows = {}
    cursor = 0
    for bit in range(width):
        candidate = next((i for i in range(cursor, len(rows)) if (rows[i] >> bit) & 1), None)
        if candidate is None:
            continue
        rows[cursor], rows[candidate] = rows[candidate], rows[cursor]
        for i in range(len(rows)):
            if i != cursor and (rows[i] >> bit) & 1:
                rows[i] ^= rows[cursor]
        pivot_rows[bit] = cursor
        cursor += 1
    basis = []
    for free in range(width):
        if free in pivot_rows:
            continue
        vector = 1 << free
        for pivot, row in pivot_rows.items():
            if (rows[row] >> free) & 1:
                vector |= 1 << pivot
        basis.append(vector)
    return tuple(basis)
