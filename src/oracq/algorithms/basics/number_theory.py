"""Small-scale modular multiplication, quantum order finding, and classical factor post-processing."""

from __future__ import annotations

import math
from fractions import Fraction
from typing import cast

from oracq.algorithms.common.estimation import phase_estimation
from oracq.algorithms.input_model.contracts import positive_integer
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import _transposition, invoke, resources_for
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError


def modular_multiply(
    multiplier: int, modulus: int, *, width: int | None = None, max_width: int = 8
) -> Operation:
    """Generate a reversible bounded-size modular multiplication permutation.

    Args:
        multiplier: Positive integer coprime with modulus.
        modulus: Modulus of at least two.
        width: Target bit width; when omitted, the smallest width that fits
            the modulus.
        max_width: Permutation synthesis budget, default 8, at most 12.

    Returns:
        Operation: Maps x<modulus to multiplier*x mod modulus; other basis
        states are left unchanged.

    The implementation enumerates a finite permutation and does not
    represent scalable Shor modular arithmetic."""
    positive_integer(modulus, "modular_multiply.modulus", minimum=2)
    positive_integer(multiplier, "modular_multiply.multiplier", minimum=1)
    positive_integer(max_width, "modular_multiply.max_width", maximum=12)
    width = (modulus - 1).bit_length() if width is None else width
    positive_integer(width, "modular_multiply.width", maximum=max_width)
    if modulus > (1 << width) or math.gcd(multiplier, modulus) != 1:
        raise ValidationError("the modulus must fit in the register and the multiplier must be coprime with the modulus")
    multiplier %= modulus
    permutation = [multiplier * x % modulus if x < modulus else x for x in range(1 << width)]
    b = Builder(
        _name("modular_multiply", multiplier, modulus, width),
        {"target": Bits(width)},
        attributes={
            "algorithm": "modular_multiply",
            "implementation_scope": "bounded permutation synthesis",
            "modulus": modulus,
        },
    )
    visited = set()
    for start in range(1 << width):
        if start in visited:
            continue
        cycle, current = [], start
        while current not in visited:
            visited.add(current)
            cycle.append(current)
            current = permutation[current]
        for other in cycle[1:]:
            _transposition(b, b["target"], start, other)
    return b.finish()


def order_finding(
    multiplier: int, modulus: int, *, precision: int = 3, max_width: int = 8
) -> Operation:
    """Run quantum order finding on the modular multiplication unitary starting from the integer one.

    Args:
        multiplier: Multiplier coprime with modulus.
        modulus: The modulus.
        precision: Phase bit width of the QPE.
        max_width: Bit width budget of the underlying modular multiplication
            permutation.

    Returns:
        Operation: target/phase interface. phase carries fractional
        information about the order and needs classical post-processing."""
    operation = modular_multiply(multiplier, modulus, max_width=max_width)
    n = operation.module.registers[0].type.width
    qpe = phase_estimation(operation, precision=precision)
    b = Builder(
        _name("order_finding", multiplier, modulus, precision),
        {"target": Bits(n), "phase": Bits(precision)},
        resources_for(("qpe", qpe)),
        attributes={
            "algorithm": "order_finding",
            "modulus": modulus,
            "multiplier": multiplier,
            "readout_register": "phase",
        },
    )
    b.x(b["target"][0])
    invoke(b, qpe, "qpe", target=b["target"], phase=b["phase"])
    return b.finish()


def factors_from_phase(
    value: int, precision: int, multiplier: int, modulus: int
) -> tuple[int, int] | None:
    """Attempt to obtain a nontrivial factor from the continued-fraction candidate order of a phase sample.

    Args:
        value: Integer readout of the phase register.
        precision: Bit width of the phase register.
        multiplier: The multiplier used in quantum order finding.
        modulus: The integer to process.

    Returns:
        tuple or None: A verified factor pair, or None when this sample
        yields no factor."""
    positive_integer(precision, "factors.precision", maximum=63)
    positive_integer(value, "factors.phase", minimum=0, maximum=(1 << precision) - 1)
    positive_integer(modulus, "factors.modulus", minimum=3)
    positive_integer(multiplier, "factors.multiplier")
    common = math.gcd(multiplier, modulus)
    if 1 < common < modulus:
        return cast("tuple[int, int]", tuple(sorted((common, modulus // common))))
    if value == 0:
        return None
    order = Fraction(value, 1 << precision).limit_denominator(modulus).denominator
    if order % 2 or pow(multiplier, order, modulus) != 1:
        return None
    half = pow(multiplier, order // 2, modulus)
    for candidate in (math.gcd(half - 1, modulus), math.gcd(half + 1, modulus)):
        if 1 < candidate < modulus:
            return cast("tuple[int, int]", tuple(sorted((candidate, modulus // candidate))))
    return None
