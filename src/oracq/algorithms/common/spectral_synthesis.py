"""Operator-level circuit synthesis optimizations (arXiv:2509.08807 appendix D).

This implements the paper's three logical-resource optimizations while
preserving the unitary semantics exactly:

1. Uniformly controlled rotations, UCR: rewrite the multi-controlled rotation
   tree of a state preparation into plain rotations plus a CNOT ladder,
   eliminating multi-controlled Toffoli gates;
2. Similar-group merging, match by signal / similar group: merge terms of a
   linear combination that differ only by a global sign or a scalar factor
   into one term, absorbing the coefficient into the preparation amplitudes;
3. Fan-out phase grid: rewrite the controlled-phase layer of an FBBE into
   single-qubit rotations plus fan-out CNOTs, eliminating controlled
   rotations, see Figure 9(a2) of the paper.

Each construction is compared elementwise against the naive construction in
tests/core/test_spectral_synthesis.py, and estimate_resources verifies the
resource reduction.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Iterable, Mapping, Sequence
from typing import cast

from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    _state_angles,
    annotate,
    gate_state_prep,
)
from oracq.algorithms.input_model.spectral import (
    _spectrum_register_width,
    normalized_spectrum,
)
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, Ref, ValidationError


def _emit_ucr(
    builder: Builder, target: Ref, controls: list[Ref], thetas: Sequence[float]
) -> None:
    """Recursive UCR: CX; UCR(diff); CX; UCR(mean), with thetas indexed by control
    prefix.

    Prefix bit order: controls[0] is the least significant bit. Equivalent to
    a multi-controlled Ry tree indexed by prefix value, the controlled layer
    of gate_state_prep, but using only single-qubit rotations and CNOTs.
    """
    if not controls:
        if thetas[0]:
            builder.ry(target, thetas[0])
        return
    control = controls[-1]
    half = len(thetas) // 2
    mean = [(thetas[j] + thetas[j + half]) / 2 for j in range(half)]
    diff = [(thetas[j] - thetas[j + half]) / 2 for j in range(half)]
    builder.xor(control, target)
    _emit_ucr(builder, target, controls[:-1], diff)
    builder.xor(control, target)
    _emit_ucr(builder, target, controls[:-1], mean)


def uniformly_controlled_prep(
    amplitudes: Iterable[complex], *, name: str | None = None
) -> StatePreparation:
    """UCR state preparation: unitarily equivalent to gate_state_prep, without
    multi-controlled rotations.

    Each layer's 2^d controlled Ry rotations are rewritten into 2^d
    single-qubit Ry rotations and 2^d CNOTs; the complex phases are still
    written by controlled gphase.

    Args:
        amplitudes: Complex amplitude sequence of the target state, with power-of-two
            length, not all zero.
        name: Name of the generated preparation module; generated automatically when
            omitted.

    Returns:
        StatePreparation: State preparation handle containing only single-qubit
        rotations and CNOTs.
    """
    values, n, nodes = _state_angles(amplitudes)
    if n == 0:
        raise ValidationError("State preparation requires at least one target qubit")
    b = Builder(
        name or _name("ucr_state", values),
        {"target": Bits(n), "work": Bits(0)},
    )
    by_level: dict[int, list[float]] = {}
    for depth, prefix, _bit, angle in nodes:
        by_level.setdefault(depth, [0.0] * (1 << depth))[prefix] = angle
    for depth in range(n):
        bit = n - depth - 1
        controls = [b["target"][q] for q in range(bit + 1, n)]
        _emit_ucr(b, b["target"][bit], controls, by_level[depth])
    for index, value in enumerate(values):
        if value and cmath.phase(value):
            with b.control(b["target"], index):
                b.global_phase(cmath.phase(value))
    return StatePreparation(
        annotate(
            b.finish(),
            "state_prep_isometry",
            zero_input=True,
            clean_work=True,
            implementation="uniformly_controlled_rotations",
        )
    )


def merge_similar(
    terms: Iterable[tuple[complex, BlockEncoding]],
) -> list[tuple[complex, BlockEncoding]]:
    """Merge linear-combination terms by block-encoding identity, match by signal /
    similar group.

    Terms with identical operations, including the case differing only by the
    global sign −1, merge into one term with summed coefficients; terms whose
    coefficient vanishes are dropped. This is an exact algebraic rewrite: the
    LCU α and the angle block are unchanged, while the SELECT branch count,
    the preparation amplitude tree and the controlled call count all decrease.

    Args:
        terms: Linear-combination term sequence of, coefficient and block encoding,
            pairs; zero-coefficient terms are dropped directly.

    Returns:
        list[tuple[complex, BlockEncoding]]: Term list merged by block-encoding identity
        with summed coefficients.
    """
    merged: dict[str, list[complex | BlockEncoding]] = {}
    order: list[str] = []
    for coefficient, operand in terms:
        if coefficient == 0:
            continue
        key = operand.operation.module.name
        if key not in merged:
            merged[key] = [0j, operand]
            order.append(key)
        merged[key][0] += complex(coefficient)  # type: ignore[operator]
    return [
        (cast("complex", merged[key][0]), cast("BlockEncoding", merged[key][1]))
        for key in order
        if merged[key][0] != 0
    ]


def fanout_spectral_diagonal(spectrum: Mapping[int, complex], width: int) -> BlockEncoding:
    """Fan-out form of an FBBE: the controlled-phase grid is rewritten as rotations
    plus fan-out CNOTs.

    Figure 9(a2) of the paper: each spectral-bit-controlled phase layer
    becomes one half-angle rotation per target bit, a fan-out CNOT, a negated
    half-angle rotation and another fan-out CNOT; CNOTs sharing the same
    control bit form a multi-target fan-out X, a native operation of surface
    code lattice surgery. In the RIR the fan-out appears as a same-source xor
    sequence. Elementwise equivalent to the sequential variant of
    spectral_diagonal.

    Args:
        spectrum: Mapping from spectral frequencies to complex coefficients; the
            frequencies are interpreted as integers mod 2^width.
        width: Target register bit width, determining the phase angle denominator
            2^width.

    Returns:
        BlockEncoding: Spectral diagonal block encoding in fan-out form, with scale
        equal to the sum of absolute coefficient values.
    """
    k_min, coefficients = normalized_spectrum(spectrum)
    s = _spectrum_register_width(coefficients)
    alpha = sum(abs(c) for c in coefficients)
    amplitudes = [
        math.sqrt(abs(c) / alpha) * (c / abs(c) if c else 0) for c in coefficients
    ]
    amplitudes += [0.0] * ((1 << s) - len(amplitudes))
    prep = gate_state_prep(amplitudes)
    b = Builder(
        _name("fanout_spectral_diagonal", coefficients, k_min, width),
        {"target": Bits(width), "signal": Bits(s)},
        attributes={
            "spectral_sparsity": sum(1 for c in coefficients if c),
            "band_low": k_min,
            "implementation": "fanout_fourier_spectral",
        },
    )
    b.call(prep.operation, target=b["signal"], work=b["signal"][:0])
    for q in range(width):
        if k_min:
            angle = 2 * math.pi * k_min * (1 << q) / (1 << width)
            if angle % (2 * math.pi):
                b.gate("phase", b["target"][q], angle)
    for bit in range(s):
        pairs: list[tuple[int, float]] = []
        for q in range(width):
            exponent = (1 << bit) * (1 << q)
            if exponent % (1 << width):
                pairs.append((q, 2 * math.pi * exponent / (1 << width)))
        if not pairs:
            continue
        control = b["signal"][bit]
        for q, angle in pairs:
            b.gate("phase", b["target"][q], angle / 2)
        for q, _ in pairs:
            b.xor(control, b["target"][q])
        for q, angle in pairs:
            b.gate("phase", b["target"][q], -angle / 2)
        for q, _ in pairs:
            b.xor(control, b["target"][q])
        # The sandwich decomposition leaves a residual e^{-i*theta/2} phase on the
        # control=1 branch, accumulating pair by pair; a single-qubit phase gate on
        # the control bit cancels it exactly.
        b.gate("phase", control, sum(angle for _, angle in pairs) / 2)
    for j, c in enumerate(coefficients):
        if c and cmath.phase(c):
            with b.control(b["signal"], j):
                b.global_phase(cmath.phase(c))
    with b.adjoint():
        b.call(prep.operation, target=b["signal"], work=b["signal"][:0])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=alpha))
