"Modular BEs, projectors, LCU, and small-matrix representation composition."

from __future__ import annotations

import cmath
import itertools
import math
from collections.abc import Iterable, Sequence

from oracq.algorithms.input_model.operators import BlockEncoding, _name, identity, scale
from oracq.algorithms.input_model.oracles import (
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, Ref, ValidationError


def reflect_zero(builder: Builder, register: Ref, *, positive: bool = False) -> None:
    """Append the sign-flip reflection ``I-2|0><0|`` about the zero state of the given register.

    Args:
        builder: Target Builder to append gates to.
        register: Register participating in the reflection; a fused view of several registers is accepted.
        positive: When True an extra global phase pi is applied, giving ``2|0><0|-I`` with the zero-state component positive.

    Shared by walk-type and amplitude-amplification-type algorithms."""
    if positive:
        builder.global_phase(math.pi)
    if register.width:
        with builder.control(register, 0):
            builder.global_phase(math.pi)
    else:
        builder.global_phase(math.pi)


def pad_signal(a: BlockEncoding, width: int) -> BlockEncoding:
    """Widen the signal register of a BE to the given width.

    Args:
        a: Input block encoding.
        width: New signal bit width; must not be smaller than the original signal width.

    Returns:
        BlockEncoding: Corner-block semantics and alpha unchanged, with the new high signal bits always zero, enabling late binding under the same signature.

    Raises:
        ValidationError: Attempted to shrink the signal space.
    """
    if width < a.signal_qubits:
        raise ValidationError("cannot shrink the BE signal space")
    b = Builder(
        _name("pad_be", a.operation, width),
        {"target": Bits(a.width), "signal": Bits(width)},
        resources_for(("a", a.operation)),
    )
    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"][: a.signal_qubits])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=a.alpha))


def tensor(a: BlockEncoding, b: BlockEncoding) -> BlockEncoding:
    """Build the tensor product of two BEs.

    Args:
        a: BE acting on the high target bits.
        b: BE acting on the low target bits.

    Returns:
        BlockEncoding: Encodes A⊗B with alpha the product of the two; signal is likewise concatenated with a low and b high.
    """
    out = Builder(
        _name("tensor", a.operation, b.operation),
        {"target": Bits(a.width + b.width), "signal": Bits(a.signal_qubits + b.signal_qubits)},
        resources_for(("a", a.operation), ("b", b.operation)),
    )
    invoke(
        out,
        b.operation,
        "b",
        target=out["target"][: b.width],
        signal=out["signal"][a.signal_qubits :],
    )
    invoke(
        out,
        a.operation,
        "a",
        target=out["target"][b.width :],
        signal=out["signal"][: a.signal_qubits],
    )
    return BlockEncoding(annotate(out.finish(), "block_encoding", be_alpha=a.alpha * b.alpha))


def adjoint_be(a: BlockEncoding) -> BlockEncoding:
    """Return the BE encoding the adjoint matrix A†.

    Args:
        a: Input block encoding.

    Returns:
        BlockEncoding: Invokes the original operation inside an adjoint context, with alpha unchanged.
    """
    out = Builder(
        _name("adjoint_be", a.operation),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        resources_for(("a", a.operation)),
    )
    with out.adjoint():
        invoke(out, a.operation, "a", target=out["target"], signal=out["signal"])
    return BlockEncoding(annotate(out.finish(), "block_encoding", be_alpha=a.alpha))


def lcu(terms: Iterable[tuple[complex, BlockEncoding]]) -> BlockEncoding:
    """Assemble a linear combination of BEs in PREPARE/SELECT structure.

    Args:
        terms: Sequence of (coefficient, BE) pairs; zero-coefficient terms are dropped.

    Returns:
        BlockEncoding: Encodes Σ c_j A_j with alpha ``Σ |c_j|*alpha_j``; degenerates to ``scale`` for a single term.

    Raises:
        ValidationError: There is no nonzero term, or the terms have inconsistent target widths.

    The phases of complex coefficients are realized by selector-controlled global
    phases; signal is the concatenation of the select bits and each branch's
    signal bits, with a trailing inverse preparation restoring the selector."""
    from oracq.algorithms.input_model.interfaces import as_block_encoding

    terms = tuple((complex(c), as_block_encoding(a)) for c, a in terms if c != 0)
    if not terms:
        raise ValidationError("LCU requires at least one nonzero term")
    if len(terms) == 1:
        return scale(terms[0][0], terms[0][1])
    width = terms[0][1].width
    if any(a.width != width for _, a in terms):
        raise ValidationError("LCU target widths do not match")
    alpha = sum(abs(c) * a.alpha for c, a in terms)
    selector_width = (len(terms) - 1).bit_length()
    work_width = max(a.signal_qubits for _, a in terms)
    weights = [math.sqrt(abs(c) * a.alpha / alpha) for c, a in terms]
    weights += [0] * ((1 << selector_width) - len(weights))
    prep = gate_state_prep(weights)
    resources = resources_for(*[(f"term{i}", a.operation) for i, (_, a) in enumerate(terms)])
    b = Builder(
        _name("lcu_many", *(a.operation for _, a in terms), tuple(c for c, _ in terms)),
        {"target": Bits(width), "signal": Bits(selector_width + work_width)},
        resources,
    )
    selector, signal = b["signal"][:selector_width], b["signal"][selector_width:]
    invoke(b, prep.operation, target=selector, work=selector[:0])
    for i, (coefficient, a) in enumerate(terms):
        with b.control(selector, i):
            b.global_phase(cmath.phase(coefficient))
            invoke(b, a.operation, f"term{i}", target=b["target"], signal=signal[: a.signal_qubits])
    with b.adjoint():
        invoke(b, prep.operation, target=selector, work=selector[:0])
    return BlockEncoding(
        annotate(b.finish(), "block_encoding", be_alpha=alpha, lcu_terms=len(terms))
    )


def kronecker_sum(a: BlockEncoding, b: BlockEncoding | None = None) -> BlockEncoding:
    """Build the Kronecker sum A⊗I+I⊗B of two BEs.

    Args:
        a: First BE.
        b: Second BE; defaults to a itself when omitted.

    Returns:
        BlockEncoding: LCU of the two tensor terms, with alpha the sum of the two alphas.
    """
    b = a if b is None else b
    return lcu([(1, tensor(a, identity(b.width))), (1, tensor(identity(a.width), b))])


def projector(width: int, accepted: Iterable[int]) -> BlockEncoding:
    """Build the projector BE onto the specified basis-state subspace.

    Args:
        width: Target bit width.
        accepted: Set of accepted basis-state integer values; duplicates are merged and sorted.

    Returns:
        BlockEncoding: The zero-signal corner block is the diagonal projector; basis states outside the set are kicked into the signal branch, with alpha 1.
    """
    accepted = tuple(sorted(set(accepted)))
    out = Builder(_name("projector", width, accepted), {"target": Bits(width), "signal": Bits(1)})
    for value in range(1 << width):
        if value not in accepted:
            with out.control(out["target"], value):
                out.x(out["signal"])
    return BlockEncoding(annotate(out.finish(), "block_encoding", be_alpha=1.0))


def direct_sum(a: BlockEncoding, b: BlockEncoding) -> BlockEncoding:
    """Build the direct sum of same-width matrices.

    Args:
        a: BE in effect when the select bit is 0.
        b: BE in effect when the select bit is 1; the width must equal a's.

    Returns:
        BlockEncoding: Block-diagonal composition using the target's most significant bit as the select bit, with alpha the sum of the two.

    Raises:
        ValidationError: The two BE widths differ.
    """
    if a.width != b.width:
        raise ValidationError("direct_sum currently requires matrices of equal width")
    return lcu([(1, tensor(projector(1, [0]), a)), (1, tensor(projector(1, [1]), b))])


def truncated_shift(width: int, last: int) -> BlockEncoding:
    """Build the truncated up-shift BE.

    Args:
        width: Target bit width.
        last: Truncation threshold, in the range 1..2**width-1.

    Returns:
        BlockEncoding: The zero-signal corner block maps basis state v to v+1 for v below last; other basis states are kicked into the signal branch, with alpha 1.

    Raises:
        ValidationError: last is outside the allowed range.
    """
    if not 1 <= last < 1 << width:
        raise ValidationError("truncated shift range is invalid")
    out = Builder(_name("shift", width, last), {"target": Bits(width), "signal": Bits(1)})
    for value in range(last, 1 << width):
        with out.control(out["target"], value):
            out.x(out["signal"])
    out.add_const(out["target"].reinterpret("uint"), 1)
    return BlockEncoding(annotate(out.finish(), "block_encoding", be_alpha=1.0))


def pauli_word(word: str) -> BlockEncoding:
    """Encode a Pauli string as a BE without signal bits.

    Args:
        word: String over I/X/Y/Z whose first character acts on the least significant bit.

    Returns:
        BlockEncoding: Explicit per-bit single-qubit gate sequence, with alpha 1.

    Raises:
        ValidationError: A character outside I/X/Y/Z appears.
    """
    b = Builder("pauli_" + word, {"target": Bits(len(word)), "signal": Bits(0)})
    for bit, letter in enumerate(word):
        if letter not in "IXYZ":
            raise ValidationError("Pauli words allow only I/X/Y/Z")
        if letter != "I":
            b.gate(letter.lower(), b["target"][bit])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))


def matrix_pauli_encoding(
    matrix: Sequence[Sequence[complex]], *, drop_tolerance: float = 1e-12
) -> BlockEncoding:
    """Explicit gate implementation for small applications; claims no quantum speedup for matrix input or classical expansion.

    Args:
        matrix: Square complex matrix of power-of-two dimension, at most 32 (i.e. 5 qubits).
        drop_tolerance: Pauli coefficient terms with magnitude not exceeding this tolerance are dropped.

    Returns:
        BlockEncoding: LCU block encoding of the Pauli expansion; the zero operator when all coefficients are dropped.
    """
    matrix = tuple(tuple(complex(v) for v in row) for row in matrix)
    d = len(matrix)
    if d < 2 or d & (d - 1) or any(len(row) != d for row in matrix):
        raise ValidationError("matrix must be a square matrix of power-of-two dimension")
    n = (d - 1).bit_length()
    if n > 5:
        raise ValidationError("explicit Pauli expansion is only for small instances of at most 5 bits; use an access oracle for larger instances")
    terms: list[tuple[complex, BlockEncoding]] = []
    for letters in itertools.product("IXYZ", repeat=n):
        coefficient = 0j
        for column in range(d):
            row, phase = column, 1 + 0j
            for bit, letter in enumerate(letters):
                value = (column >> bit) & 1
                if letter in "XY":
                    row ^= 1 << bit
                if letter == "Y":
                    phase *= -1j if value else 1j
                elif letter == "Z":
                    phase *= -1 if value else 1
            coefficient += phase.conjugate() * matrix[row][column]
        coefficient /= d
        if abs(coefficient) > drop_tolerance:
            terms.append((coefficient, pauli_word("".join(letters))))
    if not terms:
        from oracq.algorithms.input_model.operators import zero

        return zero(n)
    return lcu(terms)
