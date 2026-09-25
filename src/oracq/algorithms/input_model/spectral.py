"""Spectral input/output primitives: Fourier-basis diagonal block encoding, sparse spectral block encoding, and spectral state preparation.

Corresponds to Lemma C.8 (SS-BE) of arXiv:2509.08807 and the input side of
hierarchical spectral encoding. The diagonal matrix of Fourier basis functions
F_k(J)=exp(2*pi*1j*k*J/2**w) factors over the binary digits of J into a tensor
product of single-bit phase gates (paper Eq. C58); the controlled-product form
yields a spectral block encoding with zero Toffolis (Eq. C60). Spectral state
preparation sparsely prepares the frequency-domain amplitudes on the same
lattice register, then a QFT performs the Hilbert-space amplification from the
S-dimensional spectral space to the N-dimensional lattice space.

Convention: f(J) = sum_k c_k * exp(2*pi*1j*k*J/2**width) with integer
frequencies k (negative frequencies are understood mod 2**width); the alpha of
the spectral block encoding is sum_k abs(c_k).
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Mapping, Sequence

from oracq.algorithms.common.fourier import qft
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    _state_angles,
    annotate,
    gate_state_prep,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Ref, ValidationError


def normalized_spectrum(spectrum: Mapping[int, complex]) -> tuple[int, tuple[complex, ...]]:
    """Validate spectral coefficients and normalize them into a contiguous frequency band.

    The input ``{k: c_k}`` represents f(J) = sum_k c_k exp(2*pi*1j*k*J/2**width).
    The controlled-product decomposition requires a contiguous band; returns
    ``(k_min, coefficients)`` where coefficients covers every integer frequency
    in ``[k_min, k_max]`` (all-zero entries inside the band are allowed and
    counted as zero coefficients).

    Args:
        spectrum: Non-empty mapping from frequencies to complex coefficients; it must form a contiguous band with not all coefficients zero.

    Returns:
        tuple[int, tuple[complex, ...]]: ``(k_min, coefficients)``, i.e. the band
        lower limit and the tuple arranged in ascending frequency order,
        including zero coefficients inside the band.
    """
    if not isinstance(spectrum, dict) or not spectrum:
        raise ValidationError("spectral coefficients must be a non-empty dict mapping frequency to complex coefficient")
    ks = sorted(spectrum)
    if any(type(k) is not int for k in ks):
        raise ValidationError("spectral frequencies must be integers")
    if ks[-1] - ks[0] + 1 != len(spectrum):
        raise ValidationError("spectral frequencies must form a contiguous band; split a discontinuous band and combine via LCU")
    k_min = ks[0]
    coefficients = tuple(spectrum[k] for k in range(ks[0], ks[-1] + 1))
    if all(c == 0 for c in coefficients):
        raise ValidationError("spectral coefficients cannot all be zero")
    for c in coefficients:
        if not (math.isfinite(c.real) and math.isfinite(c.imag)):
            raise ValidationError("spectral coefficients must be finite")
    return k_min, coefficients


def _spectrum_register_width(coefficients: Sequence[complex]) -> int:
    """Spectral register width ``ceil(log2 S)`` required for ``S`` band entries, at least one bit."""
    return max(1, (len(coefficients) - 1).bit_length())


def frequency_amplitudes(
    coefficients: Sequence[complex], k_min: int, width: int, *, normalize: bool = True
) -> list[complex]:
    """Convert spectral coefficients into a 2**width-dimensional amplitude vector indexed by k mod 2**width.

    Args:
        coefficients: Complex coefficients within the band, arranged in ascending frequency order.
        k_min: Lower frequency limit of the band.
        width: Lattice register bit width, determining the dimension 2**width of the amplitude vector.
        normalize: When True, divide by the Euclidean norm of the coefficients.

    Returns:
        list[complex]: Frequency-domain amplitude vector of length 2**width.
    """
    scale = math.sqrt(sum(abs(c) ** 2 for c in coefficients)) if normalize else 1.0
    if scale == 0:
        raise ValidationError("spectral coefficient norm is zero")
    amplitudes = [0j] * (1 << width)
    for j, c in enumerate(coefficients):
        if c:
            amplitudes[(k_min + j) % (1 << width)] = c / scale
    return amplitudes


def pruned_state_prep(
    amplitudes: Sequence[complex], *, name: str | None = None
) -> StatePreparation:
    """Multiplexed-rotation state preparation with zero-angle pruning; unitarily equivalent to gate_state_prep but skipping zero rotations.

    For sparse spectral input (a 2**width-dimensional vector with only S nonzero
    amplitudes), every rotation angle in an all-zero subtree is zero, so after
    pruning the resources grow with S instead of 2**width.

    Args:
        amplitudes: Sparse complex amplitude vector; the length must be a power of two.
        name: Name of the generated operation; defaults to one derived from the amplitude content.

    Returns:
        StatePreparation: Preparation view implemented by pruned multiplexed rotations.
    """
    values, n, nodes = _state_angles(amplitudes)
    if n == 0:
        raise ValidationError("state preparation requires at least one target bit")
    b = Builder(
        name or _name("pruned_state", values),
        {"target": Bits(n), "work": Bits(0)},
    )
    for depth, prefix, bit, angle in nodes:
        if angle == 0:
            continue
        if not depth:
            b.ry(b["target"][bit], angle)
        else:
            with b.control(b["target"][bit + 1 :], prefix):
                b.ry(b["target"][bit], angle)
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
            implementation="multiplexed_rotations_pruned",
        )
    )


def fourier_phase(width: int, k: int) -> Operation:
    """Exact circuit of diag_J exp(2*pi*1j*k*J/2**width) (paper Eq. C58).

    e^{2*pi*1j*k*J/2**width} factors over the binary digits of J into a tensor
    product with one single-bit phase gate per bit; the angles are binary-power
    multiples of pi and land on the exact-angle grid. Contains only the target
    register and is invoked as a submodule.

    Args:
        width: Target register bit width, in the range 1..64.
        k: Integer frequency, understood ``mod 2**width``.

    Returns:
        Operation: Diagonal phase operation built from a tensor product of single-bit phase gates.
    """
    if type(width) is not int or not 1 <= width <= 64:
        raise ValidationError("fourier_phase.width must be within 1..64")
    if type(k) is not int:
        raise ValidationError("fourier_phase.k must be an integer")
    k %= 1 << width
    b = Builder(f"fourier_phase_{width}_{k}", {"target": Bits(width)})
    for q in range(width):
        angle = 2 * math.pi * k * (1 << q) / (1 << width)
        if angle % (2 * math.pi):
            b.gate("phase", b["target"][q], angle)
    return b.finish()


def fourier_phase_encoding(width: int, k: int) -> BlockEncoding:
    """(1, 0, 0) block encoding of fourier_phase: zero-width signal interface.

    Args:
        width: Target register bit width, in the range 1..64.
        k: Integer frequency, understood ``mod 2**width``.

    Returns:
        BlockEncoding: Block encoding with normalization constant 1 and a zero-width signal.
    """
    if type(width) is not int or not 1 <= width <= 64:
        raise ValidationError("fourier_phase_encoding.width must be within 1..64")
    if type(k) is not int:
        raise ValidationError("fourier_phase_encoding.k must be an integer")
    k %= 1 << width
    b = Builder(
        _name("fourier_phase_be", width, k),
        {"target": Bits(width), "signal": Bits(0)},
        attributes={"matrix_interpretation": f"diag(exp(2*pi*1j*{k}*J/2^{width}))"},
    )
    for q in range(width):
        angle = 2 * math.pi * k * (1 << q) / (1 << width)
        if angle % (2 * math.pi):
            b.gate("phase", b["target"][q], angle)
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))


def _phase_grid(
    builder: Builder,
    control_builder: Ref,
    k_min: int,
    coefficients: Sequence[complex],
    width: int,
    spectrum_width: int,
) -> None:
    """Spectral-register-controlled phase grid (the controlled-product part of paper Eq. C60).

    For every bit of the spectral register and every qubit of the target
    register a controlled phase 2*pi*2**bit*2**q/2**width is applied; the
    contribution of the band lower limit k_min is an unconditional phase. Terms
    whose angle is an integer multiple of 2*pi are skipped (identity).
    """
    target, spectrum = builder["target"], control_builder
    for q in range(width):
        if k_min:
            angle = 2 * math.pi * k_min * (1 << q) / (1 << width)
            if angle % (2 * math.pi):
                builder.gate("phase", target[q], angle)
    for bit in range(spectrum_width):
        for q in range(width):
            exponent = (1 << bit) * (1 << q)
            if exponent % (1 << width):
                angle = 2 * math.pi * exponent / (1 << width)
                with builder.control(spectrum[bit], 1):
                    builder.gate("phase", target[q], angle)


def spectral_diagonal(
    spectrum: Mapping[int, complex], width: int, *, variant: str = "sequential"
) -> BlockEncoding:
    """Sparse spectral block encoding: BE(diag f), with f given by a contiguous-band Fourier sum.

    variant="sequential" (paper Lemma C.8 / Eq. (7)): P_L prepares amplitudes
    sqrt(abs(c_j)/alpha) on an s=log S bit spectral register, each bit of the
    spectral register controls the phase grid on the target register, the
    per-mode phases are restored by controlled gphase, and P_R=P_L^dagger
    uncomputes. Zero Toffolis. variant="naive": every basis function is handed
    to the generic LCU as an independent unitary operator, serving as the
    baseline against which the gain of the structured construction (sequential
    form) is measured.

    Args:
        spectrum: Contiguous-band mapping from frequencies to complex coefficients.
        width: Lattice register bit width.
        variant: ``sequential`` or ``naive``, taking the controlled-product decomposition or the generic LCU respectively.

    Returns:
        BlockEncoding: ``BE(diag f)`` with alpha = sum_k abs(c_k).
    """
    if variant not in {"sequential", "naive"}:
        raise ValidationError("variant must be sequential or naive")
    k_min, coefficients = normalized_spectrum(spectrum)
    if variant == "naive":
        return lcu(
            [
                (c, fourier_phase_encoding(width, (k_min + j) % (1 << width)))
                for j, c in enumerate(coefficients)
                if c
            ]
        )
    s = _spectrum_register_width(coefficients)
    alpha = sum(abs(c) for c in coefficients)
    amplitudes: list[complex] = []
    for c in coefficients:
        amplitudes.append(math.sqrt(abs(c) / alpha) * (c / abs(c) if c else 0))
    amplitudes += [0.0] * ((1 << s) - len(amplitudes))
    prep = gate_state_prep(amplitudes)
    b = Builder(
        _name("spectral_diagonal", coefficients, k_min, width),
        {"target": Bits(width), "signal": Bits(s)},
        attributes={
            "spectral_sparsity": sum(1 for c in coefficients if c),
            "band_low": k_min,
            "implementation": "sequential_fourier_spectral",
        },
    )
    b.call(prep.operation, target=b["signal"], work=b["signal"][:0])
    _phase_grid(b, b["signal"], k_min, coefficients, width, s)
    for j, c in enumerate(coefficients):
        if c and cmath.phase(c):
            with b.control(b["signal"], j):
                b.global_phase(cmath.phase(c))
    with b.adjoint():
        b.call(prep.operation, target=b["signal"], work=b["signal"][:0])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=alpha))


def spectral_state_prep(spectrum: Mapping[int, complex], width: int) -> StatePreparation:
    """Spectral input (the input side of hierarchical spectral encoding): prepare the normalized state ``|f⟩``.

    The frequency-domain amplitudes sum_k c_k / norm(c) are sparsely prepared on
    basis states ``|k⟩`` of the same width-bit lattice register, then a
    positive-sign QFT amplifies into the lattice space: the QFT maps basis
    states to uniform states with per-bit phases, and the linear combination
    yields sum_J f(J)*ket(J)/sqrt(N*sum(abs(c)^2)), the normalized state
    (Parseval).

    Args:
        spectrum: Contiguous-band mapping from frequencies to complex coefficients; the band entry count must not exceed 2**width.
        width: Lattice register bit width.

    Returns:
        StatePreparation: Preparation view implemented by sparse frequency-domain preparation plus QFT amplification.
    """
    k_min, coefficients = normalized_spectrum(spectrum)
    if len(coefficients) > (1 << width):
        raise ValidationError("spectral band width exceeds the range representable by the lattice register")
    amplitudes = frequency_amplitudes(coefficients, k_min, width)
    prep = pruned_state_prep(amplitudes)
    b = Builder(
        _name("spectral_state_prep", coefficients, k_min, width),
        {"target": Bits(width), "work": Bits(0)},
        attributes={
            "implementation": "hierarchy_spectral_encoding",
            "spectral_sparsity": sum(1 for c in coefficients if c),
            "zero_input": True,
            "clean_work": True,
        },
    )
    b.call(prep.operation, target=b["target"], work=b["target"][:0])
    b.call(qft(width), target=b["target"])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry"))
