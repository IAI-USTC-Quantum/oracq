# Quantum Fourier Transform

**English** · <a href="../../zh/manual/algorithms/qft.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.common.fourier`](../../api/algorithms/common/fourier.rst) · Stage V1

## Overview

An $n$-bit positive-sign discrete Fourier transform with matrix elements

$$
F_{y x} = \frac{1}{\sqrt{2^n}} \exp\!\bigl(2\pi i \cdot x y / 2^n\bigr),
$$

that is, $|x\rangle \mapsto \frac{1}{\sqrt{2^n}} \sum_y \exp(2\pi i \cdot x y / 2^n) |y\rangle$. The implementation uses the standard Hadamard + controlled-phase structure with a trailing bit-order swap, with a little-endian bit-order convention. The QFT is a foundational component of modules such as [Fourier addition](fourier-addition.md) and phase estimation.

## Interface and input model

```python
qft(width)
qft_with_work(width)
inverse_qft(width)
```

API entry points: {obj}`qft <oracq.algorithms.common.fourier.qft>`, {obj}`qft_with_work <oracq.algorithms.common.fourier.qft_with_work>`, {obj}`inverse_qft <oracq.algorithms.common.fourier.inverse_qft>`

- {obj}`qft(width) <oracq.algorithms.common.fourier.qft>`: `width` ranges over 1..64; non-positive integers or out-of-range values raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time. Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with only a `target: Bits(width)` register.
- {obj}`qft_with_work(width) <oracq.algorithms.common.fourier.qft_with_work>`: an adapter preserving the earlier zero-width `work` interface; the registers are `target: Bits(width)` and `work: Bits(0)`, and it calls `qft` directly internally.
- {obj}`inverse_qft(width) <oracq.algorithms.common.fourier.inverse_qft>`: the adjoint (inverse transform) of the QFT, implemented by wrapping the module call in an adjoint context.

None of the three entry points consumes oracle input (no input model); the bit width is the entire parameter set.

## Implementation notes

The generation order of `qft`: iterate `high` from the most significant bit toward the least significant bit; first apply H to `target[high]`, then for each `low < high` apply a phase gate controlled by `target[low]` with angle $\pi / 2^{\text{high} - \text{low}}$; then $\lfloor n/2 \rfloor$ swaps exchange the mirrored bits, flipping the bit order at the end of the transform to little endian. `inverse_qft` wraps the call to `qft` in `with b.adjoint()`; the transform itself remains a module call and is not expanded during generation or serialization, and the inverse shares the same submodule definition as the forward transform.

Applicability boundary: the upper limit of 64 on `width` comes from the register bit-width contract; the output bit order is little endian, and when interfacing with an external bit order note that the trailing swaps have already been applied.

## Validation approach

Category C1 (exact discrete semantics; acceptance criteria in `../development/validation-plan.md` §2): the acting unitary must equal the Fourier matrix pointwise. Three layers of evidence:

- Structure: `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects` pins the historical import name `oracq.algorithms.elementary.qft` back to the canonical `qft` object of this module (matrix criterion).
- Numerical: `AlgorithmExpansionTests.test_qft_matches_positive_fourier_matrix` — `qft(3)` is simulated column by column over all 8 basis states, and every output amplitude is cross-checked against the positive-sign Fourier matrix element $\exp(2\pi i \cdot x y / 8)/\sqrt{8}$ to precision places = 10.
- Binding: this algorithm has no separate binding witness (matrix criterion is —; no open declaration entry point).

## Known gaps and planned stages

No known gaps (the gap column of the validation matrix is —); stage V1.

## Related links

- Same module: [Fourier addition](fourier-addition.md)
- Source: `src/oracq/algorithms/common/fourier.py`
- API reference: [Fourier transforms and arithmetic](../../api/algorithms/common/fourier.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_fourier.py` (real-backend execution, no mock substitutes). Experimental design: at the unitary level, OriginIR-ext export through UniQC `Circuit.to_matrix` obtains the full unitary (n ≤ 4; no workspace, so `effective_block` is the full matrix and leakage is 0), compared element by element against an independently constructed positive-sign DFT matrix $F_{yx} = \exp(2\pi i x y/2^n)/\sqrt{2^n}$; at the register level, rir-pysparq runs two families of pointwise cross-checks at n = 6/8/12 — all $2^n$ amplitudes of the basis state $|x\rangle$ under QFT compared against DFT rows, and a Fourier mode $\sum_x \omega^{-kx}|x\rangle/\sqrt{2^n}$ prepared independently using only H and single-qubit phase gates, checked after QFT for focusing on basis-state index $k$; the n = 4 phase-rich state additionally gets a four-path cross-check across reference / rir-pysparq / adapter-pysparq / originir-ext (including `qft_with_work` adapter equivalence); inverse-QFT round-trip identity is validated at both the unitary level (n = 3/4 composite unitary = $I$) and the register level (n = 8/12 full superposition recovers the uniform state, sampled basis-state round trips).

| Case | Scale | Backend path | Metric | Value |
|---|---|---|---|---|
| `qft-unitary-dft` | n = 1/2/3/4 | originir-ext + UniQC to_matrix | max_error | 8.7e-17 / 2.7e-16 / 1.4e-15 / 3.8e-15 |
| `inverse-qft-unitary` | n = 1/2/3/4 | originir-ext + UniQC to_matrix | max_error | 8.7e-17 / 2.7e-16 / 1.4e-15 / 3.8e-15 |
| `inverse-qft-roundtrip-unitary` | n = 3/4 | originir-ext + UniQC to_matrix | max_error | 4.4e-16 / 7.8e-16 |
| `qft-basis-row-pointwise` | n = 6/8/12 | rir-pysparq | max_error | 4.2e-15 / 9.5e-15 / 5.4e-14 |
| `qft-fourier-mode-focus` | n = 6/8/12 | rir-pysparq | min_success_probability | ≥ 1 − 4.4e-15 (leakage amplitude 0) |
| `qft-cross-path` | n = 4 | four-path cross-check | max_pairwise_deviation | 0.0 |
| `inverse-qft-roundtrip-wide` | n = 8/12 | rir-pysparq | uniform_max_error | 1.2e-16 / 4.5e-17 (basis-state round-trip failures 0) |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_fourier.py
```

Output: `out/verification/fourier.json` (39 cases, all passing; total run about 80 seconds).
