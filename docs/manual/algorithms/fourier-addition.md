# Fourier Addition

**English** · <a href="../../zh/manual/algorithms/fourier-addition.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.common.fourier`](../../api/algorithms/common/fourier.rst) · Stage V1

## Overview

A mod-$2^n$ adder implemented in the [QFT](qft.md) domain: $|a, b\rangle \mapsto |a,\ (a+b) \bmod 2^n\rangle$. The implementation follows Draper's QFT addition construction ([arXiv:quant-ph/0008033](https://arxiv.org/abs/quant-ph/0008033); see the "Implementation basis" section of the algorithm catalog). Integer addition degenerates into single-qubit phase rotations in the Fourier domain: after the QFT of $b$, the phase carried by each Fourier-basis component is linear in the addend, so it suffices to apply fixed-angle phase gates to the Fourier bits of $b$ controlled by the bits of $a$, and then transform back to the computational basis — no carry register is needed anywhere.

## Interface and input model

```python
fourier_add(width)
```

API entry point: {obj}`fourier_add <oracq.algorithms.common.fourier.fourier_add>`

- `width`: the common bit width of `a` and `b`, range 1..64; non-positive integers or out-of-range values raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time.

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `a: Bits(width)` and `b: Bits(width)`: `a` is preserved unchanged as the control, and `b` is updated to $(a+b) \bmod 2^{\text{width}}$, with no additional public workspace. This entry point consumes no oracle input (no input model).

## Implementation notes

The generation sequence is: call {obj}`qft(width) <oracq.algorithms.common.fourier.qft>` on `b` → a double loop adds controlled phases: a phase gate on `b[j]` controlled by `a[i]` with angle $2\pi / 2^{\text{width} - i - j}$ (for $0 \le i, j$ with $i + j < \text{width}$) → call {obj}`inverse_qft(width) <oracq.algorithms.common.fourier.inverse_qft>` on `b`. The QFT and inverse QFT are both kept as module calls and not expanded during generation or serialization. The angle index `width - i - j` comes from the little-endian bit weighting: bit $i$ of `a` contributes $a_i 2^i$, and the $j$-th Fourier bit of `b` needs a phase rotation by an integer multiple of $2\pi a_i 2^i / 2^{\text{width} - j}$.

Applicability boundary: the result is addition with wraparound modulo $2^{\text{width}}$ and no overflow flag is provided; when signed arithmetic or overflow detection is needed, use the `add` of [reversible fixed-point arithmetic](fixed-point-arithmetic.md) (with `status` flags).

## Validation approach

Category C1 (exact discrete semantics; acceptance criteria in `../development/validation-plan.md` §2): the acting unitary must equal the truth table pointwise. Three layers of evidence:

- Structure: `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` (construction and property assertions in the same group as QFT, matrix criterion).
- Numerical: `AlgorithmExpansionTests.test_fourier_add_all_basis_inputs` — {obj}`fourier_add(3) <oracq.algorithms.common.fourier.fourier_add>` is simulated over all $8 \times 8 = 64$ pairs of basis-state inputs $(a, b)$; the output amplitude is concentrated on $(a,\ (a+b) \bmod 8)$ and equals 1, to precision places = 10, i.e. an exhaustive cross-check over all inputs of a small instance.
- Binding: this algorithm has no separate binding witness (matrix criterion is —; no open declaration entry point).

## Known gaps and planned stages

No known gaps (the gap column of the validation matrix is —); stage V1.

## Related links

- Same module: [Quantum Fourier Transform](qft.md)
- Source: `src/oracq/algorithms/common/fourier.py`
- API reference: [Fourier transforms and arithmetic](../../api/algorithms/common/fourier.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_fourier.py` (real-backend execution, no mock substitutes). The truth-table oracle is an independent classical integer modular addition: at the unitary level, for n ≤ 4, OriginIR-ext export through UniQC `Circuit.to_matrix` obtains the full unitary, compared element by element against the permutation $|a,b\rangle \mapsto |a,(a+b)\bmod 2^n\rangle$ (`effective_block` extraction, leakage 0); at the state-vector level, for n ≤ 4, the addend $b_0$ is fixed and $a$ is superposed, iterating over all $b_0$ to cover the full $4^n$ input pairs (dual reference and originir-ext paths); at n = 8, rir-pysparq exhausts all 256 branches of $a$ in one superposition (peak intermediate sparse state $2^{16}$, several representative $b_0$); at n = 12, because the intermediate state $2^{24}$ exceeds the sparse budget, deterministic output checks on sampled basis-state pairs are used instead (noted in the case parameters).

| Case | Scale | Backend path | Metric | Value |
|---|---|---|---|---|
| `fourier-add-unitary` | n = 1/2/3/4 | originir-ext + UniQC to_matrix | max_error | 2.3e-16 / 3.8e-16 / 6.5e-16 / 1.0e-15 |
| `fourier-add-truthtable` | n = 1/2/3/4 (all $4^n$ input pairs) | reference + originir-ext | max_error | 1.2e-16 / 2.4e-16 / 2.6e-16 / 3.0e-16 |
| `fourier-add-superposed-a-n8` | n = 8, 5 groups of $b_0$ × 256 branches | rir-pysparq | max_error | 1.5e-16 |
| `fourier-add-sampled-basis-n12` | n = 12, 12 sampled groups of $(a,b)$ | rir-pysparq | failures | 0 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_fourier.py
```

Output: `out/verification/fourier.json`.
