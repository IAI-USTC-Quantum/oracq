# Sparse Matrix Block Encoding

**English** · <a href="../../../zh/manual/algorithms/sparse-block-encoding.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.input_model.sparse`](../../api/algorithms/input_model/sparse.rst) · Stage V1

## Overview

Adapts the sparse access oracle (position and entry separated, input model SO) into a block encoding. The formal entry point {obj}`real_symmetric_sparse_encoding <oracq.algorithms.input_model.sparse.real_symmetric_sparse_encoding>` targets real symmetric (Hermitian) matrices with a nonnegative diagonal and constructs the CKS-type $T^\dagger S T$ (Childs–Kothari–Somma, arXiv:1511.02306; interface conventions in [QFVM input model review](../../reference/qfvm-input-models.md)): first the amplitude-transduction preparation $T$ — prepare the uniform prefix superposition of the sparse positions on `neighbor`, query the position oracle for the row index, query the entry oracle for the entry word, and rotate the success flag by $\sqrt{|v|/a_{\max}}$; then swap the coordinates on the two sides and the failure flags on the two sides to obtain $S$. If each column enumerates $s$ structural positions and $a_{\max}$ genuinely bounds the element magnitudes, the zero-signal block of $T^\dagger S T$ is $A/(s \cdot a_{\max})$.

## Interface and input model

```python
real_symmetric_sparse_encoding(access, fmt, amax, *, diagonal_nonnegative=False, rotation=None)
```

API entry points: {obj}`real_symmetric_sparse_encoding <oracq.algorithms.input_model.sparse.real_symmetric_sparse_encoding>`

- `access`: {obj}`SparseAccess <oracq.algorithms.input_model.oracles.SparseAccess>` — the two operations `location` (column/index/work in-place permutation) and `entry` (row/column/data XOR).
- `fmt`: the {obj}`FixedFormat <oracq.algorithms.common.arithmetic.FixedFormat>` fixed-point format of the entries; its width must equal `access.value_width`; `fmt.signed` decides whether the sign phase is superimposed.
- `amax`: the element magnitude bound $a_{\max}$, a finite positive number.
- `diagonal_nonnegative`: must be passed explicitly as `True`, otherwise it errors with "the current symmetric sparse adaptation requires a nonnegative diagonal".
- `rotation`: the amplitude-transduction operation, defaulting to {obj}`magnitude_rotation(fmt, amax) <oracq.algorithms.input_model.sparse.magnitude_rotation>`.

Returns a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` with module attributes:

| Attribute | Meaning |
|---|---|
| `be_alpha` | $s \cdot a_{\max}$ ($s$ is `access.sparsity`) |
| `construction` | `"CKS_Tdag_S_T"` |
| `self_adjoint_extension` | `True` (a precondition of the Chebyshev walk) |
| `correctness` | `"pending"` |

Auxiliary entry points: `magnitude_rotation(fmt, amax)` (an explicit gate implementation for value widths ≤ 12 bits; wider values return a declaration pending binding that carries an `amplitude_contract` attribute), {obj}`prefix_state(width, count) <oracq.algorithms.input_model.sparse.prefix_state>` (uniform prefix superposition), {obj}`compare_words(width, kind) <oracq.algorithms.input_model.sparse.compare_words>` (eq/lt Boolean network, cached and reused), and {obj}`chebyshev_block(a, degree) <oracq.algorithms.input_model.sparse.chebyshev_block>` (see below). The old entry point {obj}`sparse_block_encoding <oracq.algorithms.input_model.sparse.sparse_block_encoding>` is marked legacy (`legacy_input_model=True`, `matrix_contract` declared as "legacy transduction unspecified"); it serves only the old catalog descriptions and is not a general sparse input adapter.

## Implementation notes

Register layout: `target = Bits(n)`, `signal = Bits(n+2)`; `signal[:n]` is neighbor, `signal[n]` / `signal[n+1]` are the failure flags on the two sides, and the entry word and the position work are local registers returned after being cleaned back to zero. The preparation sequence is superposition → position query → entry query → magnitude rotation (→ directional phase for signed formats); the outer `swap(target, signal[:n])` and `swap(signal[n], signal[n+1])` implement $S$, and swapping the flags on the two sides is an essential part — swapping only the indices without the flags breaks self-adjointness. Negative off-diagonal elements acquire a $\pi$ phase under the `target < neighbor` direction convention (recorded by the `sign_convention` attribute).

`chebyshev_block(a, degree)`: requires the input BE to explicitly declare `self_adjoint_extension` (an ordinary BE does not suffice to guarantee Chebyshev semantics); {obj}`Repeat(degree) <oracq.infrastructure.ir.Repeat>` alternates "reflection about the signal-zero state + a call to a" and returns a walk power with `be_alpha = 1.0` and `argument_scale = a.alpha` whose zero-signal block implements $T_k(A/\alpha)$.

Applicability boundary: real Hermitian with a nonnegative diagonal only; a general matrix must first go through an explicit Hermitian dilation. For value widths above 12 bits the amplitude transduction is an open declaration, and the program is not executable before binding.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2). Consistent with the `sparse.py` / `block_encoding.py` rows of the validation coverage matrix:

- Structure: `tests/core/test_language.py:BlockEncodingTests` (BE construction and attribute assertions).
- Numerical: the direct witness of this entry point is registered on the `qlss.py` row — `tests/core/test_qlss_input_models.py:QLSSInputTests.test_signed_sparse_encoding_and_chebyshev`: for a 2×2 matrix with negative off-diagonal entries ($s = 2$, $a_{\max} = 1$), the zero-signal block is read out column by column, multiplied back by `alpha = 2`, and cross-checked against the classical matrix (places = 11), and the zero-signal block of `chebyshev_block(be, 3)` is cross-checked element by element against $4h^3 - 3h$ ($h = A/\alpha$, i.e. $T_3$).
- Binding: `tests/core/test_language.py:BlockEncodingTests.test_alpha_survives_ir_serialization` (the JSON round trip preserves `be_alpha`).

## Known gaps and planned stages

No known gaps, stage V1 (the gap column of the `sparse.py` / `block_encoding.py` rows in the validation coverage matrix is empty). The sign phase convention and the alpha conversion already have small-instance witnesses; amplitude transduction for large word lengths is kept as an explicit declaration pending binding, and its quantization error bound has no separately listed witness.

## Related links

- Source: `src/oracq/algorithms/input_model/sparse.py`
- Same-family pages: [Block-encoding composition algebra](block-encoding-algebra.md), [CKS Chebyshev solver](cks.md), [Costa walk solver](costa-walk.md), [Sparse access](sparse-access.md), [Select-Swap QROM](select-swap.md), [VTAA-CKS variable-time solver](vtaa-cks.md)
- API reference: [Sparse access adaptation](../../api/algorithms/input_model/sparse.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical validation script: `tests/verification/verify_blockencoding.py` (real-backend execution, no mock substitutes, no skips; all 73 cases passed on 2026-09-16), artifact `out/verification/blockencoding.json`. This page corresponds to `sparse-be-*`, `chebyshev-walk-*`, `cross-tridiagonal-*`, and `cross-qram-be-*`, 14 cases in total.

Experiment design:

1. **Sparse BE across matrices and scales**: the tridiagonal matrix $\alpha I + \beta T$ (the signed format with $\beta<0$ exercises the directional phase convention; the unsigned format exercises the unsigned path), dim ∈ {2, 4, 8, 16}, $s$ structural positions per column with boundary columns padded by zero-entry positions. The zero-signal block is extracted column by column (reference over all columns + rir-pysparq / adapter-pysparq cross), multiplied back by $\alpha = s\cdot a_{\max}$, and cross-checked element by element against the classical matrix; dim = 2 (15 qubits) and unsigned dim = 4 (12 qubits) additionally go through OriginIR-ext state-vector column-by-column extraction; signed dim = 4 reaches the 24-qubit OriginIR budget and dim = 8 at 34 qubits exceeds it, so both run only the pysparq/reference paths within the budget, noted in the case parameters.
2. **QRAM data-bound sparse BE** (dim = 4): the forward/inverse position tables and the entry table are all bound via memory (never entering the IR); reference and rir agree amplitude by amplitude.
3. **Chebyshev walk**: `chebyshev_block` on a 2×2 sparse BE with degrees $k=1..4$, the zero-signal block compared against $T_k(A/\alpha)$.
4. **Independent cross-validation against pysparq's own block encodings**: the same tridiagonal matrix is encoded via the oracq sparse BE (plus a second path via the Pauli LCU BE for dim ≤ 8) and via pysparq `BlockEncodingTridiagonal`; each side's effective block is extracted, multiplied by its own normalization, and the two are cross-checked against each other (dim ∈ {2, 4, 8, 16}); additionally, pysparq `BlockEncodingViaQRAM` (C++ criterion configuration data_size=50, rational=51, exponent=15, matrices Frobenius-normalized into the table) cross-checks two dim = 4 instances, one tridiagonal and one non-tridiagonal (matched sparse graph with $s=2$).

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| `sparse-be-signed-tridiagonal-d2` | dim 2, s=2, α=3.0 | originir-ext + three paths | max_error | 4.4e-16 |
| `sparse-be-signed-tridiagonal-d4` | dim 4, s=3, α=4.5 | three paths | max_error | 1.3e-16 |
| `sparse-be-signed-tridiagonal-d8` | dim 8, s=3, α=4.5 | reference + rir | max_error | 1.3e-16 |
| `sparse-be-unsigned-unitary-d4` | dim 4, s=3, α=4.5 | originir-ext + three paths | max_error | 1.1e-16 |
| `sparse-be-unsigned-d16` | dim 16, s=3, α=4.5 | reference + rir | max_error | 1.1e-16 |
| `sparse-be-qram-access-d4` | dim 4, QRAM bound | reference + rir | max_error / backend deviation | 1.3e-16 / 0 |
| `chebyshev-walk-k1..k4` | dim 2, $\alpha=3.0$ | three paths | max_error | 1.7e-16 / 2.5e-16 / 5.6e-16 / 6.7e-16 |
| `cross-tridiagonal-*` (5 groups) | dim 2–16 | reference × pysparq | qecc / pysparq / cross error | ≤4.4e-16 / ≤4.4e-16 / ≤6.7e-16 |
| `cross-qram-be-tridiagonal-d4` | dim 4, s=3 | reference × pysparq QRAM | qecc error / pysparq quantization error | 1.3e-16 / 2.8e-5 |
| `cross-qram-be-matched-pairs-d4` | dim 4, s=2 non-tridiagonal | ditto | cross error | 2.8e-15 |

Cross-validation finding (recorded as observed, not asserted as a defect of that library): the (0,0) block normalization of pysparq `BlockEncodingTridiagonal`, measured with a 1-bit primary register (dim = 2) and $\beta \neq 0$, is $|\alpha|+2|\beta|$ rather than the Frobenius norm of its construction documentation (on a 1-bit register both +1 and −1 trigger the overflow branch, and the anc==0 corner block loses Frobenius normalization; for $\beta=0$ there is no shift branch and it still equals $A/\lVert A\rVert_F$). Its C++ correctness test domain `randint(2,5)` (dim 4–16) does not cover that size; the script cross-checks with the empirically observed normalization and records both conventions in the case parameter `psparq_alpha_effective`.

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_blockencoding.py
```

Artifacts: `out/verification/blockencoding.json`.
