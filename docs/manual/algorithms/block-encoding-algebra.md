# Block Encoding Algebra

**English** · <a href="../../../zh/manual/algorithms/block-encoding-algebra.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.input_model.block_encoding`](../../api/algorithms/input_model/block_encoding.rst) · Stage V1

## Overview

A block encoding hides the matrix $A$ in the corner block of a larger unitary: the zero-signal block of $U$ satisfies $(\langle 0| \otimes I)\, U\, (|0\rangle \otimes I) = A/\alpha$, with the normalization $\alpha$ stored alongside the operation as the RIR module attribute `be_alpha`; precision is not part of the language core. This module provides the compositional arithmetic of BEs: linear combination (LCU), tensor, adjoint, signal-space dilation, direct sum, Kronecker sum, Boolean embedding, and explicit Pauli expansion of small matrices. Upper-layer algorithms (QSVT, qubitization, QLSS, low-rank decomposition, ODE solvers) all assemble their inputs on this algebra; the lowest-layer {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` type and the binary primitives {obj}`identity <oracq.algorithms.input_model.operators.identity>` / {obj}`pauli_x <oracq.algorithms.input_model.operators.pauli_x>` / {obj}`zero <oracq.algorithms.input_model.operators.zero>` / {obj}`scale <oracq.algorithms.input_model.operators.scale>` / {obj}`linear_combination <oracq.algorithms.input_model.operators.linear_combination>` / {obj}`product <oracq.algorithms.input_model.operators.product>` live in `operators.py`.

## Interface and input model

The input model is unified as BE: all entry points receive and return `BlockEncoding` ({obj}`reflect_zero <oracq.algorithms.input_model.block_encoding.reflect_zero>` is a {obj}`Builder <oracq.infrastructure.builder.Builder>` helper).

```python
lcu(terms)
tensor(a, b)
adjoint_be(a)
pad_signal(a, width)
kronecker_sum(a, b=None)
direct_sum(a, b)
projector(width, accepted)
truncated_shift(width, last)
pauli_word(word)
matrix_pauli_encoding(matrix, *, drop_tolerance=1e-12)
reflect_zero(builder, register, *, positive=False)
```

API entry points: {obj}`lcu <oracq.algorithms.input_model.block_encoding.lcu>`, {obj}`tensor <oracq.algorithms.input_model.block_encoding.tensor>`, {obj}`adjoint_be <oracq.algorithms.input_model.block_encoding.adjoint_be>`, {obj}`pad_signal <oracq.algorithms.input_model.block_encoding.pad_signal>`

Semantics and output normalization of the main entry points:

| Entry point | Semantics | Output alpha |
|---|---|---|
| {obj}`lcu(terms) <oracq.algorithms.input_model.block_encoding.lcu>` | PREPARE/SELECT composition of $\sum_j c_j A_j$ (complex phases via a global phase; zero coefficients dropped, a single term degenerates to `scale`) | $\sum_j \lvert c_j\rvert\, \alpha_j$ |
| {obj}`tensor(a, b) <oracq.algorithms.input_model.block_encoding.tensor>` | $A \otimes B$ | $\alpha_A \alpha_B$ |
| {obj}`adjoint_be(a) <oracq.algorithms.input_model.block_encoding.adjoint_be>` | $A^\dagger$ | $\alpha_A$ |
| {obj}`kronecker_sum(a, b) <oracq.algorithms.input_model.block_encoding.kronecker_sum>` | $A \otimes I + I \otimes B$ | $\alpha_A + \alpha_B$ |
| {obj}`direct_sum(a, b) <oracq.algorithms.input_model.block_encoding.direct_sum>` | direct sum of same-width matrices (single select bit) | $\alpha_A + \alpha_B$ |
| {obj}`pauli_word(word) <oracq.algorithms.input_model.block_encoding.pauli_word>` / {obj}`matrix_pauli_encoding(matrix) <oracq.algorithms.input_model.block_encoding.matrix_pauli_encoding>` | Pauli-word BE / explicit Pauli LCU of a small matrix | 1 / LCU normalization |

`matrix_pauli_encoding` computes the expansion coefficients analytically by accumulating phases column by column according to $\operatorname{Tr}(P^\dagger M)/2^n$, drops terms with magnitude below `drop_tolerance`, and degenerates to `zero` when everything is zero; it accepts only square matrices of power-of-two dimension and at most 5 bits, is positioned as an explicit gate implementation for small-scale applications, and claims no quantum speedup for matrix input or classical expansion.

## Implementation notes

`lcu` is the core combinator: the selector bit width is $\lceil\log_2 L\rceil$, the weights $\sqrt{\lvert c_j\rvert \alpha_j / \alpha}$ are prepared via {obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>`, each term is invoked controlled on its own signal partition, and an inverse preparation closes it; the binary `linear_combination` in `operators.py` is the isomorphic special case with a single selector, and together they cover the composition semantics. {obj}`pad_signal <oracq.algorithms.input_model.block_encoding.pad_signal>` only allows enlarging the signal space (shrinking is rejected), for late binding under the same signature. {obj}`projector <oracq.algorithms.input_model.block_encoding.projector>` / {obj}`truncated_shift <oracq.algorithms.input_model.block_encoding.truncated_shift>` construct the $\alpha = 1$ BEs of Boolean embeddings, feeding structured assembly such as `direct_sum`; `reflect_zero` implements the positive/negative reflection about the signal zero state and is a shared building block for walk-type operators.

Applicability boundary: all compositions guarantee only corner-block semantics and alpha arithmetic; the mathematical claims of the input BEs are not checked; `be_alpha` must be a finite positive number — invalid values (0, negative numbers, `True`, inf) are rejected at construction time.

## Validation approach

Category C2 (approximate continuous semantics; acceptance criteria in `../development/validation-plan.md` §2). Consistent with the `sparse.py` / `block_encoding.py` rows of the validation coverage matrix; all three layers of evidence are in `tests/core/test_language.py:BlockEncodingTests`:

- Structural: in-class constructions and attribute assertions; `test_invalid_alpha` covers alpha validity (0, −1, `True`, and inf all raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`).
- Numerical: `test_unequal_normalization_and_signed_sum` ($\mathrm{lc}(-0.5,\ 2I,\ 2,\ 3X)$ with normalization 7 and matrix $[[-1, 6], [6, -1]]$), `test_complex_coefficients` (phases of complex coefficients enter the corner block), `test_product_order_and_independent_signal_spaces` (product order and signal-space partitioning, cross-checked element by element at $k = 2\cos(0.4)/\sqrt{2}$), `test_zero_and_zero_coefficient` (zero-coefficient and zero-matrix degeneration). The witnessing technique simulates each initial `target` basis state column by column, reads the zero-signal block amplitudes, and multiplies alpha back in (the same convention as `assert_block_equals` in `tests/core/witness.py`).
- Binding: `test_alpha_survives_ir_serialization` — after a JSON {obj}`dumps <oracq.infrastructure.serialization.dumps>`/{obj}`loads <oracq.infrastructure.serialization.loads>` round trip of `scale(4, identity(2))`, `be_alpha` stays 4.

## Known gaps and planned stages

No known gaps; stage V1 (the gap column of the `sparse.py` / `block_encoding.py` rows in the validation coverage matrix is empty).

## Related links

- Source: `src/oracq/algorithms/input_model/block_encoding.py` (combinators) and `src/oracq/algorithms/input_model/operators.py` (the `BlockEncoding` type and binary primitives)
- Same-family pages: [sparse matrix block encoding](sparse-block-encoding.md), [QSVT matrix inversion](qsvt-matrix-inversion.md), [truncated Taylor block encoding](taylor-block-encoding.md)
- API reference: [Block-encoding composition](../../api/algorithms/input_model/block_encoding.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-level numerical validation script: `tests/verification/verify_blockencoding.py` (executed on real backends, no mocks, no skips; as of 2026-09-16 all 73 cases pass), artifact `out/verification/blockencoding.json`. This page corresponds to 12 cases: `matrix-pauli-encoding-*`, `be-algebra-*`, `pauli-word-YZX`, `projector-w2-03`, `truncated-shift-w2-l3`.

Experiment design: random Hermitian / non-Hermitian small matrices with a fixed seed are encoded via `matrix_pauli_encoding`, and seven combinators (tensor / adjoint / pad_signal / lcu with complex coefficients / {obj}`kronecker_sum <oracq.algorithms.input_model.block_encoding.kronecker_sum>` / {obj}`direct_sum <oracq.algorithms.input_model.block_encoding.direct_sum>` / product) each assemble one; the full program unitary is exported via OriginIR-ext and extracted with UniQC `Circuit.to_matrix`; `harness.effective_block` extracts the zero-signal effective block and reports an upper bound on the leakage of failure branches (signal ≠ 0); the expected matrices are assembled independently with numpy (Pauli-word convention: `word[0]` acts on the least significant bit). `pauli_word` (no signal bit, leakage exactly 0) and the Boolean embeddings `projector` / `truncated_shift` are likewise cross-checked element by element against explicit classical matrices. Each case is additionally crossed column by column with reference / rir-pysparq / adapter-pysparq (all deviations exactly 0), and α is cross-checked against an independently computed Pauli l1 bound.

| Case | Scale | Backend paths | Metric | Value |
|---|---|---|---|---|
| `matrix-pauli-encoding-2x2-general` | 2×2 non-Hermitian | originir-ext + to_matrix, three-backend cross | max_error | 1.3e-16 |
| `matrix-pauli-encoding-4x4-hermitian` | 4×4 Hermitian | ditto | max_error | 4.6e-16 |
| same as above | — | — | α − independent Pauli l1 | 0 (exactly equal to the bound) |
| `be-algebra-tensor` | 2×2 ⊗ 2×2 | originir-ext + to_matrix | max_error | 1.3e-16 |
| `be-algebra-adjoint` | 2×2 non-Hermitian | ditto | max_error | 2.3e-16 |
| `be-algebra-pad-signal` | signal +2 bits | ditto | max_error | 1.2e-16 |
| `be-algebra-lcu-complex` | complex coefficients 0.6+0.3j / −0.8 | ditto | max_error | 2.4e-16 |
| `be-algebra-kronecker-sum` | A⊗I+I⊗B | ditto | max_error | 2.5e-16 |
| `be-algebra-direct-sum` | same-width direct sum | ditto | max_error | 1.1e-16 |
| `be-algebra-product` | A·C (non-Hermitian) | ditto | max_error | 2.2e-16 |
| `pauli-word-YZX` | 3 bits, α=1 | originir-ext + to_matrix | max_error / leakage | 0 / 0 |
| `projector-w2-03` / `truncated-shift-w2-l3` | 2 bits | ditto | max_error | 0 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_blockencoding.py
```

Artifact: `out/verification/blockencoding.json`.
