# Double Factorization

**English** · <a href="../../zh/manual/algorithms/double-factorization.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.input_model.lowrank`](../../api/algorithms/input_model/lowrank.rst) · Stage V1

## Overview

The double factorization (DF) of the quantum chemistry electronic Hamiltonian compresses the two-body integral tensor into a low-rank form

$$
H = \mathrm{scalar}\cdot I + \sum_r U_r\, \mathrm{diag}(g_r)\, U_r^\dagger ,
$$

each rank term being a Hermitian matrix whose diagonal spectrum $g_r$ is conjugated by the rotation $U_r$. The implementation follows the "low-rank tensor → LCU → BE" pipeline of Berry et al. 2019 and von Burg et al. 2021 (literature citations in the module docstring): the integral tensor is given as classical input data, this module performs no real chemistry integral computation, and the resulting block encoding can enter qubitization directly.

## Interface and input model

```python
DoubleFactorization(scalar, rotations, spectra)
DoubleFactorization.from_symmetric(scalar, rotations, factors)
diagonalize_symmetric(matrix, *, tolerance=1e-12, max_sweeps=100)
double_factorized_encoding(df, *, name=None)
```

API entry points: {obj}`DoubleFactorization <oracq.algorithms.input_model.lowrank.DoubleFactorization>`, {obj}`diagonalize_symmetric <oracq.algorithms.input_model.lowrank.diagonalize_symmetric>`, {obj}`double_factorized_encoding <oracq.algorithms.input_model.lowrank.double_factorized_encoding>`

- {obj}`DoubleFactorization <oracq.algorithms.input_model.lowrank.DoubleFactorization>`: the DF input model (input model CP — the tensor data is given directly as classical parameters, i.e. the "low-rank tensor (HAM) → BE" pipeline in algorithm-coverage). `rotations` are explicit small unitary matrices (dimensions powers of two within 2..32, unitarity checked by column orthogonality, tolerance 1e-9), and `spectra` are the corresponding real spectra; the dataclass attributes `width` gives the number of target qubits and `rank` the number of rank terms.
- `from_symmetric`: the classical preprocessing entry point of physical DF — takes explicit unitaries $U_r$ and real symmetric $G_r$, and folds the eigenvector matrix of $G_r = V_r\,\mathrm{diag}(g_r)\,V_r^{\mathsf T}$ into the rotation ($U_r V_r$).
- {obj}`diagonalize_symmetric <oracq.algorithms.input_model.lowrank.diagonalize_symmetric>`: Jacobi eigendecomposition of a real symmetric matrix, returning `(eigenvalues, eigenvector matrix)`; non-symmetric input raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time.
- {obj}`double_factorized_encoding <oracq.algorithms.input_model.lowrank.double_factorized_encoding>`: assembles the LCU block encoding (the `name` parameter is retained but currently takes no part in naming).

Returns a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` with module attributes:

| Attribute | Meaning |
|---|---|
| `be_alpha` | $\lvert\mathrm{scalar}\rvert + \sum_r \lVert g_r\rVert_1$ (sum of the spectra's 1-norms) |
| `be_form` | `"double_factorization"` |
| `df_rank` / `lcu_terms` | number of rank terms / number of nonzero terms of the outer LCU |
| `df_lambda` | a λ report under the same accounting, equal in value to `be_alpha` |

## Implementation notes

The outer LCU's PREPARE runs over the rank index $r$ with weights $\propto \lVert g_r\rVert_1$; a nonzero `scalar` joins the same LCU as `(scalar, identity)`, and all-zero coefficients (including all-zero spectra) raise `ValidationError`. Each rank term's term block encoding is $U_r^\dagger\cdot\mathrm{diag}\cdot U_r$: the diagonal part applies a controlled $R_y(\theta_t)$ on a single-bit signal for every basis state $\lvert t\rangle$, with $\cos(\theta_t/2) = g_t/\alpha_r$ and $\alpha_r = \sum_p\lvert g_p\rvert$ taken as the spectral 1-norm, so the (0,0) block is exactly $\mathrm{diag}(g_r)/\alpha_r$; the synthesized $U_r$ is applied on the two sides to complete the conjugation.

$U_r$ is synthesized by two-level decomposition: columns are eliminated into diagonal phases and replayed in reverse order; each two-level unitary is realized by a ZYZ decomposition ($e^{i\varphi}R_z(\alpha)R_y(\beta)R_z(\gamma)$) plus multi-controlled X transpositions along a Gray path. Register layout: the term block encoding is `target(n) | signal(1)`, and the outer LCU's signal is the selector (as many bits as the term count requires) concatenated with `work(1)`.

Applicability boundary: the explicit-matrix path supports only validation-sized instances of dimension 2..32 (at most 5 qubits); at larger scales the spectra and rotation angle tables of $U_r$/$G_r$ should switch to QRAM data binding (see the three-layer paradigm of `prepare_select`). In terms of normalization, the DF α of a single-bit positive definite $g$ equals $\mathrm{tr}(g)$ and introduces no off-diagonal redundancy from the Pauli expansion (tightness witness in the validation section).

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2): the block encoding's zero-signal block must agree with the target Hamiltonian within tolerance. Three layers of evidence:

- Structure: `tests/core/test_lowrank.py:DiagonalizeSymmetricTests.test_reconstructs_factor` ($V\,\mathrm{diag}(\lambda)\,V^{\mathsf T}$ reconstruction and eigenvector orthogonality, places = 10) and `test_rejects_non_symmetric`; the construction attribute assertions of `DoubleFactorizationTests` (`be_form`, `df.rank == 1` / `df.width == 1`); `ThcTests.test_invalid_inputs_fail_at_generation` covers the DF negative cases (empty rank terms, mismatched spectrum lengths, non-unitary rotations, non-numeric `scalar`, inconsistent counts in `from_symmetric`).
- Numerical: `DoubleFactorizationTests.test_single_rank_block_equals_hamiltonian` — single-rank DF cross-checked column by column against $G$ via `assert_block_equals` (places = 9); `test_scalar_and_rotated_ranks` — a scalar and a Hadamard-rotated second rank, with the expected matrix assembled independently inside the test before cross-checking, and α cross-checked against $\lvert\mathrm{scalar}\rvert + \sum_r\lVert g_r\rVert_1$ (places = 10); `test_feeds_qubitization_walk` — after the product enters `transforms.qubitization_walk` there are no unbound slots and the registers are `["target", "signal"]`.
- Binding: `DoubleFactorizationTests.test_matches_pauli_encoding_block` — the DF encoding and the Pauli LCU encoding of the same Hamiltonian are cross-checked column by column across BEs via the `block_column` helper (places = 9), covering the consistency of two independent encoding paths.

α tightness witnesses (added in V1): `test_alpha_matches_closed_form_eigenvalues` independently computes $\sum\lvert\lambda\rvert$ from the 2×2 closed-form eigenvalues $\lambda_\pm = (t \pm \sqrt{t^2-4d})/2$ ($t$ the trace, $d$ the determinant) and cross-checks it against `be.alpha` (places = 10, measured 2.0); `test_df_alpha_tighter_than_pauli` takes an off-diagonal-dominant positive definite $g$ and measures DF α = 1.4 ≤ Pauli LCU α = 1.5 (theoretical condition: on single-bit PSD matrices Pauli α − DF α = $\lvert g_{01}\rvert + \lvert\Delta/2\rvert - \mathrm{tr}/2$, positive in the off-diagonal-dominant case).

## Known gaps and planned stages

No known gaps; the stage V1 witnesses are complete (block cross-checks + α tightness against the closed-form eigenvalues and the Pauli comparison + cross-encoding binding consistency + the qubitization hookup).

## Related links

- Source: `src/oracq/algorithms/input_model/lowrank.py`
- Algorithms in the same module: [THC block encoding](thc.md)
- API reference: [Chemical low-rank decomposition block encodings](../../api/algorithms/input_model/lowrank.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical validation script: `tests/verification/verify_blockencoding.py` (real-backend execution, no mock substitutes, no skips; all 73 cases passed on 2026-09-16), artifact `out/verification/blockencoding.json`. This page corresponds to `diagonalize-symmetric-*`, `df-encoding-*`, and `diagonal-be-*`, 11 cases in total (the diagonal block encoding `_diagonal_encoding` is exactly the kernel of a DF rank term).

Experiment design:

1. **Diagonal BE across scales** (`_diagonal_encoding`, the diagonal BE named in the paper): three spectrum groups (mixed signs, containing zeros, generated from trigonometric functions) with dim ∈ {2, 4, 8} go through OriginIR-ext + UniQC `to_matrix` full-unitary extraction (the (0,0) block compared against $\mathrm{diag}(g)/\alpha$, with leakage into the failure branch reported), and dim ∈ {16, 32, 64} go through reference + rir-pysparq column by column; there is also a public assembly-consistency case — the same diagonal spectrum yields element-wise identical blocks via `double_factorized_encoding` (identity rotations, two ranks with the same spectrum) and via the direct diagonal BE.
2. **Jacobi eigendecomposition**: random real symmetric matrices with dim ∈ {2, 4, 8, 16, 32}, eigenvalues compared against `numpy.linalg.eigh`, plus checks of the reconstruction $V\mathrm{diag}(\lambda)V^{\mathsf T}$ and orthogonality.
3. **DF block encoding across ranks and scales**: 2×2 double rank + scalar (`from_symmetric` preprocessing, full-unitary extraction), 4×4 three ranks with random orthogonal rotations (full-unitary extraction), 8×8 double rank (reference + rir); the expected matrices are assembled independently with numpy, and α is compared against the closed form $|\mathrm{scalar}|+\sum_r\lVert g_r\rVert_1$.

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| `diagonal-be-unitary-d2/d4/d8` | dim 2/4/8, α=2.0/3.25/4.0 | originir-ext + to_matrix, three-backend cross | max_error | 5.6e-17 / 5.6e-17 / 1.1e-16 |
| `diagonal-be-wide-d16/d32/d64` | dim 16/32/64 | reference + rir | max_error | 1.1e-16 / 1.0e-16 / 1.1e-16 |
| `diagonal-be-public-df-d4` | dim 4 (two paths) | reference | max_error | 5.6e-17 |
| `diagonalize-symmetric-d2..d32` | 5 groups | independent classical cross-check | reconstruction / eigenvalue / orthogonality error | ≤6.2e-13 / ≤5.7e-14 / ≤7.3e-15 |
| `df-encoding-rank2-scalar-d2` | rank 2 + scalar, α=5.05 | to_matrix + three paths | max_error / α − closed form | 4.4e-16 / −8.9e-16 |
| `df-encoding-rank3-d4` | rank 3, α=8.24 | to_matrix + reference | max_error | 1.6e-15 |
| `df-encoding-rank2-d8` | rank 2 + scalar 0.7 | reference + rir | max_error | 1.2e-15 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_blockencoding.py
```

Artifacts: `out/verification/blockencoding.json`.
