# Tensor Hypercontraction

**English** · <a href="../../zh/manual/algorithms/thc.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.input_model.lowrank`](../../api/algorithms/input_model/lowrank.rst) · Stage V1

## Overview

Tensor hypercontraction (THC) compresses the two-body integral tensor of the electronic Hamiltonian into leaf-operator form

$$
H = \sum_{\mu\nu} \zeta_{\mu\nu}\, L_\mu L_\nu^\dagger ,
$$

where $\zeta$ is a real symmetric coefficient matrix and $L_\mu$ are explicit small-matrix leaf operators. The implementation follows the THC representation of Lee et al. 2021, plugged into this repository's "low-rank tensor → LCU → BE" pipeline (literature citations in the module docstring): the integral tensor is given as classical input data, this module performs no real chemistry integral computation, and the resulting block encoding can enter qubitization directly.

## Interface and input model

```python
THCDecomposition(coefficients, leaves)
thc_encoding(thc)
```

API entry points: {obj}`THCDecomposition <oracq.algorithms.input_model.lowrank.THCDecomposition>`, {obj}`thc_encoding <oracq.algorithms.input_model.lowrank.thc_encoding>`

- {obj}`THCDecomposition <oracq.algorithms.input_model.lowrank.THCDecomposition>`: the THC input model (input model CP — ζ and the leaf matrices are given directly as classical parameters, i.e. the "low-rank tensor (HAM) → BE" pipeline in algorithm-coverage). `coefficients` must be a real symmetric square matrix (symmetry tolerance 1e-12), and `leaves` are explicit small matrices of consistent dimension (powers of two within 2..32), which need not be unitary or Hermitian. The dataclass attributes `width` gives the number of target qubits and `leaf_count` the number of leaf operators.
- {obj}`thc_encoding <oracq.algorithms.input_model.lowrank.thc_encoding>`: assembles the LCU block encoding of the THC Hamiltonian.

Returns a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` with module attributes:

| Attribute | Meaning |
|---|---|
| `be_alpha` | $\sum_{\mu\nu}\lvert\zeta_{\mu\nu}\rvert\,\alpha_\mu\alpha_\nu$ ($\alpha_\mu$ is the Pauli l1 bound of a leaf operator) |
| `be_form` | `"thc"` |
| `thc_leaves` / `lcu_terms` | number of leaf operators / number of nonzero terms of the outer LCU |
| `thc_lambda` | a λ report under the same accounting, equal in value to `be_alpha` |

## Implementation notes

Each leaf operator is first expanded into a Pauli LCU block encoding via {obj}`matrix_pauli_encoding <oracq.algorithms.input_model.block_encoding.matrix_pauli_encoding>` (words with Pauli coefficients below 1e-12 are dropped; the explicit expansion supports only small instances of at most 5 bits), $L_\nu^\dagger$ takes the conjugated direction via {obj}`adjoint_be <oracq.algorithms.input_model.block_encoding.adjoint_be>`, and the $(\mu,\nu)$ term takes the product of the two BEs — signals concatenated, α multiplied. The outer LCU's PREPARE runs over $(\mu,\nu)$ pairs with weights $\propto\lvert\zeta_{\mu\nu}\rvert\alpha_\mu\alpha_\nu$, and the signal is the selector (as many bits as the term count requires) concatenated with each leaf encoding's signal; an all-zero ζ raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`.

Applicability boundary: leaf matrices of dimension 2..32; the approximation error of the THC decomposition relative to the true Hamiltonian belongs to classical preprocessing, and this module assembles the given ζ and leaves exactly. Large-scale ζ and leaf data should switch to QRAM data binding, and the product can be handed directly to `transforms.qubitization_walk`.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2): the block encoding's zero-signal block must agree with the target Hamiltonian within tolerance. Three layers of evidence:

- Structure: the construction attribute assertions of `tests/core/test_lowrank.py:ThcTests` (`be_form == "thc"`, `thc_leaves`, `leaf_count`); `ThcTests.test_invalid_inputs_fail_at_generation` covers the THC negative cases (non-symmetric ζ, non-square ζ, all-zero coefficients).
- Numerical: `ThcTests.test_thc_block_equals_hamiltonian` — two leaf matrices and a 2×2 ζ are assembled, with the expected matrix $\sum_{\mu\nu}\zeta_{\mu\nu}L_\mu L_\nu^\dagger$ assembled independently inside the test and cross-checked column by column via `assert_block_equals` (places = 9).
- Binding: this module's cross-encoding binding witness is `DoubleFactorizationTests.test_matches_pauli_encoding_block` in the same file (the DF encoding and the Pauli LCU encoding of the same Hamiltonian cross-checked column by column); THC's leaf operators are themselves the explicit Pauli LCU encoding path, so no separate abstract declaration is set up.

α tightness witness (added in V1): `ThcTests.test_thc_alpha_matches_hand_computed_bound` hand-computes the leaf bounds inside the test with an independent `pauli_l1` closed form (the coefficients of $M = c_I I + xX + yY + zZ$ expressed linearly through the matrix elements, without going through `matrix_pauli_encoding`): $\alpha_0 = 1.2$, $\alpha_1 = 1.3$, so α = 0.7·1.44 + 0.1·1.56 + 0.1·1.56 + 0.4·1.69 = 1.996, cross-checked against `be.alpha` (places = 10).

## Known gaps and planned stages

No known gaps; the stage V1 witnesses are complete (block cross-checks + an α cross-check against the independent hand computation).

## Related links

- Source: `src/oracq/algorithms/input_model/lowrank.py`
- Algorithms in the same module: [Double factorization block encoding](double-factorization.md)
- API reference: [Chemical low-rank decomposition block encodings](../../api/algorithms/input_model/lowrank.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical validation script: `tests/verification/verify_blockencoding.py` (real-backend execution, no mock substitutes, no skips; all 73 cases passed on 2026-09-16), artifact `out/verification/blockencoding.json`. This page corresponds to the two cases `thc-encoding-leaves2-d2` and `thc-encoding-leaves3-d4`.

Experiment design: 2×2 dense leaves (2 leaves, general complex matrices, neither unitary nor Hermitian) go through OriginIR-ext + UniQC `to_matrix` full-unitary extraction, with the (0,0) block multiplied by α and cross-checked element by element against the independently numpy-assembled $\sum_{\mu\nu}\zeta_{\mu\nu}L_\mu L_\nu^\dagger$, and α compared against the independently hand-computed Pauli l1 accounting $\sum_{\mu\nu}|\zeta_{\mu\nu}|\alpha_\mu\alpha_\nu$ that bypasses `matrix_pauli_encoding`; 4×4 with three leaves (diagonal leaves + a fully coupled 3×3 ζ, 7 nonzero LCU terms) goes through reference + rir-pysparq + adapter-pysparq on three backends column by column, exercising multi-leaf, multi-scale assembly.

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| `thc-encoding-leaves2-d2` | 2 leaves 2×2, α=1.381 | originir-ext + to_matrix, three-backend cross | max_error | 1.7e-16 |
| ditto | — | — | α − independent l1 accounting | 0 (exactly equal) |
| `thc-encoding-leaves3-d4` | 3 leaves 4×4 (diagonal leaves), α=1.782 | three backends | max_error / backend deviation | 5.6e-16 / 0 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_blockencoding.py
```

Artifacts: `out/verification/blockencoding.json`.
