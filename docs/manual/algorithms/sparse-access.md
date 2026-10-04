# Sparse Access

**English** · <a href="../../zh/manual/algorithms/sparse-access.html">简体中文</a>

> Category C5 · Module [`oracq.algorithms.input_model.oracles`](../../api/algorithms/input_model/oracles.rst) · Stage V4

## Overview

The sparse matrix oracle (SO in the input-model vocabulary) follows the CKS position / entry separation convention and is organized as **a bundle of two query operations** rather than a fictitious total unitary:

- Position $P_A^{(1)}$ (paradigm `sparse_location_inplace`, in place): $|\text{column}, i\rangle \mapsto |\text{column}, \nu(\text{column}, i)\rangle$, where $\nu$ enumerates the nonzero rows of the column;
- Entry $P_A^{(2)}$ (paradigm `sparse_entry_xor`): $|\text{row}, \text{column}, d\rangle \mapsto |\text{row}, \text{column}, d \oplus A[\text{row},\text{column}]\rangle$, any row and column can be queried.

The gate implementation requires a full permutation extension per column (the partial permutation over the nonzero rows is completed into a permutation on $[0, 2^w)$), while the QRAM implementation uses a forward and an inverse lookup table.

## Interface and input model

```python
abstract_sparse_access(name, width, value_width, sparsity, work_width=None)
sparse_location_gate(width, permutations, *, work_width=None, name=None)
sparse_location_qram(width)
sparse_entry(database, width)
```

API entry points: {obj}`abstract_sparse_access <oracq.algorithms.input_model.oracles.abstract_sparse_access>`, {obj}`sparse_location_gate <oracq.algorithms.input_model.oracles.sparse_location_gate>`, {obj}`sparse_location_qram <oracq.algorithms.input_model.oracles.sparse_location_qram>`, {obj}`sparse_entry <oracq.algorithms.input_model.oracles.sparse_entry>`

- {obj}`abstract_sparse_access <oracq.algorithms.input_model.oracles.abstract_sparse_access>`: an open declaration producing the two operations `name_position` / `name_entry`; the position operation carries the `sparsity` and `full_permutation_extension=True` attributes; `work_width` defaults to `width`.
- {obj}`sparse_location_gate <oracq.algorithms.input_model.oracles.sparse_location_gate>`: `permutations` gives the full $2^w \times 2^w$ permutation table per column.
- {obj}`sparse_location_qram <oracq.algorithms.input_model.oracles.sparse_location_qram>`: declares the two {obj}`QRAM(2*width, width) <oracq.infrastructure.ir.QRAM>` tables `forward` / `inverse`.
- {obj}`sparse_entry <oracq.algorithms.input_model.oracles.sparse_entry>`: adapted from {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>`, requires `address_width == 2*width` (`row|column` concatenated addressing).

{obj}`SparseAccess(location, entry, width, value_width, sparsity) <oracq.algorithms.input_model.oracles.SparseAccess>` is a frozen dataclass constraining $1 \le \text{width} \le 64$, $1 \le \text{value\_width} \le 64$, $1 \le \text{sparsity} \le 2^{\text{width}}$, with position signature `("column", "index", "work")` and entry signature `("row", "column", "data")`. Attributes:

| Attribute | Meaning |
|---|---|
| `location` / `entry` | the two underlying {obj}`Operation <oracq.infrastructure.builder.Operation>` objects |
| `width` / `value_width` / `sparsity` | dimension bits, value word width, sparsity |
| `describe()` | {obj}`OracleSpec <oracq.algorithms.input_model.contracts.OracleSpec>` (type `cks_sparse`, `anc_qubit=None`, position / entry component descriptions, sparsity parameter) |

## Implementation notes

The gate position implementation executes the permutation column by column under control: each cycle is decomposed along a single-bit-difference path into transpositions, and each transposition is further realized as a forward path of controlled X plus a return path, with the intermediate basis states restored. The gate count explodes with $2^w$; the implementation only suits witnessing on small instances.

The QRAM position implementation performs the in-place permutation in three steps: `work ^= forward[column, index]`, `swap(index, work)`, `work ^= inverse[column, index]` — the last step cleans work back to zero while index holds the new value. In the catalog cases the forward and inverse tables are bound after resource renaming (`forward` / `inverse`).

`SparseAccess` is a bundle, not a unitary: as a whole `anc_qubit=None`, which avoids disguising the two query operations as one unitary. Downstream it is consumed by `qlss.SparseSystem` (the sparse input model of the CKS Chebyshev / Costa walk solvers) and can also be turned into a block encoding via the sparse→BE adaptation in `sparse.py`; `applications/catalog.py` and `examples/ode_input_models.py` give assembly examples.

## Validation approach

Category C5 (data-access layer; acceptance criteria in `../development/validation-plan.md` §2). Three layers of evidence:

- Structure: `tests/core/test_contracts.py:OracleContractTests.test_sparse_is_a_bundle_not_a_fictitious_unitary` — `type == "cks_sparse"`, `anc_qubit` is `None`, and the component descriptions and sparsity parameter are correct.
- Numerical: `tests/core/test_qlss_input_models.py:QLSSInputTests.test_signed_sparse_encoding_and_chebyshev` — after `sparse_location_gate` + `sparse_entry(gate_database(...))` assembles a 2×2 symmetric sparse system, the (0,0) block of the sparse→BE adaptation is cross-checked element by element against the dense matrix (places = 11, α = 2); `test_matrix_probe_recovers_scalar_system_norm` reuses the same entry query as a norm probe (places = 10). The real-backend side is consumed by `tests/integration/test_qfvm_input_models.py:QfvmInputNativeTests` (the integrity of the location permutation and the raw QRAM fields).
- Binding: the paradigm capability restrictions are witnessed by `OracleContractTests.test_capability_restrictions_survive_annotation_and_wrapping`; the QRAM position implementation enters the catalog cases (the closure and serialization checks of `tests/core/test_workloads.py:WorkloadTests.test_every_catalog_case_has_a_closed_description`).

## Known gaps and planned stages

Three-binding consistency parameterization (to be rolled out in V4): a unified parameterized cross-check of the gate / QRAM position implementations against the abstract declaration is not rolled out; the existing evidence is scattered across the downstream witnesses. Consistent with the oracles.py row of `validation-coverage.md`.

## Related links

- Source: `src/oracq/algorithms/input_model/oracles.py`
- Same-group pages: [XOR database](xor-database.md)
- API reference: [Oracle declarations and implementations](../../api/algorithms/input_model/oracles.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical validation script: `tests/verification/verify_blockencoding.py` (real-backend execution, no mock substitutes, no skips; all 73 cases passed on 2026-09-16), artifact `out/verification/blockencoding.json`. This page corresponds to `sparse-access-layer` plus `sparse-lookup-helpers`, `sparse-boolean-helpers`, and `sparse-rotation-helpers`, 4 cases in total (for the end-to-end behavior of the access layer inside a block encoding see the numerical validation section of [Sparse matrix block encoding](sparse-block-encoding.md)).

Experiment design (dim = 4, tridiagonal structural positions, fixed-point format 3.1 signed):

1. Position oracle `sparse_location_gate`: all 16 (column, index) basis-state inputs are read out one by one against the classical permutation table, and work must be cleaned back to 0;
2. Entry oracle `sparse_entry`: all 16 (row, column) × two initial data values (0 and the nonzero value 5) verify the XOR semantics $d \oplus A_{rc}$;
3. QRAM position implementation `sparse_location_qram`: the forward/inverse tables are bound as memory data, the clean-up semantics are verified basis state by basis state, and the reference and rir-pysparq backends agree amplitude by amplitude;
4. Access-layer helpers: {obj}`reversible_lookup <oracq.algorithms.input_model.sparse.reversible_lookup>` with fused address/data views over 128 input groups, {obj}`batch_lookup <oracq.algorithms.input_model.sparse.batch_lookup>` with three concurrent lookups, {obj}`compare_words <oracq.algorithms.input_model.sparse.compare_words>` (eq/lt, exhaustive over all inputs for bit widths 1–4), {obj}`value_transposition <oracq.algorithms.input_model.sparse.value_transposition>` (all 512 input groups at 3 bits), {obj}`prefix_state <oracq.algorithms.input_model.sparse.prefix_state>` uniform prefix superposition, and the analytic probability semantics of {obj}`word_rotation <oracq.algorithms.input_model.sparse.word_rotation>` / {obj}`magnitude_rotation <oracq.algorithms.input_model.sparse.magnitude_rotation>`.

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| `sparse-access-layer` (gate position) | dim 4, 16 basis states | reference | permutation mismatch count | 0 |
| ditto (entry XOR) | 16 × 2 groups | reference | mismatch count | 0 |
| ditto (QRAM position) | 16 basis states | reference + rir | mismatch count / backend deviation | 0 / 0 |
| `sparse-lookup-helpers` | 128 groups of fused views + batching | reference | mismatch count | 0 |
| `sparse-boolean-helpers` | eq/lt exhaustive + transposition 512 groups | reference | mismatch count / prefix uniformity error | 0 / 1.1e-16 |
| `sparse-rotation-helpers` | bit widths 2–4 / format 4.1 | reference + rir | probability semantics error | 5.6e-16 / 1.1e-16 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_blockencoding.py
```

Artifacts: `out/verification/blockencoding.json`.
