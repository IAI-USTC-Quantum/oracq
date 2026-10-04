# Alias Sampling Preparation

**English** · <a href="../../zh/manual/algorithms/alias-preparation.html">简体中文</a>

> Category C5/C1 · Module [`oracq.algorithms.common.prepare_select`](../../api/algorithms/common/prepare_select.rst) · Stage V4

## Overview

The alias-sampling PREPARE of Babbush et al. 2018 (cited in the module docstring; CP→SP / QRAM in the input-model vocabulary): the LCU probabilities $p_i = |c_i|/\alpha$ are packaged through Vose alias preprocessing into a `(keep, alt)` table; the quantum side needs only a uniform state, one table load, one comparison, and one controlled swap, avoiding a multiplexed rotation tree. `keep` is quantized by rounding down to `precision` bits, and the total variation distance between the sampling distribution and the exact distribution does not exceed $2^{\text{selector}} \cdot 2^{-\text{precision}}$.

## Interface and input model

```python
alias_table(coefficients, *, precision=8)
alias_prepare(coefficients, *, precision=8, database=None)
```

API entry points: {obj}`alias_table <oracq.algorithms.common.prepare_select.alias_table>`, {obj}`alias_prepare <oracq.algorithms.common.prepare_select.alias_prepare>`

- {obj}`alias_table <oracq.algorithms.common.prepare_select.alias_table>`: purely classical preprocessing. The coefficients are first zero-padded to $2^{\text{selector}}$ and then split into keep / alt around the mean (zero-probability slots also work); the data word is packed as `keep | (alt << precision)`; `precision` is 2..32 and selector bit width + precision ≤ 64 (the QRAM word-width limit).
- {obj}`alias_prepare <oracq.algorithms.common.prepare_select.alias_prepare>`: quantum-side assembly. `database` defaults to {obj}`qram_database(width, precision + width) <oracq.algorithms.input_model.oracles.qram_database>` and returns the default memory (key `db__table`); {obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>` or any other {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` can also be injected (address must be the selector bit width and data the precision + selector bit widths, otherwise it is rejected at generation time).

Returned objects:

| Object / attribute | Meaning |
|---|---|
| `alias_table` → {obj}`AliasTable <oracq.algorithms.common.prepare_select.AliasTable>` | `probabilities` / `keep` / `alt` / `quantized` / `precision` / `table`; `distribution()` gives the classical sampling distribution under quantized keep and uniform drawing |
| `alias_prepare` → {obj}`AliasPreparation <oracq.algorithms.common.prepare_select.AliasPreparation>` | `.preparation` ({obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>`), `.table` (`AliasTable`), `.memory` (the default QRAM binding data); `.state_preparation()` returns the handle |

Module attributes of `AliasPreparation.preparation`: `implementation="alias_sampling"`, `alias_precision`, `prepare_terms` / `prepare_alpha` / `selector_width`, `clean_work=False`.

## Implementation notes

Register layout: `target = selector(width)`, `work = data(precision + width) | coin(precision)`; the low precision bits of `data` hold keep and the high bits hold alt. Flow: apply H to target and coin → table query `data ^= (keep|alt)[target]` → the unsigned comparison `coin < keep` yields flag → when flag = 0, a controlled `swap(target, alt)` (selecting the alias slot) → the comparison is invoked once more to uncompute flag.

data / coin are entangled with selector, forming the **purified junk** of $\sum_i \sqrt{p_i}\,|i\rangle|\mathrm{junk}_i\rangle$: `clean_work=False` honestly marks that work is not uncomputed. When downstream consumes only the marginal distribution or the (0,0) block of the block encoding, junk inner products have no effect; junk is uncomputed by the adjoint PREPARE. The comparator reuses {obj}`fixed_arithmetic("lt", FixedFormat(precision, 0, signed=False)) <oracq.algorithms.common.arithmetic.fixed_arithmetic>`.

When composed with {obj}`lcu_prepare_select <oracq.algorithms.common.prepare_select.lcu_prepare_select>`, the signal bit width of the BE is the selector bit width plus the work width (the `prepare_work` attribute); for 4 coefficients and precision = 8, target has width 2 and work has width 18 (data 10 + coin 8), and the memory key of the default QRAM route is `db__table`. The gate route expands the whole table into a gate-level sub-database, suitable for witnessing small instances. Applicability boundary: the data table does not enter the IR; the default QRAM route requires memory to be supplied alongside the circuit at execution time.

## Validation approach

Category C5/C1: the selector marginal distribution of the preparation is cross-checked against the closed-form classical sampling distribution, and the (0,0) block assembled by the block encoding is cross-checked pointwise against a dense matrix per C1 (acceptance criteria in `../development/validation-plan.md` §2). Three layers of evidence in `tests/core/test_prepare_select.py:PrepareSelectTests`:

- Structural: `test_alias_preparation_samples_target_distribution` asserts `implementation == "alias_sampling"`, `clean_work` is False, and the returned `.table` matches an independent construction; `test_bad_inputs_fail_at_generation` covers generation-time rejection of database width mismatches.
- Numerical: `test_alias_table_reproduces_normalized_distribution` — the probabilities reconstructed classically from keep / alt are cross-checked against $p_i$ (places = 12), the quantized distribution has total variation delta = 0.02, and data words do not exceed `precision + selector` bits; `test_alias_preparation_samples_target_distribution` — the selector marginal distributions of both the gate_database and the default QRAM routes are cross-checked against `table.distribution()` (places = 10); `test_alias_block_encoding_approaches_hamiltonian` — the (0,0) block of the block encoding assembled from alias PREPARE, after multiplying by α, is cross-checked against the dense Hamiltonian (delta = 0.02).
- Binding: the data table is injected through the `XorDatabase` interface, and the gate-table and QRAM-resource routes are cross-checked within the same test; three-binding consistency parameterization is pending for V4.

## Known gaps and planned stages

Three-binding consistency for alias (gate / QRAM currently lack independent scenarios): unified parameterized cross-checks are not rolled out, and the distribution cross-checks of the two routes are currently embedded in the numerical witnesses. Stage V4, consistent with the prepare_select.py row of `validation-coverage.md`.

## Related links

- Source: `src/oracq/algorithms/common/prepare_select.py`
- Same-group pages: [PREPARE–SELECT decomposition](prepare-select.md), [State preparation oracle](state-preparation.md), [XOR database](xor-database.md)
- API reference: [PREPARE-SELECT decomposition](../../api/algorithms/common/prepare_select.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

The results in this section were produced by `tests/verification/verify_stateprep.py` on real backends (the reference built-in reference executor and the rir-pysparq native RIR interpreter; the alias program workspace reaches 51 qubits, exceeding the 24-qubit budget of OriginIR-ext, so it runs the sparse simulation route).

**Experiment design** (4 cases): the selector (target) marginal distribution is cross-checked, with the oracle being two independently constructed classical distributions — the exact distribution $p_i = |c_i|/\alpha$ (computed directly from the coefficients) and the quantized alias distribution $q_i \propto \text{quantized}_i + \sum_{j:\,\text{alt}_j=i}(2^p - \text{quantized}_j)$ (expanded independently from the keep/alt table by definition, not by calling `AliasTable.distribution()`; the agreement of the latter with the closed form is reported alongside as an informational metric). Instances: dual bindings of gate_database and the default QRAM for 4 coefficients (including complex ones) with precision 8, a gate binding of 6 coefficients (including one zero-probability slot, zero-padded to 8) with precision 10, and the special case of a uniform distribution where keep is exactly all 1s. A single quantum run yields the marginal distribution over all $2^{\text{selector}+\text{precision}}$ branches.

**Key metrics**:

| Case | Scale | Binding / paths | TVD (vs quantized distribution) | TVD (vs exact distribution) | Theoretical bound $2^w \cdot 2^{-p}$ |
|---|---|---|---|---|---|
| alias-preparation-gate-w2-p8 | w=2, p=8 | gate / ref+rir | 8.3e-17 | 1.33e-3 | 1.56e-2 |
| alias-preparation-qram-w2-p8 | w=2, p=8 | QRAM / ref+rir | 8.3e-17 | 1.33e-3 | 1.56e-2 |
| alias-preparation-gate-w3-p10 | w=3, p=10 (with zero slot) | gate / ref+rir | 7.6e-17 | 1.63e-4 | 7.81e-3 |
| alias-preparation-uniform-p8 | w=2, p=8 uniform | gate / ref+rir | 5.6e-17 | 5.6e-17 | 1.56e-2 |

The quantum marginal distribution agrees with the quantized alias distribution at machine precision (the cross-path TVD of the gate and QRAM bindings is 0); the TVD against the exact distribution is dominated by keep quantization and is far below the theoretical bound; in the uniform special case quantization is lossless and the TVD drops to machine precision. `distribution_self_consistency` is 0 everywhere, showing that the module's own `distribution()` agrees with the independent closed form.

**Reproduction command**:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_stateprep.py
```

**Artifact path**: `out/verification/stateprep.json`.
