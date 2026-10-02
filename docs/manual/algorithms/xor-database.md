# XOR Database

**English** · <a href="../../../zh/manual/algorithms/xor-database.html">简体中文</a>

> Category C5 · Module [`oracq.algorithms.input_model.oracles`](../../api/algorithms/input_model/oracles.rst) · Stage V4

## Overview

The standard interface for quantum access to classical data: given a query table `memory`, it implements the paradigm `database_xor`

$$
|a,\,d\rangle \;\mapsto\; |a,\,d \oplus \mathtt{memory}[a]\rangle ,
$$

where both `address` and `data` are bits registers. This XOR semantics is self-inverse under arbitrary initial values and is the common foundation of downstream constructions such as function-table loading, matrix-element lookup, and LCU coefficient tables (QRAM in the input-model vocabulary). The module offers three binding layers — the abstract open declaration, the gate truth table, and the QRAM resource — plus a variant that splits words wider than 64 bits into multiple banks.

## Interface and input model

```python
abstract_database(name, address_width, data_width)
gate_database(address_width, data_width, table, *, name=None)
qram_database(address_width, data_width, *, name=None)
banked_database(address_width, data_width, *, name="BankedData", abstract=False)
```

API entry points: {obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>`, {obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>`, {obj}`qram_database <oracq.algorithms.input_model.oracles.qram_database>`, {obj}`banked_database <oracq.algorithms.input_model.oracles.banked_database>`

- {obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>`: an open declaration (the abstract layer for QRAM in the input model) with an empty body; the program closes after an implementation is bound via `linking.bind`.
- {obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>`: a truth-table witness. `table` is a sparse dictionary or a sequence of words (enumerated by address); addresses and words are range-checked at generation time; each entry flips the data bits controlled by its address.
- {obj}`qram_database <oracq.algorithms.input_model.oracles.qram_database>`: declares a {obj}`QRAM(address_width, data_width) <oracq.infrastructure.ir.QRAM>` resource (named `table`); at execution time the data is supplied as a memory dictionary and does not enter the IR.
- {obj}`banked_database <oracq.algorithms.input_model.oracles.banked_database>`: splits `data_width` into `data0, data1, ...` registers of at most 64 bits and the corresponding QRAM banks (little-endian); with `abstract=True` it yields an open declaration carrying `logical_data_width`.

All four entry points return {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` (an {obj}`OracleView <oracq.algorithms.input_model.contracts.OracleView>`; `xor_database()` returns itself). Attributes:

| Attribute | Meaning |
|---|---|
| `operation` | the underlying {obj}`Operation <oracq.infrastructure.builder.Operation>` (registers are exactly `address`, `data`) |
| `address_width` / `data_width` | interface bit widths |
| `spec` / `describe()` | oracle description (type `database_xor`, capabilities, open/closed state) |

Module attributes: `implementation` (`gate_truth_table` / `qram`), `implementation_status`; the banked version additionally has `logical_data_width` and `word_order="little_endian_banks"`.

## Implementation notes

The gate version applies X bit by bit under address control for each nonzero table word; the gate count grows linearly with the table size, so it suits only small-instance witnessing. The QRAM version externalizes the data as an execution-time resource: only a single query primitive remains in the IR, and under nested calls the resource name is lifted with the module prefix (e.g. `db__table`). Data tables and gate-level sub-databases are both expanded at generation time and do not enter the IR JSON.

Downstream adaptations: {obj}`sparse_entry <oracq.algorithms.input_model.oracles.sparse_entry>` concatenates `row|column` into a 2w-bit address (see [sparse access](sparse-access.md)); {obj}`diagonal_block_encoding <oracq.algorithms.input_model.oracles.diagonal_block_encoding>` controls a signal Ry by the queried word, encoding the diagonal matrix determined by the angle table; {obj}`alias_prepare <oracq.algorithms.common.prepare_select.alias_prepare>` injects the (keep, alt) table through the `database` parameter (see [Alias sampling preparation](alias-preparation.md)); the QROM family returns the same `XorDatabase` interface (see [QROM lookup](qrom-lookup.md)).

## Validation approach

Category C5 (data-access layer; acceptance criteria in `../development/validation-plan.md` §2): pointwise-correct query semantics + three-layer binding consistency + uncomputation. Three layers of evidence:

- Structural: `tests/core/test_contracts.py:OracleContractTests.test_aggregate_report_and_kernel_not_executed` — when the `abstract_database` input does not satisfy the protocol contract, the kernel is not executed (the error is located before the algorithm runs); `tests/core/test_open_ir.py:OpenIRTests.test_open_roundtrip_and_export_boundary` (the open declaration stays unresolved through a serialization round trip, and export is refused).
- Numerical: the per-address semantics of the gate version serves as the downstream cross-check baseline: `tests/core/test_data_loading.py:DataLoadingTests.test_select_swap_matches_gate_database_per_address` and `test_result_invariant_across_partitions_and_qrom_lookup_baseline` (pointwise exact equality over the 16-address table); `tests/core/test_qlss_input_models.py:QLSSInputTests.test_signed_sparse_encoding_and_chebyshev` carries matrix elements with `gate_database` and cross-checks the block encoding element by element (places = 11).
- Binding: `tests/core/test_open_ir.py:OpenIRTests.test_partial_binding_and_lifted_resources` (after `qram_database` is bound via {obj}`Binding <oracq.infrastructure.linking.Binding>`, the program closes and resource renames are lifted) and `test_shared_capture_has_one_entry_resource` (when the same implementation is bound to two slots, a single resource is shared); survival of capability restrictions at the contract layer is covered by `OracleContractTests.test_capability_restrictions_survive_annotation_and_wrapping`; wide-word splitting is covered by `tests/core/test_workloads.py:WorkloadTests.test_banked_words_do_not_violate_register_storage_width` (a 96-bit word splits into 64+32 with exactly two resources after closing).

## Known gaps and planned stages

Three-binding consistency parameterization (to be rolled out in V4): the gate / QRAM (and downstream variant) bindings of the same abstract declaration do not yet have unified parameterized cross-checks; each binding is currently witnessed independently. Consistent with the oracles.py row of `validation-coverage.md`.

## Related links

- Source: `src/oracq/algorithms/input_model/oracles.py`
- Tutorial: [Replacing an algorithm's oracle](../../tutorials/oracle-binding.md)
- Same-group pages: [QROM lookup](qrom-lookup.md), [Select-Swap QROM](select-swap.md), [sparse access](sparse-access.md), [state preparation oracle](state-preparation.md), [PREPARE–SELECT decomposition](prepare-select.md)
- API reference: [Oracle declarations and implementations](../../api/algorithms/input_model/oracles.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-level numerical validation script: `tests/verification/verify_oracles.py` (group-A cases), artifact `out/verification/oracles.json`.

Experiment design: the gate truth table exhausts the entire input domain in two modes — **pointwise over basis states** (each of the $2^{a+d}$ $(a,d)$ initial states runs once, asserting the output is exactly the single basis state $|a,\,d \oplus table[a]\rangle$) and **full superposition** (address/data both Hadamard-transformed, one run covers all branches, cross-checked amplitude by amplitude against the classical permutation). The QRAM binding runs the same superposition exhaustion with an execution-time memory dictionary, plus spot checks of basis states with nonzero data initial values. The dual-binding consistency cases bind the same `abstract_database` open declaration to the gate and QRAM implementations respectively and cross-check the three execution results pairwise — i.e. the first round of numerical evidence for the three-binding consistency parameterization (V4) in "Known gaps" above. Backend paths: reference (the built-in reference executor), rir-pysparq (the PySparQ native RIR interpreter), adapter-pysparq (the oracq event adapter), originir-ext (UniQC full-amplitude state vector).

| Case | Scale (address×data) | Paths | Metric values |
|---|---|---|---|
| xor-gate-basis-2x3 | 2×3, 32 initial states | reference, rir-pysparq, originir-ext | failures = 0 |
| xor-gate-basis-3x2 | 3×2, 32 initial states | reference, rir-pysparq | failures = 0 |
| xor-gate-superposition-2x3 | 2×3, 32 branches | all four paths | max_error = 5.6e-17 |
| xor-gate-superposition-3x4 | 3×4, 128 branches | all four paths | max_error = 4.2e-17 |
| xor-gate-superposition-4x3 | 4×3, 128 branches | all four paths | max_error = 4.2e-17 |
| xor-qram-superposition-3x4 | 3×4, 128 branches | all four paths | max_error = 4.2e-17, basis_failures = 0 |
| xor-abstract-binding-consistency-3x2 | 3×2, gate/qram dual binding | reference, rir-pysparq, originir-ext | max_error = 5.6e-17 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_oracles.py
```

Artifact path: `out/verification/oracles.json`.
