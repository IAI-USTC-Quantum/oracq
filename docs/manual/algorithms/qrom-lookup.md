# QROM Lookup

**English** · <a href="../../../zh/manual/algorithms/qrom-lookup.html">简体中文</a>

> Category C5 · Module [`oracq.algorithms.input_model.data_loading`](../../api/algorithms/input_model/data_loading.rst) · Stage V4

## Overview

QROM lookup implements the classical table `T` as an XOR database $|a, d\rangle \mapsto |a, d \oplus T[a]\rangle$ (QRAM in the input-model vocabulary). Implementation follows Low–Kliuchnikov–Schaeffer 2018 and Babbush et al. 2018 (cited in the module docstring). The {obj}`qrom_lookup <oracq.algorithms.input_model.data_loading.qrom_lookup>` on this page is the **per-address controlled-XOR unary-iteration** baseline: structurally it is Select-Swap with `partitions = 1`, serving as the λ = 1 endpoint of the query-depth-versus-ancilla trade-off curve; for the partitioned version see [Select-Swap QROM](select-swap.md).

## Interface and input model

```python
qrom_lookup(table, *, data_bits=None, name=None)
qrom_cost(n_addresses, data_bits, partitions=1)
```

API entry points: {obj}`qrom_lookup <oracq.algorithms.input_model.data_loading.qrom_lookup>`, {obj}`qrom_cost <oracq.algorithms.input_model.data_loading.qrom_cost>`

- `table`: a sequence of words or a sparse dictionary; missing addresses are treated as 0; addresses and words must be non-negative integers.
- `data_bits`: the data bit width; by default the bit width of the largest table word.
- {obj}`qrom_cost <oracq.algorithms.input_model.data_loading.qrom_cost>`: purely classical resource estimation; no circuit is generated; `n_addresses` is the number of table words N and `partitions` is the partition count λ (a power of two, no more than $2^{\lceil\log_2 N\rceil}$).

`qrom_lookup` returns {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` (`implementation="qrom_unary_iteration"`); module attributes carry the `qrom_cost` estimate: `qrom_partitions`, `qrom_address_bits`, `qrom_data_bits`, `select_toffoli`, `swap_toffoli`, `t_count`, `round_trip_t_count`, `work_qubits`, `dirty_fanout_qubits`.

`qrom_cost` returns a {obj}`QromCost <oracq.algorithms.input_model.data_loading.QromCost>` snapshot:

| Field / attribute | Meaning |
|---|---|
| `n_addresses` / `data_bits` / `partitions` / `address_bits` | echoed inputs and the address bit width |
| `select_toffoli` / `swap_toffoli` | Toffoli counts of the select and swap phases (the latter is 0 when λ = 1) |
| `t_count` | one-pass compute T count = 4 × (sum of both phases' Toffolis) |
| `t_depth` | serial depth of the select phase plus the stage-by-stage merge depth of the swap phase |
| `round_trip_t_count` | coherent uncomputation (compute + adjoint uncompute) = 2 × `t_count` |
| `work_qubits` / `fanout_qubits` / `ancilla_qubits` | data windows and low-bit fanout copies (the latter usable as dirty qubits) |
| `to_dict()` | structured dictionary for catalog and backend reports |

## Implementation notes

The baseline structure reuses `oracles.gate_database`: data bits are flipped controlled per address, with the unary-iteration tree giving `select_toffoli = 2^{k} - 1` (k is the number of address bits); at λ = 1 the T count is $4(N-1)$ (matching the formula when N is a power of two). The only differences from {obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>` are the annotation and the attached cost attributes: `implementation="qrom_unary_iteration"` lets downstream code distinguish "baseline witness" from "resource-optimized implementation".

The data table is expanded at generation time into a gate-level sub-database and does not enter the IR JSON; resource estimation is given purely classically by `qrom_cost`, and the structural attributes of the generated program agree with it (checked in both directions; see Validation approach). Applicability boundary: the gate count grows linearly with N; for large tables, switch to {obj}`select_swap_qrom <oracq.algorithms.input_model.data_loading.select_swap_qrom>` to trade ancilla qubits for depth.

## Validation approach

Category C5 (data-access layer; acceptance criteria in `../development/validation-plan.md` §2). Three layers of evidence in `tests/core/test_data_loading.py:DataLoadingTests`:

- Structural: `test_generated_operations_carry_cost_attributes` — attributes such as `implementation` / `qrom_partitions` / `t_count` agree with the `qrom_cost` formulas (baseline λ = 1).
- Numerical: `test_result_invariant_across_partitions_and_qrom_lookup_baseline` — the per-address readouts of the 16-address table equal the table values exactly and serve as the common baseline for all Select-Swap partitions; `test_cost_model_matches_formulas_and_tradeoff_curve` validates the formula values for N = 16, b = 3 (t_count = 4·15 at λ = 1) and the closed-form curve $4(N/\lambda - 1 + b(\lambda-1))$; `test_cost_report_is_structured` validates `to_dict` (N = 1024, b = 32, λ = 8: t_count = 4·(127 + 224), 277 ancilla qubits).
- Binding: the gate / QRAM bindings are covered indirectly by `test_result_invariant_across_partitions_and_qrom_lookup_baseline` (the same query semantics give invariant results across implementations / partitions); three-binding consistency parameterization is pending for V4.
- Negative cases: `test_invalid_inputs_fail_at_generation` — empty tables, negative values, out-of-range word widths, partitions not a power of two, etc. raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time.

## Known gaps and planned stages

Three-binding parameterization (V4): unified parameterized cross-checks of the gate / QRAM bindings of the same declaration are not rolled out. Consistent with the data_loading.py row of `validation-coverage.md`.

## Related links

- Source: `src/oracq/algorithms/input_model/data_loading.py`
- Same-group pages: [Select-Swap QROM](select-swap.md), [XOR database](xor-database.md)
- API reference: [Select-Swap QROM data loading](../../api/algorithms/input_model/data_loading.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

The results in this section were produced by `tests/verification/verify_stateprep.py` on real backends (four independent paths: reference, rir-pysparq, adapter-pysparq, and originir-ext + UniQC).

**Experiment design** (14 cases): truth-table exhaustion instead of per-basis-state loops — after Hadamards on the address register, a single query reads out the whole table in parallel (expected state $\frac{1}{\sqrt{N}}\sum_a |a, T[a]\rangle$), and a full address+data superposition further exhausts the XOR semantics $|a, d\rangle \mapsto |a, d \oplus T[a]\rangle$. Instances: a 16-address dense table and a sparse dictionary table with missing addresses treated as 0, plus a 64-address wide table for `qrom_lookup`; the 16-address table of `select_swap_qrom` is exhausted amplitude by amplitude over the full partition curve $\lambda \in \{1, 2, 4, 8, 16\}$, plus a 64-address table at $\lambda \in \{2, 8\}$; {obj}`qram_database <oracq.algorithms.input_model.oracles.qram_database>` (QRAM-resource binding, table data supplied as memory) runs the same two superposition exhaustions. The workspaces of $\lambda = 8/16$ and of the 64-address $\lambda = 8$ are 31/42/55 qubits respectively, exceeding the 24-qubit budget of OriginIR-ext (predicted by {obj}`workspace_table <oracq.infrastructure.layout.workspace_table>`); these cases run only the pysparq paths, as noted in the case parameters. `qrom_cost` is cross-checked bidirectionally against the paper's closed-form formulas $4(N/\lambda - 1 + b(\lambda-1))$ for T count, $N/\lambda - 1 + (\lambda - 1)$ for T depth, and $\lambda b + (\lambda-1)\log_2\lambda$ for ancilla qubits on an $(N, b, \lambda)$ grid, and the module attributes of the generated operations are validated against the estimates.

**Key metrics**:

| Case | Scale | Paths | max_error |
|---|---|---|---|
| qrom-lookup-truth-table-n16 | N=16, b=3 | four paths | 8.3e-17 |
| qrom-lookup-sparse-dict | N=8, b=3 (sparse table) | four paths | 5.6e-17 |
| qrom-lookup-wide-n64 | N=64, b=4 | four paths | 5.6e-17 |
| select-swap-truth-table-l1..l16 | N=16, b=3, λ∈{1,2,4,8,16} | ref+rir+adapter (originir added for λ≤4) | ≤ 8.4e-17 |
| select-swap-xor-superposition-l4 | N=16, b=3, 128 branches | four paths | 4.2e-17 |
| select-swap-wide-n64-l2 / l8 | N=64, b=4 | ref+rir+adapter (originir added for l2) | ≤ 5.6e-17 |
| qrom-cost-formula | 9-point (N,b,λ) grid + generated attributes | purely classical | failures = 0 |
| qram-database-truth-table / xor-superposition | N=4, b=3 (QRAM binding) | four paths | ≤ 1.2e-16 |

For all partitions the query results agree with the truth table amplitude by amplitude at machine precision, the local window registers are uncomputed (the simulator enforces this check at LocalExit), and the cost model agrees exactly with the closed-form formulas.

**Reproduction command**:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_stateprep.py
```

**Artifact path**: `out/verification/stateprep.json`.
