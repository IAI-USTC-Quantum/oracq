# Select-Swap QROM

**English** · <a href="../../../zh/manual/algorithms/select-swap.html">简体中文</a>

> Category C5 · Module [`oracq.algorithms.input_model.data_loading`](../../api/algorithms/input_model/data_loading.rst) · Stage V4

## Overview

Select-Swap QROM is the resource-optimized implementation of QROM lookup (Low–Kliuchnikov–Schaeffer 2018; Babbush et al. 2018, cited in the module docstring): the address is split into the high part `h` (k bits) and the low part `y` (l bits), $\lambda = 2^l$ window partitions share the high address bits and load the sub-tables $T[h\cdot\lambda + i]$ in parallel, then merge through controlled swaps keyed on $y = i$, trading query depth against ancilla qubits at a T count of about $4(2^k + \lambda b)$ (b is the data bit width) (QRAM in the input-model vocabulary).

## Interface and input model

```python
select_swap_qrom(table, *, partitions, data_bits=None, name=None)
qrom_cost(n_addresses, data_bits, partitions=1)
```

API entry points: {obj}`select_swap_qrom <oracq.algorithms.input_model.data_loading.select_swap_qrom>`, {obj}`qrom_cost <oracq.algorithms.input_model.data_loading.qrom_cost>`

- `table`: a sequence of words or a sparse dictionary; missing addresses are treated as 0.
- `partitions`: the partition count λ; it must be a power of two and no more than the number of addresses $2^{\lceil\log_2 N\rceil}$.
- `data_bits`: the data bit width; by default the bit width of the largest table word.

Returns {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` (`implementation="select_swap"`); the external interface is still the two bits registers `address` / `data`, with XOR semantics exactly matching the baseline. Module attributes carry the {obj}`qrom_cost <oracq.algorithms.input_model.data_loading.qrom_cost>` estimate (field table in [QROM lookup](qrom-lookup.md)), where `work_qubits = λ·b` counts the window registers and `fanout_qubits = (λ-1)·l` the low-bit fanout copies (usable as dirty qubits).

## Implementation notes

Three phases: **select** — each of the λ clean local window registers `window_i` performs its sub-table query `window_i ^= T[h·λ+i]`; all windows share the high address bits, so the unary iteration is paid only once: $2^k - 1$ (for k = 0 the sub-table is constant and it degenerates to X gates); **merge** — controlled on `y == i`, the compound swap–xor–swap is applied to `(window_i, data)`, with net effect a controlled `data ^= window_i`, preserving XOR semantics for arbitrary data initial values without destroying the window contents; **uncompute** — an adjoint replay of the sub-queries uncomputes all windows back to zero.

For the test baseline (N = 16, b = 3), the trade-off curve given by `qrom_cost`:

| λ | select | swap | t_count | t_depth | work | fanout | ancilla qubits |
|---|---|---|---|---|---|---|---|
| 1 | 15 | 0 | 60 | 15 | 3 | 0 | 3 |
| 2 | 7 | 3 | 40 | 8 | 6 | 1 | 7 |
| 4 | 3 | 9 | 48 | 6 | 12 | 6 | 18 |
| 8 | 1 | 21 | 88 | 8 | 24 | 21 | 45 |
| 16 | 0 | 45 | 180 | 15 | 48 | 60 | 108 |

t_count reaches its valley at λ = 4 (the two inequalities asserted by the test `test_cost_model_matches_formulas_and_tradeoff_curve`), and t_depth is likewise minimal at intermediate partitions; at λ = 16 the high bits are exhausted and the select phase degenerates into pure swaps. Applicability boundary: increasing λ lowers the T depth but linearly increases ancilla qubits; the optimal λ depends on the hardware's qubit / depth exchange rate.

## Validation approach

Category C5 (data-access layer; acceptance criteria in `../development/validation-plan.md` §2). Three layers of evidence in `tests/core/test_data_loading.py:DataLoadingTests` (16-address table `TABLE16`, b = 3):

- Structural: `test_generated_operations_carry_cost_attributes` — attributes such as `select_toffoli` / `swap_toffoli` / `t_count` / `work_qubits` / `dirty_fanout_qubits` agree item by item with `qrom_cost(16, 3, 4)`.
- Numerical: `test_select_swap_matches_gate_database_per_address` — for λ ∈ {1, 2, 4, 8, 16}, the readouts of all 16 addresses equal the {obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>` baseline pointwise; `test_xor_semantics_with_nonzero_data_register` (with data initial value 5 the output is `T[a] XOR 5`); `test_superposition_query_restores_clean_work` (after a query on a uniform address superposition the marginal distribution is uniform, places = 12; the windows are forced clean by LocalExit); `test_cost_model_matches_formulas_and_tradeoff_curve` — λ = 4 sits at the bottom of the curve (curve[4] < curve[1] = 4·15 and < curve[16] = 4·45), and the curve is asserted to agree pointwise with the closed form $4(N/\lambda - 1 + b(\lambda - 1))$.
- Binding: the gate / QRAM bindings are covered indirectly by `test_result_invariant_across_partitions_and_qrom_lookup_baseline` (results do not change with partitions or implementations); three-binding consistency parameterization is pending for V4.
- Negative cases: `test_invalid_inputs_fail_at_generation` — partitions not a power of two / out of range, empty tables, negative values, out-of-range word widths, etc. raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time.

## Known gaps and planned stages

Three-binding parameterization (V4): unified parameterized cross-checks of the gate / QRAM bindings of the same declaration are not rolled out; validation-plan §5 additionally suggests bringing `qrom_cost` into catalog-report attributes. Consistent with the data_loading.py row of `validation-coverage.md`.

## Related links

- Source: `src/oracq/algorithms/input_model/data_loading.py`
- Same-group pages: [QROM lookup](qrom-lookup.md), [XOR database](xor-database.md)
- API reference: [Select-Swap QROM data loading](../../api/algorithms/input_model/data_loading.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-level numerical validation script: `tests/verification/verify_blockencoding.py` (executed on real backends, no mocks, no skips; as of 2026-09-16 all 73 cases pass), artifact `out/verification/blockencoding.json`. This page corresponds to 6 cases: `select-swap-lambda{1,2,4,8,16}` and `select-swap-as-sparse-entry`.

Experiment design: 16 addresses, 3-bit random data words (fixed seed), with the partition count λ scanned over {1, 2, 4, 8, 16}; for each λ, all 16 addresses × two data initial values (0 and the nonzero 5, checking XOR semantics) are read out per basis state and cross-checked amplitude by amplitude against the `gate_database` truth-table baseline (reference path). End-to-end case: select_swap (λ = 2) serves as the element database of the CKS sparse block encoding, assembling the $T^\dagger S T$ block encoding of a 2×2 signed sparse matrix, whose zero-signal block, after multiplying by α, is cross-checked element by element against the classical matrix (see the numerical validation section of [sparse matrix block encoding](sparse-block-encoding.md)).

| Case | Scale | Backend paths | Metric | Value |
|---|---|---|---|---|
| `select-swap-lambda1` | 16 addresses × 2 initial values | reference | max deviation from baseline | 0 |
| `select-swap-lambda2` | ditto | reference | ditto | 0 |
| `select-swap-lambda4` | ditto | reference | ditto | 0 |
| `select-swap-lambda8` | ditto | reference | ditto | 0 |
| `select-swap-lambda16` | ditto | reference | ditto | 0 |
| `select-swap-as-sparse-entry` | dim 2, s=2, α=3.0 | reference | sparse BE block max_error | 4.4e-16 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_blockencoding.py
```

Artifact: `out/verification/blockencoding.json`.
