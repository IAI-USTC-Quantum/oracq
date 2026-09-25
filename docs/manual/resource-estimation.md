# Resource estimation (Toffoli+Clifford+T+QRAM)

**English** · [简体中文](../zh/manual/resource-estimation.html)

> Modules [`oracq.infrastructure.estimate`](../api/infrastructure/estimate.rst) and [`oracq.infrastructure.backends.strict`](../api/infrastructure/backends/strict.rst)

## Overview

Resource estimation answers "how much will this program cost at the fault-tolerant level": the RIR program is counted along the real OriginIR-ext → basis lowering chain down to the **Toffoli + Clifford + T + QRAM** level, outputting the qubit count, the Toffoli count, the Clifford atom count, the exact T count, the number of rotation atoms awaiting synthesis (with a Ross–Selinger-style T-cost estimate), and the **QRAM query count**. The estimator is compositional: aggregation is memoized by `(module, control count)`, and {obj}`Repeat <oracq.infrastructure.ir.Repeat>` multiplies counts symbolically — a program repeated `2^40` times can be estimated without expansion, a capability bought directly by the RIR's "structure-preserving" design.

On the same lowering chain, {obj}`export_strict(program) <oracq.infrastructure.backends.strict.export_strict>` yields **checkable compilation artifacts**: the strict netlist preserved module by module (TOFFOLI/CZ/H/S/SDG/T/TDG/X/Y/Z atoms + `RZ`/`RY` rotations awaiting synthesis + QRAM lines), with `DEF` calls, `Repeat` symbols, and QRAM declarations left unexpanded.

## Interface

Open programs use {obj}`estimate_resources(program, require_closed=False) <oracq.infrastructure.estimate.estimate_resources>`. The report lists known circuit costs and `oracle_calls` separately; each entry records the slot name, control count, adjoint marker, resource arguments, and call count. Re-estimating after binding yields the concrete implementation cost. `complete=False` means unimplemented dependencies remain, so the known gate counts are not the total program cost.

For open programs `qubits=None`; `qubits_lower_bound` is the lower bound from known registers and private workspaces, and `unknown_workspace` lists the unimplemented modules. MCX ancillas likewise cover only known implementations. Complete programs keep the original integer results. QRAM reads and writes are accounted per entry logical resource, with formal parameters substituted along call edges; distinct banks are not merged into an internal parameter of the same name. The same analysis also has a command-line entry; see [open and closed resource analysis](cli.md#open-and-closed-resource-analysis).

`rotations` is now the compact {obj}`RotationCounts <oracq.infrastructure.estimate.RotationCounts>`: `rotations.total` gives the total of arbitrary size, and `rotations.counts` stores multiplicities keyed by `(axis, angle)`. Read-only `len`, indexing, and lazy iteration are retained so small programs can still be read as a sequence; it is not a mutable list and not an execution order. Counts and JSON reports for huge Repeats do not expand the rotation list.

`epsilon` is the synthesis configuration per rotation, not an overall error guarantee for the whole scientific-computing algorithm. The resource report does not fill in default zeros for QRAM hardware, data loading, or physical error-correction costs.

```python
from oracq import estimate_resources, export_strict

estimate = estimate_resources(program)          # or operation.estimate()
estimate.qubits            # entry registers + module workspaces (excluding the MCX pool bits)
estimate.mcx_ancilla       # basis pool bits (max(0, max control count - 1))
estimate.toffoli           # total Toffolis
estimate.clifford          # total H/S/SDG/X/Y/Z/CZ atoms
estimate.t_exact           # exact T/TDG count (rotations on the π/4 grid)
estimate.rotations         # rotation atoms awaiting synthesis [(axis, angle), ...]
estimate.t_total(1e-10)    # exact T + rotations × ceil(3·log2(1/eps))
estimate.qram_queries      # queries per resource; estimate.qram_total is the total
estimate.to_dict()         # JSON-friendly ledger
```

The counting model aligns with the exporter line by line (this is not an idealized cost model but the ledger of the **actually emitted circuit**): X-type multi-controls go through the `basis.mcx` recipe (2c−3 Toffolis for c≥3); controlled single-bit gates go through the `basis.controlled_u3` recipe (a 2(c−1)-Toffoli ladder plus an intermediate network for c≥2); zero-valued control bits are flipped once on each side of every emitted line group; each {obj}`Load <oracq.infrastructure.ir.Load>` records one QRAM query. Rotations whose angle is exactly an integer multiple of π/4 are counted as exact Clifford+T (the S/T/Z families); the rest are recorded as rotations awaiting synthesis. Global phases are not counted (consistent with the strict export).

## Applicability boundaries and known gaps

- The estimator counts by **the current compiler's emission recipes**: constant addition goes through the O(n³) staircase (the literature has O(n) reversible adders), and {obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>`/{obj}`qrom_lookup <oracq.algorithms.input_model.data_loading.qrom_lookup>` go through per-branch truth-table networks (an idealized unary-iteration QROM is N−1 Toffolis; see the {obj}`qrom_cost <oracq.algorithms.input_model.data_loading.qrom_cost>` model and babbush2018encoding). The gap between real emission and the idealized model is the optimization space of the compiler lowering, and the estimator turns it into a measurable quantity.
- The T cost of rotations is a leading-order estimate (≈3·log2(1/ε), Ross–Selinger style) without a concrete synthesizer implementation; physical-layer costs (T factories, distillation, routing, surface-code cycles) are out of scope at this layer.
- The QRAM query count treats QRAM as a first-class resource: it answers "how many queries the algorithm issues to the memory system", in contrast with "how many Toffolis it would take to synthesize the same table with gates" (see the table below).

## Numerical validation

Experiment design (`tests/core/test_estimate.py` cross-checks hand-computed counts against the strict netlist line by line; `tools/build_resource_estimates.py` runs ten scaling experiments whose artifacts land in `out/resource-estimates/*.json`, registered into `tools/check_project.py`):

- Hand-computed counts: staircase addition, controlled-gate recipes, zero-valued control flips, symbolic Repeat (2^40) multiplication, and per-resource QRAM queries — all 15 cases pass; the per-atom counts of `estimate_resources` and `export_strict` agree item by item on flattened programs.
- Scaling experiments (log-log fitted slopes vs theoretical expectations):

| Group | Scale | Key reading | Fit/expected |
|---|---|---|---|
| `add_const` | n=4..64 | Toffoli 81 375 at n=64 | staircase ≈3.5 (cubic regime); literature optimum O(n) |
| `fixed_mul` | n=4..32 | Toffoli 112→6 356 | 1.94 vs O(n²) |
| {obj}`qft <oracq.algorithms.common.fourier.qft>` | n=3..12 | rotations 3→165, exact T 6→33 | ≈2.85 (exact rotations on the π/4 grid drop out at small n) |
| `qpe_add_const` | p=3..12 | ~55,000 gates at p=12 | symbolic Repeat: ~2^p × cost(add) |
| {obj}`grover <oracq.algorithms.common.search.grover>` | n=4..14 | iterations ⌊π/4·2^{n/2}⌋ | total gates ~2^{n/2}·O(n) |
| `state_prep_dense` | n=2..10 | rotations 17→5 114 | O(2^n) |
| `qram_lookup` | n=4..16 | **queries = 1**, gates = 0 | independent of size |
| `qrom_lookup_gate` | n=4..12 | Toffoli 344 064 at n=12 | real emission vs the `qrom_cost` ideal model 4 095 (≈84×) |
| `diagonal_be_qram` | n=2..12 | **queries = 2**, rotations = 2 | constant |
| `diagonal_be_gate` | n=2..8 | Toffoli 11 388 at n=8 | gate-level counterpart of the same input model |
| {obj}`roe_face <oracq.applications.roe.roe_face>` / `math_polynomial` | w=6..16 | Roe face Toffoli 446 076 at w=16 | ≈1.8, polynomial in word length |

- A documentation-level gap exposed: the `qrom_unary_iteration` cost attribute of `qrom_lookup` describes the idealized model while the current emission is a truth-table network — recorded as a lowering optimization item (alongside the two OAA and Schrödingerization records, a validation-driven finding).

Reproduce: `PYTHONPATH=src python tools/build_resource_estimates.py`; `python -m unittest tests.core.test_estimate -v`.
