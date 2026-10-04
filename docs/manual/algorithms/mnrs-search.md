# MNRS Quantum Walk Search

**English** · <a href="../../zh/manual/algorithms/mnrs-search.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.input_model.graph_walks`](../../api/algorithms/input_model/graph_walks.rst) · Stage V3

## Overview

A quantum-walk skeleton for marked-vertex search: `setup` prepares the initial state $\frac{1}{\sqrt N}\sum_v A|v\rangle$, followed by iterating `Repeat{ walk; marked phase flip }`. The implementation basis is Magniez–Nayak–Roland–Santha 2011 (Search via quantum walk; hitting time and step scaling, cited in the module docstring): the step count takes the scaling $\lceil \frac{\pi}{4}\sqrt{H_{\text{avg}}}\rceil$ of the average hitting time $H_{\text{avg}}$ of the classical random walk; on the complete graph $H_{\text{avg}} = N/|M|$, and the step count recovers the Grover scaling $\frac{\pi}{4}\sqrt{N/|M|}$.

## Interface and input model

```python
quantum_walk_search(setup, walk, marked, steps, *, name=None)
transition_matrix(neighbors)
hitting_times(transition, marked)
suggest_steps(transition, marked)
```

API entries: {obj}`quantum_walk_search <oracq.algorithms.input_model.graph_walks.quantum_walk_search>`, {obj}`transition_matrix <oracq.algorithms.input_model.graph_walks.transition_matrix>`, {obj}`hitting_times <oracq.algorithms.input_model.graph_walks.hitting_times>`, {obj}`suggest_steps <oracq.algorithms.input_model.graph_walks.suggest_steps>`

- `setup`: the initial-state preparation; it must satisfy the {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` protocol (e.g. the output of {obj}`szegedy_setup <oracq.algorithms.input_model.graph_walks.szegedy_setup>` from [Szegedy Walk](szegedy-walk.md)).
- `walk`: a walk-step handle whose registers map, in declaration order, onto contiguous slices of `setup.target`; it may be an open declaration.
- `marked`: a phase oracle acting on the first few bits of the walk space (`target[, work]` interface, constructed by {obj}`phase_marks(width, marked) <oracq.algorithms.input_model.oracles.phase_marks>`); for {obj}`szegedy_walk <oracq.algorithms.input_model.graph_walks.szegedy_walk>` it acts precisely on the `current` register.
- `steps`: a non-negative iteration count, which may be chosen with reference to {obj}`suggest_steps <oracq.algorithms.input_model.graph_walks.suggest_steps>`.
- All three handles may be open declarations; the generated program is executed after implementations are bound in batches via {obj}`bind <oracq.infrastructure.linking.bind>`. The overall input model is FO: the graph is given through an [Adjacency Oracle](adjacency-oracle.md), and the marked set through a phase oracle.

{obj}`quantum_walk_search <oracq.algorithms.input_model.graph_walks.quantum_walk_search>` returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `target: Bits(width)` (width is the total walk-space width) and `work: Bits(prep.work_width + marked_work)`. Attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"quantum_walk_search"` |
| `framework` | `"MNRS"` |
| `steps` / `walk_width` | iteration count and walk-space width |

Classical-side tools: {obj}`transition_matrix <oracq.algorithms.input_model.graph_walks.transition_matrix>` derives $P[v][u] = |\{j: N(v,j) = u\}| / D$ from the neighbor table; {obj}`hitting_times <oracq.algorithms.input_model.graph_walks.hitting_times>` solves the linear system $(I - P_{\text{free}})h = \mathbf{1}$ to give the expected number of steps for each vertex to first reach the marked set (0 for marked vertices; raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` when unreachable); `suggest_steps` returns $\max(1, \lceil \frac{\pi}{4}\sqrt{H_{\text{avg}}}\rceil)$, or 0 when the whole graph is marked.

## Implementation notes

The generated structure is a setup call followed by one {obj}`Repeat(steps) <oracq.infrastructure.ir.Repeat>` loop: the loop body first maps the walk step's registers, in declaration order, onto contiguous slices of `target` and invokes it, then applies the marked phase flip to the first `marked_width` bits of `target`. Width consistency is validated at generation time (the total walk width equals the setup width, marked does not exceed the walk space, steps is non-negative); violations raise `ValidationError`.

Applicability boundary: the step count is supplied by the caller, and `suggest_steps` is only a reference — it is a function of the classical hitting time and is unaware of the quantum spectrum; on even cycles (alternating edge-coloring tables) this skeleton does not amplify, see the known gaps section.

## Validation approach

Category C1 (exact discrete semantics; acceptance criteria in `../../development/validation-plan.md` §2): the output distribution must cross-check against an exact construction. Evidence:

- Structure: `tests/core/test_graph_walks.py:QuantumWalkSearchTests` / `HittingTimeTests`; `test_search_rejects_mismatched_handles` covers generation-time violations such as negative steps and mismatched setup widths.
- Numerical: `QuantumWalkSearchTests.test_search_amplifies_marked_vertex_on_hypercube` — Q3 hypercube, marked {0}, `suggest_steps` gives 3 steps: the baseline (0 steps) has marked probability 0.125, and 3 steps give 0.78125 (places = 9); `test_search_recovers_grover_limit_on_complete_graph` — on the K4 complete graph the probability is exactly 1.0 after 2 steps (places = 10); `HittingTimeTests.test_cycle_hitting_times_match_closed_form` — the hitting times of every vertex on an 8-cycle cross-check pointwise against the closed form $h_k = d(8-d)$, $d = \min(k, 8-k)$ (places = 9); `test_suggest_steps_matches_grover_scaling` and `test_all_marked_needs_no_steps` validate the step-count formulas.
- Binding: `QuantumWalkSearchTests.test_abstract_handles_bind_to_gate_implementations` — both the adjacency oracle and the marked oracle are open declarations ({obj}`unresolved <oracq.infrastructure.linking.unresolved>` for {G, M}); after binding gate implementations, the marked probability matches direct generation (0.78125, places = 9); `test_walk_handle_may_stay_abstract` witnesses that the setup/walk handles can likewise stay bound later.

## Known gaps and planned stages

Consistent with the validation matrix: the element distinctness witness on Johnson graphs is pending, stage V3. Another confirmed fact: on even cycles (alternating edge-coloring neighbor tables) the marked probability stays at the uniform baseline $1/N$ and does not amplify with the step count (reference simulations rechecked several step counts on the 8-cycle, all 0.125); this has been confirmed as a theoretical result rather than an implementation defect (`validation-plan.md` §3 records this position); witnesses of amplifiable behavior are carried by the hypercube and complete-graph instances.

## Related links

- Same module: [Adjacency Oracle](adjacency-oracle.md), [Szegedy Walk](szegedy-walk.md), [Coined Cycle Walk](coined-cycle-walk.md)
- Source: `src/oracq/algorithms/input_model/graph_walks.py`
- API reference: [Graph walk search](../../api/algorithms/input_model/graph_walks.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Experiment design: three end-to-end search instances, with the step count chosen from an independent Markov-chain reference (hitting times $H_{\rm avg}$ solved via `(I-P_free)h=1`) as $\lceil\frac{\pi}{4}\sqrt{H_{\rm avg}}\rceil$; the final states are compared in full amplitude against the MNRS iteration $(M W)^s|\psi_0\rangle$ assembled independently in numpy (the walk-space reflection operators are constructed directly from the neighbor table, without reusing the tested module).

- Complete graph K4 (self-loops padding to $D=4$, neighbor table $N(v,j)=j$), marked $=\{0\}$: $H_{\rm avg}=4$, 2 steps, four paths (reference, rir-pysparq, adapter-pysparq, originir-ext).
- Hypercube Q3 ($N=8, D=4$), marked $=\{0\}$: $H_{\rm avg}\approx11.048$, 3 steps, four paths.
- K4 with QRAM binding ({obj}`qram_adjacency <oracq.algorithms.input_model.graph_walks.qram_adjacency>` + memory table): OriginIR-ext circuits carry no QRAM resources, so this instance runs only the reference and rir-pysparq paths.

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| mnrs-search-complete-k4 | N=4, D=4, steps=2 | 4 paths | marked probability / max final-state error | 0.578125 / 6.3e-16 |
| mnrs-search-hypercube-q3 | N=8, D=4, steps=3 | 4 paths | marked probability / max final-state error | 0.78125 / 1.2e-15 |
| mnrs-search-qram-k4 | N=4, D=4, steps=2 | reference, rir-pysparq | marked probability / max final-state error | 0.578125 / 6.3e-16 |

Informative metrics: $\text{steps}/\sqrt{H_{\rm avg}}=1.0$ for K4 and 0.903 for the hypercube, the same order of magnitude as $\pi/4\approx0.785$, consistent with the MNRS $\sqrt{H}$ scaling; the marked probability 0.578125 of the K4 instance corresponds to the exact walk dynamics under the $N(v,j)=j$ edge labeling (the transition matrix is the same as, but the walk different from, the 1.0 instance with cyclic edge labeling in the core unit tests; each agrees amplitude by amplitude with its independent reference).

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

Artifacts: `out/verification/search_walks.json` (cases `mnrs-search-complete-k4`, `mnrs-search-hypercube-q3`, `mnrs-search-qram-k4`, plus the classical-reference routine cross-check case `markov-chain-helpers`).
