# Szegedy Walk

**English** · <a href="../../zh/manual/algorithms/szegedy-walk.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.input_model.graph_walks`](../../api/algorithms/input_model/graph_walks.rst) · Stage V3

## Overview

Constructs the Szegedy quantized walk step on the bipartite walk space `(current, peer, index)`

$$
W = R_B R_A, \qquad R_A = A\,(2|0\rangle\langle 0| - I)\,A^\dagger,
$$

where the zero reflection acts on the `peer‖index` joint register, and the neighbor state preparation $A = \mathrm{Adj} \cdot H^{\otimes g}$ maps $|v\rangle|0\rangle|0\rangle$ to $\frac{1}{\sqrt D}\sum_j |v\rangle|N(v,j)\rangle|j\rangle$; $R_A$ is the reflection about $\mathrm{span}\{A|v\rangle|0\rangle|0\rangle\}$, and $R_B$ swaps the roles of the two ends (equivalent to the same reflection conjugated by SWAP). The implementation basis is Szegedy 2004's quantized construction of the walk operator spectrum (cited in the module docstring); the marked-vertex search skeleton is covered in [MNRS Quantum Walk Search](mnrs-search.md).

## Interface and input model

```python
szegedy_walk(adjacency, *, name=None)
szegedy_setup(adjacency, *, name=None)
```

API entries: {obj}`szegedy_walk <oracq.algorithms.input_model.graph_walks.szegedy_walk>`, {obj}`szegedy_setup <oracq.algorithms.input_model.graph_walks.szegedy_setup>`

- `adjacency`: an [Adjacency Oracle](adjacency-oracle.md) ({obj}`AdjacencyOracle <oracq.algorithms.input_model.graph_walks.AdjacencyOracle>` or a bare operation adapted through {obj}`as_adjacency <oracq.algorithms.input_model.graph_walks.as_adjacency>`); the input model is FO.
- {obj}`szegedy_walk <oracq.algorithms.input_model.graph_walks.szegedy_walk>` returns the walk step {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `current: Bits(v)`, `peer: Bits(v)`, `index: Bits(g)`. Attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"szegedy_walk"` |
| `vertex_bits` / `degree_bits` | vertex and index bit widths |
| `composition` | `"R_B.R_A"` |

- {obj}`szegedy_setup <oracq.algorithms.input_model.graph_walks.szegedy_setup>` returns a {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` that prepares the initial state $\frac{1}{\sqrt N}\sum_v A|v\rangle$: `target: Bits(2v+g)` (from least to most significant: current, peer, index) and `work: Bits(0)`; labeled `state_prep_isometry`, `zero_input=True`, `clean_work=True`, `implementation="uniform_vertex_plus_neighbor"`.

## Implementation notes

Each reflection is assembled in the order Adj → H → {obj}`reflect_zero <oracq.algorithms.input_model.block_encoding.reflect_zero>` ($2|0\rangle\langle 0| - I$) → H → Adj: `Adj` is self-inverse and $A$ is overall unitary, so the sequence realizes exactly $R = A(2|0\rangle\langle 0| - I)A^\dagger$; the `index` register stays inside the walk space and is not reused as a workspace to be cleaned up. $R_B$ reuses the same subroutine with the current/peer roles swapped. The graph need not be regular: for non-regular graphs the effective transition is determined by how the neighbor table is padded (see the adjacency oracle page).

The `target` layout of `szegedy_setup` matches the register declaration order of `szegedy_walk`, and {obj}`quantum_walk_search <oracq.algorithms.input_model.graph_walks.quantum_walk_search>` uses this to map the walk registers onto contiguous slices of `target`. Applicability boundary: this module provides only the single-step operator and the initial-state preparation; detecting marked vertices requires the MNRS skeleton.

## Validation approach

Category C1 (exact discrete semantics; acceptance criteria in `../../development/validation-plan.md` §2): standard witness techniques include unitarity identities and fixed-point structure. Evidence:

- Structure: `tests/core/test_graph_walks.py:SzegedyWalkTests.test_walk_step_is_unitary` — W followed by W† restores the original basis state on three initial basis vectors (amplitude 1, places = 10), covering the C1 unitarity identity $W^\dagger W = I$.
- Numerical: `SzegedyWalkTests.test_stationary_state_is_fixed_point` — on the alternating edge-coloring neighbor table of an even cycle (satisfying the involution property $N(N(v,j),j) = v$), the uniform stationary state prepared by `szegedy_setup` is unchanged amplitude by amplitude after one walk step (places = 10).
- Binding: no independent binding witness; end-to-end abstract → gate binding evidence is covered by `QuantumWalkSearchTests.test_abstract_handles_bind_to_gate_implementations` (see [MNRS Quantum Walk Search](mnrs-search.md)).

## Known gaps and planned stages

Consistent with the validation matrix: new witnesses are needed after Johnson graphs (element distinctness), stage V3. That search on even cycles does not amplify is a theoretical result; see the discussion in [MNRS Quantum Walk Search](mnrs-search.md).

## Related links

- Same module: [Adjacency Oracle](adjacency-oracle.md), [MNRS Quantum Walk Search](mnrs-search.md), [Coined Cycle Walk](coined-cycle-walk.md)
- Source: `src/oracq/algorithms/input_model/graph_walks.py`
- API reference: [Graph walk search](../../api/algorithms/input_model/graph_walks.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Experiment design: an even cycle on 8 vertices (alternating edge-coloring neighbor table, satisfying the involution property $N(N(v,j),j)=v$, $v=3, g=1$, 7-qubit walk space). Two levels of comparison, with the classical reference being a $W=R_B R_A$ assembled independently in numpy ($R_A, R_B$ are reflections about neighbor superpositions, constructed directly from the neighbor table):

1. Unitary level: the OriginIR-ext export goes through UniQC `Circuit.to_matrix` to obtain the full unitary, compared element by element with the reference matrix;
2. State level: 1 and 3 walk steps applied to the `szegedy_setup` initial state, compared amplitude by amplitude across the four backend paths (reference, rir-pysparq, adapter-pysparq, originir-ext).

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| szegedy-walk-cycle8 | 7 qubits, steps=1,3 | 4 paths + to_matrix | max unitary-matrix error | 4.1e-16 |
| ditto | ditto | ditto | max evolved-state error | 5.3e-16 |

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

Artifacts: `out/verification/search_walks.json` (case `szegedy-walk-cycle8`; the MNRS end-to-end search is covered in the Numerical validation section of [MNRS Quantum Walk Search](mnrs-search.md)).
