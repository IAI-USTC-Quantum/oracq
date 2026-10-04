# Adjacency Oracle

**English** · <a href="../../zh/manual/algorithms/adjacency-oracle.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.input_model.graph_walks`](../../api/algorithms/input_model/graph_walks.rst) · Stage V3

## Overview

The input model of the graph-walk framework: the graph is accessed through a neighbor-table oracle with query semantics

$$
|v\rangle|j\rangle|c\rangle \mapsto |v\rangle|j\rangle|c \oplus N(v,j)\rangle,
$$

that is, given a vertex $v$ and an outgoing-edge index $j$, the index of the $j$-th neighbor is XORed into the `neighbor` register. Let $g$ be the outgoing-edge index bit width and $D = 2^g$; the outgoing-edge table of every vertex is uniformly padded to exactly $D$ columns (usually with self-loops), so that the neighbor superposition $\frac{1}{\sqrt D}\sum_j |N(v,j)\rangle|j\rangle$ normalizes without knowing the concrete degree. This input model (FO, the XOR query paradigm `database_xor`) serves both the [Szegedy quantum walk](szegedy-walk.md) and the [MNRS quantum walk search](mnrs-search.md).

## Interface and input model

The core type {obj}`AdjacencyOracle <oracq.algorithms.input_model.graph_walks.AdjacencyOracle>` (a frozen dataclass, subclass of {obj}`OracleView <oracq.algorithms.input_model.contracts.OracleView>`) wraps an operation with the `(vertex, index, neighbor)` signature and reuses the `database_xor` paradigm; the signature and the equal widths of `vertex`/`neighbor` are validated at construction. Three construction entry points:

```python
abstract_adjacency(vertex_bits, degree_bits, *, name=None)   # open declaration
gate_adjacency(neighbors, *, name=None)                      # gate-level controlled-X network
qram_adjacency(vertex_bits, degree_bits, *, name=None)       # QRAM table
```

API entry points: {obj}`abstract_adjacency <oracq.algorithms.input_model.graph_walks.abstract_adjacency>`, {obj}`gate_adjacency <oracq.algorithms.input_model.graph_walks.gate_adjacency>`, {obj}`qram_adjacency <oracq.algorithms.input_model.graph_walks.qram_adjacency>`

- {obj}`abstract_adjacency <oracq.algorithms.input_model.graph_walks.abstract_adjacency>`: an open declaration with an empty body; `vertex_bits` is at most 32 and `degree_bits` ranges over 0..32; gate/QRAM implementations can be bound via {obj}`bind <oracq.infrastructure.linking.bind>`.
- {obj}`gate_adjacency <oracq.algorithms.input_model.graph_walks.gate_adjacency>`: a gate-level implementation for small neighbor tables; the table must be a nonempty rectangle whose entries are valid vertex indices, and the bit widths are derived automatically from the table size.
- {obj}`qram_adjacency <oracq.algorithms.input_model.graph_walks.qram_adjacency>`: the QRAM implementation, with table entries addressed by `vertex | (index << vertex_bits)`.
- {obj}`as_adjacency(value) <oracq.algorithms.input_model.graph_walks.as_adjacency>` adapts a bare {obj}`Operation <oracq.infrastructure.builder.Operation>` into an `AdjacencyOracle`.

Attributes of `AdjacencyOracle`:

| Attribute | Meaning |
|---|---|
| `vertex_bits` / `degree_bits` | vertex and outgoing-edge index bit widths $v$ / $g$ |
| `vertices` / `degree` | $2^v$ and $2^g$ |
| `xor_database()` | adapts to {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>`: `address = vertex‖index` (width $v+g$), `data = neighbor` |

## Implementation notes

`gate_adjacency` enumerates all $(v, j)$ branches to build a controlled-X network, emitting gates only for nonzero target bits, so on any table the query semantics agree pointwise with the truth values. `qram_adjacency` is a single QRAM Load table lookup. Padding is the caller's responsibility: the framework only requires each row to have exactly $D$ columns; in-row repetition (self-loops) does not change the XOR semantics but determines the effective transitions of the Szegedy walk. `xor_database()` plugs graph input into the `XorDatabase` ecosystem, reusing QRAM bindings and data-loading utilities.

Applicability boundary: the gate count of `gate_adjacency` grows with the table size, so it targets small graphs and tests; larger graphs use `qram_adjacency` or stay as abstract declarations with deferred binding.

## Validation approach

Category C1 (exact discrete semantics; acceptance criteria in `../../development/validation-plan.md` §2): the query semantics must equal the neighbor table pointwise. Evidence:

- Structural: `tests/core/test_graph_walks.py:AdjacencyOracleTests`; `test_bad_tables_fail_at_generation` covers generation-time violations such as non-rectangular tables, out-of-range vertices, and divergent hitting times caused by unreachable marked sets.
- Numerical: `AdjacencyOracleTests.test_gate_adjacency_implements_neighbor_table` — for the Q3 hypercube neighbor table (8 vertices × 4 columns, self-loop-padded to $D=4$), all 32 $(v, j)$ queries are exhausted and the readout amplitude is exactly 1 (places = 10); `test_xor_database_adapter` validates the adapted address/data widths and the queried values.
- Binding: `AdjacencyOracleTests.test_qram_adjacency_matches_gate_implementation` (the gate and QRAM implementations cross-checked pointwise) and `QuantumWalkSearchTests.test_abstract_handles_bind_to_gate_implementations` (after the abstract declaration is bound to a gate implementation, the end-to-end semantics are unchanged; see [MNRS quantum walk search](mnrs-search.md)) together form the gate/qram dual-binding evidence for the adjacency oracle.

## Known gaps and planned stages

Consistent with the validation matrix: there is no witness on Johnson graphs (a dependency of element distinctness) yet; new witnesses must be added before it can be wired in. Stage V3.

## Related links

- Same module: [Szegedy quantum walk](szegedy-walk.md), [MNRS quantum walk search](mnrs-search.md), [coined cycle walk](coined-cycle-walk.md)
- Source: `src/oracq/algorithms/input_model/graph_walks.py`
- API reference: [Graph-walk search](../../api/algorithms/input_model/graph_walks.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Experiment design (end-to-end walk numerics; the adjacency query itself is already covered by unit tests): taking `gate_adjacency` on the hypercube Q3 neighbor table (8 vertices, degree 4) as an example, H gates are applied to the vertex/index registers for a single superposition call that exhausts all 32 $(v,j)$ query branches, and four backend paths (reference, rir-pysparq, adapter-pysparq, originir-ext) verify the $\lvert v\rangle\lvert j\rangle\lvert N(v,j)\rangle$ structure amplitude by amplitude. Additionally, a QRAM binding (`qram_adjacency` + memory table) runs the MNRS search end to end (K4, steps=2; OriginIR-ext has no QRAM resource, so this instance uses only the reference and rir-pysparq paths).

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| adjacency-superposition-hypercube | 8 vertices × 4 columns = 32 queries | 4 paths | max amplitude error / off-table branch weight | 5.6e-17 / 0.0 |
| mnrs-search-qram-k4 | N=4, D=4, steps=2 | reference, rir-pysparq | marked probability / final-state max error | 0.578125 / 6.3e-16 |

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

Artifact: `out/verification/search_walks.json` (cases `adjacency-superposition-hypercube`, `mnrs-search-qram-k4`).
