# QAOA for MaxCut

**English** · <a href="../../zh/manual/algorithms/qaoa-maxcut.html">简体中文</a>

> Category C4 · Module [`oracq.algorithms.optimization.variational`](../../api/algorithms/optimization/variational.rst) · Stage V3

## Overview

Given an undirected weighted graph $G=(V,E)$, MaxCut asks for a partition of the vertices into two parts that maximizes the sum of the weights of the edges whose endpoints lie in different parts. QAOA starts from the uniform superposition and alternately applies cost evolution and mixer evolution, with the layer count and angles $(\gamma,\beta)$ chosen by a classical outer loop; the cost/mixer layering follows the original algorithm of Farhi et al. ([arXiv:1411.4028](https://arxiv.org/abs/1411.4028), already cited under "Implementation basis" in the manual's algorithm catalog). The cost operator takes the value

$$
C(z) = \sum_{(u,v,w)\in E} \frac{w}{2}\,(1 - z_u z_v),
$$

on bit strings $z$ — the sum of the weights of the cut edges; the mixer is the bitwise $X$. This function generates the quantum circuit for fixed angles and reads out `target` to obtain a candidate bit string for a cut; it runs no classical optimizer.

## Interface and input model

```python
qaoa_maxcut(width, edges, gammas, betas)
```

API entry point: {obj}`qaoa_maxcut <oracq.algorithms.optimization.variational.qaoa_maxcut>`

- `width`: the number of vertices of the graph, also the target bit width, range 1..64.
- `edges`: `(u, v, weight)` triples; `u` and `v` must lie in `[0, width)`, the weights must be non-negative finite real numbers, and self-loops are rejected.
- `gammas` / `betas`: the per-layer cost / mixer evolution angles; the two lists must be non-empty and of equal length, and the angles must be finite real numbers.

The graph's edge list and the variational angles are both classical data parameterized directly, so the input model is CP. Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` whose only register is `target` (width bits), with no work bits. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"qaoa_maxcut"` |
| `layers` | the QAOA layer count, i.e. `len(gammas)` |

## Implementation notes

Initialization is the Hadamard uniform superposition over the whole `target` register. Each QAOA layer applies cost evolution first, then mixer evolution:

- the cost layer applies an XOR–`Rz(-gamma·weight)`–XOR conjugation per edge, writing the phase onto the $Z_u Z_v$ component; together with the global phase $-\gamma w/2$, the single-edge gate is exactly $e^{-i\gamma w(1-Z_uZ_v)/2}$; the per-edge ZZ gates commute with one another, and the whole layer implements $e^{-i\gamma C}$.
- the mixer layer applies `rx(2·beta)` to every qubit.

Register layout: vertex i corresponds to bit i of `target`. All parameters are validated at generation time (width cap, self-loops, negative weights, and inconsistent or empty gamma/beta lengths all raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`). Applicability boundary: only the circuit for the given angles is generated; angle initialization, expected-cut-value evaluation, and the outer optimization loop are organized by the caller; this module contains no classical optimizer.

## Validation approach

Category C4 (heuristic/optimization semantics; acceptance criteria in `../../development/validation-plan.md` §2: strictly better than the random baseline, and reaching the known optimum / theoretical fraction on small instances). Three layers of evidence:

- Structure: construction and attribute assertions in `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`; the same family's `test_bad_inputs_fail_at_generation` covers generation-time violations such as the self-loop edge `(0, 0, 1)`.
- Numerical: `AlgorithmExpansionTests.test_qaoa_single_edge_optimal_layer` — a single-layer QAOA with $\gamma=\pi/2$ and $\beta=\pi/8$ on the 2-vertex single-edge graph (weight 1.0), where the summed probability of the two optimal cut states `|01>` and `|10>` equals 1 (places = 10), i.e. the single layer already attains the optimal cut.
- Binding: this module has no separate binding witness (the inputs are classical parameters, with no abstract oracle slots).

## Known gaps and planned stages

Consistent with the `variational.py` row of `validation-coverage.md`: the strong witness for "reaching the optimal cut on small instances" is missing — currently there is only the analytically optimal layer for the single-edge graph, and a distribution-level cross-check against the exhaustively enumerated optimal cut value on multi-vertex graphs has not been done. Integration is planned for V3 (registering the new algorithm in the catalog + extending the real-backend cross-checks in tests/integration).

## Related links

- Source: `src/oracq/algorithms/optimization/variational.py`
- Same-module pages: [Hardware-Efficient Ansatz](variational-ansatz.md), [VQE Pauli measurements](vqe.md)
- API reference: [variational algorithm circuits](../../api/algorithms/optimization/variational.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_misc_algorithms.py` (misc_algorithms group), all executed on real backends.

**Experimental design**: (a) the analytic anchor $\gamma = \pi/2$, $\beta = \pi/8$ of a single-layer QAOA on the single-edge graph (2 vertices, unit weight); (b) a single-layer QAOA on the 4-vertex cycle C4 (unit weights) — the angles are selected by a grid search over an exact numpy simulation ($\gamma = 0.85$, $\beta = 0.45$; the classical outer loop takes care of itself per the module contract), and the quantum distribution is cross-checked pointwise against the exact numpy simulation (diagonal cost phase + per-bit mixer $e^{-i\beta X}$); the optimal cuts are $\lvert 0101\rangle$ and $\lvert 1010\rangle$ (random baseline 2/16 = 0.125). Backend paths: full-amplitude cross-check on `reference` and `originir-ext`.

**Key metrics**:

| Case | Scale | Path | Distribution TVD | Optimal-cut probability | Relative to random baseline |
|---|---|---|---|---|---|
| qaoa-single-edge-optimal | 2 qubits, 1 layer | reference | — | 1.0000 (analytic optimum) | — |
| qaoa-c4-distribution | 4 qubits, 1 layer | reference + originir-ext | 1.8e-16 | 0.5379 | 4.30× |

**Reproduce**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

Artifact: `out/verification/misc_algorithms.json` (all 24 cases pass; this page corresponds to the two `qaoa-*` cases).
