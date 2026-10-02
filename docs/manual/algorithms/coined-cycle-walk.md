# Coined Cycle Walk

**English** · <a href="../../../zh/manual/algorithms/coined-cycle-walk.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.common.walks`](../../api/algorithms/common/walks.rst) · Stage V1

## Overview

A minimal instance of the discrete-time coined quantum walk: on the ring lattice $\mathbb{Z}_{2^n}$ with period $2^n$, a one-bit Hadamard coin decides whether each step moves $+1$ or $-1$. The coined walk is one of the standard models of quantum walks, and this module provides its one-dimensional periodic version; for walks on arbitrary graphs driven by oracle input see [Szegedy Walk](szegedy-walk.md).

## Interface and input model

```python
cycle_walk(width, *, steps=1)
```

API entry point: {obj}`cycle_walk <oracq.algorithms.common.walks.cycle_walk>`

- `width`: the `position` bit width $n$; the period length is $2^n$, capped at 64.
- `steps`: a non-negative number of steps.
- The input model is CP: the bit width and the step count are classical parameters, with no oracle input; the walk-space initial state is prepared by the caller, with the all-zero input corresponding to position 0 and coin 0.

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `position: Bits(n)` and `coin: Bits(1)`. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"coined_cycle_walk"` |
| `steps` | echo of the call parameter |

## Implementation notes

Each step updates the coin first and then moves conditionally: H acts on `coin`, then the `coin == 0` branch adds 1 to `position` and the `coin == 1` branch adds $2^n - 1$ (i.e. $-1 \bmod 2^n$); `position` operates under unsigned interpretation, and out-of-range values wrap around naturally. The step loop is generated with the IR's {obj}`Repeat <oracq.infrastructure.ir.Repeat>` structure and is not expanded at generation time. Applicability boundary: the period is fixed to a power of two, the coin is fixed to Hadamard, and the movement is fixed to $\pm 1$; general graphs or replaceable oracles should use the `graph_walks.py` framework instead.

## Validation approach

Category C1 (exact discrete semantics; acceptance criteria in `../../development/validation-plan.md` §2): the acting unitary must equal the walk definition pointwise. Evidence:

- Structure: construction and attribute assertions in `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`.
- Numerical: the walk part of `AlgorithmExpansionTests.test_ansatz_zero_angles_and_walk_one_step` — after one step of {obj}`cycle_walk(2) <oracq.algorithms.common.walks.cycle_walk>` from the all-zero input, the state is $(|1,0\rangle + |3,1\rangle)/\sqrt{2}$ (basis states written as `(position, coin)`), with both basis-state amplitudes equal to $1/\sqrt{2}$.
- Binding: this algorithm has no independent binding witness (no open-declaration entry point).

## Known gaps and planned stages

No known gaps (the validation matrix lists "—" in the gap column), stage V1.

## Related links

- Walks on graphs: [Adjacency Oracle](adjacency-oracle.md), [Szegedy Walk](szegedy-walk.md), [MNRS Quantum Walk Search](mnrs-search.md)
- Source: `src/oracq/algorithms/common/walks.py`
- API reference: [Quantum walks](../../api/algorithms/common/walks.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Experiment design: a Hadamard coined walk on rings of period $N=4$ (width=2) and $N=8$ (width=3), starting from $|0\rangle_{\rm position}|0\rangle_{\rm coin}$ and scanning the step count $s=0\ldots8$; the four backend paths (reference, rir-pysparq, adapter-pysparq, originir-ext) are compared amplitude by amplitude against an independent numpy reference (a step-by-step simulation applying the coin H followed by a conditional $\pm1$ shift on 0/1), with the TVD of the position marginal distribution computed separately. Note: at $s=0$ the program contains no gates, and UniQC has no qubit mapping for a zero-gate circuit, so the zero-step case excludes the originir-ext path.

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| cycle-walk-w2 | N=4, s=0..8 | 4 paths | max amplitude error / max position-distribution TVD | 5.6e-17 / 1.1e-16 |
| cycle-walk-w3 | N=8, s=0..8 | 4 paths | max amplitude error / max position-distribution TVD | 5.6e-17 / 1.1e-16 |

Informative metric (mean ring distance at $s=8$, quantum ballistic transport vs classical diffusion): at $N=8$, quantum 3.0 vs 1.875 for the classical symmetric random walk; at $N=4$, quantum 0.0 vs classical 1.0 (on the small ring, 8 steps return exactly to the origin — a revival). The qualitative difference between the two matches the known behavior of discrete-time quantum walks.

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

Artifacts: `out/verification/search_walks.json` (cases `cycle-walk-w2`, `cycle-walk-w3`).
