# Quantum Counting

**English** · <a href="../../../zh/manual/algorithms/quantum-counting.html">简体中文</a>

> Category C3 · Module [`oracq.algorithms.common.estimation`](../../api/algorithms/common/estimation.rst) · Stage V2

## Overview

Estimates the number $t$ of marked elements in a search space of size $N=2^n$. Quantum counting is not a separate circuit but a direct application of amplitude estimation (same Brassard framework, [arXiv:quant-ph/0005055](https://arxiv.org/abs/quant-ph/0005055)): on a uniform superposition the good-state probability is exactly $a=t/N$; run standard amplitude estimation and multiply $N$ back:

$$
\hat t = N\cdot\hat a.
$$

This repository registers it as a gallery entry (`quantum_counting` in `applications/gallery.py`) that reuses the {obj}`amplitude_estimation <oracq.algorithms.common.estimation.amplitude_estimation>` entry point.

## Interface and input model

There is no separate entry function; the recipe is:

```python
operation = amplitude_estimation(uniform_state(n), marked, precision=p)
```

- {obj}`uniform_state(n) <oracq.algorithms.input_model.oracles.uniform_state>`: uniform state preparation (SP), making $a=t/N$.
- `marked`: the set of integer indices of the marked basis states, $t=|M|$.
- `precision`: the number of phase-register bits, range 1..63.

The returned object's registers and attributes are the same as in [amplitude estimation](qae.md) (`target`/`work`/`phase`, `decoder` is `"sin(pi*phase/2**precision)**2"`): after reading out `phase`, decode via {obj}`amplitude_from_phase <oracq.algorithms.common.estimation.amplitude_from_phase>` and multiply by $2^n$ to get $\hat t$. The gallery instance takes $n=2$, `marked=(3,)`, `precision=4` ($t=1$, $a=1/4$).

## Implementation notes

Circuit generation is identical to QAE; the difference is entirely in the classical-side decoding (multiplying by $N$). $\hat t$ need not be an integer, and rounding it into a legal range is the host's decision; the estimation precision is set by the QAE grid density (`precision`), and the point estimate is exact when $t/N$ lands exactly on the grid.

## Validation approach

Category C3 (acceptance criterion in `docs/development/validation-plan.md` §2). Quantum counting is registered in `docs/development/validation-coverage.md` under the `estimation.py` row (the algorithm column includes "quantum counting"), with no separate witness entry:

- Structure: `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`.
- Numerical: the circuit witness is QAE's `AlgorithmExpansionTests.test_amplitude_estimation_half_probability` (see [amplitude estimation](qae.md)); the gallery instance is demonstrated by the `quantum_counting` entry of `applications/gallery.py` ("multiply the amplitude estimate by 4 to estimate the number of marked items").
- Binding: the module row registers "—".

## Known gaps and planned stages

Same as QAE: the confidence-interval claim is unwitnessed (only the point estimate is verified), and the error propagation of multiplying the point estimate by $N$ is likewise not cross-checked; the `estimation.py` module is at stage V2.

## Related links

- Source: `src/oracq/algorithms/common/estimation.py` (circuit), `src/oracq/applications/gallery.py` (gallery entry)
- API reference: [Phase, amplitude, and overlap estimation](../../api/algorithms/common/estimation.rst)
- Same-group pages: [amplitude estimation](qae.md), [Grover search](grover.md), [quantum phase estimation](qpe.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Experiment design: $n=3$ ($N=8$), marked $=\{1,5,7\}$ ($t=3$, $a=3/8$), precision $=4$ (grid $M=16$). Following the recipe above, `amplitude_estimation(uniform_state(3), marked, precision=4)` is generated, and the full distribution of the phase register is read out on the four backend paths (reference, rir-pysparq, adapter-pysparq, originir-ext), compared against an independent Dirichlet-kernel reference (a QPE grid convolution over the eigenphases $\pm\theta/\pi$, $\theta=\arcsin\sqrt{3/8}$), and the count is decoded as $\hat t=N\sin^2(\pi y/M)$.

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| quantum-counting-n3-t3 | N=8, t=3, M=16 | 4 paths | phase-distribution TVD | 5.7e-15 |
| same | same | same | peak estimate $\hat t$ (error) | 2.4693 (0.531, ≤ 1 grid step) |
| same | same | same | central grid mass (QAE lower bound $8/\pi^2$) | 0.8529 (≥ 0.8106) |

The true value $\theta M/\pi\approx3.357$ lies between grid points 3 and 4; the bimodal distribution ($y=3,13$, each about 0.325) matches the theoretical positions; the point-estimate precision is set by the grid density, consistent with the stated convention that "$\hat t$ need not be an integer".

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

Artifacts: `out/verification/search_walks.json` (case `quantum-counting-n3-t3`).
