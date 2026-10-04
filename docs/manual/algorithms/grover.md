# Grover Search

**English** · <a href="../../zh/manual/algorithms/grover.html">简体中文</a>

> Category C3 · Module [`oracq.algorithms.common.search`](../../api/algorithms/common/search.rst) · Stage V1

## Overview

Searches a space of $2^{\text{width}}$ basis states for targets marked by a phase oracle: each Grover iteration rotates the state toward the marked subspace, and after a number of iterations the hit probability peaks. Implementation basis: the amplitude amplification/estimation framework of Brassard et al. ([arXiv:quant-ph/0005055](https://arxiv.org/abs/quant-ph/0005055); see the Implementation basis section of the algorithm catalog); the iteration count is given by the caller — this entry point does not choose it automatically.

## Interface and input model

```python
grover(phase_oracle, width, *, iterations=1, preparation=None)
```

API entry point: {obj}`grover <oracq.algorithms.common.search.grover>`

- `phase_oracle`: a phase {obj}`Operation <oracq.infrastructure.builder.Operation>` marking the target basis states (FO, a phase-function oracle), with interface `target` and optionally `work`.
- `width`: the target bit width of the search space.
- `iterations`: a non-negative iteration count (default 1).
- `preparation`: the initial state preparation (SP, must support adjoint calls); when omitted, {obj}`uniform_state(width) <oracq.algorithms.input_model.oracles.uniform_state>` is used. A width inconsistent with `width` raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`.

Returns a {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>`: `target` holds the search result, and `signal` keeps the workspace of the phase query and the preparation. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"grover"` |
| `iterations` | echo of the iteration count |
| `validation_stage` | `"paradigm"` |

Related entry points in the same module:

- {obj}`phase_from_database(database) <oracq.algorithms.common.search.phase_from_database>`: adapts a single-output-bit {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` (QRAM lookup) into a phase oracle — load, Z on `work`, unload; the product is annotated `phase_oracle`.
- {obj}`grover_iterate(preparation, marked) <oracq.algorithms.common.search.grover_iterate>`: returns the iteration operator $Q=A(2|0\rangle\langle 0|-I)A^\dagger S_{\mathrm{good}}$ (marks come from {obj}`phase_marks <oracq.algorithms.input_model.oracles.phase_marks>`), used by [amplitude estimation](qae.md) to assemble QPE.

## Implementation notes

A single iteration is: phase oracle → $A^\dagger$ → positive reflection about the all-zero state ({obj}`reflect_zero <oracq.algorithms.input_model.block_encoding.reflect_zero>`, `positive=True`) → $A$; the reflection about the initial state covers both `target` and the preparation workspace (per the docstring). Register layout: the leading segment of `signal` is the phase oracle's work (`phase_work` bits) and the trailing segment is the preparation's work; when the phase oracle has no work, the leading segment has width 0. The iteration count is stored as `repeat(iterations)` and is not expanded during generation or text serialization.

The product is wrapped as a `StateOracle` and can be composed further as a state oracle (for example as the input of {obj}`amplify_success <oracq.algorithms.common.search.amplify_success>`). Applicability boundary: the phase oracle's work must restore itself after its call (the reflections cover only target and the preparation workspace); too many iterations overshoot the peak and lower the hit probability, the same iteration-count constraint as in [amplitude amplification](amplitude-amplification.md).

## Validation approach

Category C3 (acceptance criterion in `docs/development/validation-plan.md` §2: the output distribution equals the closed-form expectation). The three layers of evidence match the `search.py` row of `docs/development/validation-coverage.md`:

- Structure: construction and attribute assertions in `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` (non-negative iteration count, consistent preparation width, and similar checks are validated at generation time).
- Numerical: the module row registers `AlgorithmExpansionTests.test_amplification_of_known_quarter_probability` (a direct witness of the same iteration principle, see [amplitude amplification](amplitude-amplification.md)); the full circuit of {obj}`grover() <oracq.algorithms.common.search.grover>` has no standalone numerical entry — its iteration operator `grover_iterate` is witnessed indirectly via QAE's `test_amplitude_estimation_half_probability` (the `estimation.py` row), and the uniform-initial-state special case is demonstrated by the `grover` entry of the gallery (`grover(phase_marks(2, [3]), 2)`, reading out target=3, signal=0).
- Binding: the module row registers "—".

## Known gaps and planned stages

No known gaps (the `search.py` row lists "—" in the gap column), stage V1.

## Related links

- Source: `src/oracq/algorithms/common/search.py`
- Tutorial: [Search for an element and estimate the success probability](../../tutorials/search-and-estimation.md)
- API reference: [Search and amplitude amplification](../../api/algorithms/common/search.rst)
- Same-group pages: [amplitude amplification](amplitude-amplification.md), [amplitude estimation](qae.md), [quantum counting](quantum-counting.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Experiment design: two instances with $n=3$ ($N=8$) sweep the iteration count $k$ on the four backend paths (reference, rir-pysparq, adapter-pysparq, originir-ext); the marked probability and the per-basis-state distribution are compared against the closed-form two-level expressions $\sin^2((2k+1)\theta)/t$ and $\cos^2((2k+1)\theta)/(N-t)$. Instance 1: `phase_marks(3, (5,))`, $t=1$, $\theta=\arcsin(1/\sqrt 8)$, $k=0\ldots4$ (including the oscillating decay after passing the optimal iteration count); instance 2: `phase_from_database(gate_database)` marking $\{5,6\}$, $t=2$, $\theta=\pi/6$, amplified exactly to 1 at $k=1$.

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| grover-phase-marks-n3-t1 | N=8, t=1, k=0..4 | 4 paths | max success-rate error / max distribution error | 2.3e-15 / 2.3e-15 |
| grover-xor-database-n3-t2 | N=8, t=2, k=0..3 | 4 paths | max success-rate error / max distribution error | 1.4e-15 / 7.2e-16 |

The theoretical success rates of instance 1 over the iteration counts are 0.125, 0.78125, 0.9453125, 0.330078125, and 0.012207031 ($k=2$ optimal, $k\ge3$ overshooting and falling back), reproduced pointwise on all four paths.

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

Artifacts: `out/verification/search_walks.json` (cases `grover-phase-marks-n3-t1` and `grover-xor-database-n3-t2`, with per-iteration values).
