# Amplitude Amplification

**English** · <a href="../../../zh/manual/algorithms/amplitude-amplification.html">简体中文</a>

> Category C3 · Module [`oracq.algorithms.common.search`](../../api/algorithms/common/search.rst) · Stage V1

## Overview

Coherently amplifies the success subspace of an arbitrary state oracle: given a preparation $A$ of a {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>` that produces the success state with probability $a$ (success condition signal==0), the iteration operator

$$
Q = A\,(2|0\rangle\langle 0|-I)\,A^\dagger\, S_{\mathrm{good}}
$$

rotates the success amplitude by about $2\theta$ each time ($\sin^2\theta=a$). Implementation basis: the amplitude amplification framework of Brassard et al. ([arXiv:quant-ph/0005055](https://arxiv.org/abs/quant-ph/0005055); see the Implementation basis section of the algorithm catalog). It generalizes Grover search: the initial state need not be uniform, and the marking information is carried by the signal register instead of a phase oracle.

## Interface and input model

```python
amplify_success(state, *, iterations=1)
```

API entry point: {obj}`amplify_success <oracq.algorithms.common.search.amplify_success>`

- `state`: a `StateOracle` supporting adjoint calls (a state oracle with the `target`/`signal` interface; the signal indicates the success subspace, an FO amplitude-type access). A non-`StateOracle` input raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time.
- `iterations`: a non-negative amplification count (default 1). The count must be chosen with the input success probability in mind; too many iterations may lower the success probability (per the docstring), and this entry point does not choose it automatically.

Returns a `StateOracle`: the public interface is the same as the input, and the success condition remains signal==0. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"amplitude_amplification"` |
| `success_condition` | `"signal == 0"` |

## Implementation notes

$A$ is first invoked once forward to prepare the initial state; each subsequent iteration is {obj}`reflect_zero(signal) <oracq.algorithms.input_model.block_encoding.reflect_zero>` (a phase flip on the good subspace with signal==0, acting as $S_{\mathrm{good}}$) → $A^\dagger$ → positive reflection about the all-zero state of target+signal → $A$. The operator shape matches {obj}`grover_iterate <oracq.algorithms.common.search.grover_iterate>`; the only difference is the source of the marks: here $S_{\mathrm{good}}$ is carried by the signal register itself, so it applies to success indications from any source (such as the decision bit of a QRAM lookup result, or a check bit in an algorithm's workspace). The iterations are stored as `repeat(iterations)` and are not expanded during generation or text serialization.

## Validation approach

Category C3 (acceptance criterion in `docs/development/validation-plan.md` §2: the output distribution equals the closed-form expectation). The three layers of evidence match the `search.py` row of `docs/development/validation-coverage.md`:

- Structure: construction and attribute assertions in `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` ({obj}`require_instance <oracq.algorithms.input_model.contracts.require_instance>` validates the StateOracle and the non-negative iteration count).
- Numerical: `AlgorithmExpansionTests.test_amplification_of_known_quarter_probability` — builds a state oracle with `ry(2π/3)` on the signal, success probability $P(\text{signal}=0)=\cos^2(\pi/3)=1/4$; after the default single iteration $P(\text{signal}=0)=1$ ($\theta=\pi/6$, one iteration rotates to $\sin^2(3\theta)=1$), places=10.
- Binding: the module row registers "—".

The witness technique is exact state-vector simulation, taking the signal marginal distribution and cross-checking it against the closed form, with no sampling assertions.

## Known gaps and planned stages

No known gaps (the `search.py` row lists "—" in the gap column), stage V1.

## Related links

- Source: `src/oracq/algorithms/common/search.py`
- API reference: [Search and amplitude amplification](../../api/algorithms/common/search.rst)
- Same-group pages: [Grover search](grover.md) (the special case of a uniform initial state plus a phase oracle), [amplitude estimation](qae.md) (the dual readout of the same iteration operator)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Experiment design: build a state oracle whose success subspace is `signal==0` (`target` is a two-bit uniform superposition, the signal bit is $R_y(3\pi/4)$, initial success probability $a=\cos^2(3\pi/8)=\sin^2(\pi/8)\approx0.1464$, hence $\theta_a=\pi/8$); measure the zero-signal probability of {obj}`amplify_success <oracq.algorithms.common.search.amplify_success>` for $k=0\ldots3$, compared against $\sin^2((2k+1)\theta_a)$, on the four backend paths (reference, rir-pysparq, adapter-pysparq, originir-ext). The sequence covers amplification ($k=1,2$ reaching $0.8536$) and overshoot decay ($k=3$ returning to $0.1464$).

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| amplify-success-curve | target 2 + signal 1, k=0..3 | 4 paths | max success-probability error | 1.9e-15 |

Theoretical values over the iteration counts: 0.146447, 0.853553, 0.853553, 0.146447, reproduced pointwise.

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

Artifacts: `out/verification/search_walks.json` (case `amplify-success-curve`).
