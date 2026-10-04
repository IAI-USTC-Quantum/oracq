# Hadamard Test

**English** · <a href="../../zh/manual/algorithms/hadamard-test.html">简体中文</a>

> Category C3 · Module [`oracq.algorithms.common.estimation`](../../api/algorithms/common/estimation.rst) · Stage V2

## Overview

Measures the real or imaginary part of the complex expectation value $\langle\psi|U|\psi\rangle$: after a single probe bit interferes through the controlled $U$, the Z expectation of `probe` equals the selected component. This operation performs no measurement; sampling and the classical computation of probability differences are done by the host.

## Interface and input model

```python
hadamard_test(unitary, preparation=None, *, component="real")
```

API entry point: {obj}`hadamard_test <oracq.algorithms.common.estimation.hadamard_test>`

- `unitary`: a complete unitary {obj}`Operation <oracq.infrastructure.builder.Operation>`, or a block encoding with no signal bit and `alpha=1` (UO, adapted via {obj}`as_block_encoding <oracq.algorithms.input_model.interfaces.as_block_encoding>`); inputs with a signal bit or $\alpha\ne 1$ raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time.
- `preparation`: the initial state preparation (SP, zero input + clean workspace); when omitted, the zero basis state of the target space is used ({obj}`basis_state(width) <oracq.algorithms.input_model.oracles.basis_state>`).
- `component`: `"real"` or `"imag"`.

Returns an `Operation` with registers `target`, `work`, and `probe`. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"hadamard_test"` |
| `readout_register` | `"probe"` |
| `component` | `"real"` / `"imag"` |

The preparation width must equal the unitary width; an illegal `component` value raises an error at generation time.

## Implementation notes

The circuit is preparation $\to$ H(probe) $\to$ controlled $U$ $\to$ (for the imaginary component, a $-\pi/2$ phase applied to probe) $\to$ H(probe). The imaginary readout places an $S^\dagger$-type phase after the controlled invocation, rotating the orthogonal component into the Z basis. The controlled invocation adapts to the BE interface with a zero-width signal view (`b["work"][:0]`), so BE-form inputs introduce no signal bit.

## Validation approach

Category C3 (acceptance criterion in `docs/development/validation-plan.md` §2: the output distribution equals the closed-form expectation). The three layers of evidence match the `estimation.py` row of `docs/development/validation-coverage.md`:

- Structure: construction and attribute assertions in `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`.
- Numerical: `AlgorithmExpansionTests.test_hadamard_test_real_and_imaginary` — a single-bit phase gate (angle 0.6, eigenvalue $e^{0.6i}$) acting on $|1\rangle$; the `real`/`imag` components are cross-checked against $p(0)-p(1)=\cos 0.6$ and $\sin 0.6$ respectively, with places=10.
- Binding: the module row is recorded as "—".

The witness technique is exact state-vector simulation, taking the probe marginal distribution and cross-checking it against the closed form; no sampling assertions are made.

## Known gaps and planned stages

This algorithm has no gap of its own; the only gap registered for `estimation.py` belongs to the QAE confidence-interval claim (see [amplitude estimation](qae.md)), and the module as a whole is at stage V2.

## Related links

- Source: `src/oracq/algorithms/common/estimation.py`
- API reference: [Phase, amplitude, and overlap estimation](../../api/algorithms/common/estimation.rst)
- Same-group pages: [SWAP test](swap-test.md) (the corresponding readout for the overlap of two states)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

`tests/verification/verify_estimation.py` performs expectation-value-level numerical validation of this interface on real backends (6 cases in total, all passing). Experimental design:

- A single-bit phase gate $U=\mathrm{diag}(1,e^{i\theta})$ acting on $|1\rangle$ ($\theta=0.6$ and $-1.1$): the closed-form expectation is $e^{i\theta}$.
- A two-bit diagonal unitary $U|x\rangle=e^{i\theta_x}|x\rangle$ (angle table $(0.35,-0.9,1.7,0.55)$, implemented as a controlled global phase) paired with the {obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>` complex-amplitude preparation $\psi=(0.5,\,0.5i,\,0.5,\,-0.5)$: the expectation value $\sum_x|\psi_x|^2 e^{i\theta_x}$ is computed independently with the math library.
- Backend paths: reference, rir-pysparq, adapter-pysparq, and originir-ext; the Z expectation of probe, $p(0)-p(1)$, is compared against the real/imaginary components of the expectation value.

| Case | Scale | Path | Metric value |
|---|---|---|---|
| Phase gate $\theta=0.6$ | 1 target bit | all four | real $=0.8253356$ ($\cos\theta$), imag $=0.5646425$ ($\sin\theta$), error $\le 5.6\times10^{-16}$ |
| Phase gate $\theta=-1.1$ | 1 target bit | all four | real $=0.4535961$, imag $=-0.8912074$, error $\le 4.4\times10^{-16}$ |
| Diagonal unitary + complex preparation | 2 target bits | all four | real $=0.5711657$, imag $=0.2684807$, error $\le 1.1\times10^{-16}$ |

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_estimation.py
```

Artifact path: `out/verification/estimation.json` (case names prefixed `hadamard-`).
