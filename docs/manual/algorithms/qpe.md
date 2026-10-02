# Quantum Phase Estimation

**English** · <a href="../../../zh/manual/algorithms/qpe.html">简体中文</a>

> Category C3 · Module [`oracq.algorithms.common.estimation`](../../api/algorithms/common/estimation.rst) · Stage V2

## Overview

Given an eigenstate $|\psi\rangle$ of a unitary operator $U$, $U|\psi\rangle = e^{2\pi i\varphi}|\psi\rangle$, quantum phase estimation (QPE) estimates $\varphi\in[0,1)$ with a $p$-bit phase register: the readout value approximates $2^p\varphi$, i.e. the phase is encoded on an integer grid of $2^p$ equal subdivisions. The circuit is the standard textbook construction — a Hadamard layer, controlled powers $U^{2^k}$, and the inverse QFT.

QPE is the readout core of the estimation module: the number-theory module's {obj}`order_finding <oracq.algorithms.basics.number_theory.order_finding>` assembles order finding on top of it, and this module's {obj}`amplitude_estimation <oracq.algorithms.common.estimation.amplitude_estimation>` performs amplitude estimation by running QPE on the Grover iterate.

## Interface and input model

```python
phase_estimation(operation, *, precision=2)
```

API entry point: {obj}`phase_estimation <oracq.algorithms.common.estimation.phase_estimation>`

- `operation`: a complete {obj}`Operation <oracq.infrastructure.builder.Operation>` supporting controlled invocation (input model UO, a unitary-operator oracle); the input state is prepared by the caller, and QPE itself performs no preparation.
- `precision`: the number of bits in the phase register, range 1..63.

Returns an `Operation`: all public registers of the input are preserved and `phase` is appended. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"qpe"` |

Cases that raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time: precision out of range, the callee interface occupying the `phase` parameter name, and inputs that do not support the required controlled invocation.

Readout example (the {obj}`phase_estimation <oracq.algorithms.common.estimation.phase_estimation>` entry in `applications/gallery.py`): for the phase gate $e^{i\pi/2}$ on $|1\rangle$ ($\varphi=1/4$) with `precision=3`, the `phase` readout is 2.

## Implementation notes

Register layout: the names, bit widths, and views of the input's public registers are preserved unchanged up to backend lowering, with `phase` ({obj}`Bits(precision) <oracq.infrastructure.ir.Bits>`) appended. Each power is stored as `repeat(1 << bit)` — {obj}`Repeat <oracq.infrastructure.ir.Repeat>` nodes count by power, and the gate sequence is not expanded during generation or text serialization. The inverse transform is implemented as an adjoint invocation of {obj}`qft(precision) <oracq.algorithms.common.fourier.qft>` (positive-sign Fourier convention, see `fourier.py`).

Readout and statistics are done on the host side (per the module docstring): when $\varphi$ lands exactly on a grid point, the phase distribution has a unique peak; otherwise the mass is spread over the neighborhood of the nearest grid point, and a finite-precision readout may correspond to several approximate values.

## Validation approach

Category C3 (acceptance criterion in `docs/development/validation-plan.md` §2: the output distribution equals the closed-form expectation). The three layers of evidence match the `estimation.py` row of `docs/development/validation-coverage.md`:

- Structure: construction and attribute assertions in `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`.
- Numerical: QPE has no independent numerical witness entry; it is covered indirectly through two downstream consumers — `AlgorithmExpansionTests.test_amplitude_estimation_half_probability` (QAE internally invokes `phase_estimation` on the Grover iterate, phase distribution cross-checked with places=10) and `AlgorithmExpansionTests.test_modular_multiplication_total_permutation_and_order` from the `number_theory` row (`order_finding` assembles QPE, phase distribution $p[0]=p[2]=0.5$).
- Binding: the module row is recorded as "—" (the input already requires a concrete `Operation`; there is no open slot).

The witness technique is exact state-vector simulation, taking the phase marginal distribution and cross-checking it against the closed form; no sampling assertions are made.

## Known gaps and planned stages

This algorithm has no gap of its own. The only gap registered for `estimation.py` in validation-coverage belongs to the QAE confidence-interval claim (see [amplitude estimation](qae.md)); the module as a whole is therefore at stage V2.

## Related links

- Source: `src/oracq/algorithms/common/estimation.py`
- API reference: [Phase, amplitude, and overlap estimation](../../api/algorithms/common/estimation.rst)
- Same-group pages: [amplitude estimation](qae.md), [quantum counting](quantum-counting.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

`tests/verification/verify_estimation.py` performs distribution-level numerical validation of this interface on real backends (13 cases in the QPE part, all passing). Experimental design:

- Known-phase unitaries: a single-bit diagonal phase gate $U|1\rangle=e^{2\pi i\varphi}|1\rangle$, with on-grid phases ($\varphi=0.5$ and $0.625$, $p=2..6$) and an off-grid phase ($\varphi=0.3$, $p=3..6$); the phase-register histogram is compared against the Dirichlet-kernel closed form $D(\delta)=\sin^2(\pi 2^p\delta)/(4^p\sin^2(\pi\delta))$.
- Non-diagonal unitary: the 2-bit `add_const(1)` takes its Fourier eigenstate $|\tilde\varphi_j\rangle$ ($j=1,3$) as input; the expected eigenphase $(-j/4)\bmod 1$ is obtained by independently constructing the permutation matrix and Fourier state with numpy and computing the eigenvalue, without reusing the in-library implementation.
- Unitary level: the $p=3$ circuit is exported via OriginIR-ext → UniQC `Circuit.to_matrix` to obtain the full unitary, which is compared element by element with an independent numpy assembly (H layer, controlled $U^z$, inverse DFT).
- Backend paths: reference (the built-in reference executor), rir-pysparq (PySparQ's native RIR interpreter), adapter-pysparq (the PySparQ adapter), and originir-ext (UniQC full-amplitude state vector).

| Case | Scale | Path | Metric value |
|---|---|---|---|
| On-grid phase readout | $p=2..6$, $\varphi=0.5/0.625$ | all four | peak $=2^p\varphi$, peak probability $=1$, max TVD $\le 4.4\times10^{-16}$ |
| Off-grid phase readout | $p=3..6$, $\varphi=0.3$ | all four | peak $=\mathrm{round}(2^p\varphi)$, peak probability $0.573/0.876$ ($\ge 4/\pi^2$), max TVD $\le 7.3\times10^{-16}$ |
| Fourier eigenstates | $j=1,3$, $p=3,4$ | all four | peak $=2^p\varphi$ ($\varphi=0.75/0.25$), peak probability $=1$, max TVD $\le 7.8\times10^{-16}$ |
| Full-unitary comparison | $p=3$, $16\times16$ | originir-ext+to_matrix | max_error $=1.0\times10^{-15}$ |

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_estimation.py
```

Artifact path: `out/verification/estimation.json` (case names prefixed `qpe-`).
