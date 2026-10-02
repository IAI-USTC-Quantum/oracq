# Qubitization Walk

**English** · <a href="../../../zh/manual/algorithms/qubitization-walk.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.common.transforms`](../../api/algorithms/common/transforms.rst) · Stage V1

## Overview

Given a block encoding $U$ of a Hermitian matrix $A$, construct the qubitization walk operator

$$
W = (2\Pi - I)\,U, \qquad \Pi = |0\rangle\langle 0|_{\text{signal}} \otimes I_{\text{target}},
$$

that is, $U$ acts first, followed by a positive reflection about the $|0\rangle$ state of the signal register. $W$ acts as a rotation on the two-dimensional invariant subspace associated with each spectral value $x$ of $A$ and is the elementary iteration unit that turns a block encoding into spectral-function transforms (the qubitization construction of Low–Chuang 2019, see `../../development/algorithm-coverage.md` A2; for the QSVT framework see Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)). This page covers only the assembly of the single-step walk operator; the alternating invocations of the phase sequence are covered in [QSVT phase sequence](qsvt-sequence.md).

## Interface and input model

```python
qubitization_walk(a)
```

API entry: {obj}`qubitization_walk <oracq.algorithms.common.transforms.qubitization_walk>`

- `a`: a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`, the block encoding of the walk matrix (input model BE).

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>`, with registers `target` (width `a.width`) and `signal` (width `a.signal_qubits`). Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"qubitization_walk"` |

The normalization `alpha` of the input BE does not propagate through attributes: the return value is a bare `Operation`, not a `BlockEncoding`; the spectral variable $x$ is defined relative to the input normalization $\alpha$, i.e. $x = \lambda/\alpha$.

## Implementation notes

The generation strategy is two structural calls: {obj}`invoke <oracq.algorithms.input_model.oracles.invoke>` hangs the input BE on `target | signal`, then {obj}`reflect_zero(b, signal, positive=True) <oracq.algorithms.input_model.block_encoding.reflect_zero>` realizes $2\Pi - I$ (`positive=True` makes the $|0\rangle_{\text{signal}}$ branch take $+1$ and the other branches $-1$). Repeated invocations are orchestrated by the caller (such as {obj}`qsvt_sequence <oracq.algorithms.common.transforms.qsvt_sequence>` or the low-rank decomposition pipeline); this function builds in no {obj}`Repeat <oracq.infrastructure.ir.Repeat>` structure.

Applicability boundary: the input must be a block encoding; sparse oracles or low-rank decomposition outputs must first be adapted into a BE via `sparse.py` / `lowrank.py` and the like. The walk operator itself applies no polynomial transform; its spectral properties show up only together with a phase sequence or a projective measurement.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../../development/validation-plan.md` §2). Current evidence (consistent with the `transforms.py` row of the validation coverage matrix):

- Structure: `tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` and `test_trotter_only_input_does_not_need_block_encoding` cover the protocol and register contracts of the same module's assembly chain.
- Numerical: `tests/core/test_qsvt.py:PhaseSynthesisTests.test_convention_matches_qsvt_sequence` pins, under random phases, the pointwise agreement between the circuit's zero-signal block under the same signal/reflection convention and {obj}`qsp_response <oracq.algorithms.common.qsvt.qsp_response>` (delta = 1e-10); the sequence skeleton that the walk step embeds into is thereby witnessed indirectly.
- Binding: no independent binding witness (the input already requires a concrete BE).

## Known gaps and planned stages

No known gaps; the stage V1 witnesses are complete (the sequence convention agrees pointwise to 1e-10). The standalone rotation property of the single-step walk operator (the rotation angle on the two-dimensional invariant subspace) has no direct witness yet and is covered by the sequence-level convention test.

## Related links

- Source: `src/oracq/algorithms/common/transforms.py`
- API reference: [Matrix transform sequences](../../api/algorithms/common/transforms.rst)
- Same-family pages: [QSVT phase sequence](qsvt-sequence.md), [Oblivious amplitude amplification](oblivious-amplification.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_fourier.py` (real-backend execution, no mock substitutes). The input BE is a small gate-level bound instance: {obj}`matrix_pauli_encoding <oracq.algorithms.input_model.block_encoding.matrix_pauli_encoding>` encodes the 1-bit Hermitian matrix $A = \begin{pmatrix} 0.5 & 0.3 \\ 0.3 & -0.1 \end{pmatrix}$ (an independent explicit Pauli expansion gives $\alpha = 0.8$, spectral values $\lambda/\alpha \approx 0.780, -0.280$). The full unitary of the walk operator is extracted via OriginIR-ext + UniQC `to_matrix` and compared element by element with $(2\Pi - I)U$, where $U$ is the full unitary of the same BE program (also cross-checked against the unitary assembled column by column over basis states on the reference path); the spectral properties are checked independently: the eigenangles of the walk unitary must fall in $\{\pm\arccos(\lambda/\alpha)\} \cup \{0, \pi\}$, and the zero-signal block must equal $A/\alpha$.

| Case | Scale | Backend paths | Metric | Value |
|---|---|---|---|---|
| `qubitization-walk-unitary-spectrum` | target 1 + signal 2 (8-dimensional) | originir-ext + UniQC to_matrix, reference cross-check | unitary_max_error | 1.7e-16 |
| same as above | — | — | be_path_crosscheck (BE unitary deviation between the two paths) | 1.1e-16 |
| same as above | — | — | zero_block_max_error ($\Pi W\Pi$ vs $A/\alpha$) | 1.8e-16 |
| same as above | — | — | spectrum_max_deviation (eigenangles vs $\pm\arccos(\lambda/\alpha)$) | 4.4e-16 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_fourier.py
```

Artifacts: `out/verification/fourier.json`.

## Numerical validation

Experiment design: construct a Hermitian (Householder-type) block encoding $U=R_y(2\theta)Z\otimes I_{\rm target}$ whose zero-signal block is $x=\cos(\pi/6)$ ($\theta=\pi/6$, target/signal one bit each). Three layers of comparison:

1. Rotation property (closing the gap of "no direct witness for the rotation angle on the two-dimensional invariant subspace"): the walk is repeated $k=1\ldots4$ times, and the zero-signal probabilities are compared against the Chebyshev closed form $T_k(x)^2=\cos^2(k\arccos x)$, i.e. the sequence $0.75,\ 0.25,\ 0,\ 0.25$ ($k=3$ is an exact zero), on the four backend paths (reference, rir-pysparq, adapter-pysparq, originir-ext);
2. Unitary decomposition: UniQC `Circuit.to_matrix` verifies $W=(2\Pi-I)U$;
3. Spectrum: the eigenvalues of $W$ should be $\mathrm e^{\pm i\pi/6}$ (each of the two target sectors contributes a pair, doubly degenerate).

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| qubitization-walk-chebyshev | 2 qubits, k=1..4 | 4 paths + to_matrix | max error of zero-signal probability | 5.6e-16 |
| ditto | ditto | ditto | $\|W-(2\Pi-I)U\|_{\max}$ | 2.1e-16 |
| ditto | ditto | ditto | eigenphase error (vs $\pm\pi/6$) | 3.3e-16 |

Premise: the $T_k$ rotation relation requires the BE unitary to be Hermitian; for non-Hermitian $U$, the correct qubitized iteration must alternate $U$ and $U^\dagger$ (the QSVT sequence does exactly that), and plain powers $W^k$ do not obey Chebyshev — this validation therefore uses a Hermitian construction.

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

Artifacts: `out/verification/search_walks.json` (the case `qubitization-walk-chebyshev`).
