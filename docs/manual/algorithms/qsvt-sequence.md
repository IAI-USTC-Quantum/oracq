# QSVT Phase Sequence

**English** · <a href="../../zh/manual/algorithms/qsvt-sequence.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.common.transforms`](../../api/algorithms/common/transforms.rst) · Stage V1

## Overview

Given a block encoding $U$ of a matrix $A$ and a list of real phases $\Phi = (\varphi_0, \dots, \varphi_d)$, explicitly assemble the QSVT sequence: on the two-dimensional invariant subspace associated with each singular value $x$ of $A$, the zero-signal block realizes the polynomial

$$
p(x) = \bigl[S(\varphi_0)\, W(x)\, S(\varphi_1)\, W(x) \cdots W(x)\, S(\varphi_d)\bigr]_{00},
\qquad
W(x) = \begin{pmatrix} x & s \\ s & -x \end{pmatrix},\ s = \sqrt{1 - x^2},\
S(\varphi) = \operatorname{diag}(e^{i\varphi}, e^{-i\varphi}),
$$

The phases are arranged in time order ($\varphi_0$ acts first), for $d$ BE calls and $d + 1$ phases in total. The realizable $(P, Q)$ satisfy the necessary and sufficient conditions: $\deg P \le d$, $\deg Q \le d - 1$, the parity of $P$ is $d \bmod 2$, the parity of $Q$ is $(d-1) \bmod 2$, and the polynomial identity $P\bar{P} + (1-x^2)Q\bar{Q} \equiv 1$ holds (in particular $|P(\pm 1)| = 1$ is forced). The sequence convention is taken from the QSVT framework (Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838); see the `qsvt.py` module docstring); this function is the framework's circuit-assembly layer, and valid phases are computed and self-checked by [QSP phase synthesis](qsp-phase-synthesis.md).

## Interface and input model

```python
qsvt_sequence(a, phases)
```

API entry: {obj}`qsvt_sequence <oracq.algorithms.common.transforms.qsvt_sequence>`

- `a`: a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`, the block encoding of the matrix to transform (input model BE).
- `phases`: an iterable of real phases in time order ($\varphi_0$ acts first); converted internally to a tuple of floats, with no correctness check of its own.

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` (a bare operation, not a `BlockEncoding`), with registers `target` (width `a.width`) and `signal` (width `a.signal_qubits`). Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"qsvt_sequence"` |
| `validation_stage` | `"paradigm"` (a witness-level label) |

The BE normalization does not propagate through attributes; the spectral variable $x$ is defined relative to the input normalization $\alpha$ ($x = \lambda/\alpha$), and the amplitude semantics are supplied by the caller.

## Implementation notes

Each phase is realized as a pair of global phase gates: `global_phase(-φ)` plus `global_phase(2φ)` conditioned on the zero branch of `signal`; the net effect is that the zero branch of the signal takes $e^{+i\varphi}$ and the other branches take $e^{-i\varphi}$, i.e. $S(\varphi)$. Between phases, the input BE is invoked alternately: after phases with an even index it is invoked forward with {obj}`invoke <oracq.algorithms.input_model.oracles.invoke>`, and after phases with an odd index it is invoked adjoint (`adjoint`), $d$ calls in total. The degenerate case `a.signal_qubits == 0` has a dedicated branch (the phases degenerate into unconditional global phases).

Applicability boundary: the input must be a block encoding; sparse oracles or QRAM must first be turned into a BE by an adapter layer. This function accepts arbitrary real phases, and the numerical correctness of the sequence is the caller's responsibility (the computation and round-trip self-check of valid phases are covered in [QSP phase synthesis](qsp-phase-synthesis.md)). The sequence takes the complex polynomial $P$ directly as the zero-signal block and performs no real-part extraction; the real-part extraction of real targets is done by the `qsvt.py` standard-transform family through the LCU combination $(U_\Phi + U_{-\Phi})/2$.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../../development/validation-plan.md` §2). Current evidence (consistent with the `transforms.py` row of the validation coverage matrix):

- Structure: `tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` and `test_trotter_only_input_does_not_need_block_encoding` cover the protocol and register contracts of the same module's assembly chain.
- Numerical: `tests/core/test_qsvt.py:PhaseSynthesisTests.test_convention_matches_qsvt_sequence` — with fixed-seed random phases (6 of them), the circuit's zero-signal block of a 2×2 diagonal BE agrees pointwise with {obj}`qsp_response <oracq.algorithms.common.qsvt.qsp_response>` (delta = 1e-10), pinning the circuit-assembly layer and the mathematical convention to each other.
- Binding: no independent binding witness (the input already requires a concrete BE).

## Known gaps and planned stages

No known gaps (consistent with the `transforms.py` row of the validation coverage matrix); the stage V1 witnesses are complete (the sequence convention agrees pointwise to 1e-10).

## Related links

- Source: `src/oracq/algorithms/common/transforms.py`
- API reference: [Matrix transform sequences](../../api/algorithms/common/transforms.rst)
- Same-family pages: [qubitization walk](qubitization-walk.md), [Oblivious amplitude amplification](oblivious-amplification.md), [QSP phase synthesis](qsp-phase-synthesis.md), [fixed-point search](fixed-point-search.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_fourier.py` (real-backend execution, no mock substitutes). The input BE is a small gate-level bound instance ({obj}`matrix_pauli_encoding <oracq.algorithms.input_model.block_encoding.matrix_pauli_encoding>` encoding a 1-bit Hermitian matrix $A$; an independent Pauli expansion gives $\alpha = 0.8$). Three mutually independent oracles: first, with fixed-seed random phases (7 of them, $d = 6$), the full unitary of the sequence is compared element by element with the time-ordered alternating product $S(\varphi_d)\,U/U^\dagger \cdots S(\varphi_0)$ ($U$ is the BE's measured unitary, and $S(\varphi)$ is built per this page's convention from the signal's zero branch $e^{+i\varphi}$ and the others $e^{-i\varphi}$); second, the zero-signal block is compared with $p(A/\alpha)$, where $p(x)$ is given by an independent 2×2 matrix-product implementation of this page's formula (not through the in-library `qsp_response`); third, Chebyshev end to end — with the $T_4, T_5$ phases synthesized by {obj}`qsp_phases <oracq.algorithms.common.qsvt.qsp_phases>` as input, the zero-signal block is compared with the independently recursed $T_d(A/\alpha)$. The `signal_qubits == 0` degenerate branch (a Pauli-word BE) is also covered: the phases degenerate into unconditional global phases, and the sequence operator is $e^{i\sum\varphi} X^d$.

| Case | Scale | Backend paths | Metric | Value |
|---|---|---|---|---|
| `qsvt-sequence-random-phases` | $d = 6$, 8-dimensional | originir-ext + UniQC to_matrix | unitary_max_error | 5.6e-16 |
| same as above | — | — | zero_block_max_error (vs independent 2×2 product) | 6.5e-16 |
| `qsvt-sequence-chebyshev-t4` | $d = 4$ | originir-ext + UniQC to_matrix | zero_block_max_error | 5.6e-16 |
| `qsvt-sequence-chebyshev-t5` | $d = 5$ | originir-ext + UniQC to_matrix | zero_block_max_error | 2.9e-16 |
| `qsvt-sequence-degenerate-signal0` | $d = 2, 3$ (signal width 0) | originir-ext + UniQC to_matrix | max_error | 1.6e-16 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_fourier.py
```

Artifacts: `out/verification/fourier.json`.

### Supplementary validation: multiple backend paths and the consumer-side matrix block (verify_hamiltonian.py)

`tests/verification/verify_hamiltonian.py` (all 27 cases PASS) adds two layers of evidence from the QSVT consumer's perspective: first, with fixed-seed random phases (6 of them, $d = 5$), the zero-signal amplitudes of a 1-qubit diagonal BE (spectral variables $x \in \{1.0, -0.647\}$) are read out on the **four paths** reference / rir-pysparq / adapter-pysparq / originir-ext and cross-checked pointwise against a 2×2 matrix-recursion response implemented independently in numpy per this page's convention; second, for a 2-qubit non-diagonal BE (four-term Pauli expansion, $\alpha = 0.95$), the full 4×4 zero-signal block is read out with the $T_4$ phases synthesized by `qsp_phases` and cross-checked against an independent Horner evaluation of the matrix polynomial $T_4(A/\alpha)$. It also contains a qubitization recursion witness: the zero-signal blocks of $W^n$ ($n = 1, 2, 3, 5$) of {obj}`qubitization_walk <oracq.algorithms.common.transforms.qubitization_walk>` agree with $T_n(x)$ (max_error 1.2e-15).

| Case | Scale | Backend paths | Metric | Value |
|---|---|---|---|---|
| `qsvt-sequence-convention-1q` | $d = 5$, two eigenstates | four paths | max_deviation | 1.26e-15 |
| `qsvt-sequence-matrix-block-2q` | $T_4$, 4×4 block | reference, rir-pysparq, originir-ext | max_error | 3.9e-16 |
| `qubitization-chebyshev-recurrence` | $W^n$, $n$ = 1,2,3,5 | reference, rir-pysparq, originir-ext | max_error | 1.2e-15 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_hamiltonian.py
```

Artifacts: `out/verification/hamiltonian.json`.
