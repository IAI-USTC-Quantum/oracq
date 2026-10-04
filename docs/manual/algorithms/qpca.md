# Quantum Principal Component Analysis

**English** · <a href="../../zh/manual/algorithms/qpca.html">简体中文</a>

> Category C3 · Module [`oracq.algorithms.qml.qpca`](../../api/algorithms/qml/qpca.rst) · Stage V1

## Overview

Reads out the spectrum of a density matrix $\rho$: phase estimation is applied to the unitary step $e^{-i\rho\Delta t}$, and the phase read out for the eigenstate $|\lambda\rangle$ is

$$
\varphi = -\frac{\lambda\,\Delta t}{2\pi} \pmod 1,
$$

which decodes to the eigenvalue $\lambda$. The unitary step is supplied by LMR density matrix exponentiation (see [Density Matrix Exponentiation](density-matrix-exponentiation.md)), and one call consumes $2^{\text{precision}} - 1$ copies of $\rho$ in total. The implementation basis is Lloyd–Mohseni–Rebentrost 2014 ("Quantum principal component analysis", Nature Physics 10, 631).

## Interface and input model

```python
qpca(preparation, *, precision, step_time, system=None, swap_width=None, name=None)
```

API entry point: {obj}`qpca <oracq.algorithms.qml.qpca.qpca>`

- `preparation`: the state preparation of the $\rho$ copies ({obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>`); the input model is SP + QRAM: pure states use {obj}`gate_state_prep(amplitudes) <oracq.algorithms.input_model.oracles.gate_state_prep>`, mixed states are adapted via `density.gate_purification(rho)`'s `PurificationAccess.as_state_preparation()`; the preparation must have zero work.
- `precision`: the number of phase-register bits, range 1..6; the copy count grows exponentially with the precision (already 63 copies at precision = 6, matching the copies cap of 63 in {obj}`density_matrix_exponentiation <oracq.algorithms.qml.qpca.density_matrix_exponentiation>`).
- `step_time`: the single-step evolution time $\Delta t$, which must be positive; $\lambda\Delta t \ll 2\pi$ must hold to avoid readout aliasing.
- `system`: an optional input-state preparation for the system (default $|0\rangle$) whose width must equal swap_width; feeding different eigenstates reads out the corresponding different eigenvalue peaks.
- `swap_width`: the prefix bit width participating in the swap, with the same semantics as in `density_matrix_exponentiation` (in purification scenarios take the system bit width of $\rho$).

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `system`, `copies`, and `phase`. After reading out `phase`, decode the eigenvalue with {obj}`eigenvalue_from_phase(value, precision, step_time) <oracq.algorithms.qml.qpca.eigenvalue_from_phase>`: the two's-complement phase is taken with branching on $\lambda \in [0, \pi/\Delta t)$, and when $\lambda\Delta t$ leaves this range the aliasing is the caller's responsibility. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"qpca"` |
| `readout_register` / `decoder` | `"phase"` / `"eigenvalue_from_phase"` |
| `copies` | total copy count $2^{\text{precision}} - 1$ |
| `step_time` | single-step evolution time $\Delta t$ |
| `reference` | `"Lloyd-Mohseni-Rebentrost 2014, Nature Physics 10, 631"` |

## Implementation notes

The generation chain is: optional system initial-state preparation → all copies filled from all-zero by the preparation oracle → `phase` prepared in uniform superposition → controlled-power expansion of the phase estimation: the $b$-th bit of `phase` controls $2^b$ consecutive partial swaps (copies are consumed one by one in a fixed cursor order, using up exactly $2^p - 1$ of them) → {obj}`inverse_qft(precision) <oracq.algorithms.common.fourier.inverse_qft>` applied to `phase`. Register layout `system(swap_width) | copies((2^p-1)\cdot w) | phase(p)`, where $w$ is the preparation width.

The bitwise exact decomposition of the partial swap (XX+YY+ZZ = 2·SWAP − I) is shared with the exponentiation primitive; see [Density Matrix Exponentiation](density-matrix-exponentiation.md) for details.

Design decision and applicability boundary: the LMR first-order error at large $\Delta t$ manifests as peak broadening rather than peak shifts — decoding by peak position (mode) remains exact and the circular mean of the distribution stays within tolerance (see the validation approach), so witnesses are allowed to read out directly by peak position. The decoding branch is limited to $\lambda\Delta t \in [0, \pi)$; larger eigenvalues require the caller to pre-scale $\Delta t$.

## Validation approach

Category C3 (probability-distribution semantics; acceptance criteria in `../development/validation-plan.md` §2): the peak positions of the readout distribution must cross-check against the closed-form spectrum of $\rho$. Three layers of evidence:

- Structure: `tests/core/test_qpca.py:QpcaTests.test_invalid_inputs_fail_at_generation` covers all entry violations: precision out of range (0 / 7), step_time = 0, swap_width out of range, mismatched system width, out-of-range decoded values, and so on.
- Numerical: `test_pure_state_eigenvalues` — $\rho = |+\rangle\langle+|$ (eigenvalues 1, 0): with system input $|+\rangle$ (the copy state itself; the partial swap acts on SWAP-symmetric eigenstates and is exact for any $\Delta t$), take precision = 3, $\Delta t = \pi/4$ ($\varphi = 7/8$ falls on the QPE grid); the readout is deterministic (places = 9) and decodes to $\lambda = 1$ (places = 12); with input $|-\rangle$ the peak broadens, yet the distribution's circular mean is still centered on $\lambda = 0$ (delta = 0.15). `test_mixed_state_eigenvalue_via_purification` — $\rho = \mathrm{diag}(0.75, 0.25)$ adapted via purification, with system input $|0\rangle$ reading out the dominant eigenvalue: peak-position decoding is exactly 0.75 (places = 12, peak height > 0.4), circular-mean delta = 0.15.
- Binding: the mixed-state DM input is adapted to SP and enters the circuit via {obj}`gate_purification(...).as_state_preparation() <oracq.algorithms.input_model.density.gate_purification>`; covered by `test_mixed_state_eigenvalue_via_purification`.

## Known gaps and planned stages

No known gaps; the stage V1 witnesses are complete (pure/mixed-state eigenvalue peak cross-checks + $\Delta t$ first-order-rate and copies-scaling witnesses on the exponentiation side). Per the V3 plan in validation-plan §5, QPCA still needs registered catalog cases and cross-checks on real backends.

## Related links

- Source: `src/oracq/algorithms/qml/qpca.py`
- Same-module page: [Density Matrix Exponentiation](density-matrix-exponentiation.md)
- API reference: [QPCA quantum principal component analysis](../../api/algorithms/qml/qpca.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-level numerical experiments are in `tests/verification/verify_misc_algorithms.py` (the misc_algorithms group), all executed on real backends.

**Experiment design**: three QPCA readout cases, with the eigenvalue oracle being numpy `eigvalsh` (independent of the tested implementation). (a) pure state $\rho = \lvert+\rangle\langle+\rvert$, system input $\lvert+\rangle$: precision = 3, $\Delta t = \pi/4$ ($\varphi = 7/8$ falls on the QPE grid), 11 qubits in total; (b) mixed state $\rho = \mathrm{diag}(0.75, 0.25)$ adapted via purification, system input $\lvert 0\rangle$: $\Delta t = 2\pi/6$ ($\lambda = 0.75$ exactly on the grid), 18 qubits in total; (c) the same $\rho$ with system input $\lvert 1\rangle$ reads the secondary eigenvalue $\lambda = 0.25$ (not on the 3-bit grid at this $\Delta t$; compared via the circular mean). Backend paths: `reference`, `rir-pysparq`, `adapter-pysparq`, `originir-ext` (UniQC full amplitude), with phase distributions cross-checked by TVD.

**Key metrics**:

| Case | Scale | Paths | peak-decoded λ | peak height | circular-mean λ | cross-backend TVD |
|---|---|---|---|---|---|---|
| qpca-pure-plus-deterministic | 11 qubits, 7 copies | distribution cross-check on four paths | 1 (error 0, deterministic readout) | 1.0000 | — | ≤ 1.2e-16 |
| qpca-mixed-primary | 18 qubits, 7 copies | distribution cross-check on four paths | 0.75 (error 0) | 0.4997 | 0.8306 (deviation 0.081) | ≤ 1.0e-15 |
| qpca-mixed-secondary | 18 qubits, 7 copies | distribution cross-check on four paths | 0.75 (off-grid information item) | 0.2874 | 0.3909 (deviation 0.141, ≤ 0.15) | ≤ 2.4e-15 |

The LMR first-order error at large $\Delta t$ only broadens the peaks without shifting their positions: the dominant eigenvalue is decoded exactly by the mode, and the secondary eigenvalue's circular-mean deviation (0.141) stems from broadening skew, consistent with the tolerance adopted in the core tests.

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

Artifacts: `out/verification/misc_algorithms.json` (all 24 cases pass; this page corresponds to the three `qpca-*` cases).
