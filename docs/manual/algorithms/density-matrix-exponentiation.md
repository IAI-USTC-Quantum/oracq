# Density Matrix Exponentiation

**English** · <a href="../../../zh/manual/algorithms/density-matrix-exponentiation.html">简体中文</a>

> Category C3 · Module [`oracq.algorithms.qml.qpca`](../../api/algorithms/qml/qpca.rst) · Stage V1

## Overview

Given a supply of copies of a density matrix $\rho$, approximately implements $e^{-i\rho t}$ on a system register. The implementation basis is the LMR protocol of Lloyd–Mohseni–Rebentrost 2014 ("Quantum principal component analysis", Nature Physics 10, 631): exploiting the eigenvalue structure of the SWAP operator, one $\rho$ copy and the system undergo a partial swap $e^{-i\Delta t\cdot\mathrm{SWAP}}$; after the copy is discarded, the effective channel on the system is the first-order approximation of $e^{-i\rho\Delta t}$, with single-step error $O(\Delta t^2)$; chaining `copies` copies gives the total time $t = \mathrm{copies}\cdot\Delta t$, with the first-order error $O(t\cdot\Delta t)$ decreasing linearly in the number of copies.

This primitive is the source of the unitary step for the phase estimation in [Quantum Principal Component Analysis](qpca.md).

## Interface and input model

```python
density_matrix_exponentiation(preparation, *, time, copies, swap_width=None, name=None)
```

API entry point: {obj}`density_matrix_exponentiation <oracq.algorithms.qml.qpca.density_matrix_exponentiation>`

- `preparation`: the state preparation of the $\rho$ copies ({obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>`); the input model is SP + QRAM: pure states use {obj}`gate_state_prep(amplitudes) <oracq.algorithms.input_model.oracles.gate_state_prep>`, mixed states are adapted via `density.gate_purification(rho)`'s `PurificationAccess.as_state_preparation()` (target = system ⊕ environment, zero work); the preparation must have zero work, and non-zero work raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time.
- `time`: the total evolution time $t$, which must be a positive finite real number.
- `copies`: the number of copies, range 1..63; each step $\Delta t = t/\mathrm{copies}$.
- `swap_width`: the prefix bit width participating in the swap, defaulting to the preparation's entire target; in purification scenarios take the system bit width of $\rho$, with the environment bits left in the copies and not swapped (equivalent in effect to taking the partial trace).

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `system` (swap_width bits) and `copies` (copies × preparation width). The system's input state is prepared by the caller; the copy registers are filled from all-zero by the preparation oracle. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"density_matrix_exponentiation"` |
| `copies` | number of copies |
| `step_time` | single-step duration $\Delta t = t/\mathrm{copies}$ |
| `error_scaling` | `"O(time * step_time)"` (first-order error scaling) |
| `reference` | `"Lloyd-Mohseni-Rebentrost 2014, Nature Physics 10, 631"` |

## Implementation notes

The partial swap is decomposed exactly by the identity $\mathrm{XX}+\mathrm{YY}+\mathrm{ZZ} = 2\cdot\mathrm{SWAP} - I$: the Pauli-pair rotations $e^{-i(\Delta t/2)\,P\otimes P}$ on the three axes commute mutually and are applied axis by axis, bit by bit; the single-axis rotation is parity-based (XOR–Rz–XOR), the $y$ axis switches basis via $R_x(\pi/2)$ conjugation, and the $I$ component is corrected by a global phase $-\Delta t/2$.

The generation chain is: the $k$-th copy is filled from all-zero by the preparation oracle → a partial swap of angle $\Delta t$ is applied bitwise between `system` and the first swap_width bits of that copy. Register layout `system | copies`, with the $k$-th copy occupying the segment $[k\cdot w, (k+1)\cdot w)$ ($w$ is the preparation width).

Design decision: the algorithm only requires the copies to be repeatedly preparable — multiple copies are just multiple calls to the preparation oracle, and the matrix elements of $\rho$ need not be given explicitly, which is the resource premise of QRAM-based QML. Discarding the copies amounts to taking the partial trace over the copy registers, with the purification's environment bits left entangled inside the copies. Applicability boundary: first-order approximation; for higher accuracy increase copies (smaller $\Delta t$); see the `error_scaling` attribute for the error scaling.

## Validation approach

Category C3 (probability-distribution semantics; acceptance criteria in `../development/validation-plan.md` §2): after discarding the copies, the system state must approach the closed-form expectation in trace distance. Three layers of evidence:

- Structure: `tests/core/test_qpca.py:QpcaTests.test_invalid_inputs_fail_at_generation` covers generation-time violations at the exponentiation entry such as `time = 0`, `copies = 0`, and non-`StatePreparation` inputs.
- Numerical: `test_small_step_first_order_accurate` — $\rho = |+\rangle\langle+|$, a single step $\Delta t = 0.05$ applied to $|0\rangle$; after discarding the copy, the trace distance from the exact $e^{-it\,|+\rangle\langle+|}$ action is $< 0.6\cdot\Delta t^2$ (the first-order convergence rate in $\Delta t$, measured coefficient about 0.53); `test_error_halves_with_copies` — fixed total time $t = 0.4$, copies = 1, 2, 4, the trace distance decreases level by level with each level $< 0.6\times$ the previous (doubling copies halves $\Delta t$, roughly halving the error).
- Binding: the mixed-state DM input is adapted to SP and enters the circuit via {obj}`gate_purification(...).as_state_preparation() <oracq.algorithms.input_model.density.gate_purification>`, witnessed by `QpcaTests.test_mixed_state_eigenvalue_via_purification` (see [Quantum Principal Component Analysis](qpca.md)).

## Known gaps and planned stages

No known gaps; the stage V1 witnesses are complete (the $\Delta t$ first-order rate + the error-halves-when-copies-double scaling + the purification adaptation path). Per the V3 plan in validation-plan §5, the QPCA module still needs registered catalog cases and cross-checks on real backends.

## Related links

- Source: `src/oracq/algorithms/qml/qpca.py`
- Same-module page: [Quantum Principal Component Analysis](qpca.md)
- API reference: [QPCA quantum principal component analysis](../../api/algorithms/qml/qpca.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-level numerical experiments are in `tests/verification/verify_misc_algorithms.py` (the misc_algorithms group), all executed on real backends.

**Experiment design**: two density matrix exponentiation cases — (a) pure state $\rho = \lvert+\rangle\langle+\rvert$ (`gate_state_prep`, total time $t = 0.4$); (b) mixed state $\rho = \mathrm{diag}(0.75, 0.25)$ (adapted via `gate_purification` purification, `swap_width = 1`, $t = 0.5$). The system initial state is $\lvert 0\rangle\langle 0\rvert$ in both, with copies = 1/2/4 (at most 9 qubits). The classical oracle is the exact evolution $e^{-i\rho t}\sigma e^{i\rho t}$ from a numpy eigendecomposition (fully independent of the tested implementation). Backend paths: `reference` (built-in reference executor), `rir-pysparq` (PySparQ native RIR), `adapter-pysparq` (PySparQ adapter), `originir-ext` (UniQC full-amplitude state vector), cross-checked at two levels: amplitude by amplitude, and on the density matrix after the partial trace.

**Key metrics**:

| Case | Scale | Paths | trace distance (copies = 1/2/4) | error ratio | cross-backend deviation |
|---|---|---|---|---|---|
| dm-exponentiation-pure-convergence | t = 0.4, ≤5 qubits | full-amplitude cross-check on four paths | 8.55e-2 / 4.34e-2 / 2.20e-2 | 0.508 / 0.507 | ≤ 2.1e-16 |
| dm-exponentiation-mixed-convergence | t = 0.5, ≤9 qubits | full-amplitude cross-check on four paths | 5.75e-2 / 2.97e-2 / 1.52e-2 | 0.516 / 0.512 | ≤ 3.3e-16 |

When copies double (i.e. $\Delta t$ halves), the trace-distance ratio stabilizes around 0.51, directly observing the first-order convergence scaling $O(t\cdot\Delta t)$ of the LMR protocol; the four backend paths agree to machine precision.

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

Artifacts: `out/verification/misc_algorithms.json` (all 24 cases pass; this page corresponds to the two `dm-exponentiation-*` cases).
