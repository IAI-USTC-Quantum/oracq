# VQE Pauli Measurements

**English** · <a href="../../../zh/manual/algorithms/vqe.html">简体中文</a>

> Category C4 · Module [`oracq.algorithms.optimization.variational`](../../api/algorithms/optimization/variational.rst) · Stage V3

## Overview

The variational quantum eigensolver (VQE) uses a parametrized state $\lvert\psi\rangle$ to approximate the ground-state energy of a Hamiltonian $H=\sum_k c_k P_k$, where $P_k$ are Pauli words (bitwise tensor products of I/X/Y/Z) and $c_k$ are real coefficients. The quantum side only needs to repeatedly measure $\langle P_k\rangle$ for each $P_k$; the classical side forms the weighted sum $E=\sum_k c_k\langle P_k\rangle$ for the energy and drives the parameter optimization.

This module implements the measurement-circuit family: {obj}`pauli_measurement <oracq.algorithms.optimization.variational.pauli_measurement>` rotates the prepared state into the requested Pauli measurement basis, and {obj}`vqe_measurements <oracq.algorithms.optimization.variational.vqe_measurements>` batch-generates the circuits for a whole Hamiltonian term list. Energy aggregation and parameter optimization happen on the classical side; the module itself contains no optimizer.

## Interface and input model

```python
vqe_measurements(preparation, terms)
pauli_measurement(preparation, word)
```

API entry points: {obj}`vqe_measurements <oracq.algorithms.optimization.variational.vqe_measurements>`, {obj}`pauli_measurement <oracq.algorithms.optimization.variational.pauli_measurement>`

- `preparation`: an ansatz with parameters already instantiated, or another state preparation ({obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>`); the input model is SP; it must satisfy the zero-input, clean-workspace contract (checked by {obj}`checked_state_preparation <oracq.algorithms.input_model.interfaces.checked_state_preparation>`; violations raise {obj}`ContractError <oracq.algorithms.input_model.contracts.ContractError>`).
- `terms`: a list of `(real coefficient, Pauli word)`; the coefficients must be finite real numbers, the Pauli words must be as wide as the state (I/X/Y/Z), and an empty list raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`.
- `word`: an I/X/Y/Z string of the same width; the first character corresponds to the least significant bit.

`vqe_measurements` returns a tuple whose elements are `(coefficient, measurement Operation)`. Each measurement circuit's registers are `target` (the state width) and `work` (the preparation's workspace). Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"pauli_measurement"` |
| {obj}`pauli_word <oracq.algorithms.input_model.block_encoding.pauli_word>` | the Pauli word being measured |
| `readout_register` | `"target"` |

## Implementation notes

Generation chain: the preparation oracle is first invoked to fill `target` (the work semantics are inherited from the preparation), then the bit-by-bit basis change is applied — a Y bit first receives a phase gate of angle $-\pi/2$ ($S^\dagger$) followed by a Hadamard, an X bit only a Hadamard, an I bit nothing. After the change, `target` is measured in the standard Z basis; reading out the Z parity of the non-I bits suffices to estimate the expectation of that Pauli word.

Each Pauli word generates an independent circuit; repeated measurement and statistics are organized by the caller; same-name terms can reuse samples from the same circuit, and the module itself performs no qubit-wise commuting grouping. Applicability boundary: only the measurement circuits are provided; confidence intervals for expectation values and parameter-optimization strategies are not part of this module.

## Validation approach

Category C4 (heuristic/optimization semantics; acceptance criteria in `../../development/validation-plan.md` §2). Three layers of evidence:

- Structure: construction and attribute assertions in `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`; the same family's `test_bad_inputs_fail_at_generation` covers generation-time violations of this module's entry points (Pauli-word width/character checks are executed at generation time in the source; there is no separate negative test).
- Numerical: `AlgorithmExpansionTests.test_vqe_measurement_of_y_eigenstate` — prepare the Y eigenstate $\lvert y_+\rangle$ (Hadamard followed by a $\pi/2$ phase gate), generate the measurement circuit for the term `(2, "Y")`, the coefficient is echoed back as 2, and `target` reads 0 with probability 1, i.e. $\langle Y\rangle=+1$ is recovered exactly.
- Binding: this module has no separate binding witness (the input is already a concrete state preparation, with no abstract slots).

## Known gaps and planned stages

Consistent with the `variational.py` row of `validation-coverage.md`: the only gap registered on that row is a strong witness for "reaching the optimal cut on small instances" (it lands on the QAOA side, to be wired in at V3); the VQE measurement circuits themselves have no registered gap.

## Related links

- Source: `src/oracq/algorithms/optimization/variational.py`
- Same-module pages: [MaxCut QAOA](qaoa-maxcut.md), [Hardware-Efficient Ansatz](variational-ansatz.md)
- API reference: [variational algorithm circuits](../../api/algorithms/optimization/variational.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_misc_algorithms.py` (misc_algorithms group), all executed on real backends.

**Experimental design**: a 2-qubit hardware-efficient ansatz is adapted into a preparation oracle via `StatePreparation.from_unitary`; Pauli measurement circuits are generated for a five-term Hamiltonian (ZI, IZ, XX, YY, ZZ with coefficients 0.5/−0.3/0.7/0.2/−0.1); the Z-parity expectation of the non-I bits of each measurement circuit is computed from the exact state vector (no sampling) and independently compared against $\langle\psi\lvert P\rvert\psi\rangle$ of the numpy tensor-product operators, then aggregated into the weighted total energy. Backend paths: `reference` (also matching the numpy oracle within machine precision).

**Key metrics**:

| Case | Scale | Path | Max per-term error | Total energy (measured / exact) | Energy error |
|---|---|---|---|---|---|
| vqe-pauli-expectations | 2 qubits, 5 terms | reference | 3.3e-16 | 0.4863695277137 / 0.4863695277137 | 1.1e-16 |

**Reproduce**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

Artifact: `out/verification/misc_algorithms.json` (all 24 cases pass; this page corresponds to the `vqe-pauli-expectations` case).
