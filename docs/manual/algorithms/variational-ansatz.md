# Hardware-Efficient Ansatz

**English** · <a href="../../zh/manual/algorithms/variational-ansatz.html">简体中文</a>

> Category C4 · Module [`oracq.algorithms.optimization.variational`](../../api/algorithms/optimization/variational.rst) · Stage V3

## Overview

The hardware-efficient ansatz is a family of parametrized circuits: each layer first applies independent Ry/Rz single-qubit rotations to every qubit, then produces entanglement with a CNOT chain between adjacent qubits. The angle vector is chosen by a classical optimizer so that the objective functional (energy, loss, etc.) is optimal over the states expressible by the ansatz. This module generates the circuit for given angles; parameter optimization and repeated execution are organized by the caller.

## Interface and input model

```python
hardware_efficient_ansatz(width, layers)
```

API entry point: {obj}`hardware_efficient_ansatz <oracq.algorithms.optimization.variational.hardware_efficient_ansatz>`

- `width`: the target bit width, range 1..64.
- `layers`: the angle tensor, of shape `[layer][qubit][Ry, Rz]` — non-empty, exactly width angle pairs per layer, and exactly two finite real angles per pair (radians).

The variational angles are classical data parameterized directly, so the input model is CP. Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` whose only register is `target` (width bits). Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"hardware_efficient_ansatz"` |
| `layers` | the number of layers |

## Implementation notes

Each layer is generated in two segments, "rotations first, then entanglement": every qubit receives `Ry(θ)` and `Rz(φ)` in turn, then a CNOT is applied to each adjacent qubit pair (i, i+1) one by one (the circuit primitive is XOR). The entanglement structure is fixed to a linear chain; there is no inter-layer rearrangement or entanglement-strategy switch, and hardware without chain connectivity gets swaps inserted during backend lowering. Parameter-shape violations (zero layers, a per-layer length unequal to width, an angle pair of length other than 2, non-finite angles) raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time.

With zero angles the whole circuit is the identity; this is the anchor of the current numerical witness. Applicability boundary: the ansatz is only responsible for generating the circuit; expressivity and trainability (e.g. barren plateau) are outside the validation scope, and there is no built-in parameter-initialization strategy.

## Validation approach

Category C4 (heuristic/optimization semantics; acceptance criteria in `../../development/validation-plan.md` §2). Three layers of evidence:

- Structure: `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_bad_inputs_fail_at_generation` covers parameter-layout violations (e.g. too few angle pairs per layer).
- Numerical: `AlgorithmExpansionTests.test_ansatz_zero_angles_and_walk_one_step` — after a single-layer ansatz with all-zero angles acts, the state stays $\lvert 0\ldots0\rangle$ (amplitude exactly 1), pinning down the identity of the rotations and the CNOT chain at zero parameters (the second half of the same test witnesses the cycle walk of `walks.py`, see that module).
- Binding: this module has no separate binding witness (the inputs are classical parameters).

## Known gaps and planned stages

Consistent with the `variational.py` row of `validation-coverage.md`: the only gap registered on that row is a strong witness for "reaching the optimal cut on small instances" (it lands on the QAOA side, to be wired in at V3); the ansatz circuit itself has no registered gap.

## Related links

- Source: `src/oracq/algorithms/optimization/variational.py`
- Same-module pages: [MaxCut QAOA](qaoa-maxcut.md), [VQE Pauli measurements](vqe.md)
- API reference: [variational algorithm circuits](../../api/algorithms/optimization/variational.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_misc_algorithms.py` (misc_algorithms group), all executed on real backends.

**Experimental design**: a hardware-efficient ansatz with width = 3 and two layers totaling 12 non-trivial angles (a Ry/Rz pair per qubit per layer); the classical oracle is a gate-by-gate numpy state vector (an independent implementation of Ry/Rz/CNOT with the exactly-same layer order and least-significant-bit convention). Backend paths: full-amplitude cross-check across the four paths `reference`, `rir-pysparq`, `adapter-pysparq`, `originir-ext`.

**Key metrics**:

| Case | Scale | Path | Max amplitude error | Fidelity | Cross-backend deviation |
|---|---|---|---|---|---|
| ansatz-parameter-intent | 3 qubits, 2 layers | full-amplitude cross-check on four paths | 0 | 1 − 4e-16 | 0 |

**Reproduce**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

Artifact: `out/verification/misc_algorithms.json` (all 24 cases pass; this page corresponds to the `ansatz-parameter-intent` case).
