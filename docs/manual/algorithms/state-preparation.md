# State Preparation

**English** · <a href="../../../zh/manual/algorithms/state-preparation.html">简体中文</a>

> Category C5 · Module [`oracq.algorithms.input_model.oracles`](../../api/algorithms/input_model/oracles.rst) · Stage V4

## Overview

The standard interface of the amplitude-encoding input model (SP in the input-model vocabulary): the paradigm `state_prep_isometry` defines the preparation operation $V$ as an isometry departing from the zero-state subspace — $V|0\rangle|0\rangle = |\psi\rangle|0\rangle$ — with the work register promised to be uncomputed; inverse and controlled operations require the existence of a reversible extension $U$. The amplitude vector decomposes as a binary tree: each node writes the weights of the left and right sub-blocks as

$$
\theta = 2\,\arctan\frac{\sqrt{w_{\mathrm{right}}}}{\sqrt{w_{\mathrm{left}}}} ,
$$

implemented through one of two bindings: a multiplexed Ry rotation tree or a QRAM angle table.

## Interface and input model

```python
basis_state(width, value=0, *, work_width=0)
uniform_state(width, *, work_width=0)
gate_state_prep(amplitudes, *, work_width=0, name=None)
qram_state_prep(width, angle_width=8)
qram_state_angles(amplitudes, angle_width=8)
abstract_state_prep(name, width, work_width=0, *, reversible=True)
StatePreparation.from_unitary(operation, *, target=None, work=None, clean_work=False)
```

API entry points: {obj}`basis_state <oracq.algorithms.input_model.oracles.basis_state>`, {obj}`uniform_state <oracq.algorithms.input_model.oracles.uniform_state>`, {obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>`, {obj}`qram_state_prep <oracq.algorithms.input_model.oracles.qram_state_prep>`

- {obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>`: any complex amplitude vector (length a power of two, nonzero norm); returns a gate-level rotation-tree preparation.
- {obj}`qram_state_prep <oracq.algorithms.input_model.oracles.qram_state_prep>` / {obj}`qram_state_angles <oracq.algorithms.input_model.oracles.qram_state_angles>`: the QRAM-resource version. The former takes only the width and the angle width; the latter generates the angle-table dictionary from an amplitude vector (keys are tree-node indices) and supplies it as memory at execution time; currently only non-negative real amplitudes are accepted.
- {obj}`abstract_state_prep <oracq.algorithms.input_model.oracles.abstract_state_prep>`: an open declaration; `reversible=False` declares that adjoint/controlled are unsupported.
- `from_unitary`: explicitly casts a plain unitary operation into the $U|0\rangle$ initial-state role; a nonzero-width work requires the `clean_work=True` promise.

Returns {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` (an {obj}`OracleView <oracq.algorithms.input_model.contracts.OracleView>`; `state_preparation()` returns itself). Attributes:

| Attribute | Meaning |
|---|---|
| `operation` | the underlying {obj}`Operation <oracq.infrastructure.builder.Operation>` (registers are exactly `target`, `work`) |
| `width` / `work_width` | target / work bit widths |
| `capabilities` | adjoint / controlled capabilities (from the declaration or annotation) |

Module attributes: `zero_input=True`, `clean_work`; the gate version has `implementation="multiplexed_rotations"`, the QRAM version has `implementation="qram_rotation_tree"` plus `qram_queries=2*width`.

## Implementation notes

Gate version: the root node applies an unconditional Ry, the remaining nodes are controlled on the high-bit prefix of target; complex phases are compensated pointwise for each non-zero-phase basis state through a controlled `global_phase`, so amplitudes may be complex. At least one target bit is required (single-amplitude vectors are rejected).

QRAM version: registers `target(width)`, `work(address_width + angle_width)` (`address_width = max(1, width)`). Layer by layer, the high-bit prefix of target is XORed into the address, a tree offset is added, and the angle table is queried; $R_y(2\pi k / 2^{\text{angle\_width}})$ is applied controlled on each bit of the angle word, and then the address is uncomputed by replaying in reverse; two queries per layer. The angle resolution $2\pi/2^{\text{angle\_width}}$ is the source of quantization error; distribution agreement with the gate version is cross-checked on the prepare-select side at delta = 0.02 (see the next section).

The `target=None` branch of `from_unitary` treats the entire register space as target (work width 0); {obj}`annotate <oracq.algorithms.input_model.oracles.annotate>` grants the `zero_input` / `clean_work` promises to this paradigm by default. Applicability boundary: the gate version's gate count grows exponentially with the dimension; the QRAM version's angle table must be supplied by the caller alongside the circuit.

## Validation approach

Category C5 (data-access layer; acceptance criteria in `../development/validation-plan.md` §2). Three layers of evidence:

- Structural: `tests/core/test_contracts.py:OracleContractTests.test_explicit_unitary_state_prep_adapter` (a Hadamard adapted through `from_unitary` satisfies the {obj}`StatePreparationProtocol <oracq.algorithms.input_model.interfaces.StatePreparationProtocol>` contract; the unitary itself can likewise serve as input); `test_unitary_adapter_requires_work_promise` (without the `clean_work` promise a {obj}`ContractError <oracq.algorithms.input_model.contracts.ContractError>` is raised); `tests/core/test_open_ir.py:OpenIRTests.test_isometry_requires_declared_adjoint_capability` (a `reversible=False` declaration is rejected in adjoint contexts).
- Numerical: `OracleContractTests.test_dirty_work_is_rejected_by_algorithm` — a preparation with `clean_work=False` is rejected by the {obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>` contract before the kernel runs (the numerical anchor of the oracles.py row in validation-coverage); the numerical behavior of the QRAM rotation tree is witnessed by `tests/core/test_qlss_input_models.py:QLSSInputTests.test_tree_preparation_queries_each_layer_coherently` (8 Load queries, normalized amplitudes cross-checked at delta = 0.07).
- Binding: `OracleContractTests.test_capability_restrictions_survive_annotation_and_wrapping` (capability restrictions survive annotation and wrapping); `test_legacy_state_prep_clean_work_in_bind` (a legacy preparation lacking the `clean_work` attribute can still be bound to close the program); distribution agreement of the gate / QRAM implementations is indirectly witnessed by `tests/core/test_prepare_select.py:PrepareSelectTests.test_gate_and_qram_prepare_bindings_agree` through assembly via {obj}`gate_prepare <oracq.algorithms.common.prepare_select.gate_prepare>` / {obj}`qram_prepare <oracq.algorithms.common.prepare_select.qram_prepare>` (gate version places = 11, QRAM version delta = 0.02).

## Known gaps and planned stages

Three-binding consistency parameterization (to be rolled out in V4): unified parameterized cross-checks across the abstract / gate / QRAM layers are not rolled out; the existing bindings are each witnessed independently. Consistent with the oracles.py row of `validation-coverage.md`.

## Related links

- Source: `src/oracq/algorithms/input_model/oracles.py`
- Same-group pages: [PREPARE–SELECT decomposition](prepare-select.md), [Alias sampling preparation](alias-preparation.md), [XOR database](xor-database.md)
- API reference: [Oracle declarations and implementations](../../api/algorithms/input_model/oracles.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

The results in this section were produced by `tests/verification/verify_stateprep.py` on real backends (the reference built-in reference executor, the rir-pysparq native RIR interpreter, the adapter-pysparq event adapter, and originir-ext + UniQC full-amplitude state vector — four independent paths in total). The classical oracle is independent of the implementation under test: the target vector itself, an independent numpy implementation of Pauli action and index permutations, and a classical rotation-tree expansion of the quantized angle table.

**Experiment design** (14 cases relevant to this page):

- `gate_state_prep`: dense complex amplitude vectors of widths 1–4, dense width 6, sparse width 8 (5 nonzero amplitudes), compared amplitude by amplitude against the target vector (max_error / fidelity); in addition, UniQC `Circuit.to_matrix` extracts the full unitary of the width-3 circuit and checks that its first column equals the target vector.
- Combinators (`state_preparation.py`): {obj}`extend_initial <oracq.algorithms.common.state_preparation.extend_initial>` (3+2 bits, high bits stay |0⟩); {obj}`apply_be_to_state <oracq.algorithms.common.state_preparation.apply_be_to_state>` (the diagonal block encoding $D[x,x]=\cos(\pi T[x]/4)$ acting on a 2-bit prepared state, the signal==0 block compared against $D|\psi\rangle$ and the success probability checked); {obj}`select_subspace <oracq.algorithms.common.state_preparation.select_subspace>` (a Pauli-word state oracle in two layouts, the full amplitude dictionary compared against a numpy oracle, with the signal==0 block being the postselected subvector whose high bits equal `high_value`).
- `qram_state_prep`: widths 2/3, angle widths 8/10. Implementation error (circuit vs quantized angle-tree oracle) and method error (quantized tree vs exact vector) are reported separately; without quantization the oracle agrees amplitude by amplitude with the gate-version circuit (self-consistency check deviation 0).
- Known backend issue: the PySparQ RIR interpreter (pysparq 0.1.2.dev16) executes `add_const` on a sliced reinterpret view (a Span operand in RIR text) incorrectly; the minimal reproduction case `backend-rir-sliced-add-const` has an rir deviation of 1.0, while the reference / adapter-pysparq / originir-ext paths all uncompute correctly. The address bookkeeping of `qram_state_prep` is affected by this defect, so its acceptance criteria rest on the three correct paths above and the rir deviation is recorded as informational only (the script comments give the reasoning for the workaround).

**Key metrics**:

| Case | Scale | Paths | Metric values |
|---|---|---|---|
| gate-state-prep-complex-w1..w4 | dimension 2–16 | four paths | max_error ≤ 1.6e-16, fidelity ≥ 1 − 2e-16 |
| gate-state-prep-dense-w6 | dimension 64 | four paths | max_error 8.8e-17, fidelity 1.0 |
| gate-state-prep-sparse-w8 | dimension 256 (5 nonzero) | four paths | max_error 1.1e-16, fidelity 1.0 |
| gate-state-prep-unitary-column-w3 | 8×8 unitary | to_matrix | max_error 1.1e-16, fidelity 1 − 2e-16 |
| extend-initial-w3e2 | 5 qubits | four paths | max_error 1.1e-16 |
| apply-be-diagonal-w2 | 2+3 qubits | four paths | block error 1.2e-16; success probability 0.4740002825 (differs from oracle by 3e-16) |
| select-subspace-w3-k2-hv1 / k1-hv3 | 3 qubits | ref + rir + originir | max_error ≤ 1.3e-16, postselection block error ≤ 1.3e-16 |
| qram-state-prep-w2-a8 | 12 qubits | ref + adapter + originir | implementation error 1.1e-16; method error 1.86e-3; fidelity 0.9999928 |
| qram-state-prep-w3-a10 | 16 qubits | ditto | implementation error 3.3e-16; method error 7.8e-4; fidelity 0.9999987 |
| backend-rir-sliced-add-const | 4 qubits | ref + adapter (rir informational only) | rir_deviation 1.0 (known backend issue) |

The method error is consistent with quantization at the angle resolution $2\pi/2^{\text{angle\_width}}$ (8 bits corresponds to roughly the 2e-3 scale); the implementation error is at machine precision, and both implementations are individually correct.

**Reproduction command**:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_stateprep.py
```

**Artifact path**: `out/verification/stateprep.json` (all 37 cases pass, with all metric values).
