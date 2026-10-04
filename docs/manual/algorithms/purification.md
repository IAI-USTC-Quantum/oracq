# Purification Access

**English** · <a href="../../zh/manual/algorithms/purification.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.input_model.density`](../../api/algorithms/input_model/density.rst) · Stage V1

## Overview

The core abstraction of the density matrix (DM) input model is **purification access**: access to a density matrix $\rho$ is defined as the quantum operation $U$ preparing its purification, $U|0\rangle = |\psi_\rho\rangle$, acting on the two registers system and environment, such that after tracing out the environment

$$
\operatorname{Tr}_{\mathrm{env}} |\psi_\rho\rangle\langle\psi_\rho| = \rho .
$$

Every density matrix admits a purification with environment as wide as system, so this view is a stable interface that the downstream algorithms (B2 quantum SDP, [Gibbs state preparation](gibbs-state.md)) rely on. Following the repository's three-layer paradigm, the module provides three constructions: an abstract open declaration, a gate witness for explicit small matrices, and a pure-state adapter.

## Interface and input model

```python
abstract_purification(name, width, environment_width=None, *, reversible=True)
gate_purification(rho, *, name=None)
maximally_mixed_purification(width, *, name=None)
PurificationAccess.from_state_preparation(preparation)
```

API entry points: {obj}`abstract_purification <oracq.algorithms.input_model.density.abstract_purification>`, {obj}`gate_purification <oracq.algorithms.input_model.density.gate_purification>`, {obj}`maximally_mixed_purification <oracq.algorithms.input_model.density.maximally_mixed_purification>`

- {obj}`abstract_purification <oracq.algorithms.input_model.density.abstract_purification>`: openly declares a purification access slot (the abstract layer of the DM input model); `environment_width` defaults to `width`; the program becomes closed after a witness implementation is bound via `linking.bind`.
- {obj}`gate_purification <oracq.algorithms.input_model.density.gate_purification>`: a gate witness for an explicit small density matrix. `rho` must be a Hermitian positive semidefinite matrix of unit trace whose dimension is a power of two within 2..16.
- {obj}`maximally_mixed_purification <oracq.algorithms.input_model.density.maximally_mixed_purification>`: a purification generator for the maximally mixed state $I/2^n$ ($n$ Bell pairs, width ≤ 32).
- `from_state_preparation`: a pure state is a trivial purification; this adapts a {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` that cleans its work back to zero into a purification access (work plays the role of environment).

All four entry points return a {obj}`PurificationAccess <oracq.algorithms.input_model.density.PurificationAccess>` (`oracle_kind = "purification_access"`, operation signature `("system", "environment")`). Module attributes and attributes:

| Attribute | Meaning |
|---|---|
| `width` / `environment_width` | system / environment register widths |
| `operation` | the underlying {obj}`Operation <oracq.infrastructure.builder.Operation>`; module attributes include `density_model="purification_access"`, `zero_input=True` |
| {obj}`as_state_preparation() <oracq.algorithms.input_model.interfaces.as_state_preparation>` | viewing the whole as a `StatePreparation` on system⊕environment (for B2 composition) |
| `describe()` | oracle description (`main_qubit=width`, `anc_qubit=environment_width`) |

Accompanying classical small-matrix tools: {obj}`partial_trace <oracq.algorithms.input_model.density.partial_trace>` (partial trace), {obj}`trace_distance <oracq.algorithms.input_model.density.trace_distance>` (trace distance $T(\rho,\sigma) = \lVert\rho-\sigma\rVert_1/2$), and {obj}`gibbs_state <oracq.algorithms.input_model.density.gibbs_state>` (classical reference Gibbs state).

## Implementation notes

The gate witness first performs a cyclic Jacobi eigendecomposition $\rho = \sum_j p_j |v_j\rangle\langle v_j|$ (for complex Hermitian matrices, a diagonal phase rotation precedes real Jacobi elimination; eigenvalues in descending order), then takes the purification $|\psi_\rho\rangle = \sum_j \sqrt{p_j}\,|v_j\rangle_{\mathrm{s}} |j\rangle_{\mathrm{e}}$, whose amplitude vector is prepared on the concatenated register by the multi-rotation tree of {obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>`; the basis-state index convention is system | (environment << system_width). Positive semidefiniteness is checked at generation time; negative eigenvalues below $-10^{-7}$ raise an `INPUT_PROMISE` failure.

`maximally_mixed_purification` applies H and CNOT bit by bit to produce $n$ pairs of $|\Phi^+\rangle$; `from_state_preparation` relies on the `clean_work=True` promise of the preparation contract and refuses the adaptation when it is not promised. Applicability boundary: the gate witness supports only small density matrices of dimension up to 16 (small instances for validation); large-scale DM input requires providing an implementation of the `abstract_purification` declaration yourself and binding it.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2): the reduced density matrix after the partial trace must agree with the target $\rho$ within tolerance. Three layers of evidence:

- Structure: the width and attribute assertions of the cases in `tests/core/test_density.py:PurificationTests` (e.g. `test_gate_purification_recovers_rho` checks `access.width == 1`); `GibbsTests.test_invalid_inputs_fail_at_generation` covers a non-unit trace, non-Hermitian input, and an illegal dimension for `gate_purification`, as well as `maximally_mixed_purification(0)`, mismatched `partial_trace` lengths, inconsistent `trace_distance` dimensions, and similar violations.
- Numerical: `PurificationTests.test_gate_purification_recovers_rho` traces out environment after simulation and cross-checks element by element against the 2×2 mixed state `RHO` (places = 7); `test_maximally_mixed_purification` cross-checks $I/4$ and verifies the Schmidt structure (nonzero amplitudes appear only on system == environment basis states); `test_pure_state_adapter` cross-checks the reduced matrix of the pure state $\sqrt{0.3}|0\rangle + \sqrt{0.7}|1\rangle$; `test_classical_tools` checks the analytic values of `gibbs_state` and `trace_distance` (places = 12).
- Binding: `PurificationTests.test_abstract_purification_binds` — after the abstract declaration is bound to the `gate_purification` witness via {obj}`bind <oracq.infrastructure.linking.bind>`, the partial trace still returns to `RHO`, covering the consistency of the abstract / gate layers.

## Known gaps and planned stages

No known gaps; the stage V1 witnesses are complete (partial-trace cross-check + Schmidt structure + binding consistency).

## Related links

- Source: `src/oracq/algorithms/input_model/density.py`
- Algorithms in the same module: [Gibbs state preparation](gibbs-state.md)
- API reference: [Density-matrix input model and Gibbs states](../../api/algorithms/input_model/density.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

The results in this section were produced on real backends by `tests/verification/verify_stateprep.py` (reference, rir-pysparq, adapter-pysparq, originir-ext + UniQC — four independent paths).

**Experiment design** (5 cases): for the purification prepared on each path, the system reduced density matrix is computed with an independently implemented numpy partial trace (without calling this module's `partial_trace` / `trace_distance`), and the eigenvalues of the difference matrix then give the trace distance $T(\rho_{\text{sim}}, \rho)$ and the element-wise error. Instances: a 2×2 complex Hermitian mixed state (off-diagonal entries with imaginary parts), 4×4 and 16×16 Gaussian random density matrices ($G^\dagger G$ normalized, deterministic seeds), the 3-Bell-pair purification of the maximally mixed state $I/8$ (with a Schmidt structure check: nonzero amplitudes appear only on system == environment basis states), and the `from_state_preparation` pure-state adapter (the reduced matrix should equal $|\psi\rangle\langle\psi|$).

**Key metrics**:

| Case | Scale (system+environment) | Paths | trace_distance | max_element_error |
|---|---|---|---|---|
| purification-gate-rho2 | 1+1 qubit | four paths | 1.1e-16 | 1.1e-16 |
| purification-gate-rho4 | 2+2 qubit | four paths | 3.2e-11 | 1.9e-11 |
| purification-gate-rho16 | 4+4 qubit | four paths | 8.9e-11 | 1.4e-11 |
| purification-bell-w3 | 3+3 qubit | four paths | 2.2e-16 | 5.6e-17 |
| purification-pure-adapter | 1+0 qubit | four paths | 8.3e-17 | 1.1e-16 |

The residuals (~1e-10) of the 4×4 and 16×16 cases come from the method error of the generation-time cyclic Jacobi eigendecomposition (tolerance 1e-13), not from execution error — pairwise deviations between backends are at machine-precision level; the Bell purification's `off_schmidt_probability` is 0, confirming the Schmidt structure.

**Reproduction command**:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_stateprep.py
```

**Artifact path**: `out/verification/stateprep.json`.
