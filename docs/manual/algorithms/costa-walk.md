# Costa Walk QLSS

**English** · <a href="../../../zh/manual/algorithms/costa-walk.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.qlss.qlss`](../../api/algorithms/qlss/qlss.rst) · Stage V3

## Overview

A QLSS based on a discrete parameterized quantum walk (Costa et al., arXiv:2111.08152; for the input-model discussion see [QFVM input-model review](../../reference/qfvm-input-models.md)): solving $A^{-1}b$ is recast as a sequence of walk operators of the interpolated Hamiltonian $H(s)$ (RHS zero-state reflection, $R(s)$ rotation, controlled $U_A$ / $U_A^\dagger$, and a signal zero reflection), after which a Laurent polynomial filter reads the solution off the walk powers. The walk operator ({obj}`costa_walk <oracq.algorithms.qlss.qlss.costa_walk>`), the schedule ({obj}`schedule <oracq.algorithms.qlss.qlss.schedule>`), the Dolph–Chebyshev filter ({obj}`dolph_chebyshev_plan <oracq.algorithms.qlss.qlss.dolph_chebyshev_plan>` + {obj}`lcu_filter <oracq.algorithms.qlss.qlss.lcu_filter>`), and the problem-layer entry point ({obj}`make_costa_qlss <oracq.algorithms.qlss.qlss.make_costa_qlss>`) are assembled separately. The input model is BE plus SP; sparse input is converted through an explicit adaptation ([Sparse Matrix Block Encoding](sparse-block-encoding.md)), which does not violate the paper's input assumptions.

## Interface and input model

```python
costa_walk(a, bprep, fs)                          # single-step walk operator, returns an Operation
schedule(s, kappa, power=1.5)                      # schedule point s∈[0,1] → f(s)
CostaConfig(steps=2, kappa=4.0, schedule_power=1.5,
            filter_degree=2, filter_attenuation=0.2)
dolph_chebyshev_plan(degree=2, attenuation=0.2)    # → FilterPlan
lcu_filter(walk, plan)                             # coherent LCU filter, returns an Operation
unary_weight_preparation(weights)                  # unary prefix superposition preparation
costa_qlss(a, bprep, config=None, *, filtering=None)  # kernel, returns a StateOracle
make_costa_qlss(config=None)                       # → QLSSProtocol (input_model="block_encoding")
```

API entries: {obj}`costa_walk <oracq.algorithms.qlss.qlss.costa_walk>`, {obj}`schedule <oracq.algorithms.qlss.qlss.schedule>`, {obj}`CostaConfig <oracq.algorithms.qlss.qlss.CostaConfig>`, {obj}`dolph_chebyshev_plan <oracq.algorithms.qlss.qlss.dolph_chebyshev_plan>`

- `a`: {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` (BE); `bprep`: {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` (SP); the two have the same width, `fs ∈ [0,1]`.
- `CostaConfig.kappa` is the inverse spectral bound of the **encoded matrix** $A/\alpha$ (requiring $\sigma_{\min}(A/\alpha) \ge 1/\kappa$), not cond(A) at an arbitrary scale. The problem-layer entry point derives it from the spectral declarations: `replace(config, kappa=system.inverse_norm_bound)`, keeping the alpha and σ_min declarations from drifting apart.
- Problem input: {obj}`LinearSystem(block=BlockSystem(encoding, rhs, spectrum)) <oracq.algorithms.qlss.qlss.LinearSystem>`, or sparse input adapted automatically through the CKS $T^\dagger S T$ route; returns {obj}`SolveResult <oracq.algorithms.qlss.qlss.SolveResult>` (for attributes and the norm-probe contract see [CKS solver](cks.md)).

{obj}`costa_qlss <oracq.algorithms.qlss.qlss.costa_qlss>` output attributes: `algorithm = "costa_qlss"`, `input_alpha`, `encoded_inverse_bound`, `normalization_assumption = "sigma_min(A / alpha) >= 1 / kappa"`, `steps`, `filtering`, `kernel_status = "prototype; initial walk eigenstate and readout channel unverified"`, `success_condition = "signal == 0; probability and solution accuracy unverified"`.

## Implementation notes

The `costa_walk` signal layout is `enc(a.signal_qubits) | bw(RHS work) | a1 | a2 | a3 | a4`. The circuit is assembled following the paper's structure: $U_b^\dagger$ → RHS zero-state reflection (requiring both the target and the RHS work bits to be zero — $U_b$ is a unitary extension on target+work, so the projected object must be complete) → $U_b$; the $R(s)$ rotation is written as $R_y(2\arctan\frac{f}{1-f}) \cdot Z$ acting on a2; controlled $U_A$ / $U_A^\dagger$ together with the a2 zero reflection form the body of the walk; at the end, a positive reflection is applied to all signal bits together with a $\pi/2$ global phase. The input contract is checked by {obj}`operator_state_contract <oracq.algorithms.input_model.interfaces.operator_state_contract>`.

The schedule is $f(s) = \frac{\kappa}{\kappa-1}\bigl(1 - (1 + s(\kappa^{p-1} - 1))^{1/(1-p)}\bigr)$ (degenerating to $s$ when $\kappa = 1$). `costa_qlss` generates each walk at $s_i = (i+1)/\text{steps}$ and applies the Dolph–Chebyshev filter on the last walk step: {obj}`unary_weight_preparation <oracq.algorithms.qlss.qlss.unary_weight_preparation>` prepares the weights on the clock register, the walk is repeated controlled on each clock bit (negative-offset powers via the adjoint), and the Laurent polynomial comes out as a coherent superposition. The problem layer's `solve()` then uniformly performs the physical-channel selection and the independent norm probe.

Applicability boundary: kernel_status honestly records that the initial walk eigenstate and the readout channel are unverified; a conflict between the alpha and σ_min declarations is reported as `INPUT_SPECTRUM` at the check stage instead of silently generating a circuit with a wrong schedule.

## Validation approach

Category C2 (approximately continuous semantics; acceptance criteria in `../development/validation-plan.md` §2). Consistent with the `qlss.py` row of the validation coverage matrix:

- Structure: `tests/core/test_qlss_input_models.py:QLSSInputTests` (construction and contract assertions).
- Numerical: `test_protocols_consume_different_models_and_scale_kappa` — a configuration declaring `kappa=999` is replaced at the problem layer by `encoded_inverse_bound = 8` ($\alpha = 2$, $\sigma_{\min}$ lower bound 0.25), the underlying `costa_qlss` module attributes are synchronized, and the result has the same width and alpha as the CKS route with no unbound slots in the program; `test_matrix_probe_recovers_scalar_system_norm` (the norm probe recovers $\lVert x\rVert = 4$, places = 10). The related `test_qfvm_preserves_sparse_input_for_both_solvers` witnesses that both routes output same-width solution states for the QFVM sparse problem, with `encoded_inverse_bound = 72` and only the Roe oracle slot left open.
- Binding: `test_costa_rhs_reflection_is_independent_of_unitary_extension` — two unitary extensions of the same $|b\rangle$ (direct basis-state preparation vs a swap into the work bits) yield `costa_walk(identity(1), ·, 0.3)` outputs whose amplitudes agree pointwise on all target/signal basis-state inputs (places = 11), witnessing that the RHS zero-state reflection's treatment of the work bits does not depend on the particular extension.

## Known gaps and planned stages

The end-to-end cross-check against the HHL paper's reference values is missing (stage V3, a candidate for the catalog directory); verification of the initial walk eigenstate and the readout channel is also unfinished (the prototype declaration in the source attributes). Consistent with the gap column of the `qlss.py` row in the validation coverage matrix.

## Related links

- Source: `src/oracq/algorithms/qlss/qlss.py`
- Related pages: [CKS Chebyshev solver](cks.md), [Sparse Matrix Block Encoding](sparse-block-encoding.md), [VTAA-CKS variable-time solver](vtaa-cks.md)
- API reference: [Quantum linear systems](../../api/algorithms/qlss/qlss.rst)
- Input-model review: [QFVM input-model review](../../reference/qfvm-input-models.md)
- Concepts: [Conventions owned by algorithms: starting from one gate](../contracts.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Experiment design: minimal instance (1-qubit target): `a` is a diagonal block encoding ({obj}`diagonal_block_encoding <oracq.algorithms.input_model.oracles.diagonal_block_encoding>`, angle_scale $=\pi/3$), `bprep` is the basis-state preparation {obj}`basis_state(1, 0) <oracq.algorithms.input_model.oracles.basis_state>`, schedule point $f_s=0.5$, 6 signal bits in total (enc 2 + a1..a4), 7 qubits overall. Two decidable properties are validated:

1. Unitarity: the OriginIR-ext export goes through UniQC `Circuit.to_matrix` to obtain the 128-dimensional matrix $W$, and $\|W^\dagger W - I\|_{\max}$ is computed;
2. Backend consistency: execution results on the all-zero input state are cross-checked amplitude by amplitude across the four paths reference, rir-pysparq, adapter-pysparq, and originir-ext.

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| costa-walk-unitarity | 7 qubits, fs=0.5 | to_matrix + 4 paths | $\|W^\dagger W-I\|_{\max}$ | 5.6e-16 |
| ditto | ditto | ditto | max cross-path state error | 3.1e-17 |

The walk operator is unitary and consistent across backends, showing that the circuit assembly (reflections, controlled $U_A$/$U_A^\dagger$, the rotation schedule, and the final positive reflection) is self-consistent at the operator level. Note that this section does not cover the kernel's physical channels (the initial walk eigenstate and readout); that part is still labeled `prototype; ... unverified` upstream and is an open item at the algorithm-design level rather than an assembly defect.

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

Artifacts: `out/verification/search_walks.json` (case `costa-walk-unitarity`).
