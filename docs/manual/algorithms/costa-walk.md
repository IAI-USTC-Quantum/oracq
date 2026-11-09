# Costa Walk QLSS

**English** · <a href="../../zh/manual/algorithms/costa-walk.html">简体中文</a>

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

{obj}`costa_qlss <oracq.algorithms.qlss.qlss.costa_qlss>` output attributes: `algorithm = "costa_qlss"`, `input_alpha`, `encoded_inverse_bound`, `normalization_assumption = "sigma_min(A / alpha) >= 1 / kappa"`, `steps`, `filtering`, `kernel_status = "prototype; finite-step accuracy and success probability require validation"`, `success_condition = "signal == 0; probability and solution accuracy unverified"`.

## Implementation notes

The `costa_walk` signal layout is `enc(a.signal_qubits) | bw(RHS work) | a1 | a2 | a3 | a4`. The circuit is assembled following the paper's structure: $U_b^\dagger$ → RHS zero-state reflection (requiring both the target and the RHS work bits to be zero — $U_b$ is a unitary extension on target+work, so the projected object must be complete) → $U_b$; the $R(s)$ rotation is written as $R_y(2\arctan\frac{f}{1-f}) \cdot Z$ acting on a2; controlled $U_A$ / $U_A^\dagger$ together with the a2 zero reflection form the body of the walk; at the end, a positive reflection is applied to `enc | a2 | a3` together with a $\pi/2$ global phase. The input contract is checked by {obj}`operator_state_contract <oracq.algorithms.input_model.interfaces.operator_state_contract>`.

The schedule is $f(s) = \frac{\kappa}{\kappa-1}\bigl(1 - (1 + s(\kappa^{p-1} - 1))^{1/(1-p)}\bigr)$ (degenerating to $s$ when $\kappa = 1$). `costa_qlss` generates each walk at $s_i = (i+1)/\text{steps}$ and applies the Dolph–Chebyshev filter on the last walk step: {obj}`unary_weight_preparation <oracq.algorithms.qlss.qlss.unary_weight_preparation>` prepares the weights on the clock register, the walk is repeated controlled on each clock bit (negative-offset powers via the adjoint), and the Laurent polynomial comes out as a coherent superposition. The problem layer's `solve()` then uniformly performs the physical-channel selection and the independent norm probe.

Applicability boundary: `kernel_status` records that finite-step accuracy and success probability still require validation; the default two steps are a construction example and do not certify a solve tolerance; a conflict between the alpha and σ_min declarations is reported as `INPUT_SPECTRUM` at the check stage instead of silently generating a circuit with a wrong schedule.

## Validation approach

Category C2 (approximately continuous semantics; acceptance criteria in `../development/validation-plan.md` §2). Consistent with the `qlss.py` row of the validation coverage matrix:

- Structure: `tests/core/test_qlss_input_models.py:QLSSInputTests` (construction and contract assertions).
- Numerical: `test_protocols_consume_different_models_and_scale_kappa` — a configuration declaring `kappa=999` is replaced at the problem layer by `encoded_inverse_bound = 8` ($\alpha = 2$, $\sigma_{\min}$ lower bound 0.25), the underlying `costa_qlss` module attributes are synchronized, and the result has the same width and alpha as the CKS route with no unbound slots in the program; `test_matrix_probe_recovers_scalar_system_norm` (the CKS norm probe recovers $\lVert x\rVert = 4$, places = 10). The related `test_qfvm_preserves_sparse_input_for_both_solvers` witnesses that both routes output same-width solution states for the QFVM sparse problem, with `encoded_inverse_bound = 72` and only the Roe oracle slot left open.
- Binding: `test_costa_rhs_reflection_is_independent_of_unitary_extension` — two unitary extensions of the same $|b\rangle$ (direct basis-state preparation vs a swap into the work bits) yield `costa_walk(identity(1), ·, 0.3)` outputs whose amplitudes agree pointwise on all target/signal basis-state inputs (places = 11), witnessing that the RHS zero-state reflection's treatment of the work bits does not depend on the particular extension.

## Known gaps and planned stages

A general finite-step error bound and automatic selection of steps/filter parameters for a requested accuracy are still missing. The two-dimensional witnesses below validate these specific systems, not arbitrary spectra, application success costs, tomography or hardware performance. Increasing the step count is not itself a certificate of monotonic convergence.

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

The walk operator is unitary and consistent across backends, showing that the circuit assembly (reflections, controlled $U_A$/$U_A^\dagger$, the rotation schedule, and the final positive reflection) is self-consistent at the operator level. This historical single-walk check does not test the solve channel. The separate regression below covers initialization, readout and independently computed inverse solutions.

Reproduction command:

```bash
PYTHONPATH=src python tests/verification/verify_search_walks.py
```

Artifacts: `out/verification/search_walks.json` (case `costa-walk-unitarity`).


## Initial eigenstate and physical readout (0.3.1)

For this implementation, $W(0)=i(2\Pi-I)U_H(0)$ with
$\Pi$ projecting `enc | a2 | a3` onto zero. Let $v$ be the clean-signal RHS
state. The involution gives $w=U_H(0)v$ with `a4=a3=1`, `a2=+`, `a1=0`.
Consequently $W(0)v=-iw$, $W(0)w=iv$, and $(v-iw)/\sqrt{2}$ is the
$+1$ walk eigenstate. The kernel prepares it with `H(a4)`, exact `S†(a4)`
(`Rz(-pi/2)` and global phase `-pi/4`), then `H(a2)` and `X(a3)` controlled
by `a4`. The explicit global phase is retained under coherent controls.

At the final interpolation point, the inverse solution occupies `a1=1` in
our $A_f$ convention. After filtering, `X(a1)` maps that physical branch to
the public all-zero-signal success condition. No ancilla is dropped from
postselection. This corrects the 0.3.0 kernel, which began with $v$ alone
and selected the opposite final branch.

`tests/core/test_costa_channels.py` checks Hermitian and complex nonsymmetric
matrices, even/odd schedules, a clean-work unitary extension, coherent
controlled adjoints, and the Costa problem-layer independent norm probe.
The same success-channel witnesses fail against 0.3.0; they do not use a
mock backend or a success-label-only change.

`tests/integration/test_costa_channels_native.py` runs the actual PySparQ gates
on $A=L\operatorname{diag}(0.4,0.9)R$, $L=R_y(0.61)$, with
$R=R_y(-0.61)$ or $R=R_z(0.83)$ and
$b=(\sqrt{0.3},e^{0.4i}\sqrt{0.7})$. At 512 steps, $\kappa=2.5$ and the
default filter, both have success probability about 0.499974 and conditional
phase-aligned $L_2$ error about $4.24\times10^{-4}$ against an independent
NumPy solve (required tolerance $10^{-3}$). The complex case also rejects
the $A^\dagger$ inverse solution. All gates execute without native arithmetic
shortcuts. These are numerical correctness witnesses, not a speedup claim.

```bash
PYTHONPATH=src python -m unittest discover -s tests/core -p test_costa_channels.py -v
# Requires the real optional PySparQ backend and NumPy:
PYTHONPATH=src python -m unittest discover -s tests/integration -p test_costa_channels_native.py -v
```

Derivation reference: [Costa et al., arXiv:2111.08152](https://arxiv.org/abs/2111.08152),
section IV and Appendix E; the equations above use this repository's ancilla convention.
