# CKS Chebyshev QLSS

**English** · <a href="../../../zh/manual/algorithms/cks.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.qlss.qlss`](../../api/algorithms/qlss/qlss.rst) · Stage V3

## Overview

The baseline Chebyshev/LCU solution route for sparse Hermitian linear systems $Ax = b$ (the basic construction of Childs–Kothari–Somma, arXiv:1511.02306 §4). The variable-time VTAA layer of section 5 is covered in [VTAA-CKS solver](vtaa-cks.md). The truncated Chebyshev expansion of $1/x$ is applied to the spectrum of the encoded matrix $A/\alpha$:

$$
A^{-1} \approx \sum_{j=0}^{d-1} c_j\, T_{2j+1}(A/\alpha), \qquad
c_j = 4\,(-1)^j\, 2^{-2d} \sum_{i=j+1}^{d} \binom{2d}{d+i},
$$

The BE is obtained via [Sparse Matrix Block Encoding](sparse-block-encoding.md), {obj}`chebyshev_block <oracq.algorithms.input_model.sparse.chebyshev_block>` generates the odd Chebyshev powers, the LCU combines the inverse-operator BEs with weights $c_j$, and the result is finally applied to the RHS preparation. Matrix properties (spectral bounds, Hermiticity) are declared by the caller; the language does not prove them.

## Interface and input model

```python
make_cks_qlss(config=None)          # returns QLSSProtocol (input_model="sparse")
cks_chebyshev(system, config=None)  # low-level kernel, takes a SparseSystem
CKSConfig(order=2, terms=None)      # order 1..128; terms is the truncated term count, ≤ order
```

API entries: {obj}`make_cks_qlss <oracq.algorithms.qlss.qlss.make_cks_qlss>`, {obj}`cks_chebyshev <oracq.algorithms.qlss.qlss.cks_chebyshev>`, {obj}`CKSConfig <oracq.algorithms.qlss.qlss.CKSConfig>`

The problem input is {obj}`LinearSystem(sparse=...) <oracq.algorithms.qlss.qlss.LinearSystem>`, and only one source input model may be specified:

- {obj}`SparseSystem(access, value_format, entry_bound, rhs, spectrum, diagonal_nonnegative, hermitian) <oracq.algorithms.qlss.qlss.SparseSystem>` — the input model is SO (location + entry oracles) plus SP (RHS state preparation); `hermitian` and `diagonal_nonnegative` must be explicitly declared `True`, and widths and value formats must match.
- {obj}`SpectralPromise(norm_upper, sigma_min_lower, evidence) <oracq.algorithms.qlss.qlss.SpectralPromise>` — spectral bounds declared by the caller (default `evidence = "caller_declared_unverified"`), requiring $0 < \sigma_{\min}^{\text{lower}} \le$ the norm upper bound; `inverse_bound(alpha) = max(1, alpha / sigma_min_lower)`, rejected on conflict with alpha.

`protocol(problem)` returns {obj}`SolveResult <oracq.algorithms.qlss.qlss.SolveResult>`:

| Attribute | Meaning |
|---|---|
| `state` | solution state on the physical subspace ({obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>`, with `.operation` as its RIR) |
| `norm_probe` | independent matrix norm probe (`StateOracle`) |
| `input_alpha` / `encoded_inverse_bound` | the consumed BE normalization and $\max(1, \alpha/\sigma_{\min})$ |
| `rhs_norm` / `adapter_trace` | classical RHS norm and the adaptation trace |
| `recover_norm(p_solver, p_joint)` | $\lVert x\rVert = \lVert r\rVert / (\alpha\sqrt{p_{joint}/p_{solver}})$ |

Kernel output attributes: `algorithm = "cks_chebyshev_basic"`, `input_model = "sparse_location_inplace_and_entry_xor"`, `polynomial_order` / `polynomial_terms`, `inverse_lcu_normalization`, `implementation_scope = "CKS section 4 basic LCU; no VTAA"`, `correctness = "pending"`.

## Implementation notes

Generation chain: {obj}`real_symmetric_sparse_encoding <oracq.algorithms.input_model.sparse.real_symmetric_sparse_encoding>` (CKS $T^\dagger S T$, $\alpha = s \cdot a_{\max}$) → `chebyshev_block(a, 2j+1)` odd powers → {obj}`lcu <oracq.algorithms.input_model.block_encoding.lcu>` combining the inverse-operator BEs → {obj}`apply_be_to_state <oracq.algorithms.common.state_preparation.apply_be_to_state>` applied to the RHS. The protocol layer's `solve()` uniformly performs: contract checking (`check()` does not run the kernel; it reports under `INPUT_TYPE` / `INPUT_ZERO_RHS` / `INPUT_ADAPTER` / `INPUT_PROMISE` / `INPUT_SPECTRUM`), physical-channel selection via {obj}`select_subspace <oracq.algorithms.common.state_preparation.select_subspace>`, and assembly of the independent norm probe — the success rate of Costa filtering cannot reuse the normalization factor of the CKS inverse-operator LCU, because the two success branches form in different ways; the probe attribute is `probability_contract = "joint solver success and fresh matrix-ancilla success"`. The output requires adjoint and controlled capabilities (for QFVM composition).

Applicability boundary: only the paper's §4 basic LCU, no VTAA; solution accuracy and the success channel are labeled prototype/pending and constitute no promise of solution accuracy. BE input cannot be automatically converted back to sparse oracles (the reverse adaptation does not exist); this protocol also has no two-parameter legacy call shape and only accepts `LinearSystem`.

## Validation approach

Category C2 (approximately continuous semantics; acceptance criteria in `../development/validation-plan.md` §2). Consistent with the `qlss.py` row of the validation coverage matrix:

- Structure: `tests/core/test_qlss_input_models.py:QLSSInputTests` (construction and contract assertions).
- Numerical: `test_signed_sparse_encoding_and_chebyshev` (the BE angle block of a 2×2 with negative off-diagonal entries cross-checked against alpha = 2, and `chebyshev_block(be, 3)` cross-checked against $4h^3 - 3h$, places = 11); `test_matrix_probe_recovers_scalar_system_norm` (a scalar system recovers $\lVert x\rVert = 4$ through the probe, places = 10); `test_norm_recovery_uses_conditional_matrix_probe` (`recover_norm(0.5, 0.125) = 2.0`, rejected when the probability order is invalid); `test_protocols_consume_different_models_and_scale_kappa` (the two protocols output same-width solution states with identical input_alpha for the same problem, with no unbound slots in the program).
- Binding: `test_costa_rhs_reflection_is_independent_of_unitary_extension` (independence of the U extension, registered on this row). Rejection paths are pinned by `test_no_implicit_be_to_sparse` (BE input must not go through the sparse protocol) and `test_zero_rhs_is_not_a_state_preparation_problem` (zero RHS handled on the classical side).

## Known gaps and planned stages

The end-to-end cross-check against the HHL paper's reference values is missing (stage V3, a candidate for the catalog directory): the current witnesses cover the local semantics of input adaptation, kappa conversion, and norm recovery, while solution accuracy itself remains a prototype declaration. The VTAA variable-time layer is implemented separately in [VTAA-CKS solver](vtaa-cks.md); this page keeps its position as the §4 baseline route. Consistent with the gap column of the `qlss.py` row in the validation coverage matrix.

## Related links

- Source: `src/oracq/algorithms/qlss/qlss.py`
- Related pages: [Costa walk solver](costa-walk.md), [VTAA-CKS variable-time solver](vtaa-cks.md), [Sparse Matrix Block Encoding](sparse-block-encoding.md)
- API reference: [Quantum linear systems](../../api/algorithms/qlss/qlss.rst)
- Input-model review: [QFVM input-model review](../../reference/qfvm-input-models.md)
- Concepts: [Conventions owned by algorithms: starting from one gate](../contracts.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-level numerical experiments are in `tests/verification/verify_nt_qlss_sde.py` (the nt_qlss_sde group), all executed on real backends. All classical oracles are independent: the Chebyshev matrix polynomial $P(M)=\sum_j c_j T_{2j+1}(M)$ is evaluated directly with numpy via the $T_2$ recurrence (coefficients recomputed from the `math.comb` closed form), the true solution comes from `numpy.linalg.solve`, and neither goes through the tested implementation's helper functions.

**Experiment design**: (a) kernel convergence scan — a $2\times2$ signed sparse system with $\kappa=3$ ($A=[[0.75,-0.25],[-0.25,0.75]]$, $\alpha=1.5$); {obj}`cks_chebyshev <oracq.algorithms.qlss.qlss.cks_chebyshev>` runs at orders 2/4/8/16 on reference and rir-pysparq, with the conditional solution state (signal = 0 branch normalized) and the success probability compared against the polynomial oracle (implementation error), and the method error compared against the numpy true solution; (b) protocol level — {obj}`make_cks_qlss(CKSConfig(order=8)) <oracq.algorithms.qlss.qlss.make_cks_qlss>` solves a `LinearSystem`, with $p_{\text{solver}}$ and the independent matrix norm probe's $p_{\text{joint}}$ compared against the polynomial oracle and `recover_norm` compared against $\lVert A^{-1}b\rVert$; (c) Costa components — the Dolph–Chebyshev weights against the closed-form window $\gamma\,T_d(\beta\cos\theta)$ (numpy chebval sampled at 4097 points), the {obj}`schedule <oracq.algorithms.qlss.qlss.schedule>` closed form and monotonicity, and the {obj}`unary_weight_preparation <oracq.algorithms.qlss.qlss.unary_weight_preparation>` probability distribution; the assembled {obj}`costa_qlss <oracq.algorithms.qlss.qlss.costa_qlss>` program is cross-checked across three backends (its `kernel_status` is the prototype declared in the library; solution accuracy is not a criterion).

**Key metrics**:

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| cks-kernel-order-2/4/8 | $2\times2$, $\kappa=3$ | 2 paths | implementation error / success probability error | ≤ 4.9e-16 / ≤ 1.7e-16 |
| cks-kernel-order-16 | same as above | 2 paths | implementation error / success probability error / cross-backend | 7.2e-12 / 8.7e-13 / 2.5e-7 |
| cks-method-convergence | order 2→16 | numpy oracle | method error | 0.5537 → 0.4084 → 0.2130 → 0.0662 (strictly decreasing, about $e^{-4/3}$ per order) |
| cks-protocol-norm-recovery | order 8 | 2 paths | probe probability error / recover_norm relative error / fidelity | 1.1e-16 / 0.1456 / 0.9531 |
| costa-plan-schedule-unary | degree ≤ 6 | numpy + 2 paths | DC weights / schedule / unary preparation | 2.7e-15 / 2.2e-16 / 1.7e-16 |
| costa-assembly-cross-backend | steps=1 | 3 paths | per-amplitude deviation | 0.0 |

recover_norm = 2.7018 against $\lVert A^{-1}b\rVert = 3.1623$ ($b=(2,0)$): the relative error 0.1456 is consistent with the order=8 method error 0.213 and converges geometrically with the order (kernel scan). At order=16, reference stays consistent with the polynomial oracle to ~1e-16, while rir-pysparq shows a 2.5e-7 floating-point drift on a single amplitude after about $10^6$ expansion steps (physical observables agree to ≤ 8.7e-13), recorded in the final report.

**Coverage note for structural-layer modules**: `contracts.py` and `interfaces.py` are the contract/protocol structural layer and carry no independent numerical semantics of their own — their code paths (`ProtocolContract.check`, the protocol adapters {obj}`as_sparse_access <oracq.algorithms.input_model.interfaces.as_sparse_access>` / {obj}`as_block_encoding <oracq.algorithms.input_model.interfaces.as_block_encoding>`, etc.) are genuinely executed in the protocol-level cases, and their correctness is jointly covered by the `tests/core` structural tests and this group's cross-backend cross-checks; the same applies to `algorithms/legacy.py` (LCHS/Schrödingerization factories) and `applications/legacy.py` (QFVM/QHAM assembly), whose numerical content belongs to the underlying algorithm pages and whose assembly determinism is covered by the cross-backend cross-checks of the catalog cases (`lchs`, `schrodingerisation`, `qfvm_*`, `qham_*`).

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

Artifacts: `out/verification/nt_qlss_sde.json` (7 cases in total: `cks-*`, `costa-*`).
