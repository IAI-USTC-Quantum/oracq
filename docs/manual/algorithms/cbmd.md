# CBMD

**English** · <a href="../../../zh/manual/algorithms/cbmd.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.qode.cbmd`](../../api/algorithms/qode/cbmd.rst) · Stage V2

## Overview

Matrix-function and evolution assembly based on contour decomposition: write $e^{-At}$ of $u'=-Au$ ($A=L+iH$, $L\succeq 0$) as a finite weighted sum of Hermitian branches $H+qL$ ($q$ taken on real nodes), with the weights from a residue discretization of the contour identity. The main series and the truncation follow QST 11 035027 (2026), Eq. 12–13, as cited by the source's `ContourPlan.metadata` ([arXiv:2511.10267v3](https://arxiv.org/abs/2511.10267)).

The implementation stance is "explicitly record omissions": the non-Hermitian evolution branches contributed by the auxiliary poles and the infinite-series tail terms are currently **not generated**; both are kept explicitly in the `omitted` list of `metadata()`, and the validity premise `L>=0 and norm(integral L dt) <= 2*pi*a` is likewise written into the metadata.

## Interface and input model

```python
cbmd_qode(model, time, *, plan=None, hamiltonian_function=taylor_hamiltonian)
ContourPlan(a=1.0, cutoff=2, poles=(2j, -1 + 1j, 1j, 1 + 1j))
cbmd_function(a, nodes, residue_weights, hermitian_function)
```

API entries: {obj}`cbmd_qode <oracq.algorithms.qode.cbmd.cbmd_qode>`, {obj}`ContourPlan <oracq.algorithms.qode.cbmd.ContourPlan>`, {obj}`cbmd_function <oracq.algorithms.qode.cbmd.cbmd_function>`

- `model`: {obj}`LinearODE(HermitianParts(L, H), initial) <oracq.algorithms.qode.ode_models.LinearODE>`, the input model is an ODE (the same input surface as LCHS).
- {obj}`ContourPlan <oracq.algorithms.qode.cbmd.ContourPlan>`: the contour parameter `a > 0`, the truncation `cutoff >= 0`, and finitely many distinct non-real auxiliary poles avoiding $-i$. Derived quantities: the main-series nodes $q_k=k/a$ ($k=-J..J$), the node weights, and the auxiliary-pole coefficients.
- `hamiltonian_function` and `time` of {obj}`cbmd_qode <oracq.algorithms.qode.cbmd.cbmd_qode>` follow the same conventions as LCHS.
- {obj}`cbmd_function(a, nodes, residue_weights, hermitian_function) <oracq.algorithms.qode.cbmd.cbmd_function>`: the generic $f(A)$ assembly point; `a` is an operator BE, and the nodes/residue weights are provided by the caller per the contour identity.

`cbmd_qode` returns a {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>` (`target`/`signal`) with module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"cbmd_qode"` |
| `evolution_alpha` / `branch_count` | the finite-sum LCU normalization and the node count |
| `contour_plan` | JSON of `ContourPlan.metadata()` (including `omitted`, `assumption`, `source`) |
| `remainder` | `"auxiliary pole contribution and truncation pending"` |
| `input_assumption` | `"L>=0; autonomous; homogeneous"` |

`cbmd_function` returns a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` with attributes `algorithm="cbmd_matrix_function"`, `correctness="pending"`, and `residue_sign_convention="caller supplies target-side weights from contour identity"`.

## Implementation notes

The main-series node weights are given by the closed form

$$
w_k=\frac{e^{-2\pi a}-1}{a\cdot 2\pi i\,(q_k+i)\prod_{p}\frac{q_k-p}{-i-p}}
$$

and the auxiliary-pole coefficients analogously (with the denominator replaced by $(e^{-2\pi i p a}-1)\prod_{p'\neq p}\frac{p-p'}{-i-p'}$). The generation of `cbmd_qode` shares `_lcu_dynamics` with LCHS: per-node $K_k=H+q_kL$, the replaceable `hamiltonian_function`, LCU combination, and {obj}`apply_be_to_state <oracq.algorithms.common.state_preparation.apply_be_to_state>` applied to the initial state; the only differences are the plan type and the recorded metadata.

The branch-combination order of `cbmd_function` is {obj}`lcu([(node, parts.h), (1, parts.hermitian)]) <oracq.algorithms.input_model.block_encoding.lcu>`, i.e. $qH+L$ (different from the $H+qL$ convention of `cbmd_qode`); it keeps the Hermitian-function protocol open, does not secretly substitute matrix inversion, and does not claim to have verified the residue sign direction.

Applicability boundary: beyond LCHS's dissipation/autonomous/homogeneous declarations, the contour premise $\lVert\int L\,dt\rVert\le 2\pi a$ is also required; the omitted auxiliary-pole branches and the infinite-series tail mean the current output is only the finite part of the main series, not a complete approximation of $e^{-At}$.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2). Three layers of evidence:

- Structure: the cbmd subtest of `tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` — with abstract BE/SP inputs the open slots are kept, the `algorithm` attribute is correct, and a serialization round trip holds.
- Numerical: `tests/core/test_differential.py:DifferentialStructureTests.test_cbmd_plan_records_omitted_terms` — `ContourPlan(cutoff=2)` has 5 weights, the auxiliary coefficients match the pole count, and the `omitted` list of `metadata()` contains `auxiliary_nonhermitian_evolutions`.
- Binding: same-type registration and contract checks (sharing the `_dynamics` contract path; see the contract-chain witness on [LCHS](lchs.md)).

## Known gaps and planned stages

A unified "analytically solvable ODE family" convergence benchmark is missing (stage V2); implementing the auxiliary-pole branches and the infinite-series tail is follow-up work at the algorithm layer, currently declared explicitly through the `omitted` metadata and not counted in the gap column, consistent with the cbmd.py row of `validation-coverage.md`.

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_ode.py` (the ode group; this file covers the parts of `cbmd.py`). All classical references are independent: `ContourPlan` derived quantities recomputed directly from the Eq. 12–13 closed forms, the $t=0$ contour-identity residual, a numpy per-branch truncated-Taylor sum, and the exact solution from `scipy.linalg.expm`; the quantum programs run on the real backends reference, rir-pysparq, and OriginIR-ext.

**Experiment design**: (a) plan substructure — the nodes/weights/auxiliary-pole coefficients of `ContourPlan(a=1, cutoff=2..8)` against an independent closed-form recomputation, plus the truncation-residual trend of the $t=0$ identity main+aux→1; (b) end to end — the same non-commuting $2\times2$ problem as LCHS ($L=[[1,0.3],[0.3,0.5]]$, $H=[[0.2,0.1],[0.1,-0.1]]$), with `ContourPlan(a=1.0, cutoff=2)`, Taylor degree 3, $t=0.2$; the post-selected block is compared against the independent simulation and the exact solution; (c) cross-check on the same problem — the quantum physical blocks of the two methods are reused to compare the solution direction.

**Key metrics**:

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| cbmd-contour-plan-formula-and-identity | cutoff 2/4/8 | plan object | maximum deviation of weights/auxiliary coefficients vs the closed forms | 0.0 |
| ditto | | | $t=0$ identity residual (informational) | 3.73e-2 / 3.20e-3 / 1.59e-4 |
| cbmd-parts-noncommuting | 15 qubits | reference, originir | implementation error (vs independent simulation) | 1.1e-16 |
| ditto | | | method error (omitted terms, informational) | 9.1e-3 |
| cbmd-vs-lchs-direction | same problem | reference | solution-direction error: CBMD / LCHS | 3.0e-3 / 4.5e-2 |

The identity residual decreases stably with the truncation, consistent with `infinite_series_tail` in the `omitted` metadata; the end-to-end method error is exactly the contribution of the auxiliary-pole branches declared omitted in the documentation (about 1% in this instance). On the same problem, the solution-direction error of the CBMD main series (3.0e-3) is significantly smaller than the Cauchy quadrature remainder of LCHS (4.5e-2).

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

Artifacts: `out/verification/ode.json` (3 `cbmd-*` cases in total).

## Related links

- Source: `src/oracq/algorithms/qode/cbmd.py`
- API reference: [CBMD](../../api/algorithms/qode/cbmd.rst)
- Related pages: [QODE problem and protocol](qode-problem.md) · [LCHS](lchs.md) (the same input surface and assembly skeleton) · [Carleman linearization](carleman.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)
