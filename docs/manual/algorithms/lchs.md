# Linear Combination of Hamiltonian Simulations

**English** · <a href="../../zh/manual/algorithms/lchs.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.qode.lchs`](../../api/algorithms/qode/lchs.rst) · Stage V2

## Overview

For the autonomous homogeneous equation $u'=-Au$ with $A=L+iH$ and $L\succeq 0$ ($[L,H]=0$ not required), write the non-unitary evolution as a weighted integral of Hermitian evolutions (Theorem 1 of the original LCHS paper, [arXiv:2303.01029](https://arxiv.org/html/2303.01029v2); see also [the differential equations chapter §4](../differential-equations.md)):

$$
e^{-At}=\int_{\mathbb R}\frac{e^{-it(H+kL)}}{\pi(1+k^2)}\,dk .
$$

This implementation discretizes the integral with finitely many nodes $\{k_j\}$ and complex weights $\{w_j\}$ that already include the integration kernel, obtaining the finite sum $\widetilde V=\sum_j w_j\,e^{-it(H+k_jL)}$, and then combines the Hamiltonian simulations of the branches with an LCU. The alignment with the generic interface $u'=Gu$ is $A=-G$ (done by the routing of {obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>`; see [QODE problem and protocol](qode-problem.md)).

## Interface and input model

```python
lchs_qode(model, time, *, plan=None, hamiltonian_function=taylor_hamiltonian)
QuadraturePlan(nodes, weights, kernel="user_supplied")
QuadraturePlan.cauchy(cutoff=2, spacing=1.0)
```

API entries: {obj}`lchs_qode <oracq.algorithms.qode.lchs.lchs_qode>`, {obj}`QuadraturePlan <oracq.algorithms.qode.lchs.QuadraturePlan>`

- `model`: {obj}`LinearODE(HermitianParts(L, H), initial) <oracq.algorithms.qode.ode_models.LinearODE>`, the input model is an ODE; `parts.hermitian` stores $L$ and `parts.h` stores $H$, and the initial state is an SP.
- `plan`: the nodes are finite real numbers, the weights are finite complex numbers that are not all zero, and the two are of equal length; the `cauchy` factory generates nodes $k=-J\cdot h,\dots,J\cdot h$ and weights $w=h/[\pi(1+k^2)]$, with the kernel labeled `"finite_cauchy"`.
- `hamiltonian_function`: a replaceable protocol `(K, t) -> BlockEncoding`, by default {obj}`taylor_hamiltonian <oracq.algorithms.common.hamiltonian.taylor_hamiltonian>` (a non-unitary BE of the truncated Taylor series).
- `time`: a non-negative real number.

Returns a {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>` (`target`/`signal`) with module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"lchs_qode"` |
| `evolution_alpha` | the normalization constant of the finite-sum LCU, $\sum_j|w_j|\,\alpha_{E_j}$ |
| `branch_count` | the number of quadrature nodes |
| `quadrature_kernel` / `quadrature_nodes` | the plan's kernel name and the JSON node table |
| `input_assumption` | `"L>=0; autonomous; homogeneous"` |
| `remainder` | `"finite quadrature pending"` |

## Implementation notes

For each node construct $K_j = H + k_j L$ ({obj}`lcu([(1, parts.h), (k_j, parts.hermitian)]) <oracq.algorithms.input_model.block_encoding.lcu>`, where `lcu` drops zero-coefficient terms automatically and $\alpha_{K_j}=\alpha_H+|k_j|\,\alpha_L$), hand it to `hamiltonian_function` to obtain the BE of the approximate evolution $E_j$ (which must keep the alpha), then combine $\{(w_j, E_j)\}$ with `lcu` into $\widetilde V$ and apply it to the initial state via {obj}`apply_be_to_state <oracq.algorithms.common.state_preparation.apply_be_to_state>`. The register layout is inherited from the input oracles: the target width is the width of $L/H$, and the signal is the concatenation of the LCU select bit and each branch's signal bit. The success subspace is all signals zero; what is read out is the normalized direction $\widetilde V|u_0\rangle/\alpha_V$. The physical amplitude still contains the initial-value norm, which the protocol does not recover automatically.

Applicability boundary: $L\succeq 0$, autonomous, and homogeneous are all input declarations; `remainder` states explicitly that the finite sum carries no tail-integration guarantee — do not forcibly normalize the weights and still call it the same approximation operator. The time-dependent and inhomogeneous results of the original paper are not in the current interface. Before solving, {obj}`operator_state_contract <oracq.algorithms.input_model.interfaces.operator_state_contract>` checks the capability contracts of the two BEs $L$, $H$ and of the initial state.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2): the acting operator must approximate the target continuous function within tolerance. `tests/core/test_differential.py` is positioned as a structure and open-input difference test (see its docstring); numerical cross-checks are currently borne by the input-model side. Three layers of evidence:

- Structure: the lchs subtest of `tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` — with abstract BE/SP inputs the open slots (`input_A`/`input_b`) are kept, the `algorithm` attribute is correct, and the RIR round-trips through serialization.
- Numerical: `tests/core/test_sde.py:SolverContractTests.test_qode_problem_accepted_by_lchs` — a {obj}`QODEProblem <oracq.algorithms.qode.ode.QODEProblem>` on a 4-point OU grid (`QuadraturePlan.cauchy(cutoff=0)`, Taylor degree 1) goes through the full `check`/`solve` chain of `linear_qode("lchs")`, with the output width and the dissipation-declaration pass-through correct and a serialization round trip.
- Binding: same-type registration and contract checks; `FokkerPlanckProblem.qode_problem` is precisely a consumer of this protocol (see [Fokker–Planck input model](fokker-planck.md)).

## Known gaps and planned stages

A unified "analytically solvable ODE family" convergence benchmark is missing: a batch cross-check running the same analytically solvable problem through all ODE solvers has not been established; it is assigned to the stage V2 convergence-sweep framework, consistent with the lchs.py row of `validation-coverage.md`.

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_ode.py` (the ode group; this file covers the parts of `lchs.py` and of the `_dynamics.py` assembly skeleton). All classical references are independent: a numpy per-branch truncated-Taylor sum for $K_j=H+k_jL$, the exact solution from `scipy.linalg.expm`, analytic solutions ($e^{-t}$, heat-equation Fourier eigenvalues), and `sde.matrix_exponential`, a pure-Python matrix exponential, cross-checking one another; the quantum programs run on the real backends reference, rir-pysparq, adapter-pysparq, and OriginIR-ext. Implementation error (quantum post-selected block vs independent simulation) and method error (finite quadrature + Taylor remainder vs the exact solution, marked pending in the library) are reported separately.

**Experiment design**: end-to-end cross-checks of five input-model classes for $u'=-Au$, $A=L+iH$ (all plans `QuadraturePlan.cauchy(cutoff=2, spacing=1.0)`): (1) a BE of the whole $G=-I$ (analytic solution); (2) directly given Hermitian parts (a $2\times2$ matrix with $[L,H]\neq0$, Pauli-expansion BE); (3) the diagonal spectral angle database $A=\mathrm{diag}(1,\cos\pi/4)$, the same open RIR bound to a gate table and to a QRAM table; (4) the Fokker–Planck OU 4-point zero-flux discretized generator (Pauli expansion) through `QODEProblem.solve` (Boltzmann root-amplitude initial state); (5) a structured shifted-difference BE of the periodic 4-point heat equation. There is also a quadrature convergence scan (scalar $L=I$, cutoff 2/8 quantum + a 2..32 classical-kernel trend line).

**Key metrics**:

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| lchs-given-be-scalar-decay | 7 qubits | four paths | implementation error / method error | 1.2e-16 / 1.2e-2 |
| lchs-given-parts-noncommuting | 15 qubits | three paths | implementation error / method error | 1.1e-16 / 7.3e-2 |
| lchs-diagonal-spectral-gate / -qram | 2× programs | reference | implementation error (both) / gate-QRAM amplitude difference | 5.8e-17 / 0.0 |
| lchs-fokker-planck-ou | 16 qubits | reference, originir | implementation error / method error / cross-check between classical references | 1.1e-16 / 7.9e-2 / 1.1e-16 |
| lchs-heat-structured-stencil | 15 qubits | reference, originir | implementation error / method error / Fourier-expm cross-check | 1.7e-17 / 0.19 / 1.4e-17 |
| lchs-quadrature-scan | cutoff 2 / 8 | reference | implementation error | 5.6e-17 / 1.1e-16 |
| lchs-quadrature-kernel-trend | cutoff 2→32 | classical kernel quadrature | kernel quadrature error | 7.1e-2→3.7e-3 (including the oscillating tail) |

The success probability (all-signal-zero subspace) matches the classical value (deviation < 1e-16). The method error is the Cauchy quadrature remainder: the tails of the truncation $K_{\max}$ and of the step size $h$ show, respectively, an oscillating $O(1/K_{\max})$ decay and an aliasing floor $2e^{-(2\pi/h-\lambda t)}$; on deeply nested LCU programs the two pysparq implementations have a numerical floor of about 1e-7 on the junk branches (the physical blocks are unaffected; see the group report).

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

Artifacts: `out/verification/ode.json` (12 `lchs-*` cases in total).

## Related links

- Source: `src/oracq/algorithms/qode/lchs.py`
- Tutorial: [Swapping QODE methods for the same linear problem](../../tutorials/differential-equations.md)
- API reference: [LCHS](../../api/algorithms/qode/lchs.rst)
- Related pages: [QODE problem and protocol](qode-problem.md) · [Carleman linearization](carleman.md) (lifting, then LCHS) · [Fokker–Planck input model](fokker-planck.md) · [CBMD](cbmd.md), [Schrödingerization](schrodingerization.md) (other linear routes)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)
