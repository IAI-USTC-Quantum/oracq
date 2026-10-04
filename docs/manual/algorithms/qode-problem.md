# QODE Problem and Protocol

**English** · <a href="../../zh/manual/algorithms/qode-problem.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.qode.ode`](../../api/algorithms/qode/ode.rst) · Stage V2

## Overview

This module does not implement a concrete solving algorithm; instead, for the autonomous homogeneous linear ODE

$$
u'(t) = G\,u(t),\qquad u(0)=u_0
$$

it defines three composable building blocks: the problem object {obj}`QODEProblem <oracq.algorithms.qode.ode.QODEProblem>` (generator, initial state, dissipation declaration, and evidence note), the replaceable protocol {obj}`QODEProtocol <oracq.algorithms.qode.ode.QODEProtocol>` (contract checks and declaration pass-through), and the generic entry {obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>`, which routes $G$ to a concrete method such as LCHS, CBMD, or Schrödingerization. It also provides {obj}`make_euler_history_qode <oracq.algorithms.qode.ode.make_euler_history_qode>`, which writes implicit Euler time marching as one quantum linear-system solve.

The design stance is "declaration first": `dissipative` is the mathematical declaration $\mathrm{Hermitian}(G)\le 0$, not a language proof; the `evidence` field forces the source of the declaration to be recorded. For the full user-facing workflow (input paradigm adaptation, batched binding) see [the differential equations chapter](../differential-equations.md).

## Interface and input model

```python
QODEProblem(generator, initial, dissipative=None, initial_norm=None,
            evidence="caller_declared_unverified")
linear_qode(method, *, hamiltonian_function=taylor_hamiltonian, **options)
make_euler_history_qode(qlss, *, steps=2)
```

API entries: {obj}`QODEProblem <oracq.algorithms.qode.ode.QODEProblem>`, {obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>`, {obj}`make_euler_history_qode <oracq.algorithms.qode.ode.make_euler_history_qode>`

- `generator`: a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` of $G$ (the input model is a BE); it must have adjoint/controlled capabilities.
- `initial`: a {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` (SP) with zero input and clean work, of the same width as `generator`.
- `dissipative` / `initial_norm`: a `bool | None` dissipation declaration and the non-negative initial-value norm.
- `evidence`: a non-empty string describing the source of the declaration.
- `method`: `"lchs"` / `"cbmd"` / `"schrodingerization"`; `options` accepts only `plan`, whose type must match the chosen method.
- `qlss`: a callable linear-system solver (such as the Costa route in `qlss.py`) with signature `(matrix, rhs)`.

`linear_qode` returns a `QODEProtocol` (fields `name`, `kernel`, `requires_dissipative`, which is `True` for lchs/cbmd). Protocol methods:

| Method | Behavior |
|---|---|
| `check(generator, initial=None, time=None)` | contract report: BE/SP capabilities, equal widths, non-negative time, dissipation premise |
| `__call__(generator, initial, time)` | compatibility entry: generates after check passes, returns {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>` |
| `solve(problem, time)` | problem-level entry: validates the `QODEProblem`, then generates (recommended) |

The `StateOracle` returned by `solve` (`target`/`signal`, width equal to `generator.width`) passes the problem declarations through, on top of the method's own attributes:

| Attribute | Meaning |
|---|---|
| `qode_input_evidence` | the `evidence` of the problem object |
| `qode_dissipative_promise` | the declared `dissipative` (only when not None) |
| `qode_initial_norm` | the declared initial-value norm (only when not None) |

## Implementation notes

`QODEProtocol.contract` adds per-method assumptions on top of the shared {obj}`operator_state_contract <oracq.algorithms.input_model.interfaces.operator_state_contract>`: lchs/cbmd append "Hermitian(G)<=0; problem-level solve must declare dissipative=True", and schrodingerization appends the application-layer validation responsibility for the auxiliary window, the Fourier convention, and the recovery region. `_generate` forces the output to be a `StateOracle`, of the same width as the input, that keeps adjoint/controlled capabilities; violations raise {obj}`ContractError <oracq.algorithms.input_model.contracts.ContractError>` at generation time. The lchs/cbmd routes internally first construct {obj}`LinearODE(HermitianParts.from_operator(scale(-1, generator)), initial) <oracq.algorithms.qode.ode_models.LinearODE>`, i.e. the $L+iH$ decomposition of $A=-G$; schrodingerization consumes $G$ directly. Hermitian/positive-semidefinite properties are all borne by the caller.

The history-state layout of `make_euler_history_qode`: the time register `nt = steps.bit_length()` bits sits in the high bits, the physical register in the low bits. The system matrix is the LCU symbolic sum

$$
C = I - \Delta t\,(q\otimes G) - (S\otimes I),\qquad \Delta t = T/\text{steps},
$$

where `q = projector(nt, range(1, steps + 1))` restricts $G$ to the time slots $1..\text{steps}$ and `S = truncated_shift(nt, steps)` is the inter-slot raising shift; expanded, the equation of each slot is exactly implicit Euler $(I-\Delta t\,G)\,u_k = u_{k-1}$, and slot 0 is filled with the initial state by {obj}`extend_initial <oracq.algorithms.common.state_preparation.extend_initial>`. After `qlss(c, rhs)` solves, {obj}`select_subspace <oracq.algorithms.common.state_preparation.select_subspace>` selects the snapshot of time slot `steps`, with the selection condition folded into the signal.

Not implemented: the inhomogeneous $u'=Gu+f$, time-dependent $G(t)$, and a unified physical-norm recovery object have no ready interface; `solve` only stores the source of the declarations and does not add fabricated norm-recovery capabilities.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2). This page covers the assembly and contract layers; the numerical witnesses of each concrete method are on the method pages. Three layers of evidence:

- Structure: `tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` — after generation with abstract BE/SP inputs, all three `linear_qode` methods keep the open slots (`input_A`/`input_b`) as is, the `algorithm` attribute is correct, and the RIR round-trips through serialization.
- Numerical: `tests/core/test_sde.py:SolverContractTests.test_qode_problem_accepted_by_lchs` — a `QODEProblem` on a 4-point OU grid goes through `check().require()` and `solve`; the output width, the `qode_dissipative_promise` pass-through, and the program serialization round trip are all correct.
- Binding: same-type registration and contract checks; the assembly of `make_euler_history_qode` enters the catalog-level round trip, closure, and export checks of `tests/core/test_workloads.py:WorkloadTests.test_every_catalog_case_has_a_closed_description` through the catalog's QHAM m=1 case (`applications/catalog.py`, `steps=1` connected to the Costa QLSS).

## Known gaps and planned stages

A unified "analytically solvable ODE family" convergence benchmark is missing (a batch cross-check running the same analytically solvable problem through all solvers has not been established); it is assigned to the stage V2 convergence-sweep framework, consistent with the ode.py row of `validation-coverage.md`.

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_ode.py` (the ode group; this file covers the parts of `ode.py` and `ode_models.py`). The numerical behavior of `QODEProtocol`/`linear_qode` routing and `QODEProblem.solve` is validated by the end-to-end cases of each method (see the numerical validation sections of [LCHS](lchs.md), [CBMD](cbmd.md), [Schrödingerization](schrodingerization.md), and [Carleman linearization](carleman.md)).

**Experiment design**: (a) `linear_qode("lchs")` routing — end-to-end cross-checks of five input models (the $A=-G$ alignment of `HermitianParts.from_operator(scale(-1, G))` confirmed by an independent numpy simulation); (b) `QODEProblem.solve` — a 4-point OU zero-flux Fokker–Planck problem (`dissipative=True`, evidence pass-through) goes through `check().require()` and `solve`, with the output width, the `qode_dissipative_promise` pass-through, and the numerical solution verified together; (c) declaration pass-through — the `qode_input_evidence`/`qode_dissipative_promise` attributes appear in the output module's attributes as declared.

**Key metrics**:

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| lchs-given-be-scalar-decay | 7 qubits | four paths | implementation error (routing + assembly vs independent simulation) | 1.2e-16 |
| lchs-fokker-planck-ou | 16 qubits | reference, originir | implementation error / cross-check between classical references / `qode_dissipative_promise` | 1.1e-16 / 1.1e-16 / True pass-through |
| cbmd-parts-noncommuting | 15 qubits | reference, originir | implementation error (CBMD routing) | 1.1e-16 |
| schrodingerization-* | 21 qubits | reference, originir (decay-grid) / reference (the rest) | assembly fidelity (Schrödingerization routing) | 8.7e-19 |

The implicit Euler history assembly of `make_euler_history_qode` (the LCU symbolic sum $C=I-\Delta t(q\otimes G)-(S\otimes I)$) is covered by the catalog's QHAM m=1 case structure test; this group does not run a numerical end-to-end validation of it (that requires a complete QLSS kernel, which is out of scope for this group).

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

Artifacts: `out/verification/ode.json`.

## Related links

- Source: `src/oracq/algorithms/qode/ode.py`
- Tutorial: [Swapping QODE methods for the same linear problem](../../tutorials/differential-equations.md)
- API reference: [QODE assembly interface](../../api/algorithms/qode/ode.rst)
- Pages in the same group: [LCHS](lchs.md) · [CBMD](cbmd.md) · [Schrödingerization](schrodingerization.md) · [Carleman linearization](carleman.md) · [Fokker–Planck input model](fokker-planck.md)
- Concepts: [Conventions owned by algorithms: starting from one gate](../contracts.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)
