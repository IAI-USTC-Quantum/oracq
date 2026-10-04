# Carleman Linearization

**English** · <a href="../../zh/manual/algorithms/carleman.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.qnlss.carleman`](../../api/algorithms/qnlss/carleman.rst) · Stage V2

## Overview

Lifts a polynomial-nonlinear ODE into a finite-dimensional linear system. Write

$$
u'=\sum_{p=0}^{D}F_p\,u^{\otimes p},\qquad u^{\otimes 0}=1,
$$

where $F_0$ is a constant forcing vector, $F_1$ a linear operator, and $F_p$ ($p\ge 2$) a contraction that is linear in the whole tensor input. Letting $y_k=u^{\otimes k}$, the product rule gives

$$
y_k'=\sum_{p=0}^{D}\sum_{j=0}^{k-1}\bigl(I^{\otimes j}\otimes F_p\otimes I^{\otimes(k-j-1)}\bigr)\,y_{k+p-1}.
$$

The truncation order $K$ keeps $y_0,\dots,y_K$ and drops terms sourced from orders greater than $K$, yielding the linear generator $G_K$; the finite matrix construction is deterministic, and how well it approximates the original nonlinear equation depends on the truncation and the applicability conditions (for the quantum convergence/efficiency basis see [the Carleman paper arXiv:2011.03185v4](https://arxiv.org/abs/2011.03185v4); see also [the differential equations chapter §6](../differential-equations.md)).

## Interface and input model

```python
PolynomialODE(width, coefficients, initial, initial_norm=1.0)
carleman_lift(problem, *, cutoff=2)
carleman_initial(problem, *, cutoff=2)
carleman_qode(problem, time, linear_solver, *, cutoff=2)
```

API entries: {obj}`PolynomialODE <oracq.algorithms.qnlss.carleman.PolynomialODE>`, {obj}`carleman_lift <oracq.algorithms.qnlss.carleman.carleman_lift>`, {obj}`carleman_initial <oracq.algorithms.qnlss.carleman.carleman_initial>`, {obj}`carleman_qode <oracq.algorithms.qnlss.carleman.carleman_qode>`

- `width`: the physical register width $n$ (dimension $d=2^n$, capped at 64).
- `coefficients`: `((p, the BE of F_p), ...)`, with non-negative, non-repeating degrees; the BE width of $F_p$ must be $\max(1,p)\cdot n$ (the zero-signal block contract padded to a square: the first $d$ rows are valid, the remaining rows are declared zero).
- `initial` / `initial_norm`: the {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` (SP) of the physical initial state and the raw vector norm; the norm determines the relative amplitudes of the orders in the lifted initial state — it is not a decorative attribute.
- `linear_solver`: the three-argument protocol `(generator, initial, time) -> StateOracle`, usually injected as {obj}`linear_qode(...) <oracq.algorithms.qode.ode.linear_qode>`.

The three entries return, respectively: a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` of the lifted generator $G_K$ (attributes `algorithm="carleman_lift"`, `cutoff`, `coefficient_assumption`), a `StatePreparation` of the lifted initial state, and the final {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>` (attributes `algorithm="carleman_qode"`, `cutoff`, `truncation_assumption="Carleman tail pending"`).

## Implementation notes

The register layout is the padded layout: `target = data region (K·n bits, low) + level register (K.bit_length() bits, high)`; the valid data of order $k$ occupies only the first $k\cdot n$ bits, the rest being zero — not the compact $1+d+\dots+d^K$ addressing; the resource difference is kept in reports.

For each $(k,p,\text{position})$ ($1\le k\le K$, source order $k+p-1\le K$, position $0..k-1$), {obj}`carleman_lift <oracq.algorithms.qnlss.carleman.carleman_lift>` generates one placement BE: the level register routes from the source order to the output order (XOR flips the differing bits), $F_p$ acts on the $p$ consecutive $n$-bit groups starting at position, and when $p>1$ the output stays in the lowest $n$ bits of the group with the zero rows moved to the high end of the data region; summing all terms with an LCU gives $G_K$.

{obj}`carleman_initial <oracq.algorithms.qnlss.carleman.carleman_initial>` prepares the normalized lifted vector

$$
\frac{(1,\,u_0,\,u_0^{\otimes 2},\,\dots,\,u_0^{\otimes K})}{\sqrt{\sum_{k=0}^{K}\lVert u_0\rVert^{2k}}},
$$

calling the initial-state oracle multiple times, controlled by level (work reserves $K\times$ the initial state's work width). {obj}`carleman_qode <oracq.algorithms.qnlss.carleman.carleman_qode>` hands $G_K$ and the lifted initial state to the injected linear solver, then uses {obj}`select_subspace <oracq.algorithms.common.state_preparation.select_subspace>` to select the level 1 channel (the selection value `1 << ((cutoff-1)*width)` lands exactly on the lowest bit of the level register, while also requiring the high data bits to be zero), with the selection condition folded into the signal.

Applicability boundary: $G_K$ is not guaranteed to be dissipative, even if the original PDE or $F_1$ is. Connecting Schrödingerization can pass $G_K$ directly (judge the recovery region separately); connecting LCHS must confirm $\mathrm{Hermitian}(-G_K)\succeq 0$, or apply a global shift $\mu\ge\alpha_{G_K}$ acting on the whole lifted system (including level 0), and the shift cost should be visible in the application configuration. The truncation error of the tensor tail terms is not quantified (`truncation_assumption` is marked pending).

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2). Three layers of evidence:

- Structure: the Carleman branch of `tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` — with {obj}`PolynomialODE(1, ((1, a), (2, f2)), initial) <oracq.algorithms.qnlss.carleman.PolynomialODE>` taking abstract BEs as coefficients and CBMD as the linear solver, the open slots `input_A`/`input_b`/`input_F2` are all kept, a placement module with `carleman_row_level == 1` and `carleman_column_level == 2` exists, and a serialization round trip holds.
- Numerical: `tests/core/test_differential.py` is positioned as a structure test; numerical cross-checks of truncation convergence belong to the unified ODE convergence benchmark (see known gaps).
- Binding: same-type registration and contract checks — the coefficient BE widths, degree uniqueness, and initial-state layout are validated at generation time in `PolynomialODE.__post_init__`, with violations raising {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`.

## Known gaps and planned stages

A unified "analytically solvable ODE family" convergence benchmark is missing (the same problem through all solvers, including a Carleman truncation-order sweep); it is assigned to the stage V2 convergence-sweep framework, consistent with the carleman.py row of `validation-coverage.md`.

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_ode.py` (the ode group; this file covers the parts of `carleman.py`). All classical references are independent: the lifted generator $G_K$ assembled in numpy basis vector by basis vector on the padded layout per the product rule, the lifted initial state constructed analytically from tensor powers, the truncated linear system via `scipy.linalg.expm`, and the exact nonlinear solution via `scipy.integrate.solve_ivp` (Riccati, rtol=1e-12); the quantum programs run on the real backends reference, rir-pysparq, adapter-pysparq, and OriginIR-ext (UniQC full amplitude / `to_matrix`).

**Experiment design**: (a) lifted block — the Riccati problem $u'=-u+u^2$ ($F_1=-I$, $F_2$ the contraction matrix $C e_{i_0+2i_1}=\delta_{i_0,i_1}e_{i_0}$, $d=2$) with truncation $K=2$; the zero-signal block of $G_K$ is extracted (unitary `to_matrix` + `effective_block`, cross-checked with a reference zero-state sweep) and compared against the independent assembly; (b) lifted initial state — $(1,u_0,u_0^{\otimes2})/Z$ ($u_0=(0.6,0.8)$, $r=0.5$) exhaustively enumerated at full amplitude on four paths; (c) end to end — `carleman_qode` with an injected three-argument protocol solver (a minimal truncated-Taylor $e^{Gt}$ block-encoding assembly inside the verification script, a real quantum program, degree 3), $t=0.2$, the physical channel against a full-stack numpy simulation; (d) truncation trend — the truncation errors of the $K=1,2$ quantum runs against the exact solve_ivp Riccati solution, with $K=3$ filled in classically.

**Key metrics**:

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| carleman-lift-block | $K$=2, 10 qubits | originir+to_matrix, reference | maximum element-wise deviation of the $G_K$ block | 3.4e-16 |
| carleman-initial-state | 4 qubits | four paths | maximum full-amplitude deviation | 0.0 |
| carleman-riccati-endtoend | 18 qubits | reference, rir | implementation error (vs independent simulation) | 2.8e-17 |
| ditto | | | Taylor remainder (degree 3, informational) | 1.4e-4 |
| ditto | | | Carleman truncation error (vs the exact nonlinear solution) | 1.1e-3 |
| ditto | | | quantum direction vs the exact truncated-linear solution direction | 4.8e-5 |
| carleman-cutoff-trend | $t=0.2$ | quantum (K=1,2) + classical (K=3) | truncation error K=1 / K=2 / K=3 | 9.3e-3 / 1.1e-3 / 1.05e-4 |

The truncation error drops by about one order of magnitude per added order — direct numerical evidence of Carleman truncation convergence. The success-subspace probability (all signals zero) also matches the classical value pointwise in the cross-check.

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

Artifacts: `out/verification/ode.json` (5 `carleman-*` cases in total).

## Related links

- Source: `src/oracq/algorithms/qnlss/carleman.py`
- Tutorial: [Swapping QODE methods for the same linear problem](../../tutorials/differential-equations.md)
- API reference: [Carleman linearization](../../api/algorithms/qnlss/carleman.rst)
- Related pages: [QODE problem and protocol](qode-problem.md) · [LCHS](lchs.md) · [Schrödingerization](schrodingerization.md) (typical composition: Carleman → linear solver) · [CBMD](cbmd.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)
