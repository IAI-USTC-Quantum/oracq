# Fokker–Planck Input Model

**English** · <a href="../../../zh/manual/algorithms/fokker-planck.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.qode.sde`](../../api/algorithms/qode/sde.rst) · Stage V3

## Overview

The probability density of the one-dimensional Itô stochastic differential equation $dx=a(x)\,dt+\sqrt{2D(x)}\,dW$ satisfies the Fokker–Planck equation

$$
p'=-\partial_x(a\,p)+\partial_{xx}(D\,p).
$$

This module performs a zero-flux finite-volume discretization on a uniform grid, obtaining a discrete generator $G$ whose column sums vanish (convention $p'=Gp$); after $G$ becomes a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` via an explicit Pauli expansion, it is assembled into a {obj}`QODEProblem <oracq.algorithms.qode.ode.QODEProblem>` / {obj}`LinearODE <oracq.algorithms.qode.ode_models.LinearODE>` and handed to the existing linear ODE solvers (LCHS and others). The dissipativity of the generator is an input-model declaration, not proven by the language. The module also provides pure-Python classical witness tools (matrix exponential, Euler evolution, moments, discrete and continuous stationary references) for numerical cross-checks.

## Interface and input model

```python
FokkerPlanckProblem(drift, diffusion, grid)
problem.generator_matrix() / generator_encoding() / qode_problem(initial=None)
problem.linear_ode(initial=None)
sde_state_preparation(probabilities, *, implementation="gate", angle_width=8, work_width=0)
sde_state_angles(probabilities, *, angle_width=8)
matrix_exponential(matrix, time=1.0)
evolve_distribution(matrix, initial, time, *, steps=None)
distribution_moments(points, probabilities, orders=(1, 2))
stationary_distribution(problem)
boltzmann_distribution(problem)
```

API entries: {obj}`FokkerPlanckProblem <oracq.algorithms.qode.sde.FokkerPlanckProblem>`, {obj}`sde_state_preparation <oracq.algorithms.qode.sde.sde_state_preparation>`, {obj}`sde_state_angles <oracq.algorithms.qode.sde.sde_state_angles>`, {obj}`matrix_exponential <oracq.algorithms.qode.sde.matrix_exponential>`

- {obj}`FokkerPlanckProblem <oracq.algorithms.qode.sde.FokkerPlanckProblem>`: the input model is CP (drift, diffusion, and grid coordinates parameterized directly), producing a BE + SP. `drift` / `diffusion` may be constants or pointwise vectors (diffusion requires $D\ge 0$), and `grid` must be a strictly increasing uniform grid; derived attributes `points`, `spacing`, `size`, `width` (the number of grid points must be a power of two).
- {obj}`generator_encoding() <oracq.algorithms.input_model.qham.generator_encoding>`: a small-grid BE from an explicit Pauli expansion, supporting only `width <= 5` (at most 32 points); large instances need an access oracle, and no quantum speedup is claimed.
- `qode_problem(initial=None)`: returns a `QODEProblem` (`dissipative=True`, `evidence="fokker_planck_zero_flux_finite_volume; dissipative caller-declared"`), with the default initial state a uniform superposition.
- `linear_ode(initial=None)`: a `LinearODE` view of the same problem as $u'=-Au$ ($A=-G$), with label `"fokker_planck_minus_A"`.
- {obj}`sde_state_preparation <oracq.algorithms.qode.sde.sde_state_preparation>`: encodes the discrete distribution as a {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` with amplitudes $\sqrt{p_i}$; `"gate"` is the reused-rotation implementation and `"qram"` the QRAM angle-table implementation (paired with {obj}`sde_state_angles <oracq.algorithms.qode.sde.sde_state_angles>` to generate the binding data).
- Classical witness tools: {obj}`matrix_exponential <oracq.algorithms.qode.sde.matrix_exponential>` (scaling and squaring + Taylor), {obj}`evolve_distribution <oracq.algorithms.qode.sde.evolve_distribution>` (exact via the matrix exponential or explicit Euler), {obj}`distribution_moments <oracq.algorithms.qode.sde.distribution_moments>`, {obj}`stationary_distribution <oracq.algorithms.qode.sde.stationary_distribution>` / {obj}`boltzmann_distribution <oracq.algorithms.qode.sde.boltzmann_distribution>`.

## Implementation notes

The discretization happens on the faces between adjacent grid points, with the face coefficient taken as the two-point average; the forward/backward rates are $\bigl(\tfrac{a}{2}+\tfrac{D}{h}\bigr)/h$ and $\bigl(\tfrac{a}{2}-\tfrac{D}{h}\bigr)/h$ respectively, written into $G$ in conservative form so that each column sums to exactly zero (zero-flux boundary). The discrete stationary state has an exact closed form: the ratio of adjacent points is $(1+u)/(1-u)$, where $u=ah/(2D)$ is the effective grid Péclet number on the face; when $|u|\ge 1$ the central-difference stationary state is no longer positive and an error is raised directly. The discrete adjacent ratio of the continuous stationary reference $p\propto\exp(\int a/D\,dx)$ is $e^{2u}$, differing from the discrete stationary state by $O(h^2)$.

The dissipation declaration and the initial-value norm enter the solving protocol through the `evidence` / `dissipative` fields of `QODEProblem` (see [QODE problem and protocol](qode-problem.md)); the QRAM state preparation accepts only non-negative real amplitudes, and signed or complex phases require a separate implementation.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2): compare against classical references (analytic moments / closed-form stationary states) within tolerance. Three layers of evidence, located in `tests/core/test_sde.py`:

- Structure: the construction assertions of `GeneratorDiscretizationTests` / `StatePreparationTests` / `SolverContractTests`; `SolverContractTests.test_invalid_inputs_rejected` and `StatePreparationTests.test_invalid_probabilities_rejected` cover violations such as negative diffusion, non-uniform grids, point counts that are not powers of two, invalid initial-state widths, and invalid probabilities.
- Numerical: `GeneratorDiscretizationTests.test_ou_moments_match_closed_form` — after a 16-point OU process ($\theta=1$, $D=0.5$, $t=0.4$) is evolved, the mean matches the closed form $\mu_0 e^{-\theta t}$ within tolerance 2e-3 and the variance matches $\sigma_0^2 e^{-2\theta t}+\frac{D}{\theta}(1-e^{-2\theta t})$ within tolerance 2e-2; `test_stationary_approaches_boltzmann` — after a uniform initial state relaxes for 30 time units, the pointwise difference from the discrete stationary state is 1e-6, and the discrete stationary state differs from the Boltzmann reference pointwise by 2e-2; `test_column_sums_vanish` (column sums vanish, places=12). Auxiliary witnesses: `test_evolution_conserves_probability` (exact and Euler evolutions conserve probability, places=9, pointwise delta 1e-4) and `test_matrix_exponential_matches_euler_limit` (first-order Euler convergence rate $e_{\text{fine}}<0.6\,e_{\text{coarse}}$).
- Binding: `SolverContractTests.test_qode_problem_accepted_by_lchs` — a `QODEProblem` on a 4-point OU grid with the Boltzmann distribution as the initial state goes through `check().require()` and `solve` of {obj}`linear_qode("lchs") <oracq.algorithms.qode.ode.linear_qode>` (Cauchy cutoff=0 plan, Taylor degree 1); the output width is 2, the `qode_dissipative_promise` pass-through is correct, and the program round-trips through serialization; `StatePreparationTests.test_gate_preparation_amplitudes_are_root_probabilities` and `test_qram_preparation_declares_angle_table` cover the two state-preparation implementations, gate and QRAM.

## Known gaps and planned stages

An end-to-end catalog case for SDE and LCHS (a catalog-level registration of discretization → QODEProblem → solve → readout) is missing; it is assigned to the stage V3 catalog registration plan, consistent with the sde.py row of `validation-coverage.md`.

## Related links

- Source: `src/oracq/algorithms/qode/sde.py`
- API reference: [SDE/Fokker–Planck input model](../../api/algorithms/qode/sde.rst)
- Related pages: [QODE problem and protocol](qode-problem.md) · [LCHS](lchs.md) (the consuming protocol of this input model)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_nt_qlss_sde.py` (the nt_qlss_sde group). All classical oracles are independent: the `scipy.linalg.expm` propagator, the OU closed-form moments, fixed-seed Euler–Maruyama Monte Carlo (numpy, seed=20260916), and the null-space vector of $G$ (numpy eigendecomposition); the quantum parts are validated on real backends (reference, rir-pysparq, adapter-pysparq, OriginIR-ext).

**Experiment design**: (a) input model — the column sums, evolution probability conservation, and spectral abscissa of the 16-point OU grid ($\theta=1$, $D=0.5$) generator, with the library matrix exponential against scipy.expm; (b) explicit Euler (`evolve_distribution(steps=k)`) against the time-convergence order of the scipy exact propagator; (c) three-way cross-check of OU moments — Fokker–Planck grid moments vs the closed forms $\mu_0 e^{-\theta t}$ / $\sigma_0^2 e^{-2\theta t}+\frac{D}{\theta}(1-e^{-2\theta t})$ vs a $2^{18}$-particle Monte Carlo ($dt=10^{-3}$); (d) stationary state — the discrete stationary state against the null-space vector, long-time relaxation convergence, and the $O(h^2)$ gap to the continuous Boltzmann reference; (e) quantum channel — the gate implementation of `sde_state_preparation` has amplitudes exactly $\sqrt{p_i}$ (full amplitude, four paths), the QRAM angle-table implementation is checked against an independent quantized-angle rotation tree (implementation error and method error reported separately), and the Pauli block encoding of the 4-point OU generator is driven from the zero basis state to verify that the $(\text{signal}=0)$ block amplitude times $\alpha$ is exactly $G$ (four paths).

**Key metrics**:

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| sde-generator-conservation | 16 points | scipy oracle | column sums / probability conservation / matrix-exponential deviation | 0.0 / 4.4e-16 / 6.9e-16 |
| sde-euler-convergence-order | $k$=2000→8000 | scipy oracle | error / measured order | 3.12e-6→7.80e-7 / 1.0001, 1.0001 |
| sde-ou-moments-vs-monte-carlo | $t=0.4$ | numpy MC | mean: FPE vs closed form / MC vs closed form | 7.0e-9 / 4.2e-4 |
| ditto | | | variance: FPE vs closed form / MC vs closed form | 1.94e-2 (grid truncation) / 1.1e-3 (statistical fluctuation) |
| sde-stationary-boltzmann | 16 points, $T=30$ | numpy/scipy | TVD: stationary vs null space / relaxation / Boltzmann | 2.8e-15 / 4.2e-14 / 1.2e-3 |
| sde-state-preparation-gate | 2 qubits | four paths | maximum deviation of amplitudes from $\sqrt{p}$ | 5.6e-17 |
| sde-state-preparation-qram | angle_width 8 | two paths | implementation error / quantization method error | 2.2e-16 / 2.9e-3 |
| sde-generator-encoding | 4 points, $\alpha=4$ | four paths | block amplitude × $\alpha$ vs $G$ | 3.3e-16 |

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

Artifacts: `out/verification/nt_qlss_sde.json` (7 `sde-*` cases in total).

## Numerical validation (ode group)

This section is the ode group's supplementary validation of the `QODEProblem` consumption chain (alongside the input-model validation of the nt_qlss_sde group in the previous section), in `tests/verification/verify_ode.py`. All classical references are independent: `scipy.linalg.expm` and the in-library pure-Python `sde.matrix_exponential` cross-check each other, plus a numpy per-branch LCHS simulation; the quantum programs run on the real backends reference and OriginIR-ext.

**Experiment design**: the zero-flux discretized generator of a 4-point OU grid ($\theta=1$, $D=0.5$, grid $(-1.5,-0.5,0.5,1.5)$) is assembled into a `QODEProblem` via `generator_encoding()` (explicit Pauli-expansion BE) and `qode_problem(initial)` (Boltzmann root-amplitude initial state, `dissipative=True`), then solved at $t=0.3$ by the `solve` of `linear_qode("lchs")` (Cauchy cutoff=2, Taylor degree 2); the post-selected block is checked against an independent simulation of the $L+iH$ decomposition and the exact propagator.

**Key metrics**:

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| lchs-fokker-planck-ou | 16 qubits | reference, originir | implementation error (vs independent simulation) | 1.1e-16 |
| ditto | | | method error (quadrature + Taylor remainder, informational) | 7.9e-2 |
| ditto | | | cross-check of the two classical references (scipy vs sde matrix exponential) | 1.1e-16 |
| ditto | | | `qode_dissipative_promise` pass-through | True |

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

Artifacts: `out/verification/ode.json` (the `lchs-fokker-planck-ou` case).
