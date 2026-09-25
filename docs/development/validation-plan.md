# Quantum algorithm implementation validation plan (review draft)

**English** · <a href="../../zh/development/validation-plan.html">简体中文</a>

This document systematically plans oracq's algorithm validation system: it
first classifies all algorithms by "how correctness is judged", then
prescribes a standard witness technique for each class, and finally gives the
all-algorithm validation matrix and an implementation roadmap. Review focus:
whether the classification in section 2 is complete, whether the cross-cutting
mechanisms in section 4 deserve infrastructure investment, and the stage
split and priorities of section 5.

## 1. General validation principles

- **What is validated is the semantics of RIR programs.** Algorithms are
  generators; witnesses must act on the generated program (executed through
  the reference simulator or a real backend), not on the generator's internal
  data structures.
- **At least three layers of evidence per algorithm**: structural validity
  (IR checks, signatures, contract properties), small-scale numerical
  witnesses (the reference simulator cross-checked against analytic solutions
  or classical references), and bindability (abstract declarations keep their
  semantics after being bound via bind).
- **Three-layer input-model consistency** is this repository's distinctive
  invariant: the abstract / gate / qram bindings of the same algorithm must
  produce consistent observable semantics (QRAM word-width quantization error
  is allowed, with an error bound stated).
- **Determinism first**: probabilistic assertions must operate at the
  distribution level (closed-form expectations, support, moments); single
  sample threshold assertions are forbidden; simulations must use fixed seeds
  or exact state vectors.
- The existing layering is kept: L0 static checks (ruff, layout,
  documentation warnings-as-errors) → L1 structure → L2 numerical witnesses
  (tests/core) → L3 cross-check catalog (applications/catalog.py) → L4 real
  backends (tests/integration, cross-checked on real uniqc/pysparq machines,
  no mocks) → L5 documentation doctests.

## 2. Algorithm classification (by how "correct" is judged)

| Class | Acceptance criterion | Standard witness technique | Main risk |
|---|---|---|---|
| C1 exact discrete semantics | Output distribution / action unitary **equal pointwise** to a truth table or exact construction | Exhaustive all-input cross-checks on small instances; the unitarity identity W†W=I; exact recovery of null spaces / secret strings | Essentially no flakiness |
| C2 approximate continuous semantics | Action operator ≈ the target continuous function, with an explicit error bound | Compare within tolerance against classical references (analytic solutions / matrix exponentials / high-precision numerics); **convergence scans**: fit error vs degree/κ/Δt to the theoretical rate | Tolerance choice masking systematic error; numerical ill-conditioning (e.g. high-degree QSP synthesis) |
| C3 probability-distribution semantics | The output **distribution** equals the closed-form expectation (mean, variance, support) | Derive the distribution from exact state vectors and cross-check per-outcome probabilities against closed forms (Krawtchouk/binomial etc.); for QAE-style algorithms, verify the estimate falls inside the theoretical confidence interval | Sampling jitter → sample assertions forbidden, exact distributions only |
| C4 heuristic/optimization semantics | Strictly better than a random baseline, and small instances reach known optima / theoretical scores | Distribution-level comparison against the expectation of random guessing; exhaustive-optimum cross-checks on small instances; cross-checks of the theoretical score formulas | Boundary-assertion jitter (happened once for the C3-class DQI: assertGreater failed when exactly equal to the baseline) |
| C5 data-access layer | Pointwise-correct query semantics + three-layer binding consistency + uncomputation | Per-address truth-table cross-checks; abstract/gate/qram triple-binding consistency; LocalExit-enforced uncomputation checks; marginal-distribution checks for superposed queries | QRAM quantization error needs an explicit bound |
| C6 composition skeleton | Component contracts satisfied + correct end-to-end semantics on small instances + bind invariance | End-to-end simulation + identical results before and after binding + structural-property checks (register layout, Repeat not expanded) | Components individually correct but composed semantics misaligned (register-mapping errors) |

## 3. All-algorithm validation matrix

Status: ✅ adequately witnessed; 🟡 witnessed but with gaps; ❌ no witness.
Modules are listed by layer.

### infrastructure / data-access layer (C5)

| Module | Contents | Status | Gaps and plans |
|---|---|---|---|
| oracles.py | XorDatabase, StatePreparation, SparseAccess, BlockEncoding across three layers | ✅ | Add a parametrized "triple-binding consistency" test (cross-check the same declaration binding by binding); each binding is currently witnessed independently |
| data_loading.py | QROM / Select-Swap | ✅ | Resource formulas already cross-checked; consider adding {obj}`qrom_cost <oracq.algorithms.input_model.data_loading.qrom_cost>` to catalog report attributes |
| mathfunc/ | Classical functions → reversible circuits | ✅ | Add an explicit witness for the fixed-point quantization error bound |
| prepare_select.py | PREPARE/SELECT, alias | ✅ | alias's clean_work=False is honestly documented; add alias triple-binding consistency |

### Exact discrete algorithms (C1)

| Module | Algorithms | Status | Gaps and plans |
|---|---|---|---|
| oracle_algorithms.py | DJ/BV/Simon | ✅ | — |
| fourier.py | QFT, QFT addition | ✅ | — |
| arithmetic.py | Fixed-point add/sub/mul, Boolean networks | ✅ | — |
| number_theory.py | Modular multiplication, order finding | ✅ | The small-scale witnesses already note they do not extrapolate to full Shor |
| error_correction.py | Repetition codes | ✅ | — |
| walks.py | Coined cycle walk | ✅ | — |
| graph_walks.py | Adjacency oracle, Szegedy, MNRS | ✅ | The cycle's failure has been verified as a theoretical result; adding Johnson graphs requires new witnesses (element distinctness) |

### Approximate continuous algorithms (C2)

| Module | Algorithms | Status | Gaps and plans |
|---|---|---|---|
| transforms.py | Qubitization walk, QSVT sequence, OAA | ✅ | The sequence conventions are pinned by qsvt.py (pointwise agreement 1e-10) |
| qsvt.py | Phase synthesis, inversion, filtering, HS simulation, fixed-point search | ✅ | **Convergence scans missing**: batch curves of error vs degree/κ are not automated; the ill-conditioned synthesis-degree boundary (≳16) already has negative-case coverage: tests/core/test_qsvt.py:PhaseSynthesisTests.test_input_validation and TransformWitnessTests.test_transform_input_validation (covering parity violations, \|coeff\|>1, unsaturated endpoints, imaginary-part parity, ill-conditioned κ, degree 0, t=0, Δ≥1, odd/even degree violations, and similar scenarios) |
| hamiltonian.py | Trotter, Taylor, hamsim | ✅ | Trotter order-error-rate scans not done |
| sparse.py / block_encoding.py | Sparse BE, BE algebra | ✅ | — |
| qlss.py | CKS, Costa | 🟡 | Paradigm witnesses exist; missing an end-to-end cross-check case against the HHL paper's reference values (could enter the catalog) |
| ode.py / lchs.py / cbmd.py / schrodingerization.py / carleman.py | ODE/PDE solvers | 🟡 | Each has witnesses; missing a unified "analytically solvable ODE family" convergence benchmark (the same problem through all solvers) |
| sde.py | Fokker–Planck | ✅ | OU-moment cross-checks; add an end-to-end catalog case for the SDE solver and LCHS |
| density.py | Purification, Gibbs | ✅ | Convergence scans added: tests/core/test_density.py:GibbsTests.test_error_convergence_decreases (error ∈ {0.4,0.2,0.1} monotone non-increasing + each step ≤ error, measured 5.6e-4/8.7e-5/8.7e-5) and test_error_bound_uniform_in_beta (β ∈ {0.2,0.5,1.0} grid, each point ≤ 0.1) |
| lowrank.py | DF/THC → BE | ✅ | α upper-bound tightness added: tests/core/test_lowrank.py:DoubleFactorizationTests.test_alpha_matches_closed_form_eigenvalues (independent 2×2 closed-form cross-check, places=10), test_df_alpha_tighter_than_pauli (DF α=1.4 ≤ Pauli α=1.5); ThcTests.test_thc_alpha_matches_hand_computed_bound (THC α=1.996 independent hand computation) |

### Probability-distribution algorithms (C3)

| Module | Algorithms | Status | Gaps and plans |
|---|---|---|---|
| estimation.py | QPE, QAE, quantum counting, Hadamard/SWAP test | ✅ | The QAE confidence-interval claim is unwitnessed (only the point estimate is verified) |
| search.py | Grover, amplitude amplification | ✅ | — |
| integration.py | Heinrich summation/integration | ✅ | Closed-form mean cross-check + heinrich_rate parametric validation already witnessed (tests/core/test_integration.py:SumPreparationTests / QuantumSumTests / RateTests) |
| gradient.py | Jordan gradient | ✅ | Exact on linear + perturbed convergence already witnessed (tests/core/test_gradient.py:GradientTests.test_perturbed_linear_concentrates_with_grid_bits includes the failure-probability decay rate q_{m+1} ≤ 0.34·q_m) |
| qpca.py | QPCA, density-matrix exponentiation | ✅ | Eigenvalue readout peak cross-check + Δt first-order error rate already witnessed (tests/core/test_qpca.py:DensityMatrixExponentiationTests / QpcaTests) |

### Heuristic/optimization algorithms (C4)

| Module | Algorithms | Status | Gaps and plans |
|---|---|---|---|
| variational.py | QAOA-MaxCut, VQE, ansatz | ✅ | One-sided witness; add the strong witness of "reaching the optimal cut on small instances" |
| dqi.py | DQI | ✅ | Pointwise distribution cross-check against closed forms (Krawtchouk); the identity-decoder baseline assertion was fixed to allow boundary equality |

### Application layer

| Module | Contents | Status | Gaps and plans |
|---|---|---|---|
| qfvm.py / qham/ | QFVM, QHAM | ✅ | Catalog cross-checks already exist |
| catalog.py / gallery.py | 33-case cross-check catalog | ✅ | New algorithms (qsvt inversion, graph_walks search, data_loading, integration, qpca) not yet registered |

## 4. Cross-cutting validation mechanisms (infrastructure investments)

1. **Invariant test library** (new proposal, tests/core/witness.py or similar):
   - `assert_unitary(program)`: W†W=I on the reference simulator (sampled
     basis states).
   - `assert_uncomputation(program, work_regs)`: work-qubit restoration
     (reusing the LocalExit mechanism).
   - `assert_bind_invariant(abstract_program, bindings)`: the bindings of the
     same open declaration produce consistent semantics (error bound
     parameterized).
   - `assert_block_equals(be, matrix, alpha)`: cross-check of the BE's (0,0)
     block.
   Benefit: witnesses for new algorithms go from "each hand-rolls its own" to
   declarative composition, with a unified regression yardstick.
2. **Determinism policy** (codified): every probabilistic assertion must be a
   closed-form cross-check against the exactly simulated distribution; when
   the theoretical value lands exactly on the assertion boundary, use ≥/≤
   with a numeric tolerance — strict inequalities are not allowed (the DQI
   boundary case was fixed this way).
3. **Convergence-scan framework**: provide parameter-scan tooling for C2
   algorithms (error vs degree/κ/Δt), asserting the fitted slope ≥ a lower
   bound on the theoretical rate; store results as JSON for validation.json
   to aggregate.
4. **Resource-estimation cross-checks**: for every module providing a cost
   function (qrom_cost, the QSP degree formula, Trotter step counts), the
   generated program's structural properties must agree with the cost formula
   (checked in both directions).
5. **Catalog registration rule**: once a new algorithm is finished it must
   register catalog cases (gate/qram input variants) and be wired into the
   real-backend cross-checks of tests/integration (no mocks).

## 5. Implementation roadmap

| Stage | Scope | Prerequisite |
|---|---|---|
| V1 | Invariant test library (4 assertion primitives) + witnesses for the in-flight modules (density/gradient/lowrank/integration/qpca) filled in per the matrix (done, 0.x iteration) | None |
| V2 | Convergence-scan framework + convergence benchmarks for the qsvt/hamsim/ODE solvers wired into validation.json | V1 |
| V3 | Register the new algorithms in the catalog (qsvt inversion, MNRS search, Select-Swap, Heinrich integration, QPCA) + extend the real-backend cross-checks in tests/integration | V1 |
| V4 | Roll out triple-binding consistency parametrized tests (oracles/prepare_select/graph_walks/data_loading) | V1 |

## 6. Relation to existing acceptance

`validation.md` / `validation.json` remain the per-version acceptance records;
once this plan lands, the V1–V4 deliverables are written into them by
version. The algorithm coverage workboard (algorithm-coverage.md) records
"what is implemented" while this document records "how correctness is
proven"; the two correspond line by line.
