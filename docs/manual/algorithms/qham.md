# General QHAM

**English** · <a href="../../zh/manual/algorithms/qham.html">简体中文</a>

> Category C6 · Module `oracq.applications.qham` · Stage V1

## Overview

QHAM turns finite-polynomial evolution PDEs into quantum-solvable linear systems. The input is an autonomous first-order time-evolution system

$$
u' = f + Lu + \sum_r B_r(u,\dots,u), \qquad u(0)=u_{in},
$$

allowing arbitrary finite components, finite spatial dimension, finite polynomials of the unknown fields and their spatial derivatives, known coefficients/forcing and their derivatives, and outer derivatives over complete monomials. The generator recurses out the $U_i$ of each order per the HAM (homotopy analysis) truncation order m, then applies quantum-adapted linearization (QCL) to the truncated HAM system, obtaining $Y'=GY$ on the lifted space and handing it to any QODE solver; the output selects the physical block $u_{sum}=\sum_i U_i$. The construction follows the QHAM paper ([arXiv:2411.06759](https://arxiv.org/html/2411.06759v2)); the mathematical rules are in the derivation reference.

The finite closure of the tensor words $Y_{(a_0,\dots,a_{k-1})}$ takes $\mathrm{grade}\cdot\sum a_j + \mathrm{rank} \le \mathrm{grade}\cdot m + 1$ (grade = max(1, D−1), where D is the highest polynomial degree of the PDE); nonlinear substitutions strictly decrease the HAM order sum; the closure represents the truncated HAM system exactly, with no further Carleman-style higher-order truncation.

## Interface and input model

The entry points come in four parts (`from oracq.applications.qham import ...`):

```python
PolynomialPDE.from_equations(equations, *, axes=("x",), label="polynomial_pde")   # Field/Known expressions
QHAMPlan(pde, order)                                    # lazy closure: blocks / row_terms / offset / locate
Discretization(pde, Grid(axes, shape, spacing), known)  # spatial discretization and classical reference
qham_input_model(plan, bindings, *, eta=-1.0, max_blocks=256, max_terms=4096)
```

API entry points: {obj}`QHAMPlan <oracq.applications.qham.linearization.QHAMPlan>`, {obj}`Discretization <oracq.applications.qham.reference.Discretization>`, {obj}`qham_input_model <oracq.algorithms.input_model.qham.qham_input_model>`

- PDEs are constructed with the {obj}`Field <oracq.applications.qham.pde.Field>`/{obj}`Known <oracq.applications.qham.pde.Known>` expression algebra (addition/subtraction/multiplication, non-negative integer powers, `d(axis, order)` derivatives; an unknown field in a denominator and non-integer powers raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`), and {obj}`dumps <oracq.infrastructure.serialization.dumps>`/{obj}`loads <oracq.infrastructure.serialization.loads>` round-trip through the JSON Schema.
- `bindings` is built by {obj}`gate_bindings(discretization, initial) <oracq.algorithms.input_model.qham.gate_bindings>` (small-scale materialized ports for comparison, single-port target ≤ 5 bits) or {obj}`structured_fd_bindings(discretization, initial) <oracq.applications.qham.stencils.structured_fd_bindings>` (shifts + component selection + same-point contraction, without materializing the $N^r\times N^r$ port matrix); open ports can also be declared with `QHAMBindings.declare(plan, state_width, port_specs, *, initial_norm, ...)`.
- The input model is ODE (generator + initial state + physical output window): the generator is carried by a BE, the initial value by an SP, and open slots are late-bound via {obj}`bind <oracq.infrastructure.linking.bind>`.

{obj}`qham_input_model <oracq.algorithms.input_model.qham.qham_input_model>` returns a {obj}`QHAMInputModel <oracq.algorithms.input_model.qham.QHAMInputModel>`:

| Member | Meaning |
|---|---|
| `generator` / `initial` | the BE of the lifted generator and the lifted initial-state preparation |
| `state_width` / `eta` | physical state bit width / homotopy parameter |
| `log_initial_norm` | logarithm of the initial norm (generation-time uses rescaled weights to avoid summing large powers directly) |
| `solve(qode, time)` | invokes the QODE protocol and selects the physical block, returning {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>` |
| `dissipative_shift(shift=None)` | explicit global shift G − μI (μ ≥ α_G), adapting to LCHS/CBMD, accumulating `growth_shift` |

There is also {obj}`open_qham_input(plan, bindings, *, generator_alpha, generator_signal, eta=-1.0, name="QhamGenerator") <oracq.algorithms.input_model.qham.open_qham_input>`: the whole generator stays an unresolved abstract BE (no fake empty body is fabricated for it), while the initial state is still generated structurally; {obj}`taylor_qode(generator, initial, time, *, degree=2) <oracq.algorithms.input_model.qham.taylor_qode>` provides a finite-Taylor candidate that requires no Hermitian/dissipative premises, used to wire up gate-level execution and comparison. The generator module attributes include `algorithm="qham_generated_qcl"`, `ham_order`, `pde_degree`, `tensor_rank_limit`, `raw_dimension`, `explicit_couplings`, `qham_plan`, and `correctness="linearization_witnessed; solver pending"`.

## Implementation notes

The generation chain is PDE → HAM recurrence (correction-term weight $-\eta(1+\eta)^{p}$, physical-term weight $1-(1+\eta)^{p}$) → ordered tensor-word closure → rectangular port placement and embedding → LCU composition into G. Each port is interpreted as a multilinear map (L: V→V, B_τ: r-fold tensor product of V → V, F: constants → V; the port coefficients already carry the PDE's signs and constants): {obj}`place_port <oracq.algorithms.input_model.qham.place_port>` places a port at the designated tensor position and permutes the remaining coordinates, and {obj}`embed_rectangular <oracq.algorithms.input_model.qham.embed_rectangular>` constrains both the input and output windows, preventing rectangular zero padding from contaminating neighboring blocks. {obj}`QHAMPlan <oracq.applications.qham.linearization.QHAMPlan>` is a lazy object — for quadratic without forcing at m = 20 it can describe $2^{21}$ function blocks without materializing them, and `row_terms`/`offset`/`locate` can be queried individually; explicit generation beyond the `max_blocks`/`max_terms` budget raises, while the plan itself remains valid.

The lifted initial state preserves the relative norms $r, r, r^2, \dots, r^K$ of the tensor words (r = ‖u_in‖; a one branch is added when forcing is present): branch labels are prepared from rescaled weights, a repeatable initial-value oracle is invoked multiple times on the corresponding target ranges (this is not cloning an unknown quantum state), and the labels are back-computed from the address ranges through a Boolean network; every preparation from a known zero input cleans its work, and the same work register is reused sequentially. When u_in = 0 with forcing, the initial state is just the one component; a zero initial value without forcing corresponds to the zero solution and should be returned directly on the classical side.

Applicability boundary: the target packing of a single BE is subject to the 64-bit limit (`plan.max_rank * n > 64` raises; lazy plans and row-by-row algorithms are not subject to it); `correctness="linearization_witnessed; solver pending"` indicates that what has been witnessed is the linearization algebra, while HAM/PDE convergence certification and QODE solution accuracy belong to the solver side; automatic selection of convergent h/m, time-dependent bindings, and the IQHAM outer restart are not implemented (see section 9 of the user guide for the full boundary).

## Validation approach

Category C6 (composition skeleton; acceptance criteria in `../../development/validation-plan.md` §2: component contracts satisfied + end-to-end small-instance semantic correctness + bind invariance). Three layers of evidence in `tests/core/test_qham_general.py:QhamGeneralTests`:

- Structure: construction and attribute assertions of the class; `test_rules_roundtrip_and_unsupported_nonlinearity` (PDE JSON round-trip equality, `1/u` and `u**0.5` rejected).
- Numerical: `test_lazy_closure_and_rank_unrank` — closure self-consistency for degrees 2–4 × orders 0–3 (nonlinear substitutions strictly decrease the order sum; offset/locate mutually inverse); the serialization of the $2^{21}$-block plan is < 3000 characters; exceeding the block budget on explicit generation raises. `test_quadratic_dimension_and_nontrivial_eta_weight` — for u′ = u², block_count = 2^(m+1) and raw_dimension(4) = 5^(m+1)+3; for η = −0.4 the homotopy weights are 0.24 / 0.4 / 0.4. `test_multidimensional_coupled_identity` — on a two-dimensional coupled PDE, the residual of G·lift against an independent chain-rule evaluation is < 1e-10. `test_m1_reduces_to_previous_special_case` — the m = 1 no-forcing generator is cross-checked column by column against the BE of the earlier special case `qham_lift_m1` (16×16 amplitudes times α, places = 10). There are additionally `test_chain_rule_identity_general_cases` (5 PDE families × η ∈ {−1, −0.4, 0.2}, same residual bound), `test_lazy_stencil_rows_match_explicit_linearization` (lazy row rules versus the explicit matrix, places = 11), `test_structured_initial_preserves_relative_tensor_norms`, `test_rectangular_forcing_and_contraction_do_not_spill`, `test_dissipative_adaptation_is_explicit`, and `test_zero_initial_forcing_and_reused_work`, adding details on initial values, rectangular embedding, shifts, and workspace reuse.
- Binding: `test_open_qode_inputs_keep_operator_and_initial_oracles` — the open ports declared by `QHAMBindings.declare` (`Qham_L`, `Qham_B_*`, `Qham_initial`) remain unresolved after `model.solve`; after `bind` swaps in gate implementations there are no unresolved slots and the IR serialization round-trip is unchanged; the slots kept by `open_qham_input` are exactly `{QhamGenerator, Qham_initial}`.

## Known gaps and planned stages

Consistent with the application-layer `qfvm.py`/`qham/` row of `validation-coverage.md`: no registered gaps, stage V1. Solution accuracy and convergence certification are marked pending per the `correctness` attribute; they are implementation boundaries on the QODE-solver side and are not listed as gaps of this module.

## Numerical validation

Paper-grade numerical validation executed on 2026-09-16 by `tests/verification/verify_qham_qfvm.py` (real backends: the PySparQ native RIR interpreter, the oracq built-in reference executor, the PySparQ event adapter, and the UniQC full-amplitude state vector; no mocks, no skips).

### Experiment design

- **Full generator-matrix cross-check**: an aux-register superposition (Σ|c⟩|c⟩) extracts the complete signal-0 matrix block of the BE in a single run, checked element-wise against the independent classical matrix in `applications/qham/reference.py`, and the support set is verified to be exact with no leakage into the padding subspace. Instances u′=−0.2u+0.1u² (order 2, grid 2, raw dimension 28, η=−0.4), once with structured shift ports and once with spectrally embedded Pauli ports; another instance u′=0.5·ν(x)·uₓ−0.2u (order 1, grid 4) compares the three input models stencil / spectral / QRAM angle table.
- **Homotopy step and contraction mapping**: the homotopy weights are checked against the explicit formulas; the generation row rules (`linear_action`) are checked against an independent HAM recurrence + tensor chain rule (`chain_rule`, including forced Burgers); the HAM partial sums are run against the analytic solution of the Riccati equation (Bernoulli closed form) for m=1..6, with the measured contraction factors compared against the theoretical value |1+η|.
- **Initial state and solve chain**: the lifted initial-state preparation is checked amplitude by amplitude against a direct construction from the tensor-word weights (71 qubits exceeds the UniQC budget of 24, so only the pysparq path); the `taylor_qode` end-to-end chain is checked against the classical (I+tG)Y_in, with the two input models cross-checked against each other; the explicit dissipative shift G−μI is validated as a full matrix on the complete 2ʷ space (including the −μ diagonal on the padding subspace); {obj}`make_qpde <oracq.algorithms.qpde.pde.make_qpde>` / {obj}`qpde_solver <oracq.algorithms.qpde.pde.qpde_solver>` from `algorithms/pde.py` get a numerical passthrough on the four-cycle Laplacian.

### Key metrics

| Case | Scale | Backend path | Value |
|---|---|---|---|
| generator full matrix (stencil ports) | 28×28, α=2.968 | rir-pysparq | max_error 1.46e-17, support exact |
| generator full matrix (spectral ports) | 28×28, α=3.336 | rir-pysparq | max_error 7.96e-18 |
| generator four-path cross-check | 14 qubits | reference / rir / adapter / originir-ext | pairwise deviation 8.67e-18 |
| input-model consistency | 28×28 | rir-pysparq(+QRAM) | stencil 1.40e-17, spectral 6.97e-18, QRAM 5.99e-05 (angle-quantization bound 9.20e-03) |
| QRAM coefficient angle encoding pointwise | 4 addresses, α=0.75 | rir-pysparq + originir-ext | amplitude error 2.15e-03 (bound 6.14e-03), cross-path 0.0 |
| lifted initial state | raw 28, 71 qubits | rir-pysparq | max_error 1.11e-16, log norm −1.2525710053499335 agrees, work cleaned |
| Taylor QODE end-to-end | t=0.01, degree 1 | rir-pysparq | per-model errors 1.11e-16 each, mutual cross-check 8.87e-25 |
| dissipative-shift full matrix | 32×32, μ=2.968 | rir-pysparq | max_error 5.70e-17 |
| homotopy weights/row rules | 2 PDEs × η∈{−1,−0.4,0.2} | classical independent | weights exact, row rules vs chain rule 1.39e-17 |
| HAM contraction convergence | t=0.5, m=1..6 | classical independent (RK4 vs analytic solution) | η=−0.4 mean contraction factor 0.604 (theory 0.6); η=−0.8 gives 0.208 (theory 0.2) |
| pde.py wrapper passthrough | four-cycle, t=0.05 | rir-pysparq | max_error 2.11e-18 |

### Reproduce

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_qham_qfvm.py
```

Artifact: `out/verification/qham_qfvm.json` (full metrics and criteria of the 15 cases).

## Related links

- Source: `src/oracq/applications/qham/` (pde / linearization / reference / stencils) and `src/oracq/algorithms/input_model/qham.py` (quantum assembly)
- Tutorial: [generating QHAM inputs from PDE expressions](../../tutorials/qham.md)
- User guide: [General QHAM automatic generation](../qham.md); mathematical derivation: [QHAM derivation](../../reference/qham-derivation.md)
- API reference: [QHAM](../../api/algorithms/input_model/qham.rst), [PDE models and adaptation](../../api/applications/qham/pde.rst), [structured stencil ports](../../api/applications/qham/stencils.rst)
- Same-group page: [Quantum Finite Volume Method](qfvm.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)
