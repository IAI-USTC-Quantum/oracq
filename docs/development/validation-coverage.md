# Validation coverage matrix

**English** · <a href="../zh/development/validation-coverage.html">简体中文</a>

This document records "where each algorithm's evidence lives", corresponding
line by line with the "how correctness is proven" of `validation-plan.md` and
the "what is implemented" of `algorithm-coverage.md` (see the division of
labor in `validation-plan.md` §6). `validation.md` / `validation.json` keep
carrying the per-version acceptance-snapshot duty; this table describes the
live coverage state of the current development version, and new algorithms or
new witnesses must synchronize it (maintenance rules in the final section).
Per-algorithm interface and validation details live on each algorithm page
(`../zh/manual/algorithms/index.html`).

Categories and acceptance criteria follow `validation-plan.md` §2 (C1 exact
discrete, C2 approximate continuous, C3 probability distribution, C4
heuristic/optimization, C5 data access, C6 composition skeleton); stages
(V1–V4) are in §5 of the same document. The three evidence layers mean
structural validity, small-scale numerical witnesses, and bindability (plan
general principles §1).

## Witness primitives (`tests/core/witness.py`)

`tests/core/witness.py` provides four assertion primitives shared across
algorithms; every primitive takes a `unittest.TestCase` instance as its first
argument and raises `AssertionError` on failure. The file name lacks the
`test_` prefix, so unittest discover never collects it; the self-tests of its
correctness live in `tests/core/test_witness.py` (four test classes —
`UnitaryWitnessTests`, `UncomputationWitnessTests`,
`BindInvariantWitnessTests`, `BlockEqualsWitnessTests` — 11 cases in total).

| Primitive | Signature | Purpose | Conventions |
|---|---|---|---|
| `assert_unitary` | `(case, program, *, places=9, samples=None)` | Sampled witness of W†W = I: each basis-state column normalized (Σ\|a\|²=1) and pairwise orthogonal | `samples` explicitly gives the list of basis-state indices; by default, small register spaces (≤16 basis states) take all of them, and large spaces sample 16 with the fixed seed `random.Random(0)`. Failure messages include the offending column and its \|a\|² or inner product |
| `assert_uncomputation` | `(case, program, *, initial=None, work_registers=None)` | Uncomputation witness: the simulator forces uncomputation at LocalExit, and the primitive converts {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` into `AssertionError`; the root registers listed in `work_registers` take value 0 in all basis states with nonzero amplitude | Returns the state of {obj}`simulate <oracq.infrastructure.execution.simulate>` so the caller can keep asserting. Zero is decided with the threshold `_ZERO_AMPLITUDE = 1e-9` (the simulator drops amplitudes with \|a\| < 1e-15, so the uncomputation check must be looser) |
| `assert_bind_invariant` | `(case, abstract_program, bindings, *, tolerance=0.0)` | The candidate implementations of the same abstract slot must give consistent observable distributions | `abstract_program` must have exactly one unbound slot (resolved with {obj}`unresolved <oracq.infrastructure.linking.unresolved>`); `bindings` is a `{label: Operation}` dict. With `tolerance=0` the cross-check is exact; bindings with quantization error such as QRAM pass an error bound (e.g. `0.02`), kept separate from the 1e-12 precision of the `integration.py` three-layer consistency |
| `assert_block_equals` | `(case, be, matrix, *, places=9)` | The BE's (0,0) block is cross-checked column by column against the dense matrix (multiplying `alpha` back in) | Relies only on the `block_column` helper: simulate the initial basis state first, then multiply the amplitudes of the `(row, 0)` branch by `be.alpha` to recover column `column` |

The helper `block_column(be, column)` is exported separately so callers can
reuse it in custom assertions (e.g. `test_matches_pauli_encoding_block` in
`test_lowrank.py` builds the other side's BE itself and cross-checks row by
row).

## Coverage matrix

### infrastructure / data-access layer (C5)

| Module | Algorithms | Class | Structural evidence | Numerical evidence | Binding evidence | Gap | Stage |
|---|---|---|---|---|---|---|---|
| `oracles.py` | <a href="../zh/manual/algorithms/xor-database.html">XorDatabase</a>, <a href="../zh/manual/algorithms/state-preparation.html">StatePreparation</a>, <a href="../zh/manual/algorithms/sparse-access.html">SparseAccess</a> | C5 | `tests/core/test_contracts.py:OracleContractTests` | `tests/core/test_contracts.py:OracleContractTests.test_dirty_work_is_rejected_by_algorithm` | The same class's `test_capability_restrictions_survive_annotation_and_wrapping` | Triple-binding consistency parametrization (V4 rollout) | V4 |
| `data_loading.py` | <a href="../zh/manual/algorithms/select-swap.html">Select-Swap QROM</a>, <a href="../zh/manual/algorithms/qrom-lookup.html">QROM Lookup</a> | C5 | `tests/core/test_data_loading.py:DataLoadingTests.test_generated_operations_carry_cost_attributes` | The same class's `test_select_swap_matches_gate_database_per_address` / `test_cost_model_matches_formulas_and_tradeoff_curve` | Both `gate/qram` bindings covered indirectly by `test_result_invariant_across_partitions_and_qrom_lookup_baseline`; triple-binding consistency parametrization awaits V4 | Triple-binding parametrization (V4) | V4 |
| `mathfunc/` | <a href="../zh/manual/algorithms/mathfunc.html">classical functions → reversible circuits</a> | C5 | `tests/core/test_mathfunc.py:MathFunctionTests.test_helper_module_reuse_and_roundtrip` | `test_all_elementary_families_construct` / `test_roe_is_a_compiled_pure_function` | — | Explicit witness of the fixed-point quantization error bound not yet added | V2 |
| `prepare_select.py` | <a href="../zh/manual/algorithms/prepare-select.html">PREPARE/SELECT</a>, <a href="../zh/manual/algorithms/alias-preparation.html">alias</a> | C5/C1 | `tests/core/test_prepare_select.py:PrepareSelectTests.test_gate_and_qram_prepare_bindings_agree` / `test_pauli_hamiltonian_input_and_single_term_shortcut` | The same class's `test_block_encoding_recovers_hamiltonian_over_alpha` / `test_alias_block_encoding_approaches_hamiltonian` / `test_qram_prepare_block_encoding_approaches_hamiltonian` | `test_gate_and_qram_prepare_bindings_agree` already covers both gate/qram bindings; `abstract_prepare → gate_prepare` is covered by `test_abstract_prepare_binds_inside_block_encoding` | alias triple-binding consistency (gate/qram currently lack an independent scenario) | V4 |

### Exact discrete algorithms (C1)

| Module | Algorithms | Class | Structural evidence | Numerical evidence | Binding evidence | Gap | Stage |
|---|---|---|---|---|---|---|---|
| `oracle_algorithms.py` | <a href="../zh/manual/algorithms/deutsch-jozsa.html">DJ</a>/<a href="../zh/manual/algorithms/bernstein-vazirani.html">BV</a>/<a href="../zh/manual/algorithms/simon.html">Simon</a> | C1 | `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects` | The same class's `test_bernstein_vazirani_recovers_secret_with_affine_bias` / `test_simon_sampling_and_rank_aware_postprocess` | — | — | V1 |
| `fourier.py` | <a href="../zh/manual/algorithms/qft.html">QFT</a>, <a href="../zh/manual/algorithms/fourier-addition.html">QFT addition</a> | C1 | `AlgorithmExpansionTests` | `test_qft_matches_positive_fourier_matrix` / `test_fourier_add_all_basis_inputs` | — | — | V1 |
| `arithmetic.py` | <a href="../zh/manual/algorithms/fixed-point-arithmetic.html">fixed-point add/sub/mul</a>, <a href="../zh/manual/algorithms/boolean-networks.html">Boolean networks</a> | C1 | `tests/core/test_stage2.py:Stage2StructureTests` | `test_arithmetic_construction_and_basis` / `test_gate_network_is_executable_small_example` | — | — | V1 |
| `number_theory.py` | <a href="../zh/manual/algorithms/modular-multiplication.html">modular multiplication</a>, <a href="../zh/manual/algorithms/order-finding.html">order finding</a> | C1 | `AlgorithmExpansionTests` | `test_modular_multiplication_total_permutation_and_order` (covers the factors_from_phase decoding path) | — | Already noted as not extrapolating to full Shor (original wording kept) | — |
| `error_correction.py` | <a href="../zh/manual/algorithms/repetition-codes.html">repetition codes</a> | C1 | `AlgorithmExpansionTests` | `test_repetition_codes_preserve_arbitrary_logical_amplitudes` | — | — | V1 |
| `walks.py` | <a href="../zh/manual/algorithms/coined-cycle-walk.html">coined cycle walk</a> | C1 | `AlgorithmExpansionTests` | `test_ansatz_zero_angles_and_walk_one_step` | — | — | V1 |
| `graph_walks.py` | <a href="../zh/manual/algorithms/adjacency-oracle.html">adjacency oracle</a>, <a href="../zh/manual/algorithms/szegedy-walk.html">Szegedy</a>, <a href="../zh/manual/algorithms/mnrs-search.html">MNRS</a> | C1 | `tests/core/test_graph_walks.py:AdjacencyOracleTests` / `SzegedyWalkTests` / `QuantumWalkSearchTests` / `HittingTimeTests` | `test_gate_adjacency_implements_neighbor_table` / `test_stationary_state_is_fixed_point` / `test_search_amplifies_marked_vertex_on_hypercube` / `test_cycle_hitting_times_match_closed_form` | `test_abstract_handles_bind_to_gate_implementations` (dual gate/qram binding of the adjacency oracle) | Johnson graphs require new witnesses (element distinctness) | V3 |

### Approximate continuous algorithms (C2)

| Module | Algorithms | Class | Structural evidence | Numerical evidence | Binding evidence | Gap | Stage |
|---|---|---|---|---|---|---|---|
| `transforms.py` | <a href="../zh/manual/algorithms/qubitization-walk.html">qubitization walk</a>, <a href="../zh/manual/algorithms/qsvt-sequence.html">QSVT sequence</a>, <a href="../zh/manual/algorithms/oblivious-amplification.html">OAA</a> | C2 | `tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` / `test_trotter_only_input_does_not_need_block_encoding` | `tests/core/test_qsvt.py:PhaseSynthesisTests.test_convention_matches_qsvt_sequence` (pointwise agreement 1e-10) | — | — | V1 |
| `qsvt.py` | <a href="../zh/manual/algorithms/qsp-phase-synthesis.html">phase synthesis</a>, <a href="../zh/manual/algorithms/qsvt-matrix-inversion.html">inversion</a>, <a href="../zh/manual/algorithms/eigenstate-filtering.html">filtering</a>, <a href="../zh/manual/algorithms/qsvt-hamiltonian-simulation.html">HS simulation</a>, <a href="../zh/manual/algorithms/fixed-point-search.html">fixed-point search</a> | C2 | `PhaseSynthesisTests.test_input_validation` / `TransformWitnessTests.test_transform_input_validation` | `test_roundtrip_chebyshev` / `test_roundtrip_with_explicit_imaginary_part` / `test_negated_phases_conjugate_polynomial` / `test_matrix_inversion_block` / `test_eigenstate_filter_isolates_eigenvalue` / `test_hamiltonian_simulation_block` / `test_fixed_point_search_phase_properties` / `test_fixed_point_search_circuit_amplifies` | — | Convergence scans missing; the ill-conditioned negative cases (degree ≳ 16 or too-low precision) are partially covered by the parametrized input_validation tests, and a separate batch parametrization of "ill-conditioned inputs must raise" is still to be added | V2 |
| `hamiltonian.py` | <a href="../zh/manual/algorithms/trotter.html">Trotter</a>, <a href="../zh/manual/algorithms/taylor-block-encoding.html">Taylor</a>, <a href="../zh/manual/algorithms/hamiltonian-simulation.html">hamsim</a> | C2 | `AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` / `test_trotter_validates_each_term` | Stage ordering covered indirectly by `test_trotter_protocol_keeps_phase_and_repeat` | — | Trotter order-error-rate scan | V2 |
| `sparse.py` / `block_encoding.py` | <a href="../zh/manual/algorithms/sparse-block-encoding.html">sparse BE</a>, <a href="../zh/manual/algorithms/block-encoding-algebra.html">BE algebra</a> | C2 | `tests/core/test_language.py:BlockEncodingTests` | The same class's `test_unequal_normalization_and_signed_sum` / `test_complex_coefficients` / `test_product_order_and_independent_signal_spaces` | `test_alpha_survives_ir_serialization` | — | V1 |
| `qlss.py` | <a href="../zh/manual/algorithms/cks.html">CKS</a>, <a href="../zh/manual/algorithms/costa-walk.html">Costa</a> | C2 | `tests/core/test_qlss_input_models.py:QLSSInputTests` | `test_signed_sparse_encoding_and_chebyshev` / `test_protocols_consume_different_models_and_scale_kappa` / `test_costa_rhs_reflection_is_independent_of_unitary_extension` / `test_matrix_probe_recovers_scalar_system_norm` | `test_costa_rhs_reflection_is_independent_of_unitary_extension` (independence of the U extension) | End-to-end cross-check against the HHL paper's reference values (could enter the catalog) | V3 |
| `vtaa_cks.py` | <a href="../zh/manual/algorithms/vtaa-cks.html">VTAA-CKS</a> | C2 | `tests/core/test_vtaa_cks.py:VTAAStructureTests` | `GappedPhaseEstimationTests.test_marker_circuit_matches_qsp_response_exactly` / `BandInverseTests.test_band_inverse_step_matches_chebyshev_polynomial` / `VTAAEndToEndTests.test_uniform_spectrum_single_band_is_exact` / `test_variable_time_clock_separates_bands` | The schedule-independence witnesses of `VTAAEndToEndTests` (amplification rounds do not change the conditional solution state) | Band-polynomial precision calibration and transition-band complex-amplitude cross-checks; the amplitude-estimation channel of the VTAA schedule | V3 |
| `ode.py` / `lchs.py` / `cbmd.py` / `schrodingerization.py` / `carleman.py` | <a href="../zh/manual/algorithms/qode-problem.html">QODE problem</a> and the solvers <a href="../zh/manual/algorithms/lchs.html">LCHS</a>, <a href="../zh/manual/algorithms/cbmd.html">CBMD</a>, <a href="../zh/manual/algorithms/schrodingerization.html">Schrödingerization</a>, <a href="../zh/manual/algorithms/carleman.html">Carleman</a> | C2 | `tests/core/test_differential.py:DifferentialStructureTests` | `test_four_methods_keep_input_oracles` / `test_cbmd_plan_records_omitted_terms` / `test_qfvm_uses_raw_flow_and_modular_arithmetic` / `test_local_riemann_updates_only_neighbor_faces_and_tree` | The same class's registration and contract checks | A unified "analytically solvable ODE family" convergence benchmark | V2 |
| `sde.py` | <a href="../zh/manual/algorithms/fokker-planck.html">Fokker–Planck</a> | C2 | `tests/core/test_sde.py:GeneratorDiscretizationTests` / `StatePreparationTests` / `SolverContractTests` | `test_ou_moments_match_closed_form` / `test_stationary_approaches_boltzmann` / `test_column_sums_vanish` | `SolverContractTests.test_qode_problem_accepted_by_lchs` | End-to-end catalog case for SDE and LCHS | V3 |
| `density.py` | <a href="../zh/manual/algorithms/purification.html">purification</a>, <a href="../zh/manual/algorithms/gibbs-state.html">Gibbs</a> | C2 | `tests/core/test_density.py:PurificationTests` | The same class's `test_gate_purification_recovers_rho` / `test_maximally_mixed_purification` / `test_classical_tools`; `GibbsTests.test_diagonal_hamiltonian_gibbs` / `test_non_diagonal_hamiltonian_gibbs` | `PurificationTests.test_abstract_purification_binds` | — (V1 added the `test_error_convergence_decreases` monotone non-increasing scan + the `test_error_bound_uniform_in_beta` β grid, measured distances per step 5.6e-4 / 8.7e-5 / 8.7e-5 and 2.7e-6 / 1.5e-5 / 1.9e-4) | V1 |
| `lowrank.py` | <a href="../zh/manual/algorithms/double-factorization.html">DF</a>/<a href="../zh/manual/algorithms/thc.html">THC</a> → BE | C2 | `tests/core/test_lowrank.py:DiagonalizeSymmetricTests` / `DoubleFactorizationTests` / `ThcTests` | `test_single_rank_block_equals_hamiltonian` / `test_thc_block_equals_hamiltonian` / `test_feeds_qubitization_walk` | Across BE bindings via the `block_column` helper inside `test_matches_pauli_encoding_block` | α tightness witnesses: V1 added `test_alpha_matches_closed_form_eigenvalues` (independent 2×2 closed-form cross-check, places=10) + `test_df_alpha_tighter_than_pauli` (DF α=1.4 ≤ Pauli α=1.5) + `test_thc_alpha_matches_hand_computed_bound` (THC α=1.996 independent hand computation) | V1 |

### Probability-distribution algorithms (C3)

| Module | Algorithms | Class | Structural evidence | Numerical evidence | Binding evidence | Gap | Stage |
|---|---|---|---|---|---|---|---|
| `estimation.py` | <a href="../zh/manual/algorithms/qpe.html">QPE</a>, <a href="../zh/manual/algorithms/qae.html">QAE</a>, <a href="../zh/manual/algorithms/quantum-counting.html">quantum counting</a>, <a href="../zh/manual/algorithms/hadamard-test.html">Hadamard test</a>/<a href="../zh/manual/algorithms/swap-test.html">SWAP test</a> | C3 | `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` | `test_amplitude_estimation_half_probability` / `test_hadamard_test_real_and_imaginary` / `test_swap_test_overlap` | — | The QAE confidence-interval claim is unwitnessed (only the point estimate is verified) | V2 |
| `search.py` | <a href="../zh/manual/algorithms/grover.html">Grover</a>, <a href="../zh/manual/algorithms/amplitude-amplification.html">amplitude amplification</a> | C3 | `AlgorithmExpansionTests` | `test_amplification_of_known_quarter_probability` | — | — | V1 |
| `integration.py` | <a href="../zh/manual/algorithms/heinrich-summation.html">Heinrich summation</a>/<a href="../zh/manual/algorithms/heinrich-integration.html">integration</a> | C3 | `tests/core/test_integration.py:SumPreparationTests` / `QuantumSumTests` / `RateTests` | `test_mean_on_qae_grid_is_exact` / `test_ramp_mean_within_qae_resolution` / `test_quantum_integral_trapezoid_scale` / `test_heinrich_rate_known_values` | `SumPreparationTests.test_qram_binding_matches_gate` / `test_abstract_loader_binds` (three-layer consistency) | — | V1 |
| `gradient.py` | <a href="../zh/manual/algorithms/jordan-gradient.html">Jordan gradient</a> | C3 | `tests/core/test_gradient.py:GradientTests` | `test_linear_function_exact_two_dimensions` / `test_function_phase_oracle_from_mathfunc` | `test_abstract_oracle_binds_to_gate_implementation` | — (V1 added `test_perturbed_linear_concentrates_with_grid_bits`, asserting the failure-probability decay rate q_{m+1} ≤ 0.34·q_m, measured ratios 0.292 / 0.268) | V1 |
| `qpca.py` | <a href="../zh/manual/algorithms/qpca.html">QPCA</a>, <a href="../zh/manual/algorithms/density-matrix-exponentiation.html">density-matrix exponentiation</a> | C3 | `tests/core/test_qpca.py:DensityMatrixExponentiationTests` / `QpcaTests` | `test_pure_state_eigenvalues` / `test_mixed_state_eigenvalue_via_purification` / `test_small_step_first_order_accurate` / `test_error_halves_with_copies` | Covered via the `gate_purification.as_state_preparation()` adaptation path | — | V1 |

### Heuristic/optimization algorithms (C4)

| Module | Algorithms | Class | Structural evidence | Numerical evidence | Binding evidence | Gap | Stage |
|---|---|---|---|---|---|---|---|
| `variational.py` | <a href="../zh/manual/algorithms/qaoa-maxcut.html">QAOA-MaxCut</a>, <a href="../zh/manual/algorithms/vqe.html">VQE</a>, <a href="../zh/manual/algorithms/variational-ansatz.html">ansatz</a> | C4 | `AlgorithmExpansionTests` | `test_qaoa_single_edge_optimal_layer` / `test_vqe_measurement_of_y_eigenstate` / `test_ansatz_zero_angles_and_walk_one_step` | — | The strong witness of "reaching the optimal cut on small instances" | V3 |
| `dqi.py` | <a href="../zh/manual/algorithms/dqi.html">DQI</a> | C4 | `tests/core/test_dqi.py:DqiTests` | `test_planted_instance_beats_random_guessing` / `test_identity_instance_with_weight_two` / `test_dicke_state_weight_and_uniformity` | `test_abstract_decoder_binds_to_witness` (pointwise distribution cross-check against Krawtchouk) | The identity boundary was already fixed with ≥/≤ tolerances | V1 |

### Application layer

| Module | Contents | Class | Structural evidence | Numerical evidence | Binding evidence | Gap | Stage |
|---|---|---|---|---|---|---|---|
| `qfvm.py` / `qham/` | <a href="../zh/manual/algorithms/qfvm.html">QFVM</a>, <a href="../zh/manual/algorithms/qham.html">QHAM</a> | C6 | `tests/core/test_qham_general.py:QhamGeneralTests` | `test_rules_roundtrip_and_unsupported_nonlinearity` / `test_lazy_closure_and_rank_unrank` / `test_quadratic_dimension_and_nontrivial_eta_weight` / `test_multidimensional_coupled_identity` / `test_m1_reduces_to_previous_special_case` | `test_open_qode_inputs_keep_operator_and_initial_oracles` | — | V1 |
| `catalog.py` | 33-case cross-check catalog | C6 | `tests/core/test_workloads.py:WorkloadTests` | `test_every_catalog_case_has_a_closed_description` | `test_alpha_change_requires_regeneration` | New algorithms (qsvt inversion, MNRS search, integration, qpca) not yet registered | V3 |

### Cross-cutting mechanisms

| Mechanism | Implementation | Self-tests |
|---|---|---|
| Witness primitives | `tests/core/witness.py` | `tests/core/test_witness.py:UnitaryWitnessTests` / `UncomputationWitnessTests` / `BindInvariantWitnessTests` / `BlockEqualsWitnessTests` (11 cases in total, including deliberately wrong counterexamples) |
| Core static and unit tests | `tests/core/test_*.py` | See the structural/numerical/binding evidence columns of each module above |

## Maintenance rules

- When adding an algorithm or a new witness, synchronize this table's
  module/algorithms/gap/stage columns, and the witness-primitives section if
  a new primitive is introduced.
- A new algorithm must also add its page under `docs/zh/manual/algorithms/`
  (template in the "Algorithm pages" section of
  `docs/development/writing-docs.md`).
- Table style matches `algorithm-coverage.md`: module names follow the actual
  files under `src/oracq/`; algorithm names reuse the public symbols exported
  by the corresponding module.
- Categories (C1–C6) are filled in per the acceptance criteria of
  `validation-plan.md` §2; when the same algorithm falls under different
  categories for different witnesses, take the dominant one and explain in
  the gap column.
- Evidence-layer locations are written uniformly in the
  `tests/core/<file>.py:<TestClass>` form (with a `.<test_method>` suffix
  when needed), so readers can jump straight to the test class.
- Stages (V1–V4) follow `validation-plan.md` §5; note "—" for anything that
  belongs to no stage.
