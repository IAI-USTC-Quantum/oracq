# Acceptance record for this version

**English** · <a href="../../zh/development/validation.html">简体中文</a>

## Algorithm research workflow improvements (2026-09-21)

Single-input adaptation, solver naming compatibility, partially bound reports
and resource origins, open-cost analysis, and the QRAM/arithmetic
implementation comparison are complete. Core and schema total 391 tests and
278 subtests; real integration 41; six groups of related numerical validation
with 193 cases; and 9 new implementation-comparison cases — all pass. Static,
type, and Sphinx checks pass. For detailed scope, environment, and artifact
paths see the [workflow implementation and acceptance
record](research-workflow-plan.md).

## 0.8 release acceptance snapshot

This page records the scope of the 0.8 repository reorganization, algorithm
expansion, and documentation build. Final test counts and artifact summaries
are kept in `validation.json` in the same directory.

Validation falls into four kinds: regression of the existing core and schema;
independent mathematical witnesses for the new algorithms;
complex-amplitude cross-checks on real backends; and Sphinx HTML, tutorial,
and package checks.

The mathematical witnesses cover the positive-sign QFT, modular addition over
all three-bit inputs, BV secret recovery, the Simon null space, amplitude
amplification and estimation, Hadamard/Swap tests, one-sided QAOA, Pauli
measurements, the periodic walk, modular multiplication for order finding,
and repetition-code recovery. The existing validation boundaries of
QLSS/QODE/QHAM were not thereby extended into full numerical certification.

Documentation builds use warnings-as-errors. The API is generated from the
canonical source paths; historical documents are kept out of the build.
Legacy import paths have their own compatibility witnesses.

## Acceptance results

- Core and schema: 144 tests and 135 subtests pass.
- Real backend: 29 tests pass, of which the algorithm gallery covers 22 cases;
  the maximum amplitude difference across the three execution results is
  `3.89e-16`.
- Documentation: 56 API module pages and 7 tutorials; the HTML build has zero
  warnings, and 10 tutorial tests and 2 search tests pass.
- Browser: home, tutorials, API navigation, English and Chinese search,
  formulas, and Mermaid diagrams have all been checked.
- Packaging: sdist and wheel builds pass; canonical imports, legacy-import
  compatibility, and basic algorithm execution in an isolated environment
  pass.

Machine-readable results are in [validation.json](validation.json).

## V1 validation coverage progress (2026-09-10)

The progress below landed as stage V1 of [`validation-plan.md`](validation-plan.md)
§5; the previous section remains the acceptance snapshot as of the 0.8
release, while this section reflects iteration after 0.8 without rewriting
history. For the complete "where each algorithm's evidence lives" view, see
[`validation-coverage.md`](validation-coverage.md).

- Invariant test library: `tests/core/witness.py` provides the four assertion
  primitives `assert_unitary` / `assert_uncomputation` / `assert_bind_invariant`
  / `assert_block_equals`; the self-tests live in `tests/core/test_witness.py`
  (4 test classes, 11 cases, covering both the correct and the deliberately
  wrong paths).
- In-flight module witnesses were filled in as the matrix requires:
  - density: Gibbs error scan (`GibbsTests.test_error_convergence_decreases`,
    measured 5.6e-4 / 8.7e-5 / 8.7e-5, monotone non-increasing + each step
    ≤ error) and the β grid (`test_error_bound_uniform_in_beta`).
  - gradient: failure-probability decay rate
    (`GradientTests.test_perturbed_linear_concentrates_with_grid_bits` asserts
    q_{m+1} ≤ 0.34·q_m, measured ratios 0.292 / 0.268).
  - lowrank: independent 2×2 closed-form cross-check places=10
    (`DoubleFactorizationTests.test_alpha_matches_closed_form_eigenvalues`),
    DF α=1.4 ≤ Pauli α=1.5 (`test_df_alpha_tighter_than_pauli`), THC α=1.996
    independent hand computation
    (`ThcTests.test_thc_alpha_matches_hand_computed_bound`).
  - integration: closed-form mean cross-check + heinrich_rate validation
    (`SumPreparationTests` / `QuantumSumTests` / `RateTests`).
  - qpca: eigenvalue readout peak cross-check + Δt first-order error rate
    (`DensityMatrixExponentiationTests` / `QpcaTests`).
- Total core tests: 237 → 253 (`python -m unittest discover -s tests/core`
  fully green). `src/` was not changed in this round.
