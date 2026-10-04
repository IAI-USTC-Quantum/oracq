# VTAA-CKS QLSS

**English** · <a href="../../zh/manual/algorithms/vtaa-cks.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.qlss.vtaa_cks`](../../api/algorithms/qlss/vtaa_cks.rst) · Stage V3

## Overview

The variable-time amplitude amplification (VTAA) linear-system solver of section 5 of Childs–Kothari–Somma (arXiv:1511.02306, SIAM J. Comput. 2017): a variable-time layer stacked on the $\kappa^2$ route of the [CKS section 4 baseline solver](cks.md), reducing the query complexity to $O(d\,\kappa\,\mathrm{polylog}(d\kappa/\epsilon))$. As of the time of implementation, no public implementation of VTAA or CKS §5 existed anywhere (Qiskit/PennyLane/Qrisp and the like offer only fixed-time amplitude amplification or classical cost models), making this module the first circuit-level implementation of this algorithmic route.

The variable-time layer is assembled following the paper's structure:

- **Clock bits $C_j$** (Lemma 22 GPE): at each step, a decision polynomial is applied to the qubitization walk of the block encoding, coherently writing "whether the eigenvalue is large enough to be inverted in this band" into the clock bit; the amplitude of a 1 verdict is $|P(\lambda/\alpha)|$.
- **Per-band inverse LCU** (Lemma 23 $W(\lambda,\delta)$): controlled on $C_j=1$, the band's truncated Chebyshev inverse polynomial is applied; the rotation of Eq. (98) uniformly compresses each band's success amplitude to $1/\alpha_{\max}$.
- **Clock-aware nested amplification**: Ambainis's (arXiv:1010.4458) VTAA cascade, with the operator form taken from Low–Su (arXiv:2410.18178) Eqs. (47)–(53); correctness does not depend on the amplification schedule, which only affects the success rate.
- **$A'$ uncomputation** (Eqs. (99)–(112)): inverting an $A'$ that replaces $W_j$ with a pure flag flip erases the GPE clock and garbage, leaving the solution state on target, with success condition `signal == 0`.

The GPE verdict follows the deterministic route of Low–Su Prop 23 (the same $O((\alpha/\theta_j)\log(1/\epsilon))$ query order as CKS Lemma 22's PEA + majority vote): it reuses the verified Yoder–Low–Chuang fixed-point polynomial synthesis of `qsvt.fixed_point_search_phases` and applies {obj}`qsvt_sequence <oracq.algorithms.common.transforms.qsvt_sequence>` to the walk. The fire side ($|x| \geq \theta_j$) has a hard bound $\epsilon$; the transition band at the low end of the spectrum (a band CKS does not commit to) behaves deterministically and can be computed exactly, but when neither side falls in the verified region only mixed-semantics guarantees hold.

## Interface and input model

```python
make_vtaa_cks_qlss(config=None)   # returns QLSSProtocol (input_model="sparse")
vtaa_cks(system, config=None)     # low-level kernel, takes a SparseSystem
VTAAConfig(order=2, terms=None, clock_steps=None, marker_epsilon=0.02,
           degree_cap=40, rounds=None)
gapped_phase_estimation(a, threshold, x_edge, *, epsilon=0.02, degree_cap=40)
band_inverse_step(a, coefficients, alpha_max)
tunable_rounds(stage_amplitudes, thresholds=None)  # Low–Su style (52)–(53) schedule
```

API entries: {obj}`make_vtaa_cks_qlss <oracq.algorithms.qlss.vtaa_cks.make_vtaa_cks_qlss>`, {obj}`vtaa_cks <oracq.algorithms.qlss.vtaa_cks.vtaa_cks>`, {obj}`VTAAConfig <oracq.algorithms.qlss.vtaa_cks.VTAAConfig>`, {obj}`gapped_phase_estimation <oracq.algorithms.qlss.vtaa_cks.gapped_phase_estimation>`

The problem input is the same {obj}`LinearSystem(sparse=...) <oracq.algorithms.qlss.qlss.LinearSystem>` as for {obj}`make_cks_qlss <oracq.algorithms.qlss.qlss.make_cks_qlss>` (with the Hermitian and non-negative diagonal declarations). The default `clock_steps` is derived from the declared physical condition number $\kappa_{\mathrm{phys}} = $ `norm_upper / sigma_min_lower` ($\lceil\log_2\kappa\rceil+1$), requiring $2^{m-1} \geq \kappa_{\mathrm{phys}}$ to cover the finest band. `rounds` holds the amplification rounds per stage (all zero by default, i.e. the pure variable-time layer plus post-selection); production deployments should choose them via the amplitude estimation of Ambainis Algorithm 2 or the deterministic schedule of {obj}`tunable_rounds <oracq.algorithms.qlss.vtaa_cks.tunable_rounds>`.

The encoding normalization must leave slack relative to the spectral upper bound (`norm_upper / alpha < 1`): the decision geometry of the walk phase $\arccos(\lambda/\alpha)$ is determined by it, and an over-tight encoding (e.g. a diagonal spectrum paired with `entry_bound = max|A|`) is rejected — widen `entry_bound` in that case.

Kernel output attributes: `algorithm = "vtaa_cks"`, `clock_steps`, `fire_thresholds`, `band_orders`, `band_lcu_normalizations`, `alpha_max`, `marker_degrees`, `marker_epsilon`, `rounds`, `implementation_scope` (citing the three paper threads), `kernel_status = "prototype; band polynomial accuracy and VTAA schedule pending"`, `success_condition = "signal == 0"`. The protocol layer reuses the physical-subspace selection and independent norm probe of `QLSSProtocol.solve()`.

## Implementation notes

Generation chain: {obj}`real_symmetric_sparse_encoding <oracq.algorithms.input_model.sparse.real_symmetric_sparse_encoding>` ($\alpha = s\cdot a_{\max}$) → per band {obj}`gapped_phase_estimation <oracq.algorithms.qlss.vtaa_cks.gapped_phase_estimation>` (`qsvt_sequence` verdict + the clock bit flipped controlled on `signal==0`, with P_j garbage kept as in the paper) and {obj}`band_inverse_step <oracq.algorithms.qlss.vtaa_cks.band_inverse_step>` ({obj}`chebyshev_block <oracq.algorithms.input_model.sparse.chebyshev_block>` odd-power LCU + the Eq. (98) uniformizing rotation) → `vtaa_variable_step` (prefix-all-0 controlled GPE, $C_j=1$ controlled $W_j$) nested level by level into `vtaa_prefix` / `vtaa_amplified_stage` (the reflection $R_f$ flips the phase of the stopped∧failed branch, and $R_s$ is built from the prefix inverse + the all-zero reflection) → after the top level invokes the amplification chain, everything is erased by the inverse of `vtaa_uncompute_step` (GPE replay + pure flag flip). Module calls and {obj}`Repeat <oracq.infrastructure.ir.Repeat>` are all kept symbolic, and the signal registers of each verdict and inversion step are allocated independently.

Two deviations from the paper (both recorded in `implementation_scope`): the GPE replaces PEA + majority vote with a deterministic QSP verdict (Low–Su Prop 23, a modern implementation of the same lemma, avoiding a probabilistic verdict distribution and making reference-simulation cost linear in the polynomial degree); the order schedule of the band inverse polynomials is a prototype geometric ramp-up ($2^{j-1}$), whose approximation accuracy for $1/x$ has not been theoretically calibrated.

## Validation approach

Category C2 (approximately continuous semantics). Consistent with the `vtaa_cks.py` row of the validation coverage matrix:

- Structure: `tests/core/test_vtaa_cks.py:VTAAStructureTests` (protocol output attributes, serialization round-trip, Repeat preservation of amplification rounds, symbolic resource estimation, configuration rejection paths).
- Numerical: `GappedPhaseEstimationTests.test_marker_circuit_matches_qsp_response_exactly` (the verdict circuit and the QSP response agree pointwise on fixed-point-grid eigenvalues, places = 11); `BandInverseTests.test_band_inverse_step_matches_chebyshev_polynomial` (the success branch of the band inverse LCU cross-checked against the Chebyshev polynomial, places = 11); `VTAAEndToEndTests.test_uniform_spectrum_single_band_is_exact` (a single band carries the $|b\rangle$ direction exactly end-to-end); `test_variable_time_clock_separates_bands` (a two-band variable-time structure: the per-eigenvalue path amplitudes are given exactly by the verdict response, the end-to-end amplitudes fall in the normalized coupling interval, and the direct-inverse component of band1 dominates).
- Binding: `test_amplification_schedule_preserves_conditional_solution` (the same conditional solution state is unchanged when the amplification rounds change — a witness of schedule independence).

## Known gaps and planned stages

The accuracy of the band inverse polynomials for $1/x$ is uncalibrated (stage V3): the current coefficients reuse the {obj}`CKSConfig <oracq.algorithms.qlss.qlss.CKSConfig>` closed form with a geometric order ramp-up, and the correspondence with the paper's per-band $\widetilde O(2^j)$ degree schedule remains to be verified; the end-to-end cross-check of transition-band mixed semantics reaches only amplitude intervals. The VTAA amplification schedule is not connected to an amplitude-estimation channel (`tunable_rounds` already gives the Low–Su formulas; the per-stage norm-estimation circuits are missing). The degree cap of 40 (inherited from the `qsvt.py` synthesis cap) limits the range of declarable $\kappa$. The end-to-end cross-check against HHL paper reference values is the same gap as the `qlss.py` row.

## Related links

- Source: `src/oracq/algorithms/qlss/vtaa_cks.py`
- Papers: [CKS arXiv:1511.02306](https://arxiv.org/abs/1511.02306) §5, [Ambainis arXiv:1010.4458](https://arxiv.org/abs/1010.4458), [Low–Su arXiv:2410.18178](https://arxiv.org/abs/2410.18178)
- Related pages: [CKS Chebyshev solver](cks.md), [Costa walk solver](costa-walk.md), [Sparse Matrix Block Encoding](sparse-block-encoding.md), [Fixed-Point Search](fixed-point-search.md)
- API reference: [VTAA-CKS variable-time linear-system solver](../../api/algorithms/qlss/vtaa_cks.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)
