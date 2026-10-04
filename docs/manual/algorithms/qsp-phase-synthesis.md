# QSP Phase Synthesis

**English** · <a href="../../zh/manual/algorithms/qsp-phase-synthesis.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.common.qsvt`](../../api/algorithms/common/qsvt.rst) · Stage V2

## Overview

Quantum signal processing (QSP) phase synthesis: given a real-coefficient target polynomial, find the phase sequence $\Phi$ under which [QSVT phase sequence](qsvt-sequence.md) realizes that polynomial. The module exports a pair of mutually inverse routines: {obj}`qsp_response(x, phases) <oracq.algorithms.common.qsvt.qsp_response>` evaluates forward, returning the top-left block $p(x)$ of the phase sequence under the reflection convention ($x \in [-1, 1]$); {obj}`qsp_phases(coeffs, imag=None) <oracq.algorithms.common.qsvt.qsp_phases>` synthesizes in reverse, returning phases in time order (length $d + 1$, where $d$ is the target degree).

By default, synthesis uses complementary-polynomial root finding plus layer stripping from Gilyén et al. 2019 ([arXiv:1806.01838](https://arxiv.org/abs/1806.01838), cited in the module docstring). Real targets that are not saturated at the endpoints (such as a truncated approximation of $1/x$) must rely on a nonzero imaginary-part completion $h$: after synthesizing $P = f + i h$, the real part is extracted with the LCU combination $(U_\Phi + U_{-\Phi})/2$ — $-\Phi$ realizes exactly the conjugate polynomial $\bar{P}$ — and the block encoding is then $f(A/\alpha)$.

Phase synthesis is a replaceable component: `qsp_phases` and the upstream standard-transform family accept an optional `synthesizer` parameter (the {obj}`PhaseSynthesizer <oracq.algorithms.common.qsvt.PhaseSynthesizer>` protocol). The built-in route remains the default and keeps the degree guard; a replacement synthesizer is not subject to the degree guard, but its output must pass the same convention check — the pinned completion mode (with `imag` given) requires pointwise reproduction of $P = f + i h$, while the free completion mode (`imag=None`) pins down only the real part $\mathrm{Re}\,p = f$, leaving the complementary imaginary part to the synthesizer's own choice. All consumers in the library depend only on the real part (the real-part-extraction LCU), so the free completion mode suffices; the YLC target $P$ of fixed-point search is generally complex-valued and goes through the pinned completion mode.

## Interface and input model

```python
qsp_response(x, phases)
qsp_phases(coeffs, imag=None, *, synthesizer=None)
```

API entries: {obj}`qsp_response <oracq.algorithms.common.qsvt.qsp_response>`, {obj}`qsp_phases <oracq.algorithms.common.qsvt.qsp_phases>`, {obj}`PhaseSynthesizer <oracq.algorithms.common.qsvt.PhaseSynthesizer>`

- `x`: a scalar spectral variable ($x \in [-1, 1]$); `phases`: the phase sequence. Returns the complex number $p(x)$.
- `coeffs`: real coefficients of the target polynomial $f$ in ascending powers (constant term first); `imag`: the optional imaginary-part completion $h$, likewise real coefficients in ascending powers.
- `synthesizer`: an optional replacement phase synthesizer implementing the `PhaseSynthesizer` call contract (input: real coefficients in ascending powers plus an optional completion; output: a reflection-convention phase sequence in time order of length $d+1$). With `imag` given, the synthesizer must reproduce $P = f + i h$ exactly; with `imag=None`, only $\mathrm{Re}\,p = f$ is required. The repository ships an optional pyqsp-based adapter `oracq.algorithms.common.qsp_pyqsp.PyqspSynthesizer` (symmetric-QSP Newton method, free completion mode only; requires the `pyqsp` extra, no new hard dependency for the core).

Both entries are purely numerical routines that generate no quantum program and therefore have no input model; their output (a phase tuple) is consumed by {obj}`qsvt_sequence <oracq.algorithms.common.transforms.qsvt_sequence>` and the standard-transform family. Realizability conditions (violations raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`): the parities of both $f$ and $h$ equal $d \bmod 2$; $f^2 + h^2 \le 1$ on $[-1,1]$; endpoint saturation $f(\pm 1)^2 + h(\pm 1)^2 = 1$; the complementary polynomial $R = (1 - f^2 - h^2)/(1 - x^2)$ is non-negative and satisfies the root-multiplicity condition of its spectral decomposition. A $d = 0$ target must be a unit-modulus constant, in which case a single phase is returned.

## Implementation notes

Synthesis pipeline: input validation (real coefficients, parity, upper bound, endpoint saturation, $R$ non-negative) → root finding of $R$ (pure-Python Durand–Kerner iteration, repeated roots clustered by relative distance) → root-structure classification (real roots in pairs, purely imaginary roots in pairs, complex roots in conjugate sign-flipped quadruples) to construct $Q$ with $Q\bar{Q} = R$ → layer stripping recovers the phases in reverse (double-precision complex numbers) → a round-trip self-check with `qsp_response` on a 512-point grid; an error above 1e-5 raises `ValidationError`.

The evaluation of `qsp_response` is a 2×2 matrix recursion: starting from the identity, each step left-multiplies by $W(x)$ (the action that precedes this step's phase in time, except for the first phase) and then by $S(\varphi)$, returning the $(0,0)$ entry of the accumulated matrix. Besides standalone calls, it is also the self-check tool inside the synthesizer: after stripping finishes, both `qsp_phases` and {obj}`fixed_point_search_phases <oracq.algorithms.common.qsvt.fixed_point_search_phases>` use it for the grid round-trip check.

There are two numerical boundaries, which must be distinguished. **Synthesis side**: the built-in route has a degree cap of 40; for degrees within a few tens and a complementary polynomial with well-separated roots, the round-trip error is typically on the order of 1e-9. At higher degrees (≳ 16) or with nearly repeated roots of the complementary polynomial, root finding and stripping are numerically unstable, and the self-check rejects ill-conditioned inputs instead of emitting wrong phases. **Representation side**: the library's coefficient-level pipeline is a monomial (ascending-power) representation, and double-precision grid evaluation distorts overall above roughly degree 40 (the leading coefficient of $T_{41}$ is $2^{40}$, and Horner evaluation suffers catastrophic cancellation at $|x| \approx 1$); this boundary is independent of the synthesis route — a replacement synthesizer crosses the synthesis degree guard, not the representation limit. High-degree targets (such as the degree-585 polynomial of a κ=8, ε=1e-2 inversion) must be constructed and checked in the Chebyshev basis: `qsp_pyqsp.chebyshev_phases_via_pyqsp` accepts Chebyshev coefficients directly. The adapted convention mapping (pyqsp's Wx rotation convention, θ₀ applied last, the target landing in $\mathrm{Im}\,U_{00}$, mapped to reflection-convention phases in time order with the endpoint phases corrected by $d \bmod 4$) is independently validated by the `qsp_response` grid round trip. Measured on the κ=8 case: degree 585, round-trip residual 1.0e-13, approximation error 1.8e-3 ≤ 1e-2 against $c/x$ on $|x| \ge 1/8$ (`tools/expressiveness/t1_qsp_phases.py`; execution-level direction error 5.6e-13 on a 2×2 block encoding, see `t2_qsvt_inversion.py`); the library-level coefficient entry `qsvt_matrix_inversion` is still explicitly rejected by the round-trip guard under the same parameters, with no silent behavior. Applicability boundary: only real-coefficient targets are accepted — complex targets must be split into a real part plus an imaginary-part completion; the time-order convention of the returned phases is aligned step by step with the circuit convention of `qsvt_sequence` (see [QSVT phase sequence](qsvt-sequence.md)). The phases of the standard-transform family ([matrix inversion](qsvt-matrix-inversion.md), [eigenstate filtering](eigenstate-filtering.md), [Hamiltonian simulation](qsvt-hamiltonian-simulation.md), [fixed-point search](fixed-point-search.md)) are synthesized through this pipeline by default.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../../development/validation-plan.md` §2). Three layers of evidence:

- Structure: `tests/core/test_qsvt.py:PhaseSynthesisTests.test_input_validation` — five classes of negative cases (parity violation, upper-bound violation, unsaturated endpoints, non-real coefficients, imaginary-part parity) all raise `ValidationError` at synthesis time.
- Numerical: `test_roundtrip_chebyshev` (Chebyshev polynomials $T_n$, $n = 1..6$; the 51-point grid round trip is consistent with delta = 1e-8, phase length $n + 1$); `test_roundtrip_with_explicit_imaginary_part` ($P = (0.5 + i\sqrt{0.75})x$, an exactly realizable case with $|P|^2 = x^2$ and $Q = 1$, delta = 1e-9); `test_negated_phases_conjugate_polynomial` ($-\Phi$ realizes $\bar{P}$, the property the real-part extraction relies on, delta = 1e-12); `test_convention_matches_qsvt_sequence` (`qsp_response` agrees pointwise with the circuit's zero-signal block, delta = 1e-10).
- Replaceability: `tests/core/test_qsvt.py:PhaseSynthesizerTests` — acceptance and self-check rejection of injected synthesizers (wrong phases, wrong length, non-finite values), the exactness of the pinned completion mode, the degree guard applying only to the built-in route, upstream consumers' outputs being bit-identical after a delegating synthesizer is injected (`dumps` cross-check), and the explicit rejection of κ=8/ε=1e-2 under the built-in route. With the `pyqsp` extra installed, `tests/core/test_qsp_pyqsp.py` covers: the convention mapping agreeing pointwise with an independently evaluated target at several parities and degrees (delta 1e-8), rejection in pinned mode, an end-to-end block-encoding cross-check for κ=3/error=0.05 (degree 51, beyond the built-in guard), the explicit failure of the κ=8 library-level entry, and κ=8 phase-level target attainment (residual ≤ 1e-2 and approximation error ≤ 1e-2 on $|x| \ge 1/8$).
- Binding: purely numerical routines, no binding witness.

## Known gaps and planned stages

Consistent with the `qsvt.py` row of the validation coverage matrix: convergence scans are missing (batch curves of the error versus degree are not automated, pending the stage V2 convergence-scan framework); the ill-conditioned negative cases (degree ≳ 16 or too-low precision) are partially covered by the parametrized `test_input_validation`, and a separate batch parametrization of "ill-conditioned inputs must raise" is still to be added.

## Related links

- Source: `src/oracq/algorithms/common/qsvt.py`; pyqsp adapter: `src/oracq/algorithms/common/qsp_pyqsp.py`
- API reference: [QSVT standard transforms](../../api/algorithms/common/qsvt.rst)
- Same-family pages: [QSVT phase sequence](qsvt-sequence.md) (the phase consumer), [QSVT matrix inversion](qsvt-matrix-inversion.md) (a typical use of imaginary-part completion)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_hamiltonian.py` (real-backend execution, no mock substitutes). Experiment design: with the target polynomials taken as the Chebyshev $T_1,\dots,T_6$ and the exactly realizable case $P = (0.5 + i\sqrt{0.75})x$, the phases synthesized by `qsp_phases` are evaluated on a uniform 401-point grid with a 2×2 matrix-recursion response **implemented independently** in numpy following this page's $W(x)$ / $S(\varphi)$ convention (the in-library `qsp_response` is not called; the two are independent implementations of each other), and cross-checked pointwise against an independent Horner evaluation of the target polynomial.

| Case | Scale | Backend paths | Metric | Value |
|---|---|---|---|---|
| `qsp-phase-synthesis-roundtrip` | 7 targets, 401-point grid | numpy independent response (the synthesizer generates no quantum program) | max_error | 3.8e-15 |
| same as above (Chebyshev $T_1$–$T_6$) | — | — | per-target error | 0 – 3.8e-15 |
| same as above (explicit-imaginary target) | — | — | error | 2.2e-16 |

For the circuit-level consumer cross-check of the synthesized phases (four backend paths, random phases, matrix polynomial block), see the Numerical validation section of [QSVT phase sequence](qsvt-sequence.md).

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_hamiltonian.py
```

Artifacts: `out/verification/hamiltonian.json`.
