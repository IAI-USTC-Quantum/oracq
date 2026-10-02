# QSVT Matrix Inversion

**English** · <a href="../../../zh/manual/algorithms/qsvt-matrix-inversion.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.common.qsvt`](../../api/algorithms/common/qsvt.rst) · Stage V2

## Overview

Given a block encoding of a Hermitian matrix $A$, construct a block encoding that approximates $A^{-1}$. The implementation follows the quantum singular value transformation framework of Gilyén et al. (Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)): a polynomial transform is applied to the singular values $x$, with the target polynomial being the odd-extension inversion polynomial

$$
J_b(x) = \frac{1-(1-x^2)^b}{x}, \qquad f(x) = c \cdot J_b(x),
$$

where the odd-power coefficients $(-1)^m \binom{b}{m+1}$ are given analytically and the scaling constant $c$ is chosen so that $\lVert f \rVert_\infty \le 1/3$. On $|x| \ge 1/\kappa$, $f(x)$ approximates $c/x$ with relative error at most `error`; the truncation parameter $b$ is determined by the condition number $\kappa$ and `error`, and the polynomial degree is $d = 2b - 1$.

## Interface and input model

```python
qsvt_matrix_inversion(a, kappa, *, error=0.05)
```

API entry: {obj}`qsvt_matrix_inversion <oracq.algorithms.common.qsvt.qsvt_matrix_inversion>`

- `a`: a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`, the block encoding of the matrix to invert (input model BE).
- `kappa`: the condition number $\kappa$, which must be a finite number not less than 1.
- `error`: the relative approximation error, which must lie in $(0, 1)$.

Returns a `BlockEncoding` whose zero-signal block is approximately `inverse_scale · A⁻¹`. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"qsvt_matrix_inversion"` |
| `be_alpha` | normalization of the output BE (1.0 for this algorithm) |
| `kappa` / `error` | echo of the call parameters |
| `qsp_degree` | the synthesized polynomial degree $d = 2b - 1$ |
| `inverse_scale` | the output scaling $c \cdot \alpha$ ($\alpha$ is the input BE normalization) |

## Implementation notes

The degree is derived in closed form from $\kappa$ and `error`: $b = \max\!\bigl(1,\ \lceil \log(1/\varepsilon) / -\log(1 - 1/\kappa^2) \rceil\bigr)$, and $b = 1$ when $\kappa = 1$; if the required degree exceeds the synthesis cap of 40 and no replacement synthesizer is provided, a {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` is raised, and the caller must relax `error` or plug in an external synthesizer through the `synthesizer` parameter (such as `qsp_pyqsp.PyqspSynthesizer`). Note that crossing the degree guard is not the same as crossing the representation limit: the coefficient-level pipeline distorts above roughly degree 40 (see the numerical boundaries in [QSP phase synthesis](qsp-phase-synthesis.md)), and high-degree targets should be constructed in the Chebyshev basis.

The target polynomial $f$ is not saturated at the endpoints ($|f(\pm 1)| < 1$), so phases can be synthesized only with a nonzero imaginary-part completion $h$ (see the realizability conditions in the module docstring). After synthesis, the real part is extracted with the LCU combination $(U_\Phi + U_{-\Phi})/2$ — $-Φ$ realizes exactly the conjugate polynomial $\bar{P}$ — and the resulting block encoding is $f(A/\alpha)$. Phase synthesis uses complementary-polynomial root finding plus layer stripping; after each synthesis a round-trip self-check runs through {obj}`qsp_response <oracq.algorithms.common.qsvt.qsp_response>`, and ill-conditioned inputs are rejected outright.

Applicability boundary: the input must be a block encoding, not a sparse oracle or QRAM; a matrix not yet block-encoded must first be adapted through `block_encoding.py` / `sparse.py` / `lowrank.py` and the like. The synthesis degree of the default route is subject to the cap of 40, so combinations of large $\kappa$ and small `error` fail at generation time rather than silently degrading; with a replacement synthesizer plugged in, the coefficient-level representation limit still applies, and out-of-range inputs are explicitly rejected by the round-trip guard (the degree-585 target at κ=8 and error=1e-2 is one example; for the phase-level route that meets the target, see `tools/expressiveness/t1_qsp_phases.py`).

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2): the acting operator must approximate the target continuous function within tolerance. Three layers of evidence:

- Structure: `tests/core/test_qsvt.py:PhaseSynthesisTests.test_input_validation` (parity violation, upper-bound violation, unsaturated endpoints, non-real coefficients, imaginary-part parity) and `TransformWitnessTests.test_transform_input_validation` (κ < 1, out-of-range error, excessive required degree, etc.).
- Numerical: `TransformWitnessTests.test_matrix_inversion_block` — for a 2×2 matrix (eigenvalues 0.4 and 0.8, κ = 2, error = 0.15), the zero-signal block amplitudes are read out column by column on the reference simulator and cross-checked against `inverse_scale` times the exact inverse $A^{-1} = \begin{pmatrix} 1.875 & 0.625 \\ 0.625 & 1.875 \end{pmatrix}$, with a relative tolerance of 35% of the expected value.
- Binding: this algorithm has no independent binding witness (the input already requires a concrete BE).

The underlying sequence convention is pinned by `PhaseSynthesisTests.test_convention_matches_qsvt_sequence` (with random phases, the circuit's zero-signal block agrees pointwise with `qsp_response`, delta = 1e-10); the conjugation property that the real-part extraction relies on is witnessed by `test_negated_phases_conjugate_polynomial` (delta = 1e-12).

## Known gaps and planned stages

Convergence scans are missing: batch curves of the error versus degree / $\kappa$ are not automated (to be wired in once the stage V2 convergence-scan framework lands). The ill-conditioned negative cases (degree ≳ 16 or too-low precision) are partially covered by `test_input_validation` / `test_transform_input_validation`; a separate batch parametrized test of "ill-conditioned inputs must raise" is still to be added.

## Related links

- Source: `src/oracq/algorithms/common/qsvt.py`
- Same-family pages: [eigenstate filtering](eigenstate-filtering.md), [QSP phase synthesis](qsp-phase-synthesis.md), [QSVT Hamiltonian simulation](qsvt-hamiltonian-simulation.md), [block encoding algebra](block-encoding-algebra.md) (BE composition entry points)
- API reference: [QSVT standard transforms](../../api/algorithms/common/qsvt.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_hamiltonian.py` (real-backend execution, no mock substitutes). Experiment design: three groups of matrices — 2×2 (eigenvalues 0.4/0.8, $\kappa = 2$, `error` ∈ {0.15, 0.10}, corresponding degrees $d = 13, 17$), 2×2 (0.2/0.6, $\kappa = 3$, `error = 0.4`, $d = 15$), and a 4×4 diagonal (0.3/0.4/0.6/0.9, $\kappa = 3$, `error = 0.4`, $d = 15$); on the reference / rir-pysparq / originir-ext paths, the full zero-signal block is read out and divided by `inverse_scale`, then cross-checked against the exact inverse from `numpy.linalg.inv`, reporting the relative spectral error, the condition number of the recovered inverse against the theoretical interval $\kappa\cdot(1\pm e)/(1\mp e)$ (each singular value carries a relative deviation ≤ `error`), and the success probability of each column.

| Case | Scale | Backend paths | Metric | Value |
|---|---|---|---|---|
| `qsvt-inversion-kappa2-error0.15` | 2×2, $d = 13$ | reference, rir-pysparq | relative spectral error | 0.1335 (≤ 0.15) |
| same as above | — | — | recovered $\kappa$ vs exact (theoretical interval) | 1.733 vs 2.0 ([1.478, 2.706]) |
| `qsvt-inversion-kappa2-error0.1` | 2×2, $d = 17$ | reference, rir-pysparq | relative spectral error | 0.0751 (≤ 0.10) |
| same as above | — | — | recovered $\kappa$ vs exact | 1.850 vs 2.0 ([1.636, 2.444]) |
| `qsvt-inversion-kappa3-1q` | 2×2, $d = 15$ | reference, rir-pysparq, originir-ext | relative spectral error | 0.3897 (≤ 0.40) |
| same as above | — | — | recovered $\kappa$ vs exact | 1.831 vs 3.0 ([1.286, 7.0]) |
| `qsvt-inversion-kappa3-2q` | 4×4, $d = 15$ | reference, rir-pysparq, originir-ext | relative spectral error | 0.3897 (≤ 0.40) |
| same as above | — | — | recovered $\kappa$ vs exact; per-column success probability | 1.863 vs 3.0; 0.031–0.109 |

The relative spectral errors all fall within the requested `error`, the recovered condition numbers all fall within the theoretical intervals, and they improve as `error` tightens (0.1335 → 0.0751). Measured numerical boundary: the inversion polynomial is rejected by the phase-synthesis self-check at $d \ge 19$ (e.g. $\kappa = 3$ with `error` ≤ 0.3) (`ValidationError`, consistent with the documented promise of "rejecting ill-conditioned inputs rather than silently degrading"); the synthesizable region $d \le 17$ passes entirely.

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_hamiltonian.py
```

Artifacts: `out/verification/hamiltonian.json`.
