# Eigenstate Filtering

**English** · <a href="../../zh/manual/algorithms/eigenstate-filtering.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.common.qsvt`](../../api/algorithms/common/qsvt.rst) · Stage V2

## Overview

Eigenstate filtering: construct a "spike" polynomial block encoding from a block encoding of a matrix $A$ — components whose spectral variables fall inside the passband around the filter center are preserved ($|f|$ close to 1), while components outside the passband are suppressed below $1/T_d(r)$. The target is a Lin–Tong-type filter polynomial (the construction basis cited in the source docstring):

$$
f(x) = \frac{T_d\!\bigl(g(x^2)\bigr)}{T_d(r)}, \qquad
g(y) = \frac{2(y - \Delta^2)}{1 - \Delta^2} - 1, \qquad
r = \frac{1 + \Delta^2}{1 - \Delta^2},
$$

where $\Delta$ is the filter width and $d$ the Chebyshev degree. $f$ saturates at $x = 0$ ($|f(0)| = 1$), and $|f(x)| \le 1/T_d(r)$ on $|x| \ge \Delta$, with the suppression rate decaying exponentially in $d \cdot \operatorname{arcosh}(r)$. The phase-synthesis framework is the same as [QSP phase synthesis](qsp-phase-synthesis.md) (Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)).

## Interface and input model

```python
eigenstate_filter(a, gap, degree, *, center=0.0)
```

API entry: {obj}`eigenstate_filter <oracq.algorithms.common.qsvt.eigenstate_filter>`

- `a`: a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`, the block encoding of the filtered matrix (input model BE; the spectral variable is $x = \lambda/\alpha$, which for the Hermitian case is just the normalized eigenvalue).
- `gap`: the filter width $\Delta$, which must lie in $(0, 1)$.
- `degree`: the Chebyshev degree $d$, which must be a positive integer (floats are rejected).
- `center`: the filter center, which must lie in $(-1, 1)$; when nonzero, the spectrum is shifted first (see below).

Returns a `BlockEncoding` whose zero-signal block is suppressed below `suppression` outside the passband. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"eigenstate_filter"` |
| `be_alpha` | normalization of the output BE (the real-part-extraction LCU in equal halves; 1.0) |
| `gap` / `filter_degree` / `center` | echo of the call parameters |
| `qsp_degree` | the synthesized polynomial degree $2d$ |
| `suppression` | the out-of-passband upper bound $1/T_d(r)$ |

## Implementation notes

The synthesized degree is $2d$ ($f$ is a degree-$d$ Chebyshev combination of $x^2$, monomial degree $2d$), subject to the synthesis cap of 40 (i.e. $d \le 20$); exceeding it raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`. When `center` is nonzero, the spectrum is first shifted by the BE linear combination {obj}`linear_combination(1.0, a, -center, identity(a.width)) <oracq.algorithms.input_model.operators.linear_combination>`: a block encoding of $A - \text{center}\cdot I$ is constructed (the normalization becomes $\alpha + |\text{center}|$), the filtering acts on the shifted spectral variable $x' = (\lambda - \text{center})/(\alpha + |\text{center}|)$, and `gap` is measured in that variable.

$f$ saturates at $x = 0$ but $|f(\pm 1)| = 1/T_d(r) < 1$ leaves the endpoints unsaturated, so phases can be synthesized only with the imaginary-part completion $h = \text{sat}\cdot x^2$ ($\text{sat} = \sqrt{1 - f(1)^2}$); after synthesis the real part is extracted via $(U_\Phi + U_{-\Phi})/2$ (the mechanism is described in [QSP phase synthesis](qsp-phase-synthesis.md)). Applicability boundary: the input must be a block encoding; the shifted spectral values should lie in $[-1, 1]$ — whatever exceeds it loses the upper-bound guarantee.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../../development/validation-plan.md` §2): the acting operator must approximate the target continuous function within tolerance. Three layers of evidence:

- Structure: `tests/core/test_qsvt.py:TransformWitnessTests.test_transform_input_validation` (negative cases such as an out-of-range gap and a non-positive-integer degree) and `PhaseSynthesisTests.test_input_validation` (the input validation of the underlying synthesis).
- Numerical: `tests/core/test_qsvt.py:TransformWitnessTests.test_eigenstate_filter_isolates_eigenvalue` — a 2×2 diagonal BE (eigenvalues 0.05 and 0.5, $\alpha = 0.5$), gap = 0.2, $d = 8$, center = 0.1: `suppression` < 0.08; the passband eigenvalue after shifting is $\approx -0.083$ ($|x'| < \Delta$) with zero-signal amplitude > 0.7; the stopband eigenvalue after shifting is $\approx 0.667$, with amplitude < 0.1.
- Binding: this algorithm has no independent binding witness (the input already requires a concrete BE).

## Known gaps and planned stages

Consistent with the `qsvt.py` row of the validation coverage matrix: convergence scans are missing (the exponential-decay curve of the suppression rate versus degree is not automated, pending the stage V2 convergence-scan framework); the ill-conditioned negative cases (degree ≳ 16 or too-low precision) are partially covered by the parametrized input_validation, and a separate batch parametrization of "ill-conditioned inputs must raise" is still to be added.

## Related links

- Source: `src/oracq/algorithms/common/qsvt.py`
- API reference: [QSVT standard transforms](../../api/algorithms/common/qsvt.rst)
- Same-family pages: [QSP phase synthesis](qsp-phase-synthesis.md), [QSVT matrix inversion](qsvt-matrix-inversion.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_hamiltonian.py` (real-backend execution, no mock substitutes). Experiment design: first, a 1-qubit diagonal BE with center shift (eigenvalues 0.05 and 0.5, $\alpha = 0.5$; `gap = 0.2`, `d = 8`, `center = 0.1`; shifted spectral variables $x' \approx -0.083$ in the passband and $0.667$ in the stopband); second, a 2-qubit diagonal BE (eigenvalues 0.02, −0.03, 0.4, 0.55; `gap = 0.2`, `d = 6`, two passbands and two stopbands). On the three paths reference / rir-pysparq / originir-ext, the diagonal amplitudes of the zero-signal block are read out and cross-checked pointwise against the independent analytic evaluation of the Lin–Tong filter polynomial $f(x) = T_d(g(x^2))/T_d(r)$ (`math.cosh`/`acosh`, not through the in-library helpers).

| Case | Scale | Backend paths | Metric | Value |
|---|---|---|---|---|
| `eigenstate-filter-centered-1q` | $2d = 16$, center=0.1 | reference, rir-pysparq, originir-ext | max_error | 1.32e-7 |
| same as above | — | — | passband amplitude / stopband amplitude (suppression attribute 0.0779) | 0.748 / 0.0235 |
| `eigenstate-filter-2q` | $2d = 12$, four eigenvalues | reference, rir-pysparq, originir-ext | max_error | 2.17e-10 |
| same as above | — | — | passband minimum / stopband maximum (suppression attribute 0.174) | 0.914 / 0.174 |

In both cases the stopband amplitudes do not exceed the `suppression` attribute bound; the 1.3e-7 deviation of the centered case is numerical noise of layer stripping at degree $2d = 16$, within the module's self-check tolerance (1e-5). Measured numerical boundary: $d = 10$ ($2d = 20$) is rejected by the phase-synthesis self-check, and the synthesizable region $d \le 8$ passes entirely.

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_hamiltonian.py
```

Artifacts: `out/verification/hamiltonian.json`.
