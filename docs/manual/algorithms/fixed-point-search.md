# Fixed-Point Search

**English** · <a href="../../../zh/manual/algorithms/fixed-point-search.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.common.qsvt`](../../api/algorithms/common/qsvt.rst) · Stage V2

## Overview

Fixed-point amplitude amplification: a search amplification that requires no advance knowledge of the initial amplitude and never overshoots from too many iterations. The construction follows the closed-form complementary polynomial of Yoder–Low–Chuang (cited in the source docstring): let $L$ be the number of BE calls (it must be odd) and $c = T_{1/L}(1/\delta)$; the complementary polynomial is $Q(x) = \delta c\, R(c^2(1 - x^2))$ ($R(u) = T_L(\sqrt{u})/\sqrt{u}$, an analytic polynomial), and $P$ is obtained from the root-finding spectral decomposition of $1 - (1 - x^2)Q^2$. The success probability achieved by the phase-sequence implementation is exactly

$$
P_S(x) = 1 - \delta^2\, T_L^2\!\bigl(c\sqrt{1 - x^2}\bigr),
$$

so $P_S \ge 1 - \delta^2$ whenever $|x| \ge \sqrt{1 - 1/c^2}$, and the threshold decreases monotonically toward 0 as $L$ grows (the fixed-point property).

## Interface and input model

```python
fixed_point_search_phases(delta, degree)
fixed_point_search(a, delta, degree)
```

API entry points: {obj}`fixed_point_search_phases <oracq.algorithms.common.qsvt.fixed_point_search_phases>`, {obj}`fixed_point_search <oracq.algorithms.common.qsvt.fixed_point_search>`

- {obj}`fixed_point_search_phases <oracq.algorithms.common.qsvt.fixed_point_search_phases>`: a purely numerical routine with no input model. `delta` is the error $\delta \in (0, 1)$; `degree` is $L$, which must be a positive odd number no greater than 20. Returns a tuple of phases in time order (length $L + 1$).
- {obj}`fixed_point_search <oracq.algorithms.common.qsvt.fixed_point_search>`: `a` is a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` (input model BE); `delta` / `degree` as above. Returns a `BlockEncoding` whose zero-signal block is the complex polynomial $P(A/\alpha)$, with success probability $|P(x)|^2$ satisfying the YLC guarantee above.

Module attributes (`fixed_point_search`):

| Attribute | Meaning |
|---|---|
| `algorithm` | `"fixed_point_search"` |
| `be_alpha` | normalization of the output BE (1.0 for this algorithm) |
| `delta` / `qsp_degree` | echo of the error parameter / $L$ |
| `threshold` | amplification threshold $\sqrt{1 - 1/c^2}$: for $\|x\| \ge$ this value, the success probability is $\ge 1 - \delta^2$ |

## Implementation notes

$Q$ is expanded analytically from the odd Chebyshev coefficients plus powers of $(1 - x^2)$; the spectral decomposition of $P$ uses Durand–Kerner root finding, clustering of repeated roots, and conjugate sign-flipped quadruple factors, with a factor-count self-check $1 + 2\times(\text{factor count}) = L$; layer stripping then recovers the phases (same pipeline as [QSP phase synthesis](qsp-phase-synthesis.md)), with a self-check that $|P(x)|^2$ and $1 - (1 - x^2)Q(x)^2$ agree pointwise on a 512-point grid (an error above 1e-5 raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`).

This algorithm performs no real-part extraction: $P$ itself is the target ($|P| \le 1$ keeps the block encoding valid). The degree cap is 20 (half of the synthesis pipeline's cap of 40). Applicability boundary: unlike Grover / amplitude amplification in `search.py`, which act on marked oracles, this entry point acts on the singular-value structure of a block encoding, and the spectral variable must lie in $[-1, 1]$.

## Validation approach

Category C2 (approximate continuous semantics; acceptance criterion in `../../development/validation-plan.md` §2). Three layers of evidence:

- Structure: `tests/core/test_qsvt.py:TransformWitnessTests.test_transform_input_validation` (negative cases such as $\delta \ge 1$ and even degrees).
- Numerical: `tests/core/test_qsvt.py:TransformWitnessTests.test_fixed_point_search_phase_properties` — $\delta = 0.4$, $L \in \{3, 5, 7\}$: phase length $L + 1$; a 400-point sweep gives $|P(x)|^2 \le 1 + 10^{-9}$ and agrees pointwise with the YLC closed form $1 - \delta^2 T_L^2(c\sqrt{1-x^2})$ (delta = 1e-5); above the threshold, $P_S \ge 1 - \delta^2$; the threshold decreases monotonically with $L$. `test_fixed_point_search_circuit_amplifies` — a 2×2 diagonal BE ($\alpha = 0.6$, spectral variables 1.0 and $1/6$), $\delta = 0.4$, $L = 5$: the `threshold` attribute is cross-checked against $\sqrt{1 - 1/c^2}$ (places = 3); the squared zero-signal amplitude in the circuit's $x = 1$ column is $\ge 1 - \delta^2 - 10^{-9}$.
- Binding: this algorithm has no independent binding witness (the input already requires a concrete BE).

## Known gaps and planned stages

Consistent with the `qsvt.py` row of the validation coverage matrix: convergence scans are missing (batch curves of the threshold–error trade-off versus $L$ are not automated, pending the convergence-scan framework of stage V2); the ill-conditioned negative cases (degree ≳ 16 or too-low precision) are partially covered by the parametrized input_validation tests, and a separate batch parametrization of "ill-conditioned inputs must raise" is still to be added.

## Related links

- Source: `src/oracq/algorithms/common/qsvt.py`
- API reference: [QSVT standard transforms](../../api/algorithms/common/qsvt.rst)
- Same-family pages: [QSP phase synthesis](qsp-phase-synthesis.md), [QSVT phase sequences](qsvt-sequence.md), [VTAA-CKS variable-time solver](vtaa-cks.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Experiment design: with $\delta=0.3$ and $L=5$ (threshold $\sqrt{1-1/c^2}\approx0.3582$), a diagonal block encoding encodes the scalar $x$ into the zero-signal block ({obj}`diagonal_block_encoding <oracq.algorithms.input_model.oracles.diagonal_block_encoding>`, both basis-state diagonal entries equal to $x$); the zero-signal success probability of the output BE of `fixed_point_search` is measured on the four backend paths (reference, rir-pysparq, adapter-pysparq, originir-ext) and compared against the YLC closed form $P_S(x)=1-\delta^2 T_L^2(c\sqrt{1-x^2})$. Two spectral points sit on either side of the threshold: $x=\cos(\pi/6)\approx0.8660$ (above the threshold, where $P_S\ge1-\delta^2=0.91$ must hold) and $x=0.3$ (below the threshold).

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| fixed-point-search-above-threshold | x=0.8660, δ=0.3, L=5 | 4 paths | success probability (closed form 0.991311) / error | 0.991311 / 8.9e-16 |
| fixed-point-search-below-threshold | x=0.3, δ=0.3, L=5 | 4 paths | success probability (closed form 0.772035) / error | 0.772035 / 1.1e-16 |

The above-threshold instance satisfies the $1-\delta^2$ guarantee (0.991311 ≥ 0.91), and the below-threshold instance still agrees pointwise with the closed form, showing that the success-probability curve of the phase-sequence implementation matches the theory across the whole spectral interval.

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

Artifacts: `out/verification/search_walks.json` (cases `fixed-point-search-above-threshold` and `fixed-point-search-below-threshold`).
