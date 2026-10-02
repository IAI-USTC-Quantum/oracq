# Gibbs State Preparation

**English** · <a href="../../../zh/manual/algorithms/gibbs-state.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.input_model.density`](../../api/algorithms/input_model/density.rst) · Stage V1

## Overview

Given a block encoding of a Hamiltonian $H$ and an inverse temperature $\beta \ge 0$, prepares an approximate purification of the Gibbs state $\rho = e^{-\beta H}/Z$ ($Z = \operatorname{Tr} e^{-\beta H}$). The implementation takes the QSVT purification route (Chowdhury–Somma 2017, van Apeldoorn–Gilyén 2019, Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)): first prepare the purification of the maximally mixed state on system and environment ($n$ Bell pairs, see [Purification access](purification.md)), then apply to system a QSVT block encoding whose zero-signal block is proportional to $g(H/\alpha)$, where on the spectral variable $x = \lambda/\alpha \in [-1,1]$

$$
g(x) = \exp\!\bigl(-\tfrac{\beta\alpha}{2}(x+1)\bigr) \in (0, 1], \qquad c = \beta\alpha/2 .
$$

$g$ is split by parity into the two branches $e^{-c}\cosh(cx)$ and $-e^{-c}\sinh(cx)$, each approximated by a modified Bessel truncation; after post-selecting signal == 0, the reduced density matrix of system is proportional to $g(H/\alpha)^2 = e^{-\beta H}$ (up to normalization).

## Interface and input model

```python
gibbs_purification(hamiltonian, beta, *, error=0.01)
```

API entry points: {obj}`gibbs_purification <oracq.algorithms.input_model.density.gibbs_purification>`

- `hamiltonian`: a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`; by convention the spectrum of $H$ lies in $[-\alpha, \alpha]$ (input model BE + DM; $\alpha$ = `be_alpha`).
- `beta`: the inverse temperature $\beta \ge 0$; at $\beta = 0$ it degenerates to the maximally mixed purification (`qsp_degree = 0`, `gibbs_scale = 1.0`).
- `error`: the uniform polynomial approximation error, which must lie in $(0, 1)$; each branch's truncation tail is controlled to $\le$ `error`/8.

Returns an {obj}`ApproximatePurification <oracq.algorithms.input_model.density.ApproximatePurification>` (`oracle_kind = "approximate_purification"`, operation signature `("system", "environment", "signal")`). Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"gibbs_purification"` |
| `beta` / `error` | echo of the call parameters |
| `qsp_degree` | the larger of the even/odd branch truncation degrees |
| `gibbs_scale` | the zero-signal block normalization $2s$ (independent of the partition function; does not affect the reduced state) |
| `success_condition` | `"signal == 0"` |
| `width` / `environment_width` / `signal_qubits` | widths of the three registers |

The classical reference {obj}`gibbs_state(hamiltonian, beta) <oracq.algorithms.input_model.density.gibbs_state>` (small-matrix $e^{-\beta H}/\operatorname{Tr}$) and {obj}`trace_distance <oracq.algorithms.input_model.density.trace_distance>` serve the witnesses and downstream reuse.

## Implementation notes

The generation chain is Bell-pair preparation → even/odd branch truncation → phase synthesis → LCU addition. Both branches are convex; the endpoint-matched imaginary-part completion always satisfies the unit-disk constraint ($f^2(x)$ does not exceed the chord joining the endpoints), so phase synthesis is always feasible; each branch, after its phases are synthesized via imaginary-part completion, has its real part extracted with $(U_\Phi + U_{-\Phi})/2$, and the two branches are then added via {obj}`linear_combination <oracq.algorithms.input_model.operators.linear_combination>`. The scaling is $s = 1.5\max(\lVert g_{\mathrm{even}}\rVert_\infty, \lVert g_{\mathrm{odd}}\rVert_\infty, 10^{-3})$, with the two branches' combined uniform error at most `error`/4.

Register layout: `system(n) | environment(n) | signal(gibbs_be.signal_qubits)`, where environment carries the other half of the Bell pairs. Applicability boundary: when $\beta\alpha$ is too large the polynomial degree exceeds the synthesis cap (40), and a {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` is raised at generation time (suggesting a smaller $\beta$ or first shrinking the spectral scale of $H$), never silently degrading; the input must be a block encoding, not a sparse oracle.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2): the post-selected reduced density matrix must approximate the classical reference Gibbs state within tolerance. Three layers of evidence:

- Structure: `tests/core/test_density.py:GibbsTests.test_beta_zero_is_maximally_mixed` (attributes of the $\beta=0$ degenerate path and a cross-check against $I/2$) and `test_invalid_inputs_fail_at_generation` ($\beta < 0$, `error` equal to 0 or $\ge 1$, a hamiltonian that is not a BlockEncoding, etc.).
- Numerical: `GibbsTests.test_diagonal_hamiltonian_gibbs` and `test_non_diagonal_hamiltonian_gibbs` — for 2×2 diagonal / non-diagonal Hamiltonians ($\beta = 0.6 / 0.8$, `error = 0.05`), the signal == 0 branch is taken on the reference simulator and, after the partial trace and renormalization, cross-checked against the `gibbs_state` classical reference with trace distance < 0.1.
- Binding: this algorithm has no separate binding witness (the input already requires a concrete BE).

Convergence witnesses (added in V1): `test_error_convergence_decreases` fixes $\beta = 0.8$ and sweeps `error` ∈ {0.4, 0.2, 0.1}, asserting that the trace distance is monotonically non-increasing across levels and that each level is ≤ `error` — measured 5.6e-4 / 8.7e-5 / 8.7e-5 (`error` 0.2 and 0.1 land on the same truncation degree and give equal distances, hence the assertion is non-increasing rather than strictly decreasing); `test_error_bound_uniform_in_beta` fixes `error = 0.1` and sweeps $\beta$ ∈ {0.2, 0.5, 1.0}, with each point's trace distance ≤ `error` — measured 2.7e-6 / 1.5e-5 / 1.9e-4, growing with $c = \beta\alpha/2$ but always far below `error`.

## Known gaps and planned stages

No known gaps; the stage V1 witnesses are complete (diagonal/non-diagonal cross-checks + a monotonically non-increasing sweep over three error levels + β-grid consistency).

## Related links

- Source: `src/oracq/algorithms/input_model/density.py`
- Algorithms in the same module: [Purification access](purification.md)
- API reference: [Density-matrix input model and Gibbs states](../../api/algorithms/input_model/density.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_hamiltonian.py` (real-backend execution, no mock substitutes). Experiment design: the non-diagonal Hamiltonian $H = \begin{pmatrix} 0.5 & 0.2 \\ 0.2 & -0.3 \end{pmatrix}$ ($\alpha = 0.7$, $\beta = 1.2$) and the diagonal Hamiltonian $H = \operatorname{diag}(0.8, -0.4)$ ($\beta = 0.6$), with `error = 0.02`; plus a $\beta = 0$ degenerate case. The purification program is executed on the reference / rir-pysparq / originir-ext paths, the signal == 0 branch is taken, the partial trace is computed with an independent numpy implementation and renormalized, the trace distance against the classical Gibbs state $e^{-\beta H}/Z$ computed by `scipy.linalg.expm` is compared (`numpy.linalg.eigvalsh`), and the post-selection success probability is reported.

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| `gibbs-trace-distance-nondiag` | 13 qubits, $\beta = 1.2$ | reference, rir-pysparq, originir-ext | trace distance | 5.10e-5 |
| ditto | — | — | success probability | 0.0953 |
| `gibbs-trace-distance-diag` | diagonal $H$, $\beta = 0.6$ | reference, rir-pysparq | trace distance | 1.28e-5 |
| ditto | — | — | success probability | 0.0992 |
| `gibbs-beta0-maximally-mixed` | $\beta = 0$ degenerate path | reference, rir-pysparq, originir-ext | trace distance | 0.0 (exactly $I/2$) |

The trace distances of all three cases are far below the `3×error = 0.06` criterion (highest measured 5.1e-5), consistent with the convergence witnesses recorded on this page (the trace distance is monotonically non-increasing in `error` and far below `error`).

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_hamiltonian.py
```

Artifacts: `out/verification/hamiltonian.json`.
