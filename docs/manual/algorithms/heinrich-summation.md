# Heinrich Quantum Summation

**English** · <a href="../../../zh/manual/algorithms/heinrich-summation.html">简体中文</a>

> Category C3 · Module [`oracq.algorithms.common.integration`](../../api/algorithms/common/integration.rst) · Stage V1

## Overview

Estimates the mean $E[f] = \frac{1}{N}\sum_i f(i)$ of a non-negative integer-valued function $f$ on a uniform grid, and assembles one-dimensional numerical integration from it. Implementation basis: Heinrich 2002 ("Quantum Summation with an Application to Integration", J. Complexity 18(1); see also Novak 2001 for quantum quadrature rates over function classes).

The readout uses a comparator construction: the threshold register is put in a uniform superposition and compared against the function value, so the probability of a good state (flag = 1) is **exactly** $E[v]/2^w$ — strictly linear in $v$, with no small-angle approximation needed. Standard amplitude estimation on this marking gives query complexity $O(1/\varepsilon)$, a quadratic improvement over classical Monte Carlo's $O(1/\varepsilon^2)$.

## Interface and input model

```python
quantum_sum(database, *, precision=4, name=None)
```

API entry point: {obj}`quantum_sum <oracq.algorithms.common.integration.quantum_sum>`

- `database`: a function-value loader ({obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>`, address = index, data = value); the input model is FO + QRAM, directly compatible with the repository's three-layer binding (abstract / gate / qram). It is usually constructed by {obj}`table_loader(values, data_width=None, backend="gate"|"qram") <oracq.algorithms.common.integration.table_loader>`.
- `precision`: the number of bits in the phase register, range 1..63; the estimation error is of order $O(1/2^{\text{precision}})$.

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `target`, `work`, and `phase`. After reading out `phase`, use {obj}`mean_from_phase(value, precision, data_width) <oracq.algorithms.common.integration.mean_from_phase>` to decode the mean estimate. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"quantum_sum"` |
| `readout_register` / `decoder` | `"phase"` / `"mean_from_phase"` |
| `value_bits` / `index_bits` | value word width $w$ and index bit count $n$ |
| `query_complexity` / `classical_query_complexity` | `O(1/epsilon)` / `O(1/epsilon**2)` |

Related entry points: {obj}`quantum_integral <oracq.algorithms.common.integration.quantum_integral>` (one-dimensional integration, the mean multiplied by the interval length, decoded with {obj}`integral_from_phase <oracq.algorithms.common.integration.integral_from_phase>`) and {obj}`heinrich_rate(smoothness, dimension) <oracq.algorithms.common.integration.heinrich_rate>` (optimal convergence rates over function classes: deterministic $s/d$, randomized $s/d + 1/2$, quantum $s/d + 1$).

## Implementation notes

The generation chain is {obj}`sum_preparation <oracq.algorithms.common.integration.sum_preparation>` (uniform index + function-value loading + threshold comparison) → {obj}`sum_iterate <oracq.algorithms.common.integration.sum_iterate>` (a Grover iterate whose marking is driven by the flag bit of the prepared target) → {obj}`phase_estimation <oracq.algorithms.common.estimation.phase_estimation>`. Register layout: `target = index(n) | threshold(w) | flag(1)`, `work = value(w)`; the value word stays entangled with index in work (the comparator readout does not need it uncomputed, and the preparation is annotated `clean_work=False`), so after marking by flag the caller should invoke the preparation in reverse to restore it.

Applicability boundary: function values must be quantized as $w$-bit non-negative integers (`table_loader` validates the word width value by value); the precision of the mean estimate is determined by the QAE grid, and each additional `precision` bit doubles the grid density. The total error of the integration route (`quantum_integral`) = discretization error (determined by the grid and smoothness, see `heinrich_rate`) + QAE estimation error.

## Validation approach

Category C3 (probability-distribution semantics; acceptance criterion in `../development/validation-plan.md` §2): the output distribution must equal the closed-form expectation. Witnesses are in place; three layers of evidence:

- Structure: construction and attribute assertions in `tests/core/test_integration.py:SumPreparationTests` / `QuantumSumTests` / `RateTests`; `RateTests.test_invalid_inputs_fail_at_generation` covers parameter violations at every entry point (empty tables, negative values, out-of-range word widths, invalid backend, out-of-range precision / interval / smoothness parameters, etc.).
- Numerical: `SumPreparationTests.test_flag_probability_matches_mean` and `test_constant_table_exact` cross-check the linear identity $P(\text{flag}=1) = E[v]/2^w$ (places = 12); `QuantumSumTests.test_mean_on_qae_grid_is_exact` requires, when the mean lands exactly on the QAE grid, that all nonzero-probability readouts decode exactly to the true value (places = 9); `test_ramp_mean_within_qae_resolution` (delta = 0.5) and `test_quantum_integral_trapezoid_scale` (cross-checked against the quantized mean with delta = 0.06 and against the integral's true value 0.5 with delta = 0.09) cover the off-grid resolution bound; `RateTests.test_heinrich_rate_known_values` validates the closed-form convergence rates.
- Binding: `SumPreparationTests.test_qram_binding_matches_gate` (the flag probabilities of the gate / qram bindings cross-checked pointwise, places = 12) and `test_abstract_loader_binds` (probabilities unchanged after the abstract declaration is bound via {obj}`bind <oracq.infrastructure.linking.bind>`) together cover three-layer consistency.

## Known gaps and planned stages

No known gaps; stage V1 witnesses are in place (closed-form mean cross-checks + `heinrich_rate` validation + three-layer consistency). Per the V3 plan in validation-plan §5, Heinrich integration is pending catalog registration and wiring into real-backend cross-checks.

## Related links

- Source: `src/oracq/algorithms/common/integration.py`
- Same-family pages: [Heinrich quantum integration](heinrich-integration.md)
- API reference: [Quantum summation and integration](../../api/algorithms/common/integration.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_nt_qlss_sde.py` (the nt_qlss_sde group), all executed on real backends; the classical oracles (the classical sum $E[v]$ and the bimodal Dirichlet theoretical distribution of QAE) are constructed independently and do not go through the helper functions of the implementation under test.

**Experimental design**: (a) comparator identity — three instances (a constant table, a ramp table, and a pseudorandom 16-value table) verify on reference, rir-pysparq, adapter-pysparq, and OriginIR-ext that the good-state probability is **exactly** $E[v]/2^w$; (b) a constant table whose mean lands exactly on the QAE grid ($E[v]/2^w = 1/2$), requiring all nonzero-probability readouts to decode to the true value; (c) the QAE readout distribution compared against the independent theory $\frac{1}{2}\left[D^2(y - y_\theta) + D^2(y + y_\theta)\right]$ ($y_\theta = 2^p\theta/\pi$, $\sin^2\theta = E[v]/2^w$), plus an expected-error convergence sweep over precision 3→5; (d) the gate and qram loader bindings cross-checked pointwise on the preparation's full amplitudes and the summation's phase distribution (the qram program's resources carry the nested prefixes `prep__db__table` and `qpe__u__prep__db__table`, binding the same data table by entry-point resource name); (e) the closed-form convergence rates of `heinrich_rate`.

**Key metrics**:

| Case | Scale | Path | Metric | Value |
|---|---|---|---|---|
| sum-preparation-flag (3 tables) | $n \le 4$, $w \le 4$ | four paths | $\lvert P(\text{flag}) - E[v]/2^w\rvert$ | 1.1e-16 |
| quantum-sum-on-grid | $w=2$, precision 4 | three paths | max readout-decode deviation / cross-backend TVD | 4.4e-16 / 0.0 |
| quantum-sum-qae-p3/p4/p5 | $w=2$ | two paths | TVD vs QAE theory | 3.5e-15 / 6.2e-15 / 1.3e-14 |
| quantum-sum-qae-convergence | precision 3→5 | reference | expected absolute error $E\lvert\hat\mu - 1.5\rvert$ | 0.7205 → 0.4752 → 0.2267 (strictly decreasing) |
| table-loader-qram-vs-gate | $w=3$, precision 3 | two paths | preparation amplitudes / flag probability / summation-distribution deviations | consistent to < 1e-9 |
| heinrich-rate-closed-form | 4 $(s,d)$ pairs | classical oracle | deviation from $s/d$, $+1/2$, $+1$ | < 1e-15 |

The QAE expected error shrinks as $\Theta(\log M / M)$ (the first moment of the Dirichlet kernel's heavy tails); the monotone decrease verifies the practical behavior of the $O(1/\varepsilon)$ query complexity, and the mode estimate falls within the first-order resolution bound $2^w\pi/2^p$ at every precision.

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

Artifact: `out/verification/nt_qlss_sde.json` (`sum-preparation-*`, `quantum-sum-*`, `table-loader-*`, `heinrich-rate-*` — 10 cases in total).
