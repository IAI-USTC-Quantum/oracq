# Order Finding

**English** · <a href="../../zh/manual/algorithms/order-finding.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.basics.number_theory`](../../api/algorithms/basics/number_theory.rst) · Stage —

## Overview

For a multiplier $a$ coprime with the modulus $m$, find the smallest positive integer $r$ (the order) satisfying $a^r \equiv 1 \pmod m$. The implementation follows Shor's construction ([arXiv:quant-ph/9508027](https://arxiv.org/abs/quant-ph/9508027)): apply the [modular multiplication permutation](modular-multiplication.md) $U_{a,m}$ to the superposed components of the eigenstate $|1\rangle$, run quantum phase estimation on the phase register, read out an approximation of $s/r$, and let classical continued fractions recover candidate orders and attempt to factor $m$. When the order is even and $a^{r/2} \not\equiv -1 \pmod m$, $\gcd(a^{r/2} \pm 1, m)$ yields a non-trivial factor.

## Interface and input model

```python
order_finding(multiplier, modulus, *, precision=3, max_width=8)
factors_from_phase(value, precision, multiplier, modulus)
```

API entry points: {obj}`order_finding <oracq.algorithms.basics.number_theory.order_finding>`, {obj}`factors_from_phase <oracq.algorithms.basics.number_theory.factors_from_phase>`

- `multiplier` / `modulus`: the same CP parameters as modular multiplication (classical integers parameterized directly, no oracle input); the underlying modular-multiplication bit-width budget is controlled by `max_width`.
- `precision`: the bit width of the QPE phase register.
- {obj}`factors_from_phase <oracq.algorithms.basics.number_theory.factors_from_phase>` is a purely classical function: `value` is the integer readout of the phase register ($0 \le \text{value} < 2^{\text{precision}}$); it returns a verified factor pair (ascending tuple), or `None` when the sample fails to produce factors.

{obj}`order_finding <oracq.algorithms.basics.number_theory.order_finding>` returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `target: Bits(n)` ($n$ is the modular-multiplication bit width) and `phase: Bits(precision)`. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"order_finding"` |
| `modulus` / `multiplier` | echo of the call parameters |
| `readout_register` | `"phase"` (classical post-processing required after readout) |

## Implementation notes

The generation chain is {obj}`modular_multiply <oracq.algorithms.basics.number_theory.modular_multiply>` → {obj}`phase_estimation <oracq.algorithms.common.estimation.phase_estimation>`: first prepare $|1\rangle$ by applying X on the least significant bit of `target`, then call QPE. $|1\rangle$ decomposes uniformly over the cyclic eigenstates of the order, so the phase readout approximates $s/r$ ($s$ uniform), and a single sample is not guaranteed to yield the order — `factors_from_phase` uses `Fraction(value, 2**precision).limit_denominator(modulus)` to obtain candidate orders, checks evenness and $a^r \equiv 1$ for each, then verifies $\gcd(a^{r/2}\pm1, m)$, returning `None` when everything fails. When the multiplier and modulus are not coprime and the common factor is non-trivial ($1 < \gcd(a,m) < m$), the quantum path is skipped and the ascending factor pair split off from that common factor is returned directly.

Applicability boundary: the underlying modular multiplication is bounded permutation synthesis capped at 8 bits by default (see [modular multiplication permutation](modular-multiplication.md)), so order finding likewise only targets interface and circuit checks on small instances and **does not represent scalable Shor modular arithmetic**. With insufficient `precision`, continued fractions may fail to recover the correct order; callers must choose the phase bit width based on $r < m$ and accept multiple samples.

## Validation approach

Category C1 (exact discrete semantics; acceptance criteria in `../../development/validation-plan.md` §2): the output distribution and the classical post-processing results must equal the exact construction pointwise. Evidence:

- Structure: construction and property assertions of `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` (register layout `target/phase`, attribute echo).
- Numerical: `AlgorithmExpansionTests.test_modular_multiplication_and_order` and the same-name witness `test_modular_multiplication_total_permutation_and_order` — the phase distribution of `order_finding(2, 3, precision=2)` is cross-checked pointwise: $p(0) = p(2) = 0.5$ (order $r=2$, $s \in \{0,1\}$ split evenly); the post-processing path is covered by `factors_from_phase(1, 2, 2, 15) == (3, 5)` and `factors_from_phase(0, 2, 2, 15) is None`.
- Binding: this algorithm has no separate binding witness (no open declaration entry point).

## Known gaps and planned stages

Gap per the validation-matrix wording: the small-scale witnesses **are explicitly noted as not extrapolable to full Shor** (preserving the original wording); quantitative sweeps of the multi-sample success rate and of phase bit width versus recovery success rate have not been done. The stage is listed as — (not in any of the V1–V4 tiers).

## Related links

- Same module: [modular multiplication permutation](modular-multiplication.md)
- Source: `src/oracq/algorithms/basics/number_theory.py`
- API reference: [Modular multiplication and order finding](../../api/algorithms/basics/number_theory.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_nt_qlss_sde.py` (nt_qlss_sde group), all executed on real backends.

**Experimental design**: (a) For 7 $(a, m, p)$ instances ($m \in \{3,5,7,15\}$, order $r \in \{2,3,4\}$, precision 3–5), run the full order-finding routine on reference and rir-pysparq; the phase marginal distribution is compared against independent theory: $|1\rangle$ is uniform over the $r$ cyclic eigenstates, each eigenphase $s/r$ produces a Dirichlet-kernel peak through QPE, and the total distribution is $\frac{1}{r}\sum_{s=0}^{r-1} D^2(y - 2^p s/r)$ (computed independently with numpy); the $(2,3)$ instance additionally runs the adapter-pysparq and full-amplitude OriginIR-ext paths. The probability of continued fractions recovering the full order is reported as well. (b) `factors_from_phase` exhaustive sweep: $m \in \{15,21,33,35\}$ × all $a \in [2, m-2]$ × all $2^8$ phase readouts (precision=8, 23552 evaluations in total); any returned factor pair must multiply to $m$ with both factors non-trivial; the success rate is weighted by the theoretical QPE distribution and compared per modulus against an independent textbook prediction ($s$ uniform on $[1,r)$, candidate order $d = r/\gcd(s,r)$, requiring $d$ even and $a^d \equiv 1$ and $\gcd(a^{d/2}\pm1, m)$ non-trivial) — the measured shortfall below the prediction is the resolution loss of continued-fraction recovery at precision=8 for larger orders (real physics, not an implementation defect).

**Key metrics**:

| Case | Scale | Path | Metric | Value |
|---|---|---|---|---|
| order-finding-2-3-p3 | $r=2$ | four paths | TVD vs Dirichlet theory / order-recovery probability | 4.4e-16 / 0.500 |
| order-finding-2-5-p4, 3-5-p4 | $r=4$ | two paths | TVD / order-recovery probability | 6.7e-16 / 0.500 |
| order-finding-2-7-p4, 4-7-p4 | $r=3$ | two paths | TVD / order-recovery probability | 8.6e-16 / 0.459 |
| order-finding-2-15-p5, 7-15-p5 | $r=4$ | two paths | TVD / order-recovery probability | 7.8e-16 / 0.500 |
| factors-from-phase-exhaustive | 4 moduli × 56 units | classical oracle | invalid factor pairs / weighted mean success rate / minimum ratio vs textbook prediction | 0 / 0.2925 / 0.791 |

Per-modulus success rates (measured vs textbook prediction): $m=15$: 0.500 vs 0.500; $m=21$: 0.214 vs 0.233; $m=33$: 0.185 vs 0.233; $m=35$: 0.271 vs 0.318. Units that inherently fail due to odd order or $a^{r/2}\equiv-1$ account for 16/56, consistent with the theoretical expectation of Shor's construction.

**Reproduce**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

Output: `out/verification/nt_qlss_sde.json` (`order-finding-*`, `factors-from-phase-exhaustive`, 8 cases in total).
