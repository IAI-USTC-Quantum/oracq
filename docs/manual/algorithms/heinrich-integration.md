# Heinrich Quantum Integration

**English** · <a href="../../../zh/manual/algorithms/heinrich-integration.html">简体中文</a>

> Category C3 · Module [`oracq.algorithms.common.integration`](../../api/algorithms/common/integration.rst) · Stage V1

## Overview

Estimates the one-dimensional definite integral $\int_0^{L} f(x)\,dx$: the interval is uniformly divided into $N$ grid points, the function values are quantized as $w$-bit non-negative integers, the grid mean $E[v]$ is first estimated by Heinrich quantum summation, and multiplying by the interval length gives the integral (composite rectangle rule). Implementation basis: Heinrich 2002 ("Quantum Summation with an Application to Integration", J. Complexity 18(1); see also Novak 2001 for quantum quadrature rates over function classes).

The integral estimate is a direct assembly of quantum summation: the readout side still uses the comparator construction to keep the good-state probability strictly linear in the function value, standard amplitude estimation is run on the marking, and the query complexity is $O(1/\varepsilon)$, a quadratic improvement over classical Monte Carlo's $O(1/\varepsilon^2)$. For the details of the summation primitive itself, see [Heinrich quantum summation](heinrich-summation.md).

## Interface and input model

```python
quantum_integral(database, *, precision=4, interval=1.0, name=None)
```

API entry point: {obj}`quantum_integral <oracq.algorithms.common.integration.quantum_integral>`

- `database`: a function-value loader ({obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>`, address = index, data = value); the input model is FO + QRAM. It is usually constructed by {obj}`table_loader(values, data_width=None, backend="gate"|"qram") <oracq.algorithms.common.integration.table_loader>`, with function values quantized as `v/full_scale` (`full_scale` defaults to $2^w - 1$).
- `precision`: the number of bits in the phase register, range 1..63; the QAE estimation error is of order $O(1/2^{\text{precision}})$.
- `interval`: the interval length $L$, which must be a positive finite real number.

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `target`, `work`, and `phase` (the same as {obj}`quantum_sum <oracq.algorithms.common.integration.quantum_sum>`). After reading out `phase`, decode with {obj}`integral_from_phase(value, precision, data_width, interval=1.0, full_scale=None) <oracq.algorithms.common.integration.integral_from_phase>`: the integral estimate is $= E[v]/\text{full\_scale} \times L$. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"quantum_integral"` |
| `readout_register` / `decoder` | `"phase"` / `"integral_from_phase"` |
| `interval` | the interval length $L$ echoed back |
| `value_bits` / `index_bits` | value word width $w$ and index bit count $n$ |
| `query_complexity` / `classical_query_complexity` | `O(1/epsilon)` / `O(1/epsilon**2)` |

Related entry point: {obj}`heinrich_rate(smoothness, dimension) <oracq.algorithms.common.integration.heinrich_rate>` gives the optimal convergence rates over function classes (deterministic $s/d$, randomized $s/d + 1/2$, quantum $s/d + 1$), used to estimate the discretization error from the grid size.

## Implementation notes

{obj}`quantum_integral <oracq.algorithms.common.integration.quantum_integral>` generates no new circuit: it calls `quantum_sum` to obtain the same module and merely rewrites `algorithm` and `decoder` to the integral versions and appends the `interval` attribute (the module attribute table is rebuilt via `dataclasses.replace`). The generation chain, register layout (`target = index(n) | threshold(w) | flag(1)`, `work = value(w)`), and Grover-iterate structure are therefore completely identical to the summation.

Design decision: the interval length never enters the circuit, only the decoder, avoiding the introduction of floating-point scaling into the IR; the quantization scale `full_scale` is likewise left to `integral_from_phase`, so the same circuit can serve reinterpretations with different intervals and ranges. Total error = discretization error (determined by the grid density and smoothness; see `heinrich_rate` for the convergence rate) + QAE estimation error (determined by `precision`).

Applicability boundary: one dimension only, uniform grids, non-negative quantized function values; multidimensional integration and adaptive grids are not implemented. `interval` must be positive; non-positive values raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time.

## Validation approach

Category C3 (probability-distribution semantics; acceptance criterion in `../development/validation-plan.md` §2): the output distribution must equal the closed-form expectation. Three layers of evidence:

- Structure: `tests/core/test_integration.py:QuantumSumTests.test_quantum_integral_trapezoid_scale` asserts the module attribute `algorithm == "quantum_integral"`; `RateTests.test_invalid_inputs_fail_at_generation` covers parameter violations such as `interval = 0.0` raising at generation time.
- Numerical: `test_quantum_integral_trapezoid_scale` — $f(x) = x$ on $[0,1]$ with an 8-point midpoint grid ($w = 4$, precision = 4); after readout decoding it is cross-checked against the closed-form quantized mean (delta = 0.06) and against the integral's true value $0.5$ (delta = 0.09); `RateTests.test_heinrich_rate_known_values` validates the closed-form convergence rates (the $s/d$ family and the quadratic spacing of the quantum rates).
- Binding: three-layer consistency of the underlying loader is covered by the summation-side witnesses (`SumPreparationTests.test_qram_binding_matches_gate` and `test_abstract_loader_binds`, gate / qram binding flag probabilities cross-checked pointwise with places = 12); the integral reuses the same circuit and needs no independent binding witness.

## Known gaps and planned stages

No known gaps; stage V1 witnesses are in place (integral closed-form cross-checks + `heinrich_rate` validation + three-layer consistency reused from the summation side). Per the V3 plan in validation-plan §5, Heinrich integration is pending catalog registration and wiring into real-backend cross-checks.

## Related links

- Source: `src/oracq/algorithms/common/integration.py`
- Same-module pages: [Heinrich quantum summation](heinrich-summation.md)
- API reference: [Quantum summation and integration](../../api/algorithms/common/integration.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_nt_qlss_sde.py` (the nt_qlss_sde group), all executed on real backends; the classical oracles are an independent classical sum and the bimodal Dirichlet theoretical distribution of QAE.

**Experimental design**: (a) a 4-point midpoint-grid instance — $f(x)=x$ on $[0,1]$ with 4 midpoints, $w=3$ (the quantized mean is exactly 0.5), precision 4; on reference, rir-pysparq, and adapter-pysparq the full phase distribution is verified against the QAE theory, and re-decoding the same distribution with `interval=2.0` verifies the interval-scaling identity (the interval length does not enter the circuit, only the decoder); (b) the canonical 8-point midpoint instance (the same one as the core test, $w=4$, precision 4, a 14-bit register) runs on rir-pysparq and is compared against the QAE theory (the reference crossover is covered by the three-backend cross-check of the structurally identical 4-point instance). Both instances' criteria include the deviation from the integral's true value $0.5$.

**Key metrics**:

| Case | Scale | Path | Metric | Value |
|---|---|---|---|---|
| quantum-integral-midpoint4 | 4 points, $w=3$, p 4 | three paths | TVD vs QAE theory / cross-backend TVD | 8.0e-15 / 7.1e-17 |
| same as above | | | mode estimate vs quantized true value 0.5 (resolution bound 0.2244) | 0.5714 (deviation 0.0714) |
| same as above | | | interval=2 scaling deviation | 0.0 |
| quantum-integral-midpoint8 | 8 points, $w=4$, p 4 | rir-pysparq | TVD vs QAE theory | 9.9e-15 |
| same as above | | | mode estimate vs quantized true value / integral true value 0.5 | 0.5333 (deviation 0.0333 / 0.0333) |

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

Artifact: `out/verification/nt_qlss_sde.json` (`quantum-integral-*` — 2 cases in total).
