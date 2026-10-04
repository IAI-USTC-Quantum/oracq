# Modular Multiplication

**English** · <a href="../../zh/manual/algorithms/modular-multiplication.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.basics.number_theory`](../../api/algorithms/basics/number_theory.rst) · Stage —

## Overview

Given a multiplier $a$ coprime with the modulus $m$, implement the reversible map on a $w$-bit register

$$
U_{a,m}\,|x\rangle =
\begin{cases}
|a \cdot x \bmod m\rangle, & x < m,\\
|x\rangle, & x \ge m,
\end{cases}
$$

that is, a finite-size version of the modular-multiplication unitary required by Shor order finding. It is the underlying building block of [order finding](order-finding.md); see that page for classical factor post-processing.

## Interface and input model

```python
modular_multiply(multiplier, modulus, *, width=None, max_width=8)
```

API entry point: {obj}`modular_multiply <oracq.algorithms.basics.number_theory.modular_multiply>`

- `multiplier`: a positive integer coprime with `modulus` (input model CP; multiplier and modulus are given directly as classical parameters, no oracle input).
- `modulus`: the modulus, at least 2.
- `width`: the target bit width; when omitted, the smallest bit width that can hold `modulus - 1`.
- `max_width`: the permutation synthesis budget, default 8, at most 12.

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with register `target: Bits(width)`. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"modular_multiply"` |
| `implementation_scope` | `"bounded permutation synthesis"` |
| `modulus` | echo of the call parameters |

## Implementation notes

The generation strategy enumerates the full permutation on the $2^w$ basis states, decomposes it into disjoint cycles, and synthesizes pairwise transpositions (`_transposition`), so the action on any basis state agrees pointwise with the definition. The multiplier is reduced modulo the modulus at generation time; `modulus > 2^width` or $\gcd(a, m) \neq 1$ raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` directly.

Applicability boundary: the number of enumerated permutations grows exponentially with the bit width; the `max_width` default of 8 and cap of 12 is a deliberate synthesis budget. **The current modular multiplication uses capped permutation synthesis, at most 8 bits by default; it is used to check the order-finding interface and circuits and does not represent implemented scalable Shor modular arithmetic** (following the established wording of the algorithm catalog page). Scenarios that need large-scale modular arithmetic should switch to scalable constructions, which this module does not provide.

## Validation approach

Category C1 (exact discrete semantics; acceptance criteria in `../../development/validation-plan.md` §2): the acting unitary must equal the truth table pointwise. Evidence:

- Structure: construction and property assertions of `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`; `test_bad_inputs_fail_at_generation` covers generation-time violations such as a multiplier not coprime with the modulus ({obj}`modular_multiply(2, 4) <oracq.algorithms.basics.number_theory.modular_multiply>`).
- Numerical: `AlgorithmExpansionTests.test_modular_multiplication_total_permutation_and_order` — for `modular_multiply(2, 5)` (bit width 3) all 8 basis states are enumerated, verifying pointwise that $x<5$ maps to $2x\bmod 5$ and $x\ge 5$ is unchanged, with the target basis-state amplitude exactly 1.
- Binding: this algorithm has no separate binding witness (no open declaration entry point).

## Known gaps and planned stages

Gap per the validation-matrix wording: the small-scale witnesses **are explicitly noted as not extrapolable to full Shor** (preserving the original wording); scalable modular-arithmetic constructions are not implemented. The stage is listed as — (not in any of the V1–V4 tiers).

## Related links

- Same module: [order finding and factor post-processing](order-finding.md)
- Source: `src/oracq/algorithms/basics/number_theory.py`
- API reference: [Modular multiplication and order finding](../../api/algorithms/basics/number_theory.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_nt_qlss_sde.py` (nt_qlss_sde group), all executed on real backends; the classical oracle is independent number-theoretic permutation ground truth ($x \mapsto ax \bmod m$ for $x<m$, other basis states unchanged), not going through any helper of the implementation under test.

**Experimental design**: (a) For 5 instances $(a, m) \in \{(2,5), (3,5), (2,7), (2,9), (7,100)\}$ (bit widths 3–7), a superposition state (H on target, then the modular multiplication) **exhausts all $2^w$ basis-state inputs in a single run**, and the output is cross-checked amplitude by amplitude against the uniform superposition relabeled by the classical permutation; the paths are reference, rir-pysparq, adapter-pysparq, and OriginIR-ext (UniQC full-amplitude state vector, workspace 0). (b) For $(2,5)$, an OriginIR-ext `Circuit.to_matrix` full-unitary check: the effective block after workspace uncomputation must equal the classical permutation matrix, and the leakage into the workspace is zero.

**Key metrics**:

| Case | Scale | Path | max_error |
|---|---|---|---|
| modmul-superposition-2-5 / 3-5 / 2-7 | $2^3$ basis states, exhaustive | four paths | 5.6e-17 |
| modmul-superposition-2-9 | $2^4$ basis states, exhaustive | four paths | 8.3e-17 |
| modmul-superposition-7-100 | $2^7$ basis states, exhaustive | four paths | 4.2e-17 |
| modmul-unitary-2-5 | $8 \times 8$ full unitary | originir-ext + to_matrix | 0.0 (leakage 0.0) |

**Reproduce**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

Output: `out/verification/nt_qlss_sde.json` (`modmul-*`, 6 cases in total).
