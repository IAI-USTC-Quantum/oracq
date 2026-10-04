# Bernstein–Vazirani

**English** · <a href="../../zh/manual/algorithms/bernstein-vazirani.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.basics.oracle_algorithms`](../../api/algorithms/basics/oracle_algorithms.rst) · Stage V1

## Overview

Given a Boolean function promised to be affine, $f(x) = s \cdot x \oplus c$ ($s \in \{0,1\}^n$ is the secret string, $c \in \{0,1\}$ is the bias, with the dot product and XOR over GF(2)), read out $s$ in full with a single oracle query. On the same phase-kickback circuit as [Deutsch–Jozsa](deutsch-jozsa.md), the linear part of the phase $(-1)^{f(x)} = (-1)^c (-1)^{s \cdot x}$ collapses input to $|s\rangle$ after Hadamard interference, and the constant factor $(-1)^c$ does not affect the readout.

## Interface and input model

```python
bernstein_vazirani(function)
affine_boolean_oracle(width, secret, *, bias=0)
```

API entry points: {obj}`bernstein_vazirani <oracq.algorithms.basics.oracle_algorithms.bernstein_vazirani>`, {obj}`affine_boolean_oracle <oracq.algorithms.basics.oracle_algorithms.affine_boolean_oracle>`

- {obj}`bernstein_vazirani(function) <oracq.algorithms.basics.oracle_algorithms.bernstein_vazirani>`: `function` is a single-result-bit XOR database (input model FO: the affine Boolean function's truth table); internally it first validates `data_width == 1` through {obj}`deutsch_jozsa <oracq.algorithms.basics.oracle_algorithms.deutsch_jozsa>`. The affine promise is the caller's responsibility. The oracle may stay an open declaration ({obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>`), bound to a gate or QRAM implementation after generation.
- {obj}`affine_boolean_oracle(width, secret, *, bias=0) <oracq.algorithms.basics.oracle_algorithms.affine_boolean_oracle>`: constructs a gate-level oracle for $f(x) = s \cdot x \oplus c$ and returns an {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>`. `width` ranges over 1..64, `secret` over $0 .. 2^{\text{width}}-1$, and `bias` takes 0 or 1; bit $i$ of `secret` corresponds to bit $i$ of the address (little endian).

`bernstein_vazirani` returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with the same registers as DJ: `input: Bits(address_width)` and `answer: Bits(1)`. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"bernstein_vazirani"` |
| `input_promise` | `"f(x)=dot(s,x) XOR c over GF(2)"` (promise declaration) |
| `readout_register` | `"input"` (reading input yields $s$ directly) |
| `validation_stage` | `"paradigm"` (inherited from the DJ circuit) |

## Implementation notes

The implementation reuses the DJ circuit: after calling `deutsch_jozsa(function)`, only the module name and attributes (`algorithm` / `input_promise`) are relabeled via `dataclasses.replace`; the dependencies and circuit structure are completely unchanged. When the promise holds, the readout `input` is exactly $s$; no guarantee is made for non-affine functions. The oracle is kept as a module call; the open-declaration route (`abstract_database` → {obj}`bind <oracq.infrastructure.linking.bind>`) has an end-to-end example in the runtime gallery: first generate a program with unbound slots via `bernstein_vazirani(abstract_database(...))`, then bind with a concrete oracle's `operation`.

## Validation approach

Category C1 (exact discrete semantics; acceptance criterion in `../development/validation-plan.md` §2): the readout distribution must equal $|s\rangle$ pointwise. Three layers of evidence:

- Structure: `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects` (matrix criterion).
- Numerical: `AlgorithmExpansionTests.test_bernstein_vazirani_recovers_secret_with_affine_bias` — all 8 secret strings at `width = 3` × bias $\{0, 1\}$, 16 combinations in total; the input distribution concentrates on $s$ with assertion precision places = 10.
- Binding: tests/core has no independent binding witness (matrix criterion is —); the open-declaration → gate-level binding usage is demonstrated by the runtime gallery and the `dj_gate` / `dj_qram` cases in `applications/catalog.py` (the same XOR database binding mechanism).

## Known gaps and planned stages

No known gaps (the gap column of the validation matrix is —); stage V1.

## Related links

- Same module: [Deutsch–Jozsa query](deutsch-jozsa.md), [Simon sampling](simon.md)
- Source: `src/oracq/algorithms/basics/oracle_algorithms.py`
- Tutorial: [Replacing an algorithm's oracle](../../tutorials/oracle-binding.md)
- API reference: [Oracle query algorithms](../../api/algorithms/basics/oracle_algorithms.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical validation script: `tests/verification/verify_oracles.py` (group B cases), artifact `out/verification/oracles.json`.

Experimental design: this algorithm is named explicitly in the paper, so per the task requirements it is exercised across **multiple widths × multiple secret strings × both biases × multiple backend paths** — widths 3–8 bit, with the secret-string set $\{1,\ 2^n-1,\ 0101\ldots,\ \texttt{0x9E37\ldots} \bmod 2^n\}$ per width (3–4 after deduplication) × bias $\{0,1\}$, 6–8 instances per width in total; each instance runs once on each of the four paths reference, rir-pysparq, adapter-pysparq, and originir-ext, with metrics the probability that the readout string equals the secret string (the input marginal distribution), the TVD from the $\delta_s$ distribution, and the amplitude-by-amplitude deviation from the closed-form final state $(-1)^c|s\rangle|-\rangle$. `affine_boolean_oracle` itself additionally gets an exhaustive superposition truth-table enumeration ($f(x)=s\cdot x\oplus c$ cross-checked branch by branch). The open-declaration route generates a program with unbound slots via `bernstein_vazirani(abstract_database(...))`, then binds the gate truth table and the QRAM implementation separately for end-to-end recovery.

| Case | Scale | Path | Metric value |
|---|---|---|---|
| bv-recovery-w3 | 3 bit, 6 instances | all four paths | success_probability = 1 - 1.1e-15, tvd = 5.6e-16, max_error = 3.3e-16 |
| bv-recovery-w4 | 4 bit, 6 instances | all four paths | success_probability = 1 - 1.4e-15, tvd = 7.2e-16, max_error = 4.4e-16 |
| bv-recovery-w5 | 5 bit, 6 instances | all four paths | success_probability = 1 - 1.8e-15, tvd = 8.9e-16, max_error = 5.6e-16 |
| bv-recovery-w6 | 6 bit, 6 instances | all four paths | success_probability = 1 - 2.1e-15, tvd = 1.1e-15, max_error = 6.7e-16 |
| bv-recovery-w7 | 7 bit, 8 instances | all four paths | success_probability = 1 - 2.3e-15, tvd = 1.2e-15, max_error = 7.8e-16 |
| bv-recovery-w8 | 8 bit, 8 instances | all four paths | success_probability = 1 - 2.7e-15, tvd = 1.3e-15, max_error = 8.9e-16 |
| bv-oracle-truth-table-w4 | 4 bit, s=0b1011, c=1 | reference, rir-pysparq, originir-ext | max_error = 5.6e-17 |
| bv-oracle-truth-table-w8 | 8 bit, s=0xA5, c=0 | reference, rir-pysparq, originir-ext | max_error = 2.8e-17 |
| bv-open-bind-w4-s11b1 | 4 bit, gate/qram dual binding | reference, rir-pysparq, originir-ext | success_probability = 1 - 1.4e-15, max_error = 4.4e-16 |
| bv-open-bind-w4-s6b0 | 4 bit, gate/qram dual binding | reference, rir-pysparq, originir-ext | success_probability = 1 - 1.4e-15, max_error = 4.4e-16 |

The deviations of success_probability from 1 all come from backend floating-point accumulation ($\le 3\times10^{-15}$); the recovered string is exactly the secret string in every instance and on every path.

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_oracles.py
```

Artifact path: `out/verification/oracles.json`.
