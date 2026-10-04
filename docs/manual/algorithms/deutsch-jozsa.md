# Deutsch–Jozsa

**English** · <a href="../../zh/manual/algorithms/deutsch-jozsa.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.basics.oracle_algorithms`](../../api/algorithms/basics/oracle_algorithms.rst) · Stage V1

## Overview

Given a Boolean function $f: \{0,1\}^n \to \{0,1\}$ promised to be constant or balanced (constant means $f$ is identically 0 or identically 1; balanced means it takes 0 on half the inputs and 1 on the other half), decide which case holds with a single oracle query. The circuit exploits phase kickback: the single-bit answer is prepared in $|-\rangle$, making the XOR oracle's action equivalent to the phase $(-1)^{f(x)}$; after Hadamard interference on the uniformly superposed input, the input of a constant function collapses entirely back to $|0^n\rangle$, while the readout of a balanced function is necessarily nonzero. Classical deterministic decision requires $2^{n-1}+1$ queries in the worst case.

## Interface and input model

```python
deutsch_jozsa(function: XorDatabase)
```

API entry point: {obj}`deutsch_jozsa <oracq.algorithms.basics.oracle_algorithms.deutsch_jozsa>`

- `function`: a single-result-bit XOR database (input model FO: the Boolean function's truth table); `data_width` must be 1, otherwise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` is raised at generation time. {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` supports three-layer construction: {obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>` (gate-level truth table), {obj}`qram_database <oracq.algorithms.input_model.oracles.qram_database>` (QRAM lookup), and {obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>` (an open declaration, bound later).

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `input: Bits(address_width)` and `answer: Bits(1)`. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"deutsch_jozsa"` |
| `readout_register` | `"input"` (the decision reads only input; answer is not read) |
| `validation_stage` | `"paradigm"` (witness-level annotation) |

## Implementation notes

The generation sequence is: X then H on answer (preparing $|-\rangle$) → H on all of input (uniform superposition) → the oracle resource invoked with the `function` prefix via {obj}`invoke <oracq.algorithms.input_model.oracles.invoke>` → H on all of input. The oracle is kept as a module call (a resource binding with the `function__` prefix) and is not expanded during generation or serialization. The constant/balanced promise is the caller's responsibility: an all-zero input readout indicates constant, a nonzero one indicates balanced; no guarantee is made about the output for functions that do not satisfy the promise.

Applicability boundary: only XOR databases with `data_width == 1` are accepted; the phase-kickback route requires the answer's initial state to be preparable as $|-\rangle$, which is handled inside the circuit — the caller need not (and should not) preprocess answer.

## Validation approach

Category C1 (exact discrete semantics; acceptance criterion in `../development/validation-plan.md` §2): the output distribution must equal the exact construction of the promised case pointwise. Three layers of evidence:

- Structure: `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects` (consistency of the module's canonical objects with the historical import names, matrix criterion).
- Numerical: this algorithm has no dedicated DJ numerical test (matrix criterion), but the [Bernstein–Vazirani](bernstein-vazirani.md) witness `AlgorithmExpansionTests.test_bernstein_vazirani_recovers_secret_with_affine_bias` runs exactly the same circuit — {obj}`bernstein_vazirani <oracq.algorithms.basics.oracle_algorithms.bernstein_vazirani>` directly wraps the product of {obj}`deutsch_jozsa <oracq.algorithms.basics.oracle_algorithms.deutsch_jozsa>`, cross-checked distribution by distribution on all 3-bit secret strings with places = 10, so DJ's interference structure is covered indirectly.
- Binding: tests/core has no independent binding witness (matrix criterion is —). The `dj_gate` / `dj_qram` cases in `applications/catalog.py` bind the open declaration `BooleanFunction` to a gate-level truth table and to QRAM respectively at the L3 layer, as a catalog demonstration of the binding path.

## Known gaps and planned stages

No known gaps (the gap column of the validation matrix is —); stage V1.

## Related links

- Same module: [Bernstein–Vazirani secret-string readout](bernstein-vazirani.md), [Simon sampling](simon.md)
- Source: `src/oracq/algorithms/basics/oracle_algorithms.py`
- API reference: [Oracle query algorithms](../../api/algorithms/basics/oracle_algorithms.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical validation script: `tests/verification/verify_oracles.py` (group C cases), artifact `out/verification/oracles.json`.

Experimental design: widths 2–5 bit, with five promised functions per width — constant 0, constant 1, balanced parity ($s\cdot x$ with all-ones $s$), balanced MSB, and a balanced threshold function (nonlinear); all oracles are `gate_database` truth tables. The decision metric is $P(\mathrm{input}=0)$: it must be exactly 1 for the constant cases and exactly 0 for the balanced cases; the constant and linear-balanced families have closed-form final states $(-1)^c|0^n\rangle|-\rangle$ / $|s\rangle|-\rangle$, which additionally get a phase-sensitive amplitude-by-amplitude cross-check. There are also end-to-end cases with a {obj}`BooleanNetwork <oracq.algorithms.common.arithmetic.BooleanNetwork>`-compiled oracle (parity3 / const3): after the network is compiled via `operation()` and wrapped in the XOR database interface, the 16 branches of `net.evaluate` are first exhaustively enumerated amplitude by amplitude under superposition, then the DJ decision is run (see the numerical-validation section of [Boolean networks](boolean-networks.md)). Backend paths: reference, rir-pysparq, originir-ext.

| Case | Scale | Path | Metric value |
|---|---|---|---|
| dj-decision-w2 | 2 bit, 5 functions | reference, rir-pysparq, originir-ext | p_zero_max_error = 7.8e-16, decision_errors = 0, closed-form max_error = 2.2e-16 |
| dj-decision-w3 | 3 bit, 5 functions | same as above | p_zero_max_error = 1.1e-15, decision_errors = 0, closed-form max_error = 3.3e-16 |
| dj-decision-w4 | 4 bit, 5 functions | same as above | p_zero_max_error = 1.4e-15, decision_errors = 0, closed-form max_error = 4.4e-16 |
| dj-decision-w5 | 5 bit, 5 functions | same as above | p_zero_max_error = 1.8e-15, decision_errors = 0, closed-form max_error = 5.6e-16 |
| dj-boolean-network-parity3 | 3 bit balanced (network-compiled) | rir-pysparq, originir-ext | max_error = 8.3e-17, p_zero_error = 0.0 |
| dj-boolean-network-const3 | 3 bit constant (network-compiled) | rir-pysparq, originir-ext | max_error = 8.3e-17, p_zero_error = 1.1e-15 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_oracles.py
```

Artifact path: `out/verification/oracles.json`.
