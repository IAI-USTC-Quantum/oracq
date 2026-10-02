# Simon Sampling

**English** · <a href="../../../zh/manual/algorithms/simon.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.basics.oracle_algorithms`](../../api/algorithms/basics/oracle_algorithms.rst) · Stage V1

## Overview

Given a function $f: \{0,1\}^n \to \{0,1\}^m$ promised to be two-to-one with a nonzero XOR period $s$ (i.e. $f(x) = f(y)$ if and only if $y = x \oplus s$), find the period $s$. One quantum-side sample yields a uniformly random vector $y$ satisfying $y \cdot s = 0$ (GF(2) dot product): the circuit applies Hadamard to input, queries the oracle, entangles away the function value on output, then applies Hadamard to input again, collapsing input onto a constraint orthogonal to $s$. After collecting about $n - 1$ linearly independent samples, GF(2) elimination is done on the classical side — consistent with the repository's and the algorithm catalog's convention, elimination lives on the classical side and the quantum circuit is responsible only for sampling.

## Interface and input model

```python
simon_sample(function)
simon_nullspace(samples, width)
```

API entry points: {obj}`simon_sample <oracq.algorithms.basics.oracle_algorithms.simon_sample>`, {obj}`simon_nullspace <oracq.algorithms.basics.oracle_algorithms.simon_nullspace>`

- {obj}`simon_sample(function) <oracq.algorithms.basics.oracle_algorithms.simon_sample>`: `function` must be an {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` instance (input model FO: the two-to-one function's truth table); other types raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time. The two-to-one promise is declared by the `input_promise` attribute and borne by the caller.
- {obj}`simon_nullspace(samples, width) <oracq.algorithms.basics.oracle_algorithms.simon_nullspace>`: purely classical post-processing. `samples` is a list of integer samples encoded little endian, and `width` (1..64) is the secret string's bit width; it returns a tuple of basis vectors of the GF(2) null space of the samples' row space. When the samples exactly span $s^\perp$, the null space is one-dimensional and its nonzero element is the period $s$; when the samples are insufficient, multiple basis vectors are kept and the caller keeps sampling.

`simon_sample` returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `input: Bits(address_width)` and `output: Bits(data_width)`. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"simon_sample"` |
| `input_promise` | `"two-to-one XOR period"` (the two-to-one period promise) |
| `readout_register` | `"input"` (reading input alone yields one orthogonality constraint) |

## Implementation notes

The generation sequence is: H on all of input → the oracle invoked with the `function` prefix via {obj}`invoke <oracq.algorithms.input_model.oracles.invoke>` → H on all of input. The output register need not be measured early inside the circuit — reading out input is already a complete sample, and the function-value entanglement on output does not enter the decision. The oracle is kept as a module call. The elimination is implemented as bit-by-bit Gaussian elimination with pivot selection: for each bit, a row among the current rows with a 1 in that bit is found and swapped in as the pivot, then the bit is cleared in the other rows by XOR; the free variables are back-substituted one by one to obtain the null-space basis.

Applicability boundary: for functions that do not satisfy the two-to-one promise, neither the readout distribution nor the elimination result is guaranteed; out-of-range samples passed to `simon_nullspace` (not in $0 .. 2^{\text{width}}-1$) raise `ValidationError` at generation time.

## Validation approach

Category C1 (exact discrete semantics; acceptance criterion in `../development/validation-plan.md` §2): the sampling distribution and the elimination result must equal the closed-form construction pointwise. Three layers of evidence:

- Structure: `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects` (matrix criterion).
- Numerical: `AlgorithmExpansionTests.test_simon_sampling_and_rank_aware_postprocess` — for the two-to-one function with period $s = 3$, {obj}`gate_database(2, 1, [0, 1, 1, 0]) <oracq.algorithms.input_model.oracles.gate_database>`, the input distribution is cross-checked pointwise against $p(0) = p(3) = 0.5$; the elimination gets three closed-form cross-checks: `simon_nullspace([3], 2) == (3,)` (samples spanning the null space), `simon_nullspace([], 3)` yielding 3 basis vectors (no samples), and `simon_nullspace([1, 2, 4], 3) == ()` (samples already full rank, the null space contains only the zero vector).
- Binding: this algorithm has no independent binding witness (matrix criterion is —).

## Known gaps and planned stages

No known gaps (the gap column of the validation matrix is —); stage V1.

## Related links

- Same module: [Deutsch–Jozsa query](deutsch-jozsa.md), [Bernstein–Vazirani secret-string readout](bernstein-vazirani.md)
- Source: `src/oracq/algorithms/basics/oracle_algorithms.py`
- API reference: [Oracle query algorithms](../../api/algorithms/basics/oracle_algorithms.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical validation script: `tests/verification/verify_oracles.py` (group D cases), artifact `out/verification/oracles.json`.

Experimental design: the two-to-one linear function with period exactly $s$, $f(x)=\mathrm{delete}_p(x\oplus x_p s)$ ($p$ is the lowest set bit of $s$; after construction, two-to-oneness and the period are self-verified by classical exhaustive enumeration). Widths 3–6 bit use a `gate_database` truth-table oracle, while 7–8 bit use a CNOT linear implementation of the same function (the truth table's gate count grows with $2^n$, so wide instances switch to the register-level path; the case parameters are annotated `oracle: cnot_linear`); each width takes $s\in\{1,\ 2^n-1,\ (\texttt{0x9E37\ldots}\bmod 2^n)\mathbin{|}1\}$. The sampling distribution is compared **phase-sensitively**, amplitude by amplitude and by TVD, against the exact joint distribution (a support set of $2^{2n-2}$ equal-amplitude basis states, with the phase $(-1)^{x_0\cdot y}$ given in closed form); then $n-1$ linearly independent vectors are selected deterministically from each path's own sampled support set and handed to `simon_nullspace` for elimination, requiring the null space to be exactly $\mathrm{span}(s)$. There is also an insufficient-samples case ($n-2$ independent samples): the null space stays two-dimensional and $s$ lies in its span. Backend paths: reference and rir-pysparq, with originir-ext added for the dense small 3–5 bit instances.

| Case | Scale | Path | Metric value |
|---|---|---|---|
| simon-sampling-recovery-w3 | 3 bit, support set 16 | reference, rir-pysparq, originir-ext | max_error = 1.1e-16, tvd = 4.4e-16, recovered 6/6 |
| simon-sampling-recovery-w4 | 4 bit, support set 64 | same as above | max_error = 8.3e-17, tvd = 6.7e-16, recovered 6/6 |
| simon-sampling-recovery-w5 | 5 bit, support set 256 | same as above | max_error = 4.9e-17, tvd = 7.8e-16, recovered 6/6 |
| simon-sampling-recovery-w6 | 6 bit, support set 1024 | reference, rir-pysparq | max_error = 3.1e-17, tvd = 1.0e-15, recovered 6/6 |
| simon-sampling-recovery-w7 | 7 bit, support set 4096 | rir-pysparq | max_error = 1.7e-17, tvd = 1.1e-15, recovered 3/3 |
| simon-sampling-recovery-w8 | 8 bit, support set 16384 | rir-pysparq | max_error = 9.5e-18, tvd = 1.2e-15, recovered 3/3 |
| simon-rank-deficiency-w4 | 4 bit, 2 independent samples | rir-pysparq | nullspace_dim = 2, secret_in_span = True |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_oracles.py
```

Artifact path: `out/verification/oracles.json`.
