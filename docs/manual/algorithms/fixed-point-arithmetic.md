# Fixed-Point Arithmetic

**English** · <a href="../../zh/manual/algorithms/fixed-point-arithmetic.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.common.arithmetic`](../../api/algorithms/common/arithmetic.rst) · Stage V1

## Overview

A family of reversible arithmetic operations on the fixed-point format: add, subtract, multiply (plus negation, absolute value, division, reciprocal, square root, comparison, and bitwise logic). Fixed-point numbers are described by {obj}`FixedFormat(width, fraction, signed) <oracq.algorithms.common.arithmetic.FixedFormat>` — a `width`-bit word whose low `fraction` bits are fractional, and, when `signed`, whose most significant bit is the two's-complement sign bit. Outputs use XOR semantics with a two-bit `status` flag: `status[0]` marks domain failures (division by zero, square root of a negative), and `status[1]` marks results outside the word length (rather than rounding within the precision bound); the rounding policy is round toward zero with modular wraparound (`rounding = "toward_zero; modular_wrap"`). Circuits are generated with the compute/XOR/uncompute pattern of [Boolean networks](boolean-networks.md), with no dirty workspace at any point.

## Interface and input model

```python
FixedFormat(width=8, fraction=3, signed=True)
fixed_arithmetic(kind, fmt=DEFAULT_FIXED_FORMAT)
```

API entry points: {obj}`FixedFormat <oracq.algorithms.common.arithmetic.FixedFormat>`, {obj}`fixed_arithmetic <oracq.algorithms.common.arithmetic.fixed_arithmetic>`

- `FixedFormat`: requires $2 \le \text{width} \le 64$ and $0 \le \text{fraction} < \text{width} - \text{signed}$ (the fraction bits do not occupy the sign bit); violations raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`; {obj}`encode(value) <oracq.infrastructure.serialization.encode>` / {obj}`decode(value) <oracq.infrastructure.serialization.decode>` convert between real numbers and two's-complement integers. `DEFAULT_FIXED_FORMAT` is the default (8, 3, signed).
- {obj}`fixed_arithmetic(kind, fmt) <oracq.algorithms.common.arithmetic.fixed_arithmetic>`: `kind` is one of `add` / `sub` / `neg` / `abs` / `mul` / `div` / `reciprocal` / `sqrt` / `lt` / `eq` / `select` / `and` / `or` / `xor`; unknown values raise `ValidationError`; results are cached and reused per `(kind, fmt)` with `lru_cache`.

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with input `a` (a unary `kind` has only `a`; `select` adds a 1-bit `select`), plus output `out` (1 bit for `lt` / `eq`, otherwise `width` bits) and `status: Bits(2)`. This entry point consumes no oracle input (CP: the fixed-point format and operation kind are parameterized directly). Module attributes:

| Attribute | Meaning |
|---|---|
| `arithmetic_kind` / `fixed_width` / `fixed_fraction` / `fixed_signed` | echo of the call parameters |
| `rounding` | `"toward_zero; modular_wrap"` |
| `arithmetic_network` | Boolean network payload (used by the native execution entry) |

## Implementation notes

The generation strategy first builds the classical combinational logic on Boolean SSA (see [Boolean networks](boolean-networks.md)) and then compiles it as a whole into a reversible circuit. Signed multiplication and division take the magnitude-domain route: absolute values on both sides, the multiplication ($2n$-bit shift-and-accumulate, then take the `[fraction:]` segment) or restoring division completed in a widened magnitude domain, then the sign reapplied from the XOR-computed sign bit with negation; `status[0]` is driven by whether the denominator is all zero (or by the sign bit of a square-root input), and on failure the output is muxed to all zeros. Overflow criterion `status[1]`: unsigned add/sub uses carry/borrow, signed add/sub uses sign-consistency criteria (same-sign addition flipping the sign / opposite-sign subtraction flipping the sign), and multiplication, division, and square root check whether the magnitude intrudes into the sign bit. Division and reciprocal share the `1 / b` division core (the numerator of the reciprocal is fixed at $2^{2f}$).

Applicability boundary: results are truncated within the word length and wrap around modulo it; the `status` bits are the only way for the caller to detect overflow; the precision is fixed in advance by `fraction`, and there is no mechanism to grow the word length inside the circuit.

## Validation approach

Category C1 (exact discrete semantics; acceptance criteria in `../development/validation-plan.md` §2): the acting operator must equal classical evaluation pointwise. Three layers of evidence:

- Structure: `tests/core/test_stage2.py:Stage2StructureTests` — `test_arithmetic_construction_and_basis` checks, for all 14 `kind` values × `FixedFormat(4, 1)`: the module has a private workspace, the program JSON round-trips equal, and the Toffoli/U3/CZ basis export contains only the three basis gates and no `controlled_by`.
- Numerical: `Stage2StructureTests.test_gate_network_is_executable_small_example` — the unsigned 2-bit `add` is simulated at `a = 1, b = 2`, with the amplitude exactly concentrated on `(a, b, out, status) = (1, 2, 3, 0)`; the `arithmetic_network` payload of the same module is restored via `from_payload` and classically evaluates `evaluate(a=1, b=2) == {"out": 3, "status": 0}`, cross-checking the circuit against the classical semantics. The other `kind` values are covered by the structural witness and payload classical-evaluation consistency; circuit-level numerical cross-checks concentrate on the small `add` instance.
- Binding: this algorithm has no separate binding witness (matrix criterion is —). Native execution on real backends ({obj}`arithmetic_native_registry <oracq.algorithms.common.arithmetic.arithmetic_native_registry>`) is separately cross-checked against the reference simulator on PySparQ by `tests/integration/test_stage2_native.py:Stage2NativeTests` (L4 smoke, not the matrix criterion).

## Numerical validation

Experimental design (`tests/verification/verify_arithmetic.py`, output `out/verification/arithmetic.json`, 20 cases): primitive arithmetic runs three paths — "OriginIR-ext full amplitude (1–8 bit) + UniQC `Circuit.to_matrix` unitary (3–4 bit) + PySparQ RIR wide registers (16/32/64 bit)"; the compiled fixed-point functions (the Boolean SSA lowering of {obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>`) exhaust all input encodings in one superposition, checked branch by branch against independent classical semantics; elementary functions report "implementation error (against the Chebyshev coefficients in the module attribute `math_approximation`)" and "method error (coefficient polynomial vs the true function)" separately. Reproduce: `PYTHONPATH=src /path/to/backend/python tests/verification/verify_arithmetic.py`.

| Case | Scale | Path | Metric |
|---|---|---|---|
| `add_const` superposition exhaustive | w=1..8 × 4 constants | originir-ext + reference/rir/adapter cross-check | max_error ≤ 1.2e-16 |
| `add_const` unitary | w=3/4 | originir-ext + to_matrix | element-wise difference vs the classical permutation matrix = 0.0 |
| `add_const` wide register | w=16/32/64 sampled basis states | rir-pysparq | failures = 0; w=12 full superposition max_error = 1.6e-17 |
| xor/swap/controlled arithmetic | 4+4+2 bit structured program | pairwise four-path cross-check | max_pairwise_deviation = 0.0 |
| compiled polynomial `x*x+0.5` | FixedFormat(4,1)/(6,2)/(8,3) full domain | rir-pysparq | unflagged output error ≤ 1 quantum (criterion 2 quanta); the status out-of-range set matches the prediction |
| compiled `sin(x)` (degree=3) | FixedFormat(8,4), interval [−1,1] | rir-pysparq | implementation error 0.055 ≤ 2 quanta; method error 0.0453 (informational); misflagged = 0 |
| compiled-function superposition cross-check | FixedFormat(4,1) all inputs | four paths | max_pairwise_deviation = 0.0 |

## Known gaps and planned stages

No known gaps (the gap column of the validation matrix is —); stage V1.

## Related links

- Same module: [Boolean networks](boolean-networks.md), [Fourier addition](fourier-addition.md)
- Source: `src/oracq/algorithms/common/arithmetic.py`
- API reference: [Reversible arithmetic](../../api/algorithms/common/arithmetic.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)
