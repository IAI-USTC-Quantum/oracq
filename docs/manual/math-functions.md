# Automatically generating reversible quantum modules from ordinary math functions

**English** · <a href="../../zh/manual/math-functions.html">简体中文</a>

First write a pure Python math function, then call {obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>` to generate a reversible quantum module. The function remains usable for classical computation; on the quantum side the compiler handles temporary registers, alias copying, result XOR, and uncomputation. The [QFVM](qfvm.md) Roe face already uses this path.

```python
from oracq import FixedFormat, compile_function, export_toffoli_u3_cz

def pressure(rho, momentum, energy, gamma=1.4):
    velocity = momentum / rho
    return (gamma - 1) * (energy - 0.5 * momentum * velocity)

compiled = compile_function(
    pressure,
    fmt=FixedFormat(width=12, fraction=6),
    constants={"gamma": 1.4},
)
operation = compiled.operation
program = compiled.program()
originir = export_toffoli_u3_cz(program).text
```

The module's interface is rho, momentum, energy, out, status; gamma is a generation-time parameter. It is invoked like any other {obj}`Operation <oracq.infrastructure.builder.Operation>`:

```python
from oracq import Builder, Bits

b = Builder("flow_pressure", {
    "rho": Bits(12), "momentum": Bits(12), "energy": Bits(12),
    "value": Bits(12), "flags": Bits(2),
})
# Inputs usually come from QRAM first. The Call below is not automatically expanded in the RIR.
b.call(operation, rho=b["rho"], momentum=b["momentum"],
       energy=b["energy"], out=b["value"], status=b["flags"])
application = b.finish().program()
```

Inputs stay unchanged and the result is XORed into the target register — that is the public contract of this feature. It computes a finite-word-length implementation F_tilde of the function; it never claims that a finite circuit exactly represents an arbitrary continuous mathematical function.

## The implemented mathematical scope

- Real numbers: the four arithmetic operations, negation, absolute value, comparisons, conditional selection, integer constant powers, and general powers; reuses the existing restoring division and digit-by-digit square-root circuits.
- Real functions: sqrt, exp, log, log10, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh.
- Complex numbers: addition, subtraction, multiplication, division, absolute value, real/imaginary parts, conjugate, and complex decompositions of the cmath functions above; additionally phase, polar, rect.
- Supplementary decompositions: magnitude-ratio and quadrant selection for atan2, hypot, the polynomial Clenshaw recursion, complex algebraic combinations, and status circuits selected by branch.
- Source composition: pure helpers, statically bounded ranges, multiple return values, and static tuples. Helpers remain separate MIR/RIR modules.

```python
import cmath
from oracq import MathConfig

def response(z: complex):
    return cmath.exp(1j * z) / (1 + z * z)

compiled = compile_function(
    response,
    fmt=FixedFormat(16, 8),
    config=MathConfig(
        degree=8,
        intervals=(("exp", -2.0, 2.0), ("sin", -3.14159, 3.14159)),
    ),
)
print(compiled.input_layout)   # z_real, z_imag
print(compiled.output_layout)  # out_real, out_imag
```

Ordinary unannotated parameters default to real; inputs can explicitly specify real/complex/bool. {obj}`Index <oracq.infrastructure.mathfunc.graph.Index>`(width) provides shorter public registers for unsigned integers such as QFVM row/column indices, automatically converted to the fixed-point representation on entry. Parameters with default values are treated as generation-time parameters; placed explicitly into inputs they can also become quantum inputs. Multiple results can be named with output_names, for example Roe's left/right.

Python's cmath distinguishes signed zeros on branch cuts; the current fixed-point encoding does not carry that information. This is why a family of functions and replaceable approximate implementations are provided here, without promising full float/branch compatibility. <a href="https://docs.python.org/3/library/cmath.html">See the official documentation for the corresponding cmath behavior</a>.

## Approximation and status

{obj}`MathConfig <oracq.infrastructure.mathfunc.numeric.MathConfig>`.degree controls the Chebyshev polynomial degree of real functions; intervals can replace the interval of every math kernel. The default intervals are the BOUNDS in [numeric.py](../api/infrastructure/mathfunc/numeric.rst). The configuration must cover every **intermediate math-kernel input**, not just the outermost function's input. The status flags paths that leave the interval or hit domain/word-length problems. There is currently no automatic interval derivation, error proof, or optimal arithmetic-circuit selection.

Function coefficients are generated from degree+1 classical sample points, followed by a quantum multiply-accumulate recursion; the function value for every input is never enumerated in advance. Word lengths, approximation degrees, intervals, and coefficients go into the generation configuration/module attributes; eps still does not enter the language core.

```python
def safe_inverse(x):
    return 0.0 if x == 0 else 1.0 / x
```

At x=0 this function does not set the final domain flag due to a division by zero on the unselected branch. Both branch circuits may be computed, but their results and status are selected by the quantum condition, and all intermediate bits are uncomputed at the end.

## Intermediate representations and backends

The full chain is Python pure function → [MIR 0.1](../reference/math-ir.md) → RIR 0.1 → modular OriginIR-ext / PySparQ. The MIR round-trips through JSON independently and can then be lowered by {obj}`lower_math_ir <oracq.infrastructure.mathfunc.lower_math_ir>` under a different configuration. The final RIR contains no Python callbacks; calls carrying math kernels and helpers remain {obj}`Module <oracq.infrastructure.ir.Module>`/{obj}`Call <oracq.infrastructure.ir.Call>`.

The existing {obj}`arithmetic_native_registry <oracq.algorithms.common.arithmetic.arithmetic_native_registry>` recognizes the arithmetic/Boolean implementations inside generated modules. PySparQ executes real custom C++ operators at these module boundaries, skipping the internal Boolean workspaces. The native path and the gate-level path use the same arithmetic networks and have been checked by actual execution. The direct evaluation of the original Python function is never passed off as quantum simulation here.

```python
from oracq import arithmetic_native_registry, run_pysparq
state = run_pysparq(
    application,
    native_registry=arithmetic_native_registry(application),
    memory=memory,  # if the entry declared a QRAM
)
```

QFVM still queries the raw conserved quantities; [roe_formulas.py](../api/applications/roe_formulas.rst) holds ordinary classical formulas, and {obj}`roe_face <oracq.applications.roe.roe_face>` generates a module with the original ABI through compile_function. The classical Riemann update and the QRAM data structures remain separated from the quantum matrix-element computation.

## Command line and cases

```bash
oracq compile-function examples/math_functions.py \
  --function pressure --width 12 --fraction 6 \
  --constants '{"gamma": 1.4}' \
  --mir-output out/pressure.mir.json -o out/pressure.rir.yaml

oracq emit out/pressure.rir.yaml --basis toffoli-u3-cz -o out/pressure.originir

PYTHONPATH=src .venv/bin/python tools/build_math_functions.py
```

--inputs takes a JSON type mapping such as '{"z":"complex"}' or '{"row":{"index":2},"x":"real"}'. Unsupported source statements raise {obj}`FunctionCompileError <oracq.infrastructure.mathfunc.frontend.FunctionCompileError>`; this is not a tool that can compile arbitrary Python programs.

out/math-functions/ contains pressure, roe_speed, phase_response, guarded_reciprocal, polynomial, the automatic Roe face, and the wired-up QFVM — seven groups of artifacts. See the [implementation board](../development/contributing.md) and the [validation records](../archive/function-compiler-validation.json).

Mathematical accuracy, complicated branch cuts, Roe numerical results, and resource optimization are still to be verified. Existing tests cover compilation, typing, reversible updates, module reuse, and real-backend consumption.

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_mathfunc.py` (real-backend execution, no mock substitutes). Experiment design: the 5 functions of this page (pressure, roe_speed, phase_response, guarded_reciprocal, polynomial) plus `roe_formulas.frozen_roe_face` are compiled by `compile_function` into both {obj}`FixedFormat <oracq.algorithms.common.arithmetic.FixedFormat>`(6,2)/(8,3) formats (default parameters such as gamma, order, and entropy_delta act as generation-time constants; roe_face's row/col are Index(2) with left/right dual output registers), then exhaust the input domain in one superposition on rir-pysparq — single-input functions cover all 2^6/2^8 branches of the full domain; multi-input functions are exhausted axis by axis in fibers plus a joint cube; roe_face additionally runs the full 16-combination row×col joint (including the out-of-range index 3) and a 16-point full-magnitude grid on six real axes. Expected values come from two layers: an independent bit-exact simulation of the fixed-point semantics (_Fx: mul/div/sqrt magnitudes truncated toward zero, add/sub modular wrap, status bit 0 = domain failure, bit 1 = range/word-length overflow, helper calls merge all parameter flags) and the original float64 expression. phase_response's elementary kernels reproduce the Chebyshev coefficients of the module attribute `math_approximation` bit for bit (degree=3), with the method error of the coefficient recipe versus the true function reported separately. Status flags are checked branch by branch; reference / adapter-pysparq are cross-checked three ways at the amplitude level on representative programs. Compiled function workspaces measure 358–2178 qubits, far exceeding the 24-qubit budget of OriginIR-ext, so the state-vector path is not used.

| Case | Scale | Backend path | Metric | Value |
|---|---|---|---|---|
| `polynomial-exhaustive-6.2/8.3` | 64 / 256 branches, full domain | rir-pysparq | max_error | 0 / 0 (bit-exact; flags on 44 / 197 branches all match the range-overflow prediction) |
| `guarded-reciprocal-exhaustive-6.2/8.3` | 64 / 256 branches, full domain | rir-pysparq | max_error | 0.235 / 0.123 (0.94 / 0.98 quanta, i.e. within the 1-quantum division-truncation bound; no flags on the full domain, the x=0 guard works) |
| `pressure-fibers-6.2/8.3` | 3 axes × 64 / 256 branches | rir-pysparq | max_error / method_error | 0 / 0 (constant-quantization deviation 1.4 / 0.578; rho=0 division-by-zero flags and overflow flags match branch by branch) |
| `roe-speed-fibers-6.2/8.3` | 4 axes × 64 / 256 branches | rir-pysparq | max_error / method_error | 0 / 0 (including sqrt-truncation simulation; rho≤0 domain flags match branch by branch) |
| `phase-response-fibers-6.2/8.3` | 2–3 axes × 64 / 256 branches | rir-pysparq | max_error / method_error | 0 / 0 (against the coefficient recipe; method errors 1.770 / 1.349 are the inherent approximation error of the degree=3 kernels; z=±i division-by-zero and kernel-interval overflow flags match branch by branch) |
| `phase-response-kernel-method` | 4001-point dense grid | classical comparison | method_error (exp/sin/cos) | 0.148 / 0.195 / 0.364 |
| `roe-face-index-joint-6.2/8.3` | row×col, 16 combinations | rir-pysparq | max_error | 0 / 0 (including out-of-range index 3 → 0.0; dual outputs left/right) |
| `roe-face-fibers-6.2/8.3` | 6 axes × 16-point full-magnitude grid | rir-pysparq | max_error | 0 / 0 (rho≤0, c2≤0 domain flags match branch by branch) |
| `cross-backend-*` (6 cases) | 16–64 branches | rir-pysparq / reference / adapter-pysparq | max_pairwise_deviation | 0.0 (three-way amplitudes fully identical) |
| `basis-determinism-6.2` | 1 physical point for each of 6 functions | rir-pysparq | error / failures | 0 / 0 (single basis state in → single basis state out, inputs preserved) |

Reproduction command:

```bash
PYTHONPATH=src <interpreter with pysparq+uniqc> tests/verification/verify_mathfunc.py
```

Artifacts: `out/verification/mathfunc.json` (all 28 cases pass; total runtime about 185 seconds; VERIFY_WORKERS controls the parallelism).
