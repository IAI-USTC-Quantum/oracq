# Math Function Compilation

**English** · <a href="../../zh/manual/algorithms/mathfunc.html">简体中文</a>

> Category C5 · Module [`oracq.infrastructure.mathfunc`](../../api/infrastructure/mathfunc.rst) · Stage V2

## Overview

Compiles pure Python math functions (a `math` / `cmath` subset) into fixed-point reversible quantum circuits (FO in the input-model vocabulary, function oracles). The public contract is the `reversible_function` paradigm: inputs are preserved, outputs and status are updated by bitwise XOR, and all temporary registers are automatically restored to zero — the circuit implements the finite-word-length function $\tilde F$ and makes no claim of exactly representing the continuous mathematical function. The compilation chain is Python AST → MIR 0.1 (a typed SSA graph) → modular RIR; elementary functions are approximated by Chebyshev–Clenshaw polynomials on configured intervals. A usage tutorial is in [automatic generation from ordinary math functions](../math-functions.md); the MIR format is in [math-ir](../../reference/math-ir.md).

## Interface and input model

```python
compile_function(function, *, fmt=None, inputs=None, constants=None, helpers=None,
                 output_names=None, config=None, max_unroll=128, entry=None)
lower_math_ir(program, *, fmt=None, config=None, output_names=None)
```

API entry points: {obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>`, {obj}`lower_math_ir <oracq.infrastructure.mathfunc.lower_math_ir>`

- `function`: an ordinary function object or a `def` source string; the source allows only `math` / `cmath` imports and pure function definitions (module docstrings and `__future__` imports are ignored); `entry` selects the entry point of a multi-function source (default the last one).
- `inputs`: `{parameter name: "real" | "complex" | "bool" | Index(width)}`; by default inferred from annotations (`float` / `int` → real, `complex` → complex, `bool` → bool).
- `constants`: generation-time constants (including default arguments); `helpers`: a namespace of helper functions (recursion unsupported); `max_unroll`: the unroll limit for static `range` loops.
- `fmt`: {obj}`FixedFormat <oracq.algorithms.common.arithmetic.FixedFormat>` (default `(12, 6)`, must be signed); `config`: {obj}`MathConfig(degree=6, intervals=()) <oracq.infrastructure.mathfunc.numeric.MathConfig>`, where degree is the polynomial degree of the math kernels (1..32) and intervals overrides the approximation interval per function.
- {obj}`lower_math_ir <oracq.infrastructure.mathfunc.lower_math_ir>`: lowers from the MIR alone (fmt / config can be swapped), for rebuilding after a JSON round-trip.

Returns a {obj}`CompiledFunction <oracq.infrastructure.mathfunc.lowering.CompiledFunction>` (re-exported at the top level as `oracq.compile_function`). Attributes:

| Attribute | Meaning |
|---|---|
| `operation` / `program()` | the generated {obj}`Operation <oracq.infrastructure.builder.Operation>` / the complete RIR program |
| `math_ir` | the entry {obj}`MathProgram <oracq.infrastructure.mathfunc.graph.MathProgram>` (including all helper functions) |
| `fmt` | echo of the fixed-point format |
| `input_layout` / `output_layout` | the expansion of parameters to physical registers (complex split into `_real` / `_imag`, bool 1 bit, {obj}`Index <oracq.infrastructure.mathfunc.graph.Index>` `width` bits) |

The registers are each input and output (default `out` / `out_i`, nameable via `output_names`) plus `status` (2 bits). Module attributes: `oracle_paradigm="reversible_function"`, `math_function`, `math_ir_version="0.1"`, `fixed_width` / `fixed_fraction`, `math_config`, `update_semantics`, `correctness="pending"`; the entry module carries the full `math_ir` JSON.

## Implementation notes

The frontend does **not execute** the source: side-effect statements, `while`, recursion, and lambda are rejected at compile time; static `range` loops are unrolled into straight-line code; `if` / ternary expressions / Boolean short-circuits compile to `select` nodes — both branch circuits are computed, and arithmetic errors on the untaken branch are masked through status (e.g. `0 if x == 0 else 1/x` sets no flag at x = 0). Numeric global names captured from the enclosing scope, together with `constants`, become generation-time constants.

The MIR is a sequential SSA graph (the indices referenced by a node are strictly smaller than its own); `MathProgram.validate()` checks operation arities and type constraints, and {obj}`dumps <oracq.infrastructure.serialization.dumps>` / {obj}`loads <oracq.infrastructure.serialization.loads>` round-trip as JSON. The lowering emitter {obj}`NumericEmitter <oracq.infrastructure.mathfunc.numeric.NumericEmitter>` uses only local registers and, at the exit, XORs the outputs and then replays the entire frame, realizing automatic compute / XOR / uncompute; overlapping read-only parameters are copied following an alias-free call ABI. Complex arithmetic and elementary functions decompose into real kernels via identities (e.g. $e^{x+iy} = e^x(\cos y + i\sin y)$; the $\arctan$ family goes through log / sqrt identities); integer-exponent powers use square-and-multiply (cap 128); general powers go through $\exp(b\log a)$; `atan2` is synthesized from a magnitude ratio plus quadrant correction.

The elementary kernel {obj}`elementary_kernel <oracq.infrastructure.mathfunc.numeric.elementary_kernel>`: samples degree+1 Chebyshev nodes on the interval classically and evaluates on the quantum side with the Clenshaw recurrence; the recipe (interval, degree, coefficients) is written into the module attribute `math_approximation`; leaving the approximation interval sets status bit 1, and domain violations of `log` / `asin` / `atanh` and the like set bit 0 (bit 0 = domain failure, bit 1 = range / word-length overflow). Kernels are cached and reused by `(function name, fmt, config)`. {obj}`Index <oracq.infrastructure.mathfunc.graph.Index>` inputs are embedded as unsigned integers at the fraction offset of the fixed-point word (requiring `width + fraction < fmt.width`). Applicability boundary: there is no automatic interval inference, error proof, or optimal arithmetic-circuit selection; the approximation intervals must cover the inputs of every intermediate math kernel.

## Validation approach

Category C5 (data-access layer; acceptance criteria in `../development/validation-plan.md` §2): pointwise-correct query semantics + input preservation / XOR update / zero restoration. Three layers of evidence in `tests/core/test_mathfunc.py:MathFunctionTests`:

- Structure: `test_helper_module_reuse_and_roundtrip` — helper deduplication (only one call node per call site), MIR JSON round-trip, and the RIR obtained by lowering the rebuilt MIR is byte-identical to direct compilation; `test_source_is_not_executed_and_effects_rejected` (side effects / recursion / `while` / file operations rejected, and the source is not executed).
- Numerical: `test_all_elementary_families_construct` — all 16 real elementary families of `math` / `cmath` construct closed programs under real / complex inputs, with the multi-return paths of `polar` / `rect` / `conjugate` also covered; `test_roe_is_a_compiled_pure_function` — the whole Euler ROE face-flux chain compiles, with complete helper family tags. There are additionally `test_nonzero_output_xor_superposition_and_inverse` (XOR semantics under superposed inputs and nonzero initial outputs, self-inverse on double invocation) and `test_unused_branch_error_is_masked` (division by zero on an untaken branch masked through status).
- Binding: no separate binding witness — the compiled artifact is a closed module with no abstract slots, consistent with the mathfunc row of `validation-coverage.md`.

## Known gaps and planned stages

An explicit witness for the fixed-point quantization error bound is not yet provided: word length, degree, intervals, and coefficients are already in the module attributes, but there is no automated assertion of the form "given fmt / degree, the output deviation is at most such-and-such". Registered as stage V2, consistent with the mathfunc row of `validation-coverage.md`.

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_mathfunc.py` (real-backend execution, no mock substitutes), directly closing the error-bound gap of the previous section. Experiment design: covers all 5 functions of `examples/math_functions.py` and `applications/roe_formulas.frozen_roe_face` (6 real inputs + 2 Index(2) + 2 generation-time constants + dual outputs), in both formats FixedFormat(6,2)/(8,3); rir-pysparq exhaustive superposition enumeration (full-domain 2^6/2^8 branches for single inputs, per-axis fibers plus a joint cube for multiple inputs, and the full 16-combination Index joint). The oracle is `_Fx`, a fixed-point simulator implemented independently per the semantics of {obj}`fixed_arithmetic <oracq.algorithms.common.arithmetic.fixed_arithmetic>` (mul/div/sqrt truncate magnitudes toward zero, add/sub wrap modularly, status bit 0 = domain failure, bit 1 = range/word-length overflow, `select` masks the unselected branch while helper calls merge all argument flags, and helper formals are not folded), and the elementary kernels additionally reproduce the Clenshaw recurrence bit for bit from the `math_approximation` coefficients — so the "implementation error" is measured against the **quantized coefficient recipe**, while the "method error" (recipe vs the true function) is reported separately. Results: every case except guarded_reciprocal is **bit-identical** to the _Fx simulation (max_error = 0, covering thousands of flag-free branches); guarded_reciprocal has full-domain error ≤ 1 quantum (the division-truncation bound) and no flags anywhere in the domain; the status flags (division by zero at rho=0, square root at rho<0, c2≤0, z=±i, kernel-interval overflow, and word-length overflow of intermediates including division by zero after wrap) match the predictions branch by branch (misflagged = 0); the three backends (rir-pysparq / reference / adapter-pysparq) agree pairwise with amplitude deviation 0.0; basis-state determinism (inputs preserved, outputs XOR-ed, locals restored to zero) passes for all 6 functions. The measured workspace is 358–2178 qubits, so the 24-qubit OriginIR-ext budget does not apply.

| Case | Scale | Backend path | Metric | Value |
|---|---|---|---|---|
| polynomial full-domain enumeration | 6.2 / 8.3, 64 / 256 branches | rir-pysparq | max_error | 0 / 0 (exact_fraction 1.0) |
| guarded_reciprocal full-domain enumeration | 6.2 / 8.3, 64 / 256 branches | rir-pysparq | max_error | 0.94 / 0.98 quanta (status = 0 on the full domain) |
| pressure fiber + cube | both formats, 192+64 / 768+64 branches | rir-pysparq | max_error / constant-quantization deviation | 0 / 1.4 (6.2), 0 / 0.578 (8.3) |
| roe_speed fiber + cube | both formats, 256+16 / 1024+16 branches | rir-pysparq | max_error | 0 / 0 |
| phase_response fiber + cube | both formats, 192+16 / 512+16 branches | rir-pysparq | max_error / method error | 0 / 1.770 (6.2), 0 / 1.349 (8.3) |
| frozen_roe_face joint index + fiber | 16 + 96 branches × both formats | rir-pysparq | max_error | 0 / 0 (left/right dual outputs) |
| kernel method error (informational) | 4001-point grid, degree=3 | classical comparison | exp / sin / cos | 0.148 / 0.195 / 0.364 |
| cross-backend check (6 cases) | 16–64 branches | rir-pysparq / reference / adapter-pysparq | max_pairwise_deviation | 0.0 |
| basis-state determinism | 6 functions | rir-pysparq | failures | 0 |

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_mathfunc.py
```

Artifact: `out/verification/mathfunc.json` (all 28 cases pass, total runtime about 185 seconds).

## Related links

- Source: `src/oracq/infrastructure/mathfunc/` (`frontend.py` / `graph.py` / `numeric.py` / `lowering.py`)
- Usage tutorial: [automatic generation of reversible quantum modules from ordinary math functions](../math-functions.md)
- API reference: [math function compilation entry points](../../api/infrastructure/mathfunc.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)
