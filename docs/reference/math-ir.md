# MIR 0.1: the pure mathematical function graph

**English** · <a href="../../zh/reference/math-ir.html">简体中文</a>

MIR is the generation-layer representation between the Python pure-function front end and RIR. It round-trips through JSON independently; the full format is in the [Schema](schemas/math-ir.schema.json). It neither replaces RIR 0.1 nor writes Python callbacks into the quantum IR. For how to write functions, the implemented mathematical coverage, and numerical validation see the manual chapter [Automatically generating reversible quantum modules from ordinary math functions](../manual/math-functions.md).

## Objects and types

A {obj}`MathProgram <oracq.infrastructure.mathfunc.graph.MathProgram>` contains version="0.1", entry, and functions. Each {obj}`MathFunction <oracq.infrastructure.mathfunc.graph.MathFunction>` contains a unique name, a source-code label, parameters, ordered nodes, and returns. The function graph has no recursion; helper calls are kept as call nodes.

Parameter types are real, complex, bool, or index. index additionally carries an unsigned bit width in 1..64; other parameters have width=0, and the concrete quantum word length is decided by the lowering's {obj}`FixedFormat <oracq.algorithms.common.arithmetic.FixedFormat>`. An {obj}`Index <oracq.infrastructure.mathfunc.graph.Index>` converts to the current fixed-point representation when entering mathematical computation, and the format must accommodate its full range. Value nodes come in only three kinds: real, complex, and bool; complex lowers to two independent real registers.

A {obj}`MathNode <oracq.infrastructure.mathfunc.graph.MathNode>` contains op, kind, args, and data. A node's number is its position in nodes, and args may only reference nodes with smaller numbers. returns is a non-empty array of node numbers; a returned tuple here means several independent results — recursively nested tuple outputs are not allowed. Function parameters, node types, arities, call interfaces, and the DAG are all checked by MathProgram.validate.

| op | args / data | result |
|---|---|---|
| input | data=[parameter name] | the parameter's mathematical type; index converts to real |
| const | data=[real or boolean]; complex as [real part, imaginary part] | the declared finite numeric type |
| add/sub/mul/div/pow | two numeric nodes | real, or promoted to complex |
| neg/conj | one numeric node | same type as the input |
| real/imag/abs | one numeric node | real |
| lt/eq | two comparable nodes | bool; lt is unsupported for complex |
| and/or/not | boolean nodes | bool |
| complex | two real nodes | complex |
| select | bool, true value, false value | the common type of the two branches |
| intrinsic | data=[mathematical function name] | per the built-in signature |
| call | data=[function symbol, return index], args in formal-parameter order | the callee's corresponding return type |

Calls to the same helper with the same arguments and common subexpressions can be reused. A multi-result call has distinct return indices in MIR, and lowering merges them into a single RIR Call. Python's static range loops are unrolled into SSA during the bounded generation stage; this does not change the retention rules for helper modules or for existing RIR Repeat.

## Finite numerics and generation configuration

MIR expresses function structure; the generation configuration comprises FixedFormat and {obj}`MathConfig <oracq.infrastructure.mathfunc.numeric.MathConfig>`. The former specifies the two's-complement word length/fraction bits; the latter specifies the math-kernel order and approximation interval. These parameters and the sampled coefficients are recorded in RIR attributes. The same mathematical graph under different configurations yields different RIR module symbols that can be assembled side by side.

Lowering uses the existing Boolean circuits for fixed-point add/subtract/multiply/divide/square-root and the like. Non-polynomial real functions use sampled Chebyshev coefficients with the Clenshaw recursion; the sample count is decided by the order, is independent of the number of input bit patterns, and never generates a truth table of the whole function. Complex functions decompose into real kernels and algebraic relations.

The concrete quantum circuit implements exactly the reversible extension of its finite Boolean function F_tilde; how close F_tilde is to the ideal mathematical function F is not proven by this language. Integer constant powers use square-and-multiply with an absolute exponent bound of 128; general real powers use exp(y log(x)), and the real path requires a positive base. The complex path uses a chosen branch of the complex logarithm.

## Outputs and status

The public interface is input ports, output ports, and status[2]:

```text
|inputs, outputs, status, 0_private>
  ↦ |inputs, outputs XOR F_tilde(inputs),
     status XOR flags(inputs), 0_private>
```

Outputs are not overwriting assignments; a non-zero output is equally legal. The compiler computes the needed intermediates into private locals, copies out the results and flags, and then uncomputes, leaving the input ports unchanged. Structural checks do not constitute a proof of uncomputation; the current implementation checks this contract with small-scale references and real backends.

status[0] records domain failures, such as division by zero, square roots of negative reals, and logarithms of non-positive reals. status[1] aggregates underlying overflows/out-of-range constants and math kernels leaving the configured interval. They are neither precision bounds nor IEEE NaN/Inf indicators.

The result status of select is the condition status OR the selected branch's status. Although the unselected branch is computed and uncomputed in the reversible implementation, its domain flags do not pollute the final status. A helper's return status is the merge of the statuses of all its return values, and the caller then merges the statuses of the already-evaluated arguments. Dead computations with no return dependency can be eliminated; Python floating-point exception timing is not simulated.

There is no IEEE signed zero, and the boundary conventions of complex functions on branch cuts cannot reproduce Python cmath bit for bit; the current convention adopts the non-negative side for zero imaginary parts. Complex branch cuts, approximation precision, and whole-domain status behavior remain to be verified.

## Front-end boundary

{obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>` reads the source of an ordinary Python function and interprets a restricted AST; it does not execute the function being compiled, nor does it use dynamic eval/exec. The source string may designate entry; when none is given, the last def is used. Supported are numeric constants, local assignment/unpacking, arithmetic, comparisons, conditional expressions, structured if, static bounded range, tuple returns, pure helpers, and whitelisted math/cmath calls.

Rejected are I/O, object mutation, arbitrary object methods, dynamic loops, recursion, exception handling, generators, lambdas, and callables whose source cannot be read. Closures and explicit constants capture only finite numerics. When source is unavailable, a def string can be passed instead. math and cmath aliases and from-imports are recognized; math intrinsics currently take positional arguments.

The mathematical function-name catalog follows the <a href="https://docs.python.org/3/library/cmath.html">Python cmath documentation</a>, and AST nodes follow the <a href="https://docs.python.org/3/library/ast.html">Python AST documentation</a>. The supported syntax/numeric subset and the fixed-point behavior are delimited by this specification and must not be taken as an arbitrary-Python or complete-cmath compatible implementation.
