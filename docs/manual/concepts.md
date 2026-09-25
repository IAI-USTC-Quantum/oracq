# Operations, registers, and the generation process

**English** · [简体中文](../zh/manual/concepts.html)

An oracq program goes through three distinct stages.

1. Python executes the generation function, selecting algorithms, fixing
   parameters, and constructing modules.
2. The RIR stores module definitions, register signatures, call
   relationships, and unfinished oracles.
3. A backend exports or executes this description.

{obj}`Builder <oracq.infrastructure.builder.Builder>` constructs modules,
{obj}`Operation <oracq.infrastructure.builder.Operation>` represents a module
together with its dependencies, and {obj}`Program <oracq.infrastructure.ir.Program>`
holds the complete set of modules available to the entry point. Calling
`operation.program()` checks the structure; it does not run the quantum
algorithm.

## Registers and views

Registers are defined by name and bit width. Index 0 is the least significant
bit, and slicing follows Python's half-open convention.
{obj}`fuse <oracq.infrastructure.ir.fuse>` can join fragments of different
registers into a single logical view; `reinterpret` only changes the numeric
interpretation of the bit pattern (full semantics in the
[RIR specification](../reference/rir.md#23-register-references-and-views)).

A register or fused view is at most 64 bits wide — a limit kept consistent
with the PySparQ storage model. A program may contain multiple registers; the
total number of qubits is not restricted to 64.

Slicing, fusing, and reinterpretation generate no gates. When calling a
module, arguments must have the same width and storage type as the parameters
and must not overlap each other.

## Modules and unfinished implementations

A {obj}`Call <oracq.infrastructure.ir.Call>` in the RIR references a shared
module definition; {obj}`Repeat <oracq.infrastructure.ir.Repeat>` stores a
repetition count and a body. Neither is expanded when saved as YAML/JSON text.

`body=None` on an open oracle means no implementation has been supplied yet.
An empty instruction list, in contrast, denotes an already implemented identity
operation. This distinction runs through validation, binding, and backend
export.

{obj}`bind <oracq.infrastructure.linking.bind>` binds implementations by slot
signature and propagates newly added QRAM resources along the call chain. When
the register layout or the alpha of a block encoding changes, the generation
function should be run again. The full binding and replacement workflow is in
[Tutorial: replacing oracles](../tutorials/oracle-binding.md).

## Precision and mathematical assumptions

The language checks structure and the capabilities available for calls.
Numeric formats, quadrature nodes, Taylor orders, and linearization truncation
are decided by the algorithm generators. Whether a matrix is Hermitian,
whether a state satisfies the physical model, and whether approximation errors
are small enough all require declarations and validation at the algorithm
layer.

This division of responsibility lets new algorithms extend their own Python
protocols incrementally without changing the RIR. See
[Algorithm inputs and outputs](contracts.md) for the detailed conventions.

For a first tutorial see [The first register
program](../tutorials/first-program.md); complete listings of the core objects
are in [RIR objects](../api/infrastructure/ir.rst), [Module
builders](../api/infrastructure/builder.rst), and [Binding and capability
analysis](../api/infrastructure/linking.rst).
