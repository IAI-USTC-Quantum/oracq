# Language and algorithm library conventions

**English** · [简体中文](../zh/reference/language.html)

This specification corresponds to oracq 0.8. It defines the boundary between the Python generation layer and RIR. For RIR objects and instruction semantics see [RIR 0.3](rir.md); for the mathematical function graph see [MIR 0.1](math-ir.md). For an introduction to operations during the generation stage see the manual chapter [Operations, registers, and the generation process](../manual/concepts.md).

## The Python generation layer

Python is the authoring environment for algorithms. Generation functions may use ordinary parameters, functions, closures, and classes to select implementations. Once generation finishes, the quantum program is represented by {obj}`Operation <oracq.infrastructure.builder.Operation>` and {obj}`Program <oracq.infrastructure.ir.Program>`; arbitrary Python callbacks must not remain in them as undefined quantum instructions.

The ordinary {obj}`Builder <oracq.infrastructure.builder.Builder>` flow does not parse the Python syntax tree. {obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>` is a separate restricted mathematical-function front end that lowers pure mathematical computation into reversible XOR operations. Its accepted syntax, numeric formats, and approximation configuration are specified in the math-functions documentation.

Generation parameters fix the concrete structure, for example register bit widths, repetition counts, quadrature nodes, truncation orders, and rotation angles. RIR has no general symbolic shape-resolution mechanism.

## Operations and modules

An `Operation` contains the entry Module and its dependency definitions. `Operation.program()` collects modules, checks for name conflicts, and forms a Program. Module definitions are shared; a call does not copy the callee body.

Registers and resources are a module's public interface. Private workspaces are declared through `Module.locals`, borrowed from the zero state and uncomputed before return. Structural validation checks references and widths; uncomputation remains an implementation obligation, and the simulator can check the actual returned state.

All calls must reference an existing module record. Unfinished implementations use `body=None`; an empty instruction body is an implemented identity operation. See [open declarations and binding](open-ir.md) for details.

## Algorithm protocols

Algorithm libraries may use Python structural protocols to describe the required methods. One object can satisfy several protocols at once; an algorithm must not infer all capabilities from a single exclusive string tag.

Common methods include `unitary()`, `state_preparation()`, {obj}`block_encoding() <oracq.algorithms.input_model.operators.block_encoding>`, `sparse_access()`, and `trotter_list()`. New algorithms can define protocols in their own modules without extending RIR.

{obj}`requires <oracq.algorithms.input_model.contracts.requires>` checks whether an object provides a protocol interface. After a concrete algorithm invokes the corresponding method, it should also check the returned value, the layout, the normalization constant, and the required call capabilities. The existence of a method is not proof of mathematical correctness.

Results of algorithm generators can likewise satisfy protocols. For example, a QLSS return object provides `state_oracle()`; a Hamiltonian evolution implementation can return a BE. A state oracle with a success signal is not automatically the same as a clean state preparation without post-selection.

## Matrix representation and normalization

The standard quantum interface of a BlockEncoding is `target` and `signal`, with a finite positive `alpha`. Its matrix interpretation is that the zero-signal corner block equals the target matrix divided by alpha.

Combinators must propagate normalization constants. The alpha of a product is the product of the alphas; the alpha of a weighted sum is the sum of each term's absolute coefficient times its own alpha. The phases of complex coefficients and the global phases inside controlled operations must be preserved.

A complete unitary can be viewed as a BE with alpha=1 on the zero-signal space, and it can also define a state preparation through its action on the zero state. When these views are used on a raw Operation, all public quantum registers belong to the complete space. If a work register must be excluded from the target of the state, the caller must explicitly declare its zero input and uncomputation conventions.

## Data access

An XOR database keeps the address unchanged and XORs the queried word into data; a query must not be treated as an overwriting assignment. QRAM resource declarations are separate from memory contents; the same read-only resource can be reused across multiple modules.

Sparse position access, element value queries, row state preparation, and block encoding are distinct interfaces. An adapter must state explicitly the data structures, auxiliary registers, and normalization it requires. There is no general guarantee of recovering an efficient sparse oracle from an arbitrary BE.

## Call capabilities

`supports_adjoint` and `supports_controlled` are existing RIR capability fields. The effective capability is the conservative conjunction of the current module and its dependencies; a call inside a controlled or adjoint context must possess the corresponding capability.

Mathematical properties such as Hermiticity, dissipation, spectral bounds, sparsity, and uncomputation are declared by algorithms and implementers. They do not automatically become new language types, nor are they proven true by a field name.

## Precision and readout

The language core does not define a universal eps type. Algorithms can use eps, word lengths, orders, or other parameters to choose a concrete generation method and record the configuration in host results or module attributes.

Approximate matrices, normalized solution states, success probabilities, and physical norms must be kept distinct. The language does not automatically perform measurement, post-selection, amplitude estimation, or classical parameter optimization. Host readout tools can carry out these steps explicitly.

## Serialization and backends

RIR text (YAML or JSON) may contain only the prescribed records and scalar attributes. Serialization preserves modules, calls, Repeat, Control, and Adjoint; an always-expanded gate list is not used as the intermediate representation.

OriginIR-ext export preserves DEF and QRAMDECL; strict-gate-set lowering produces Toffoli, U3, and CZ inside modules. PySparQ can execute at register and module boundaries. Optional native implementations do not enter RIR, nor do they automatically give a module without a gate-level body a gate-level export capability.

Backends may impose stricter execution budgets and must raise an error when a limit is exceeded. No budget limitation may be handled by silently truncating the program.

## Versions and compatibility

Package versions and intermediate-representation versions are independent. The current package is 0.8, RIR is 0.3, and MIR, PDE, and the QCL plan are all 0.1. Old RIR 0.1/0.2 continue to be read as specified; the new source layout does not change these formats.

Old Python import paths are forwarded to the canonical implementations through the compatibility layer. New applications should use the canonical paths in `infrastructure`, `algorithms`, and `applications`; see the [migration notes](../manual/compatibility.md).
