# Open declarations, capabilities, and binding

**English** · <a href="../../zh/reference/open-ir.html">简体中文</a>

This chapter defines the behavior of open declarations in RIR 0.1. For background from the algorithm-protocol and contract perspective see the manual chapter [Conventions owned by algorithms: starting from one gate](../manual/contracts.md).

## An open declaration is not an identity operation

`Module.body` has two forms:

| body | meaning |
|---|---|
| null | An oracle declaration with a signature but not yet implemented. |
| instruction array, including the empty array | A defined implementation. The empty array explicitly denotes the identity operation. |

An open declaration must carry the string attribute oracle_paradigm. Registers, resources, and scalar attributes must still be explicit. Declarations may appear inside module calls, Control, Adjoint, and Repeat, and they participate in the normal type and aliasing checks.

A Program may hold a partial implementation. Generation, validation, and text serialization do not require all declarations to be bound. Execution and backend export require that every oracle reachable from the entry be closed; otherwise the slot name, paradigm, and call path are reported. Open declarations not invoked by the entry do not prevent exporting that entry.

All call targets must still exist. An open call therefore references an explicit declaration, not a typo or an unknown symbol.

## Capabilities and isometric state preparation

The attributes supports_adjoint and supports_controlled are booleans. Their defaults remain as in 0.1. A concrete module's capability is the conservative conjunction of its own declaration and the capabilities of the modules it calls.

If a call sits in a controlled or adjoint context, the validator checks the corresponding capability. An isometric state preparation may declare that it provides only the preparation, or it may be required to provide a complete unitary extension together with its inverse/controlled versions. Mathematical isometric preparation is no longer defined as "physically irreversible preparation".

zero_input is a library contract. At the current stage it neither proves that the quantum state at the call site is zero nor that the auxiliary space has been uncomputed. All physical qubits remain explicitly present in the interface.

## Binding

{obj}`bind(program, mapping) <oracq.infrastructure.linking.bind>` returns a new Program and does not modify the input. The keys of mapping must be the names of open declarations; the values are {obj}`Operation <oracq.infrastructure.builder.Operation>` or {obj}`Binding <oracq.infrastructure.linking.Binding>`.

An implementation's register parameters match the declared kind and width positionally. The implementation may use different local parameter names. Explicitly declared paradigm, capabilities, and be_alpha must be compatible. The bound implementation can itself contain new open dependencies, so a binding does not mean the whole program is immediately closed.

The implementation is attached to the original slot through a wrapper module of the same name. The original calling module and the interface intent remain visible; the actual implementation is not inlined by copying.

`oracle_paradigm` is an extensible interface role identifier; `PARADIGMS` lists the built-in conventions, and applications may declare new roles with custom identifiers. A concrete implementation carrying a role marker must be compatible with the declaration; historical operations without a role marker are still checked by ABI and capability. A role name does not automatically prove any mathematical property of the matrix or function.

When a declaration carries `fixed_width`, `fixed_fraction`, `fixed_signed`, or `rounding`, binding requires the corresponding implementation attributes to agree; `be_alpha` must likewise agree. Binding is rejected when the implementation explicitly negates a `zero_input` or `clean_work` required by the declaration. Historical implementations lacking these promises can still be bound, so a successful binding indicates structural compatibility, not that zero input or uncomputation has been verified.

{obj}`bind_with_report(program, mapping) <oracq.infrastructure.linking.bind_with_report>` performs the link once and returns a {obj}`BindingResult <oracq.infrastructure.linking.BindingResult>`. `result.require()` obtains the program; `result.report.to_dict()` emits the SHA-256 of the input and output programs, the explicit resource mapping, the remaining dependencies, and structured errors. Errors include the slot call path and the expected and actual values. A partial binding can succeed while dependencies remain. The report stores no Python callbacks and contains no QRAM data contents; run experiments must record memory-snapshot fingerprints separately.

Same name but different definition, cyclic dependencies, wrong argument counts, and incompatible types are all rejected. These are assembly-structure errors, not verdicts on the numerical correctness of an algorithm. For an end-to-end example of binding an operation see the tutorial [replacing an algorithm's oracle](../tutorials/oracle-binding.md).

## QRAM capture

An abstract slot may declare no QRAM resource while its concrete implementation needs one. Binding.resources maps the implementation's resource parameters to logical resource names at the entry:

```python
bound = bind(open_program, {
    "A": Binding(implementation, {"table": "matrix_values"})
})
```

The linker explicitly adds capture parameters to the intermediate modules, updates the resource arguments along the call chain, and places matrix_values in the entry's resource table. The same logical resource can be shared by several read-only oracles. Mismatched types or entry-name collisions raise errors.

Resource capture binds only handle identity; it does not embed the data table into the IR. Memory contents are still supplied separately at run time. When binding in batches, resources already promoted are retained.

Newly generated captures record their origin in the `binding_captures` attribute; when independent bindings happen in a different order, the newly added resources and the corresponding call arguments are ordered by logical name. Continuing to bind the same bank after a serialization round trip reuses the existing capture. For the attribute format and the cross-node checks see the [RIR specification](rir.md).

## The boundary of openness

What is open is the implementation body; the widths in the signature and the constants used for assembly remain concrete. Changing the M or alpha of a BE requires rerunning the same Python generator to obtain the corresponding new IR. The linker will not secretly substitute alpha inside LCU angles that were already generated.

This matches the current two-stage model: Python decides structure and configuration, and RIR keeps the quantum implementations that have not been supplied. Arbitrary symbolic widths, deferred constant-expression solving, and general dependent-type inference are not part of this version.

RIR allows multiple registers and multiple QRAM banks. banked_database demonstrates representing a 96-bit logical data word in 64/32-bit registers. Some algorithm convenience interfaces still pack the signal into a single register of at most 64 bits; larger instances should split the signal interface — that is a library-interface extension and should not widen the PySparQ single-register bit width.

## Migration from the old rules

coverage.md documents 61 old positive cases and 16 old negative cases item by item. It records paradigm expression paths and does not claim verbatim migration of old .qec files or equivalence with old CLIR goldens. Python higher-order functions, closures, and explicit resource binding replace the old require and the prohibition of partial application.

Terminal measurement/reset uses the ReadoutAction host plan. Normal RIR and OriginIR exports still preserve modules. The current UnifiedQuantum dynamic parser does not accept DEF; only the explicit host-readout adapter flattens closed modules in the final step before handing them to the dynamic parser. The directory stores both the modular program and this execution-adapted text.
