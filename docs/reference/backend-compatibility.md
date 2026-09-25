# OriginIR-ext and PySparQ review

**English** · <a href="../zh/reference/backend-compatibility.html">简体中文</a>

The review date is 2026-09-08. The actual revisions are recorded in backend-revisions.json. The conclusions below come from local implementations, tests, and real compatibility experiments; neither upstream repository was modified.

## OriginIR-ext strengths and boundaries

OriginIR-ext supports QRAMDECL, DEF/ENDDEF, scalar parameters, named registers, controlled, and adjoint structures. It is suitable as oracq's first circuit-exchange backend. The oracq-side adapter is documented in the API reference under [OriginIR-ext backend](../api/infrastructure/backends/originir.rst).

However, DEF and named registers are, in the current implementation, mostly a textual authoring surface. OriginIR_BaseParser._expand_def_call inlines calls, and Circuit's usual internal form is opcode_list. Named registers are also mapped to global physical qubits after parsing. It therefore cannot directly carry this language's modular register-level IR.

| check | observed behavior | oracq handling |
|---|---|---|
| DEF parameters | Multiple fixed-width registers can be declared. | One or more DEFs are generated per RIR module. |
| Nested DEF | Can be expanded, but inner calls must remap arguments correctly. | All calls use explicit bit-by-bit arguments. |
| Whole-register arguments inside DEF | In the smallest experiment among the recorded validations, the nested call load(a,d) was rejected by the parser. | Export load(a[0],d[0],...), keeping the module call. |
| CONTROL blocks | The documentation requires ENDCONTROL to carry a list, but the current parse_line returns empty qubits for it and parsing still fails. | Use per-gate controlled_by; controlled modules are specialized with an extra control parameter while keeping DEF. |
| QRAMDECL | Declares address and data widths; memory values are supplied at run time. | Type declarations go into the exported text; data is passed separately to the execution adapter. |
| QRAM query | The address is unchanged and the data target is XORed. The least significant bit is interpreted as the first list element. | The interface and cross-checks cover non-zero targets and superposed addresses. |
| Module formal parameters for QRAM | DEF has no resource-handle formal parameters. | Module definitions are specialized on the actual global QRAM name. |
| Repetition structures | The current static subset has no execution IR that preserves arbitrary Repeat. | Export a logarithmic-size auxiliary DEF call graph. |
| Round-trip read and re-export | Re-exporting a Circuit emits a flat circuit. | The oracq serialization text (YAML/JSON) is the authoritative structural round trip. |

Reviewed locations include:
- UnifiedQuantum/uniqc/circuit_builder/originir_ext_spec.py.
- _expand_def_call, _process_statement, and _apply_op in UnifiedQuantum/uniqc/compile/originir/originir_base_parser.py.
- UnifiedQuantum/uniqc/test/core/test_originir_def.py, test_originir_named_registers.py, and test_qram.py.
- _register_qrams in UnifiedQuantum/uniqc/simulator/base_simulator.py.
- UnifiedQuantum/docs/source/1_basic_usage/originir.md.

### Data injection and execution budgets

Simulator.simulate_preprocess creates the QRAM objects. The adapter fills in the data for qram_objects and then calls simulate_statevector. The current implementation retains the data bound to the same declaration.

The adapter turns off least_qubit_remapping so that state-vector indices agree with the register mapping of the exported artifact. By default at most 24 qubits are simulated; real projects should usually pick smaller instances. The QRAM container additionally limits the sum of address and data bit widths to 30; this is not an RIR bit-width limit.

Before handing text to the parser, the adapter computes a conservative expansion budget, by default at most one million estimated steps. A Repeat(2^40) can therefore be stored and exported compactly without being accidentally expanded for execution in UnifiedQuantum.

Parse errors encountered in the recorded validations get wrapped by Simulator's format auto-fallback into OpenQASM register errors. When debugging, inspect OriginIR_BaseParser directly instead of judging the source-text format from the final exception alone.

## PySparQ's register model

PySparQ uses named registers and a global registry; each basis state stores the integer value of each register. The C++ constructor of AddRegister rejects registers wider than 64 bits. This limit does not constrain the total number of system qubits, nor may it be read as at most 64 registers. The oracq-side adapter is documented in the API reference under [PySparQ backend](../api/infrastructure/backends/pysparq.rst).

RIR's Bits, UInt, SInt, and Rational map to the General, UnsignedInteger, SignedInteger, and Rational storage interpretations respectively. Rational means an unsigned word divided by 2^width and is not the arbitrary-precision QFixed type. Current bit operations work on raw bit patterns, and add_const applies only to UInt with modulo-2^width addition semantics.

PySparQ's SplitRegister extracts the parent register's low bits and changes the parent's storage and width. CombineRegister glues them back under the corresponding low-bit rule. An RIR slice is only an immutable view and does not call these stateful structural transformations directly. The two are compatible in bit ordering but differ in lifetime semantics.

A single same-kind condition call in the native control API may overwrite the previous condition. The adapter passes all positive control bits at once to the list overload of conditioned_by_bit; zero controls are realized by flipping the relevant bits before and after the call.

### The global registry and QRAM views

The System registry is process-global. The adapter requires the registry to be empty before takeover and refuses to clear registers of an existing simulation. The adapter uses its own mutex and cleans up the state it created in a finally block; this does not mean the registry may be modified concurrently with other PySparQ code that does not use this lock.

For an arbitrary QRAM view, the adapter XOR-copies the address in low-bit order into a temporary integer register, performs the native QRAMLoad, XORs the result into the original target under the control condition, and then performs the inverse query and the inverse copy. The temporary space is removed only after being uncomputed. It supports non-zero targets and cross-register views and does not rely on copying unknown quantum states.

This version of the adapter limits materialized QRAM data to at most 2^20 entries. RIR itself stores only declarations and can represent external memory with a sparse dictionary; this execution limit does not change the IR in return.

Reviewed locations include:
- QRAM-Simulator/PySparQ/pysparq/_core.pyi.
- QRAM-Simulator/PySparQ/pysparq/operators/condition_mixin.py.
- QRAM-Simulator/SparQ/src/system_operations.cpp.
- QRAM-Simulator/PySparQ/test/test_register_split.py, test_register_capacity.py, and test_controlled_dagger.py.

## Compatibility validation in the recorded validations

Real tests cover nested DEF, resource-binding specialization, superposed addresses, non-zero data targets, QRAM views, zero controls, register modular addition, Repeat, Adjoint, and complex-coefficient LCU. The comparison target is the full complex amplitude, not just measurement probabilities.

Supporting modularity in syntax export does not mean the downstream executor reuses resources at module level. If large-scale native module interpretation is needed later, extend the downstream execution architecture or use direct PySparQ interpretation instead of bypassing the problem through one full flattening.

### The reproduced control-parsing defect

In the reviewed revision, the qubit field returned by OriginIR_LineParser.parse_line("ENDCONTROL q[0]") is None, while OriginIR_BaseParser._apply_op subsequently iterates over that field, triggering a TypeError. Instead of modifying upstream, the adapter emits inline controlled_by. For controlled module calls the adapter generates a DEF with extra control formal parameters and passes them to per-gate controls inside. Repeat's auxiliary modules work the same way.

The PySparQ execution adapter limits sparse basis states to 65536 by default and checks after native gate calls. Both normal completion and exceptional exit clean up the global registers owned by the adapter.


## 2026-09-09 addendum: application-description acceptance

In the recorded validations, 33 binding cases pass real OriginIR-ext parsing. Numerical results were not used as an acceptance criterion for the existing implementations.

BaseParser does not recognize RESET, and the dynamic parser supports RESET but does not accept DEF. Terminal measurement/reset is therefore described by a separate ReadoutAction; the normal program.originir keeps DEF. Only the explicit host-readout adapter flattens closed modules, producing execution-with-readout.originir and handing it to the real dynamic parser. This local adaptation changes neither RIR nor the modular export of ordinary applications.

## 2026-09-09 custom-operator path

Additional review of dynamic_operator/{compiler,operator_wrapper}.py and SparQ/include/basic_components.h. Dynamic-operator control is wired through split_systems/combine_systems; the generated C++ accesses the instance registers storage directly, avoiding the static name registry of the dynamic shared library being separated from the Python core. It is not a fixed-capacity problem of CACHED_REGISTER_SIZE. Real controlled calls, inverses, views, private workspaces, and Roe arithmetic have been run. The QRAM logic patch is managed by oracq, and native QRAM bank updates currently use re-materialization. Details are in the [phase-two notes](../manual/backends.md).
