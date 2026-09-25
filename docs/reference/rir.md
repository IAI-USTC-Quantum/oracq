# RIR: the modular register-level intermediate representation

**English** · <a href="../../zh/reference/rir.html">简体中文</a>

This specification defines RIR 0.1. Detailed rules for open declarations, capabilities, and binding are in [the open IR](open-ir.md); for an introduction to working with registers and views in the Python generation layer see the manual chapter [Operations, registers, and the generation process](../manual/concepts.md#registers-and-views).

Specification date: 2026-09-19. The Python API, serializers, reference executor, and backend adapters must comply with this document. The JSON Schema covers object shapes only; cross-node constraints are checked by validate.

This specification organizes RIR as a complete language-layer specification in five parts:

- **Part 1: design positioning** — what RIR expresses and what it does not;
- **Part 2: data model** — the object model, storage types and bit order, register references and views;
- **Part 3: instruction set** — the semantics and constraints of each of the seven instructions;
- **Part 4: program-level rules** — the call graph, text encoding, formal grammar, backend lowering, and private workspaces;
- **Part 5: complete examples** — six directly runnable end-to-end programs checked field by field against the rules of the first four parts.

## Part 1: Design positioning

### 1.1 Design scope

RIR represents quantum operations whose compile-time parameter evaluation has already completed. All register widths, integer constants, and rotation angles are concrete. RIR keeps module definitions, module calls, static repetition, control, and adjoint blocks, and does not require expansion into qubit-level circuits.

Implemented bodies in RIR consist of unitary operations; unimplemented modules are represented as open declarations. RIR contains no measurement, reset, runtime classical feedback, or opaque Python callback. For fixed external QRAM contents, every legal instruction has the unitary semantics prescribed below on the full quantum space; the sole exception is Store — a QRAM random write is modeled on classical memory cells and evolves the memory mapping itself rather than the quantum state (see section 3.2) — so it may appear only in non-controlled, non-adjoint positions of a module body.

Public workspaces and signal registers appear in the call interface; a module may additionally declare private zero-input, zero-output workspaces through locals. Structural validation provides no proof of auxiliary-qubit uncomputation in general.

## Part 2: Data model

### 2.1 Object model

Program has three fields: entry, modules, and version. The current version is the string "0.1". entry must reference a defined module.

Module has name, registers, resources, body, attributes, and locals. registers is an ordered list of quantum parameters, resources is an ordered list of QRAM parameters, and body is an ordered instruction body, or null for an open declaration. attributes is an ordered list of key/value pairs. Attribute values may only be strings, integers, finite floats, or booleans.

Module names are unique within a Program. Register names and resource names are unique within their module and must not conflict with each other. Identifiers match [A-Za-z_][A-Za-z0-9_]*. A module has at least one non-empty quantum interface; the register parameters themselves may contain zero-width items.

Attributes do not change instruction semantics. Libraries may agree on mathematical interpretations based on attributes. For example, be_alpha together with the library-prescribed zero-projection layout defines the scale of a block encoding, but the executor does not additionally scale the quantum state based on attributes.

The linker preserves the attribute `binding_captures`, whose value is a JSON object string mapping logical resource names to this module's local resource names. The mapping must reference declared resources and the local names must not repeat; validate checks these constraints. It lets batched bindings reuse already-captured resources and canonicalizes, by logical name, the parameter order of newly added resources and of all call arguments. The original explicit resource parameter order is kept; the attribute does not change execution semantics. Historical programs without this attribute remain readable. Binding reports and the open-cost ledger are stored separately, adding no RIR nodes.

### 2.2 Storage types and bit order

A Register consists of name and RegType. A RegType consists of kind and width. width must be a strict Python/JSON/YAML integer; a boolean must not impersonate a width.

| kind | bit-pattern interpretation |
|---|---|
| bits | a bit string with no numeric meaning assigned; corresponds to PySparQ General. |
| uint | unsigned integer in the range 0 to 2^width−1. |
| sint | two's-complement integer; bit operations still act on the raw word. |
| rational | unsigned word divided by 2^width; corresponds to PySparQ Rational. |

The width of a register or view must lie in 0..64. Zero width denotes an empty interface; native adapters create no physical register for it. The total qubit count and the number of registers are not capped at 64.

Register index zero is the least significant bit. Suppose the entry signature declares registers r0, r1, and so on in order; then the exported state-vector index satisfies:

```text
index = value(r0) + (value(r1) << width(r0)) + ...
```

This convention defines vectorization and backend physical mapping only. RIR itself stores no global physical qubit numbers.

### 2.3 Register references and views

A Ref consists of parts and type. parts is an ordered list of Spans. A Span has register, start, and width, denoting a contiguous range of one root register in the current module.

Ref.parts are concatenated from low bits to high bits. The sum of the segment widths must equal Ref.type.width. Every segment must lie within its root register's range. Segments within one Ref must not reference duplicate qubits.

Slicing, fuse, and reinterpret form Refs at the Python generation stage and generate no quantum instructions. Slicing produces a bits interpretation; reinterpret explicitly changes the interpretation without changing the bit width or the quantum state. fuse may cross root registers, but the merged view still may not exceed 64 bits.

If f is a formal register of some module and the call argument is a Ref spanning root registers, then a local slice of f must map to the sub-view of that Ref with the same low-bit offsets. The mapping must not assume that the actual argument is a contiguous run of physical qubits.

RIR references use lexical name binding. They are declarative references, not proofs about Python variable lifetimes or unique object references.

## Part 3: Instruction set

Every instruction has fixed parameter fields. Unknown instructions are not a legal extension point; the specification must be revised first and backend behavior provided.

### 3.1 Primitive

Primitive has op, operands, angle, and value. operands is a list of Refs. An unused angle or value must be null.

| op | operands | parameters | semantics |
|---|---|---|---|
| h, x, y, z, s, t | one register view | none | Broadcast the standard gate onto each bit of the view in low-to-high order. |
| rx, ry, rz | one register view | angle | Apply exp(−i angle P/2) to each bit. |
| phase | one register view | angle | Apply diag(1, exp(i angle)) to each bit. |
| gphase | none | angle | Multiply the current quantum state by exp(i angle); the relative phase under control must be preserved. |
| xor | two views of equal width | none | Keep the first input; the second input is XORed bit-wise with the first input. |
| swap | two views of equal width | none | Exchange the raw bit patterns of the two words. |
| add_const | one uint view | value | Perform addition modulo 2^width on the word. |

All angles are finite real numbers in radians. The value of add_const must be an integer within 0..2^width−1. The two views of a binary operation must not overlap at all.

A broadcast gate remains a single register-level operation in the IR; it does not become multiple RIR nodes just because it contains several physical gates. Broadcast, xor, swap, and adding zero on zero width are all no-ops.

### 3.2 Load and Store

Load has resource, address, and data. resource references a QRAM formal parameter of the current module. The widths of address and data must match the declaration, and the two views must not overlap.

The QRAM type has address_width and data_width, both in 1..64. Memory is a mapping M from addresses to unsigned data words, with unspecified cells equal to zero. Memory data is not written into the program serialization text but is bound separately as an execution input (the binding file format is defined in [QRAM memory definition](qram-memory.md)). In a program without Store, M is fixed; in a program with Store, M evolves in instruction order.

```text
|address>|data>  ↦  |address>|data XOR M[address]>
```

Load holds for an arbitrary data target; the target is not required to start at zero. Load is self-inverse and modifies neither the address nor the classical memory.

Store has the same resource, address, and data fields as Load and the same width and overlap constraints. Its semantics is a random write: at the moment the instruction executes, the address and data views must be in a definite basis state (their values are unique on the full quantum state); then the classical cell is assigned:

```text
M[address] := data
```

Store changes no qubit and is not counted toward backend gate cost. Storage cells are modeled as classical cells; a write under a superposed address or superposed data has no linear semantics, and an executor must raise an error upon encountering one. Structurally, Store may appear only in a module body or a Repeat body; it is forbidden inside Control and Adjoint bodies, and a module containing Store (including one reachable through calls) has neither the supports_adjoint nor the supports_controlled capability. Store is legal in RIR 0.1.

### 3.3 Call

Call has module, arguments, and resources. module references the called Module. Quantum arguments and resource arguments are both bound in the signature order of the called module.

The kind and width of every quantum argument must be exactly the same as the formal parameter's. When the bit interpretation must change, the caller must reinterpret explicitly. All quantum arguments must be pairwise non-overlapping.

A resource argument references a resource name of the current module, and its QRAM type must exactly match the formal parameter. Multiple read-only resource formal parameters may bind to the same actual QRAM; since resource bindings point at the same memory mapping, a Store executed through one formal name is visible to the other aliases. Resource binding binds neither the quantum address nor the data register.

The call semantics replaces all local references of the called module with the actual views and executes the body on the same quantum state. This definition does not require inlining the body at storage time.

### 3.4 Repeat

Repeat has count and body. count is an integer in 0..2^63−1. The semantics applies body count times in sequence. A count of zero is the identity operation.

Even when count is zero, the body must be structurally legal. Serialization and the generation stage must not unconditionally copy the instruction body according to count.

### 3.5 Control

Control has register, value, and body. register must be a non-empty view, and value must lie within its unsigned bit-pattern range.

The body applies when the control view equals value; otherwise the identity applies. The control register participates as a quantum condition in coherent control and is neither measured nor turned into a Python conditional.

Control bits are protected throughout the body. Quantum operands of primitives, QRAM, and module calls must not overlap them. Module calls take the conservative rule: even if the called module does not actually modify a formal parameter, control bits still must not be passed as that parameter. Store must not appear inside a Control body.

Views of nested controls must not overlap each other. Distinct control conditions combine by logical conjunction. gphase has no operands, so a control is allowed to cover all qubits of a module; the semantics is still a controlled phase.

### 3.6 Adjoint

Adjoint has body. Its semantics reverses the order of the body's instructions and takes the adjoint of each operation. The angles of rx, ry, rz, phase, and gphase are negated; add_const becomes modular subtraction; Load, xor, swap, h, x, y, z are self-inverse; s and t take the corresponding inverse phases.

The adjoint of a Call refers to the inverse of the called module, the adjoint of a Repeat repeats its inverse body, and the adjoint of a Control keeps the control condition and inverts its body. A double adjoint restores the original operation. Store is a non-unitary side effect and must not appear inside an Adjoint body; a module containing Store has no adjoint capability.

These rules are instruction semantics; RIR is not required to rewrite or expand the body immediately when an Adjoint is created.

## Part 4: Program-level rules

### 4.1 Call graph and module reuse

A Program's call graph must be acyclic. All declared modules are checked, including unreachable ones. All call targets must exist. The current implementation limits module call depth and structural-block nesting depth to 127.

Modules share one definition across call sites; a call does not copy the signature or body into Program.modules. A Python generator using one name for two different definitions must raise an error.

Every call target must have an explicit Module record. The record may be an open declaration; an undefined name is not the same as an unfinished implementation.

### 4.2 Text encoding (YAML and JSON)

RIR text uses YAML by default; {obj}`dumps(program, format="json") <oracq.infrastructure.serialization.dumps>` emits the equivalent canonical JSON, and {obj}`loads <oracq.infrastructure.serialization.loads>` accepts both formats (YAML is a superset of JSON, and the object trees parsed from the two texts agree field by field). Every data class uses an object with a tag field. The tag values correspond to the following record names:

```text
Program, Module, Register, RegType, Span, Ref, QRAM, Resource,
Primitive, Load, Store, Call, Repeat, Control, Adjoint
```

All fields must be written out, including null, empty lists, and empty attributes. The decoder rejects extra fields, missing fields, unknown tags, duplicate keys, non-finite floats, and unknown versions; the YAML path additionally rejects implicit scalar types such as dates, which must be quoted when they need to be treated as strings.

Immutable Python tuples are encoded as arrays. Deserialization restores them as tuples; deserializing arbitrary Python objects into callable code is not allowed.

The following example is an H broadcast instruction acting on a two-bit integer register:

```yaml
angle: null
op: h
operands:
  - parts:
      - register: address
        start: 0
        tag: Span
        width: 2
    tag: Ref
    type:
      kind: uint
      tag: RegType
      width: 2
tag: Primitive
value: null
```

Canonical YAML output uses UTF-8, block style, two-space indentation, key sorting, nested sequences indented relative to their owning key, and one trailing newline, without anchors or aliases. Canonical JSON output uses UTF-8, two-space indentation, key sorting, and one trailing newline. In both formats module definitions are sorted by module name, while the list order of signature parameters and instructions is preserved. Floats use Python's finite float representation; NaN and Infinity are forbidden.

The order of module attributes is also preserved. The Builder sorts attributes by key; externally constructed IR should adopt the same order to obtain the same canonical output. The current specification guarantees deterministic output for one and the same IR; it does not require all semantically equivalent circuits to have identical text.

The schema file is [rir.schema.json](schemas/rir.schema.json). It describes the encoded object structure and local ranges shared by YAML and JSON; reference resolution, aliasing, arity, cross-node types, control protection, and call-graph rules still require running the semantic validator. A formal production summary of the object shapes is in section 4.3.

### 4.3 Formal grammar

This section summarizes the object shapes defined in Parts 2 and 3 and in sections 4.1 and 4.2 as productions, for independent implementations to check against. The grammar covers structure and fields only; semantic rules such as aliasing, control protection, the call graph, and numeric ranges are governed by the body text and the semantic validator.

Lexical conventions: name matches `[A-Za-z_][A-Za-z0-9_]*`; integer is a strict integer, and a boolean must not impersonate one; float is a finite float; string is an arbitrary string scalar. `∅` denotes an absent field (the empty body of an open declaration, or an unused angle/value). `X*` denotes an ordered immutable tuple, possibly empty.

```text
program     = Program { entry: name;
                        modules: module*;
                        version: "0.1" } .

module      = Module { name: name;
                       registers: register*;
                       resources: resource*;
                       body: instruction* | ∅;
                       attributes: attribute*;
                       locals: register* } .

register    = Register { name: name; type: regtype } .
regtype     = RegType { kind: "bits" | "uint" | "sint" | "rational";
                        width: integer } .
resource    = Resource { name: name; type: qram } .
qram        = QRAM { address_width: integer; data_width: integer } .
attribute   = name "×" (string | integer | float | boolean) .
span        = Span { register: name; start: integer; width: integer } .
ref         = Ref { parts: span*; type: regtype } .

instruction = primitive | load | store | call | repeat | control | adjoint .

primitive   = Primitive { op: gate-op;
                          operands: ref*;
                          angle: float | ∅;
                          value: integer | ∅ } .
gate-op     = "h" | "x" | "y" | "z" | "s" | "t"
            | "rx" | "ry" | "rz" | "phase" | "gphase"
            | "xor" | "swap" | "add_const" .
load        = Load { resource: name; address: ref; data: ref } .
store       = Store { resource: name; address: ref; data: ref } .
call        = Call { module: name; arguments: ref*; resources: name* } .
repeat      = Repeat { count: integer; body: instruction* } .
control     = Control { register: ref; value: integer; body: instruction* } .
adjoint     = Adjoint { body: instruction* } .
```

Whether angle and value appear is decided by op (see section 3.1); unused ones serialize as null. An open module's body is ∅ and it must not declare locals (see section 4.5). Modules of versions "0.1" and "0.2" carry no locals field.

Main numeric ranges: RegType.width lies in 0..64; the two QRAM widths lie in 1..64; Repeat.count lies in 0..2^63−1; Control.value lies within the unsigned range of the control view; the value of add_const lies in 0..2^width−1; module call depth and structural-block nesting depth do not exceed 127.

The text encoding maps each record to an object carrying "tag", with field names identical to the record fields; tuples map to arrays, and scalars and null pass through unchanged:

```text
obj(T, f1: v1, …, fn: vn) = { "tag": T, "f1": enc(v1), …, "fn": enc(vn) }
enc((e1, …, ek)) = [ enc(e1), …, enc(ek) ]
enc(∅)           = null
enc(scalar)      = scalar
```

Decoding requires the field set to match the tag's record exactly and rejects unknown tags, extra or missing fields, duplicate keys, non-finite numbers, and unknown versions (see section 4.2). Canonical output is UTF-8, key-sorted, two-space indented, with one trailing newline; YAML uses block style and JSON uses braced objects.

### 4.4 Backend lowering rules

The OriginIR-ext backend may lower register operations into physical bit operations, but it must preserve module definitions and calls. QRAM resources specialize modules according to the actual binding. Repeat can be converted into a shared auxiliary module graph. Load lowers to a QRAM query line with the resource name as the operand word; Store lowers to the extension line `QRAMWRITE <resource name> <address bits>, <data bits>` and does not participate in gate-level counting. Downstream text executors (UnifiedQuantum, PySparQ) do not yet accept runtime writes; when they meet a program containing Store they must raise an error at the execution entry. Text export is unaffected.

Modules may be traversed or expanded according to an executor's capabilities only while actually running a concrete downstream executor. The downstream side may flatten on its own, but such processing must not change RIR in return or replace RIR's structural serialization.

The PySparQ backend maps non-empty entry registers to native named integer registers. Temporary rearrangements or copies performed for a view must restore their temporary space and must not modify RIR's root-register interpretation.

Backends may impose execution budgets stricter than the IR, for example state-vector qubit counts, QRAM materialization lengths, and expansion counts. When a limit is hit they must raise an error; silently truncating Repeat, qubits, or memory contents is not allowed.

### 4.5 Module-private work registers

Module.locals is an ordered array of Registers that is not part of the public call signature. Each call borrows them from the zero state and must uncompute them before returning; the IR checks only widths and references, uncomputation is an implementation obligation, and the simulator provides a runtime check. Open modules must not declare locals. Adjoint and Control include the complete module behavior, and workspaces cannot escape across calls. The OriginIR export renders them as module work parameters, and sequential calls reuse one stretch of physical workspace; PySparQ can intercept native implementations at module boundaries, skipping their internal workspaces and decompositions. The native registry is not part of the IR, and a native executable must not be misreported as gate-level closure.


## Part 5: Complete examples

The six examples below are complete programs that run on their own. Each example first gives the full generation code (with line-by-line comments), then the canonical YAML output of `dumps()` — **byte for byte from the real serializer**, unabridged — and finally a field-by-field explanation against the rules of the first four parts. Read example 1 first for the overall shape, then skip around by feature.

### Example 1: minimal unitary program (Bell pair)

```python
from oracq import Bits, Builder, dumps

# Declare the module: the public interface is one two-bit bits register named pair.
b = Builder("bell_pair", {"pair": Bits(2)})
# Broadcast an H gate over the least significant bit (index 0).
b.h(b["pair"][0])
# Bit-wise CNOT: source pair[0], target pair[1], producing a Bell state.
b.xor(b["pair"][0], b["pair"][1])
# Emit canonical YAML; the YAML below is the byte-for-byte result of this call.
print(dumps(b.finish().program()))
```

```yaml
entry: bell_pair
modules:
  - attributes: []
    body:
      - angle: null
        op: h
        operands:
          - parts:
              - register: pair
                start: 0
                tag: Span
                width: 1
            tag: Ref
            type:
              kind: bits
              tag: RegType
              width: 1
        tag: Primitive
        value: null
      - angle: null
        op: xor
        operands:
          - parts:
              - register: pair
                start: 0
                tag: Span
                width: 1
            tag: Ref
            type:
              kind: bits
              tag: RegType
              width: 1
          - parts:
              - register: pair
                start: 1
                tag: Span
                width: 1
            tag: Ref
            type:
              kind: bits
              tag: RegType
              width: 1
        tag: Primitive
        value: null
    locals: []
    name: bell_pair
    registers:
      - name: pair
        tag: Register
        type:
          kind: bits
          tag: RegType
          width: 2
    resources: []
    tag: Module
tag: Program
version: '0.1'

```

Explanation:

- {obj}`Program <oracq.infrastructure.ir.Program>` has only three fields: `entry` points at the sole module; `modules` is output sorted by module name (section 4.2); `version` is `"0.1"`.
- The `registers` of the {obj}`Module <oracq.infrastructure.ir.Module>` are the public interface (`pair: bits/2`); `resources`, `locals`, and `attributes` must be written out even when empty.
- The first {obj}`Primitive <oracq.infrastructure.ir.Primitive>`: `op=h`, the operand is a single {obj}`Ref <oracq.infrastructure.ir.Ref>` whose {obj}`Span(pair, 0, 1) <oracq.infrastructure.ir.Span>` is the lowest bit of the root register. Slicing produces a bits interpretation (section 2.3), so `Ref.type.kind` is `bits`.
- The second `Primitive`: `op=xor` with two equal-width operands — the source is bit 0 and the target is bit 1 of `pair`; the semantics is target XOR-equals source (section 3.1).
- The unused `angle` and `value` of both instructions are explicitly `null`: canonical output requires all fields to appear (section 4.2).

### Example 2: slicing, fusing, and type reinterpretation

```python
from oracq import Bits, Builder, UInt, fuse, dumps

b = Builder("views", {"x": UInt(4), "y": Bits(2)})
# Slicing x[1:3] yields a bits view; add_const needs a uint,
# so reinterpret explicitly. The addition carries modulo 4 on the 2-bit
# view and leaves the high two bits of x untouched.
b.add_const(b["x"][1:3].reinterpret("uint"), 3)
# Fuse across root registers: y occupies the low bits and the low two bits
# of x the high bits, then reinterpret as sint.
word = fuse(b["y"], b["x"][:2]).reinterpret("sint")
# Broadcast H over the 4-bit view.
b.h(word)
print(dumps(b.finish().program()))
```

```yaml
entry: views
modules:
  - attributes: []
    body:
      - angle: null
        op: add_const
        operands:
          - parts:
              - register: x
                start: 1
                tag: Span
                width: 2
            tag: Ref
            type:
              kind: uint
              tag: RegType
              width: 2
        tag: Primitive
        value: 3
      - angle: null
        op: h
        operands:
          - parts:
              - register: y
                start: 0
                tag: Span
                width: 2
              - register: x
                start: 0
                tag: Span
                width: 2
            tag: Ref
            type:
              kind: sint
              tag: RegType
              width: 4
        tag: Primitive
        value: null
    locals: []
    name: views
    registers:
      - name: x
        tag: Register
        type:
          kind: uint
          tag: RegType
          width: 4
      - name: y
        tag: Register
        type:
          kind: bits
          tag: RegType
          width: 2
    resources: []
    tag: Module
tag: Program
version: '0.1'

```

Explanation:

- The `add_const` operand `Ref` has a single `Span(x, start=1, width=2)` — the middle two bits of `x`. `Ref.type.kind` is `uint`: slicing naturally yields `bits`, and `uint` is the result of an explicit `reinterpret`. The addition carries modulo 4 on the 2-bit view, and the high two bits of `x` are unaffected; `value=3` lies within `0..2^2-1` (section 3.1).
- The `h` operand `Ref` holds two `Span`s: first `y(0..2)`, then `x(0..2)`. `parts` are concatenated from low bits to high bits (section 2.3), so `y` occupies the low bits and the low two bits of `x` the high bits; the sum of the segment widths, 2+2, must equal `type.width=4`.
- `kind: sint` likewise comes from `reinterpret`: it changes only the numeric interpretation, no qubit and no bit width.

### Example 3: module call and QRAM resource binding

```python
from oracq import Bits, Builder, QRAM, dumps

# The called module: declares the QRAM formal resource table(2,3);
# its body is one Load.
lookup = Builder("lookup", {"address": Bits(2), "data": Bits(3)},
                 resources={"table": QRAM(2, 3)})
lookup.qram("table", lookup["address"], lookup["data"])
# The calling module: declares its own actual resource mem(2,3); after
# pushing the address into superposition it calls lookup,
# binding the formal resource table to the actual resource mem.
demo = Builder("demo", {"address": Bits(2), "data": Bits(3)},
               resources={"mem": QRAM(2, 3)})
demo.h(demo["address"])
demo.call(lookup.finish(), address=demo["address"], data=demo["data"],
          resources={"table": "mem"})
print(dumps(demo.finish().program()))
```

```yaml
entry: demo
modules:
  - attributes: []
    body:
      - angle: null
        op: h
        operands:
          - parts:
              - register: address
                start: 0
                tag: Span
                width: 2
            tag: Ref
            type:
              kind: bits
              tag: RegType
              width: 2
        tag: Primitive
        value: null
      - arguments:
          - parts:
              - register: address
                start: 0
                tag: Span
                width: 2
            tag: Ref
            type:
              kind: bits
              tag: RegType
              width: 2
          - parts:
              - register: data
                start: 0
                tag: Span
                width: 3
            tag: Ref
            type:
              kind: bits
              tag: RegType
              width: 3
        module: lookup
        resources:
          - mem
        tag: Call
    locals: []
    name: demo
    registers:
      - name: address
        tag: Register
        type:
          kind: bits
          tag: RegType
          width: 2
      - name: data
        tag: Register
        type:
          kind: bits
          tag: RegType
          width: 3
    resources:
      - name: mem
        tag: Resource
        type:
          address_width: 2
          data_width: 3
          tag: QRAM
    tag: Module
  - attributes: []
    body:
      - address:
          parts:
            - register: address
              start: 0
              tag: Span
              width: 2
          tag: Ref
          type:
            kind: bits
            tag: RegType
            width: 2
        data:
          parts:
            - register: data
              start: 0
              tag: Span
              width: 3
          tag: Ref
          type:
            kind: bits
            tag: RegType
            width: 3
        resource: table
        tag: Load
    locals: []
    name: lookup
    registers:
      - name: address
        tag: Register
        type:
          kind: bits
          tag: RegType
          width: 2
      - name: data
        tag: Register
        type:
          kind: bits
          tag: RegType
          width: 3
    resources:
      - name: table
        tag: Resource
        type:
          address_width: 2
          data_width: 3
          tag: QRAM
    tag: Module
tag: Program
version: '0.1'

```

Explanation:

- The `Program` holds two modules sorted by name (`demo` before `lookup`); a {obj}`Call <oracq.infrastructure.ir.Call>` references the called module by name, and the definition is not copied into the call site (section 4.1).
- `lookup` declares the formal resource `table: QRAM(2,3)`; its body is a single {obj}`Load <oracq.infrastructure.ir.Load>`: `|address>|data> ↦ |address>|data XOR M[address]>` (section 3.2).
- `demo` declares its own actual resource `mem: QRAM(2,3)`. The `Call` node's `arguments` bind the argument views one by one in the callee's signature order (`address`, `data`); `resources: ["mem"]` binds `lookup`'s formal resource `table` to `demo`'s actual resource `mem` — the two QRAM types must match exactly (section 3.3).
- The `Load` inside `lookup` references the formal name `table`; the name replacement happens at the call site, while the data table itself is still left to be bound separately at execution time (section 3.2).
- The call is not expanded: `demo`'s body is just the broadcast `h` and one `Call`; `lookup`'s definition stays in the `modules` table as-is, shared by multiple call sites.

### Example 4: structured control and a private workspace

```python
from oracq import Bits, Builder, UInt, dumps

b = Builder("structured", {"word": UInt(4), "flag": Bits(1)})
# Private workspace: zero in, zero out; not part of the public call signature.
scratch = b.local("scratch", Bits(2))
b.h(b["word"])
# When flag equals 1, controlled execution of the three-step
# borrow–use–uncompute XOR.
with b.control(b["flag"], 1):
    b.xor(b["word"][:2], scratch)   # scratch ^= low two bits of word (borrow)
    b.xor(scratch, b["word"][2:])   # high two bits of word ^= scratch (use)
    b.xor(b["word"][:2], scratch)   # scratch ^= low two bits of word (uncompute back to zero)
# Statically repeat 3 times; each pass is the adjoint of rz(0.5).
with b.repeat(3):
    with b.adjoint():
        b.rz(b["word"][0], 0.5)
print(dumps(b.finish().program()))
```

```yaml
entry: structured
modules:
  - attributes: []
    body:
      - angle: null
        op: h
        operands:
          - parts:
              - register: word
                start: 0
                tag: Span
                width: 4
            tag: Ref
            type:
              kind: uint
              tag: RegType
              width: 4
        tag: Primitive
        value: null
      - body:
          - angle: null
            op: xor
            operands:
              - parts:
                  - register: word
                    start: 0
                    tag: Span
                    width: 2
                tag: Ref
                type:
                  kind: bits
                  tag: RegType
                  width: 2
              - parts:
                  - register: scratch
                    start: 0
                    tag: Span
                    width: 2
                tag: Ref
                type:
                  kind: bits
                  tag: RegType
                  width: 2
            tag: Primitive
            value: null
          - angle: null
            op: xor
            operands:
              - parts:
                  - register: scratch
                    start: 0
                    tag: Span
                    width: 2
                tag: Ref
                type:
                  kind: bits
                  tag: RegType
                  width: 2
              - parts:
                  - register: word
                    start: 2
                    tag: Span
                    width: 2
                tag: Ref
                type:
                  kind: bits
                  tag: RegType
                  width: 2
            tag: Primitive
            value: null
          - angle: null
            op: xor
            operands:
              - parts:
                  - register: word
                    start: 0
                    tag: Span
                    width: 2
                tag: Ref
                type:
                  kind: bits
                  tag: RegType
                  width: 2
              - parts:
                  - register: scratch
                    start: 0
                    tag: Span
                    width: 2
                tag: Ref
                type:
                  kind: bits
                  tag: RegType
                  width: 2
            tag: Primitive
            value: null
        register:
          parts:
            - register: flag
              start: 0
              tag: Span
              width: 1
          tag: Ref
          type:
            kind: bits
            tag: RegType
            width: 1
        tag: Control
        value: 1
      - body:
          - body:
              - angle: 0.5
                op: rz
                operands:
                  - parts:
                      - register: word
                        start: 0
                        tag: Span
                        width: 1
                    tag: Ref
                    type:
                      kind: bits
                      tag: RegType
                      width: 1
                tag: Primitive
                value: null
            tag: Adjoint
        count: 3
        tag: Repeat
    locals:
      - name: scratch
        tag: Register
        type:
          kind: bits
          tag: RegType
          width: 2
    name: structured
    registers:
      - name: word
        tag: Register
        type:
          kind: uint
          tag: RegType
          width: 4
      - name: flag
        tag: Register
        type:
          kind: bits
          tag: RegType
          width: 1
    resources: []
    tag: Module
tag: Program
version: '0.1'

```

Explanation:

- `locals` contains `scratch` (`bits/2`): a module-private workspace outside the public signature; each call borrows it from the zero state and must uncompute it before returning, and the executor checks this when the module returns (section 4.5).
- The {obj}`Control <oracq.infrastructure.ir.Control>` node: `register` is `flag` (`bits/1`) and `value=1` — the body executes when `flag` equals 1 and is the identity otherwise. The control bit is protected throughout the body: no operand in the body overlaps `flag` (section 3.5).
- The body's three `xor`s form the borrow–use–uncompute pattern: `scratch ^= word[:2]`; `word[2:] ^= scratch`; `scratch ^= word[:2]`. After the third one, `scratch` is back to zero.
- The {obj}`Repeat <oracq.infrastructure.ir.Repeat>` node has `count=3`; its body holds a single {obj}`Adjoint <oracq.infrastructure.ir.Adjoint>` whose body is one `rz(0.5)`. Structural blocks are stored nested as data, and serialization never copies the instruction body by count (sections 3.4 and 4.2); the `Adjoint` is semantically equivalent to `rz(-0.5)`, but the IR does not rewrite the body (section 3.6).
- The root register `word` is `uint/4`, yet the `xor` operand slices are `bits` — the slicing-yields-bits rule once more.

### Example 5: open declaration (an oracle slot)

```python
from oracq import Bits, dumps
from oracq.algorithms.input_model.oracles import declare

# An open declaration: an oracle slot with body null; paradigm database_xor,
# adjoint and controlled capabilities declared by default,
# implementation status unresolved.
op = declare("BooleanFunction", {"address": Bits(3), "data": Bits(1)},
             paradigm="database_xor")
print(dumps(op.program()))
```

```yaml
entry: BooleanFunction
modules:
  - attributes:
      - - implementation_status
        - unresolved
      - - oracle_paradigm
        - database_xor
      - - supports_adjoint
        - true
      - - supports_controlled
        - true
    body: null
    locals: []
    name: BooleanFunction
    registers:
      - name: address
        tag: Register
        type:
          kind: bits
          tag: RegType
          width: 3
      - name: data
        tag: Register
        type:
          kind: bits
          tag: RegType
          width: 1
    resources: []
    tag: Module
tag: Program
version: '0.1'

```

Explanation:

- `body` is `null`: this is an open module representing an oracle slot rather than an implemented circuit (sections 2.1 and 4.1).
- `attributes` are key-sorted: `oracle_paradigm="database_xor"` declares one of the nine named paradigms; `supports_adjoint` / `supports_controlled` are capability declarations that {obj}`bind <oracq.infrastructure.linking.bind>` verifies against the implementation when linking; `implementation_status="unresolved"` records the implementation status. Attributes do not change instruction semantics (section 2.1); the full rules for capabilities and binding are in [the open IR](open-ir.md).
- The register interface is declared as usual (`address: bits/3`, `data: bits/1`): callers use this slot through its signature even though no implementation exists yet.
- An open module must not declare `locals`; here `locals` is an empty list (section 4.5).

### Example 6: QRAM random write and read-back (pointer-style access)

```python
from oracq import Builder, QRAM, QMem, UInt, dumps, simulate

b = Builder("store_load", {"addr": UInt(2), "val": UInt(4)}, {"ram": QRAM(2, 4)})
mem = QMem(b, "ram")            # bind the resource ram as an array view
mem[b["addr"]].store(b["val"])  # random write: M[addr] := val
mem[b["addr"]].load(b["val"])   # XOR-Load: val ^= M[addr]
print(dumps(b.finish().program()))
print(simulate(b.finish().program(), {"ram": [0, 0, 0, 0]}, initial={"addr": 2, "val": 13}).amplitudes)
```

When the address expression is exactly a single full-width register with no constant component, {obj}`QMem <oracq.infrastructure.qmem.QMem>` uses that register directly as the Load/Store address and introduces no addressing arithmetic. `dumps()` outputs:

```yaml
entry: store_load
modules:
  - attributes: []
    body:
      - address:
          parts:
            - register: addr
              start: 0
              tag: Span
              width: 2
          tag: Ref
          type:
            kind: uint
            tag: RegType
            width: 2
        data:
          parts:
            - register: val
              start: 0
              tag: Span
              width: 4
          tag: Ref
          type:
            kind: uint
            tag: RegType
            width: 4
        resource: ram
        tag: Store
      - address:
          parts:
            - register: addr
              start: 0
              tag: Span
              width: 2
          tag: Ref
          type:
            kind: uint
            tag: RegType
            width: 2
        data:
          parts:
            - register: val
              start: 0
              tag: Span
              width: 4
          tag: Ref
          type:
            kind: uint
            tag: RegType
            width: 4
        resource: ram
        tag: Load
    locals: []
    name: store_load
    registers:
      - name: addr
        tag: Register
        type:
          kind: uint
          tag: RegType
          width: 2
      - name: val
        tag: Register
        type:
          kind: uint
          tag: RegType
          width: 4
    resources:
      - name: ram
        tag: Resource
        type:
          address_width: 2
          data_width: 4
          tag: QRAM
    tag: Module
tag: Program
version: '0.1'

```

Explanation:

- The `body` is a {obj}`Store <oracq.infrastructure.ir.Store>` followed by a `Load`: the Store first assigns the classical cell `M[addr]` the value `val` (section 3.2), and the subsequent XOR-Load reads the new value back; the simulator outputs `(2, 13)`.
- The `address`/`data` widths of both instructions match the {obj}`QRAM(2, 4) <oracq.infrastructure.ir.QRAM>` declaration one by one; the Store sits at the top level of the module body, satisfying the structural constraint of never appearing inside a Control/Adjoint body (sections 3.5 and 3.6).
- A module containing Store has neither the `supports_adjoint` nor the `supports_controlled` capability; this example has no controlled or adjoint calls, so validation passes (section 4.1 and the open IR).
- `QMem`'s pointers, offsets, and multidimensional views are addressing sugar at the Python generation stage: once the address expressions are materialized as register arithmetic, what lands in the IR is still only Primitive, Load, and Store (Part 1, design positioning).
