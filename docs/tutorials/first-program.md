# Your first register program

**English** · [简体中文](../zh/tutorials/first-program.html)

This tutorial generates a Bell state. You will use one register, two gate operations, and the reference executor. Afterwards, the same program is saved as RIR and exported as OriginIR-ext.

## Defining the operations

{obj}`Bits(2) <oracq.infrastructure.ir.Bits>` defines a two-bit register. We first put the low bit into a superposition, then use it as the control of an XOR applied to the high bit.

```{testcode}
from oracq import Bits, Builder, simulate

# Declare a module named bell_pair whose public interface is a two-bit bits register named pair.
b = Builder("bell_pair", {"pair": Bits(2)})
# Broadcast an H gate over bit 0 (the least significant bit) of pair: |0> -> (|0>+|1>)/sqrt(2).
b.h(b["pair"][0])
# XOR (bitwise CNOT) with bit 0 as source and bit 1 as target: yields (|00>+|11>)/sqrt(2).
b.xor(b["pair"][0], b["pair"][1])
# Finish building and obtain an immutable Operation (module plus instruction body).
operation = b.finish()
# Take out its RIR Program: entry points to bell_pair, and modules holds exactly this one module.
program = operation.program()
# Simulate with the dependency-free reference executor; it returns sparse amplitudes.
state = simulate(program)
# Print the execution result; the output block below reproduces this print verbatim.
print(state.amplitudes)

# Amplitude keys are integer tuples formed by concatenating the entry registers in declaration order; (0,) is 00 and (3,) is 11.
assert set(state.amplitudes) == {(0,), (3,)}
# Both basis states should have magnitude 1/sqrt(2).
assert abs(state.amplitudes[(0,)] - 2**-0.5) < 1e-12
assert abs(state.amplitudes[(3,)] - 2**-0.5) < 1e-12
```

```{testoutput}
{(0,): (0.7071067811865475+0j), (3,): (0.7071067811865475+0j)}
```

{obj}`Builder <oracq.infrastructure.builder.Builder>` collects the calls above one by one, and `finish()` returns the immutable {obj}`Operation <oracq.infrastructure.builder.Operation>`. The printed dictionary is the execution result: keys are tuples of the entry registers' integer values, with `0` meaning `00` and `3` meaning `11`; `0.7071067811865475` is the double-precision rendering of `1/√2`.

## Saving and exporting

```{testcode}
from oracq import dumps, loads, export_originir

# Serialize to canonical YAML: keys sorted, two-space indentation, module and instruction structure preserved as is.
text = dumps(program)
# Print the whole text; the output block abbreviates the repeated middle section with ...
print(text)
# Deserialize back into a Program object; decoding is strict and re-runs semantic validation.
restored = loads(text)
# The round trip must be equal: register names, views, and module structure are all preserved.
assert restored == program
# Export as modular OriginIR-ext text.
artifact = export_originir(restored)
# The text should contain a DEF module definition: export neither inlines calls nor flattens them into a bare gate list.
assert "DEF" in artifact.text
# Print the OriginIR-ext text; it already ends with a trailing newline, hence end="".
print(artifact.text, end="")
```

```{testoutput}
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
                ...
              ...
        tag: Primitive
        value: null
      ...
    locals: []
    name: bell_pair
    ...
tag: Program
version: '0.3'

QINIT 2
CREG 0
DEF m_bell_pair_d029b74cecc4332d414ece8a(v_pair[2])
H v_pair[0]
CNOT v_pair[0], v_pair[1]
ENDDEF
m_bell_pair_d029b74cecc4332d414ece8a(q[0], q[1])
```

{obj}`dumps <oracq.infrastructure.serialization.dumps>` produces the text in the first half, and {obj}`loads <oracq.infrastructure.serialization.loads>` re-runs the semantic validation of {obj}`validate <oracq.infrastructure.validation.validate>` while decoding; the {obj}`OriginIRArtifact <oracq.infrastructure.backends.originir.OriginIRArtifact>` returned by {obj}`export_originir <oracq.infrastructure.backends.originir.export_originir>` carries the text in the second half. The first half is the skeleton of canonical YAML: the {obj}`Program <oracq.infrastructure.ir.Program>` holds the module table, a module's `body` holds the instructions, and an instruction's operands keep the register names and views in fields such as `register` and `start`; everything outside the `...` matches the actual print verbatim. The second half is OriginIR-ext text: `QINIT 2` declares a two-qubit register and `CREG 0` says there is no classical register; `DEF m_bell_pair_<fingerprint>` defines the module — the suffix is a fingerprint determined by the module's content, and a given module is exported with only one definition — and the last line binds the entry qubits `q[0], q[1]` to the module's formal parameters. The export does not first copy all calls into a flat gate list.

{obj}`simulate <oracq.infrastructure.execution.simulate>` is well suited to small-scale checks like this one. When a real backend is needed, keep using the same {obj}`Program <oracq.infrastructure.ir.Program>` and call {obj}`run_pysparq <oracq.infrastructure.backends.pysparq.run_pysparq>` or {obj}`run_originir <oracq.infrastructure.backends.originir.run_originir>` instead.

## Related pages

- Manual: [Operations, registers, and the generation process](../manual/concepts.md), [Source code structure](../manual/architecture.md), [Export and execution backends](../manual/backends.md)
- Specification: [RIR specification](../reference/rir.md) (JSON structure and the OriginIR-ext grammar)
- API reference: [RIR objects](../api/infrastructure/ir.rst), [module builder](../api/infrastructure/builder.rst), [register reference executor](../api/infrastructure/execution.rst), [RIR serialization](../api/infrastructure/serialization.rst)
- Continue with: [Replacing an algorithm's oracle](oracle-binding.md), [Running and modifying the algorithm gallery](gallery.md)
