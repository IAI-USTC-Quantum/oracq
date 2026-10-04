# QRAM data structures: pointers, offsets, and random writes

**English** · <a href="../zh/manual/qmem.html">简体中文</a>

{obj}`QMem <oracq.infrastructure.qmem.QMem>` abstracts QRAM resources into C-style array access: base addresses, constant and quantum offsets, multidimensional views, and random reads and writes. All of the addressing is Python generation-stage sugar — what lands in the RIR is only register arithmetic, {obj}`Load <oracq.infrastructure.ir.Load>`, and {obj}`Store <oracq.infrastructure.ir.Store>` (see [the RIR specification](../reference/rir.md), section 3.2). The API is documented in [QRAM pointer-style reads and writes](../api/infrastructure/qmem.rst).

## Pointers and reads

{obj}`QMem <oracq.infrastructure.qmem.QMem>` binds a QRAM resource already declared in the builder. `ptr()` returns a pointer to cell 0; integer addition and subtraction produce constant offsets; passing a register view yields a pointer holding a quantum address. Dereferencing with `load` is the existing XOR-Load: `|addr⟩|d⟩ ↦ |addr⟩|d ⊕ M[addr]⟩`, with superposed addresses supported for free.

```{doctest}
>>> from oracq import Builder, QRAM, QMem, UInt, simulate
>>> b = Builder("demo", {"idx": UInt(2), "out": UInt(4)}, {"rom": QRAM(2, 4)})
>>> mem = QMem(b, "rom")
>>> (mem.ptr() + 1).load(b["out"])
>>> state = simulate(b.finish().program(), {"rom": [10, 11, 12, 13]})
>>> sorted(state.amplitudes.items())
[((0, 11), (1+0j))]
```

Quantum offsets (`p + b["idx"]`) and multidimensional flattening require register addition: `QMem` automatically synthesizes bit-by-bit carry adders and shift-concatenation; the temporary address registers are cleaned by a whole-block Adjoint after dereferencing, keeping the module unitary. The gate cost of address arithmetic is counted faithfully in [resource estimation](resource-estimation.md). When the address expression is exactly a single full-width register with no constant component, no addressing arithmetic is introduced.

## Multidimensional views

`shape` declares row-major flattening; subscripts may be integers, register views, or slices, and the width of a quantum subscript must satisfy 2^width not exceeding that dimension's length. A slice returns a re-addressed subarray view.

```{doctest}
>>> b = Builder("grid", {"row": UInt(2), "out": UInt(4)}, {"rom": QRAM(4, 4)})
>>> grid = QMem(b, "rom", shape=(4, 4))
>>> grid[b["row"], 2].load(b["out"])
>>> state = simulate(b.finish().program(), {"rom": list(range(16))}, initial={"row": 3})
>>> sorted(state.amplitudes.items())
[((3, 14), (1+0j))]
```

## Random writes

`store` is the RIR `Store` instruction: `M[addr] := data`. Memory cells are modeled classically and random writes carry no gate cost (counted separately as `qram_writes` in resource estimation). At execution the address and data registers must be in a definite basis state; a write under a superposed address has no linear semantics and the executor raises immediately. Structurally, `Store` cannot appear inside a {obj}`Control <oracq.infrastructure.ir.Control>` or {obj}`Adjoint <oracq.infrastructure.ir.Adjoint>` body, and modules containing `Store` have neither the controlled nor the adjoint capability.

```{doctest}
>>> b = Builder("write", {"addr": UInt(2), "val": UInt(4), "out": UInt(4)}, {"ram": QRAM(2, 4)})
>>> mem = QMem(b, "ram")
>>> mem[b["addr"]].store(b["val"])
>>> mem.ptr(3).load(b["out"])
>>> state = simulate(b.finish().program(), {"ram": [0, 0, 0, 0]}, initial={"addr": 1, "val": 9})
>>> sorted(state.amplitudes.items())
[((1, 9, 0), (1+0j))]
```

A subsequent `Load` reads the evolved memory: in the example above cell 1 now holds 9 while cell 3 is still 0. The [OriginIR-ext](backends.md#originir-ext) export lowers `Store` to the `QRAMWRITE` extension line; the text executors (UnifiedQuantum, PySparQ) do not yet accept runtime writes and raise at the execution entry when a program contains `Store`; text export is unaffected.

## Capability overview

| Capability | Form | Lowering result |
|---|---|---|
| Constant offset | `mem.ptr(2) + 3`, `p - 1` | `add_const` (modulo 2^address width) |
| Quantum pointer | `mem.ptr(b["addr"])` | direct addressing or shift-concatenation |
| Quantum offset | `p + b["idx"]` | when overlapping the pointer range, a carry adder is synthesized (Toffoli/CNOT) |
| Multidimensional indexing | `grid[b["row"], 3]` | shift-concatenation; no adder needed when the dimensions occupy disjoint bit segments |
| Subarray view | `grid[2:]` | constant re-addressing, no instructions |

The full semantic tests are in `tests/core/test_qmem.py`.
