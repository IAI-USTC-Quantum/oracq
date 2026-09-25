"""Demo of pointer-style QRAM reads and writes: offset addressing, quantum pointers, multi-dimensional views, and random writes."""

from __future__ import annotations

from typing import cast

from oracq import (
    QRAM,
    Builder,
    Operation,
    QMem,
    QPtr,
    UInt,
    dumps,
    estimate_resources,
    simulate,
)


def pointer_loads() -> Operation:
    """Constant offset and quantum pointer: read M[base + off]."""
    b = Builder(
        "pointer_loads",
        {"idx": UInt(2), "out": UInt(4), "out2": UInt(4)},
        {"rom": QRAM(2, 4)},
    )
    mem = QMem(b, "rom")
    (mem.ptr() + 1).load(b["out"])          # constant offset: address = 0 + 1
    (mem.ptr(b["idx"]) + 1).load(b["out2"])  # quantum pointer: address = idx + 1, adder synthesized automatically
    return b.finish()


def grid_load() -> Operation:
    """Multi-dimensional view: row-major flattening of grid[row, 3]."""
    b = Builder("grid_load", {"row": UInt(2), "out": UInt(4)}, {"rom": QRAM(4, 4)})
    grid = QMem(b, "rom", shape=(4, 4))
    cast(QPtr, grid[b["row"], 3]).load(b["out"])  # exact subscript in all dimensions; always a QPtr at runtime
    return b.finish()


def store_roundtrip() -> Operation:
    """Random write then read-back: M[addr] := val, then XOR-Load into an independent register."""
    b = Builder(
        "store_roundtrip", {"addr": UInt(2), "val": UInt(4), "out": UInt(4)}, {"ram": QRAM(2, 4)}
    )
    mem = QMem(b, "ram")
    cast(QPtr, mem[b["addr"]]).store(b["val"])  # exact subscript in all dimensions; always a QPtr at runtime
    cast(QPtr, mem[b["addr"]]).load(b["out"])
    return b.finish()


def main() -> None:
    """Run the three QRAM pointer examples and print amplitudes and the resource ledger."""
    loads = pointer_loads()
    state = simulate(loads.program(), {"rom": [10, 11, 12, 13]}, initial={"idx": 0})
    print("constant offset reads M[1], quantum pointer reads M[idx+1] (idx=0):", sorted(state.amplitudes.items()))

    grid = grid_load()
    state = simulate(grid.program(), {"rom": list(range(16))}, initial={"row": 2})
    print("grid[2, 3]:", sorted(state.amplitudes.items()))

    writer = store_roundtrip()
    state = simulate(writer.program(), {"ram": [0, 0, 0, 0]}, initial={"addr": 2, "val": 13})
    print("read-back after store (M[2]=13):", sorted(state.amplitudes.items()))

    estimate = estimate_resources(writer.program())
    print(
        "resource ledger: qram_queries =", dict(estimate.qram_queries),
        "qram_writes =", dict(estimate.qram_writes),
        "gate_total =", estimate.gate_total,
    )
    print(dumps(writer.program()))


if __name__ == "__main__":
    main()
