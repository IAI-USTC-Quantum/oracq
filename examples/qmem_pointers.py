"""QRAM 指针式读写演示：偏移寻址、量子指针、多维视图与随机写。"""

from __future__ import annotations

from typing import cast

from pyqecclang import (
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
    """常量偏移与量子指针：读 M[base + off]。"""
    b = Builder(
        "pointer_loads",
        {"idx": UInt(2), "out": UInt(4), "out2": UInt(4)},
        {"rom": QRAM(2, 4)},
    )
    mem = QMem(b, "rom")
    (mem.ptr() + 1).load(b["out"])          # 常量偏移：地址 = 0 + 1
    (mem.ptr(b["idx"]) + 1).load(b["out2"])  # 量子指针：地址 = idx + 1，加法器自动合成
    return b.finish()


def grid_load() -> Operation:
    """多维视图：row-major 展平 grid[row, 3]。"""
    b = Builder("grid_load", {"row": UInt(2), "out": UInt(4)}, {"rom": QRAM(4, 4)})
    grid = QMem(b, "rom", shape=(4, 4))
    cast(QPtr, grid[b["row"], 3]).load(b["out"])  # 全维精确下标，运行时必为 QPtr
    return b.finish()


def store_roundtrip() -> Operation:
    """随机写后读回：M[addr] := val，再 XOR-Load 到独立寄存器。"""
    b = Builder(
        "store_roundtrip", {"addr": UInt(2), "val": UInt(4), "out": UInt(4)}, {"ram": QRAM(2, 4)}
    )
    mem = QMem(b, "ram")
    cast(QPtr, mem[b["addr"]]).store(b["val"])  # 全维精确下标，运行时必为 QPtr
    cast(QPtr, mem[b["addr"]]).load(b["out"])
    return b.finish()


def main() -> None:
    """运行三段 QRAM 指针示例并打印振幅与资源台账。"""
    loads = pointer_loads()
    state = simulate(loads.program(), {"rom": [10, 11, 12, 13]}, initial={"idx": 0})
    print("常量偏移读 M[1]，量子指针读 M[idx+1]（idx=0）：", sorted(state.amplitudes.items()))

    grid = grid_load()
    state = simulate(grid.program(), {"rom": list(range(16))}, initial={"row": 2})
    print("grid[2, 3]：", sorted(state.amplitudes.items()))

    writer = store_roundtrip()
    state = simulate(writer.program(), {"ram": [0, 0, 0, 0]}, initial={"addr": 2, "val": 13})
    print("store 后读回（M[2]=13）：", sorted(state.amplitudes.items()))

    estimate = estimate_resources(writer.program())
    print(
        "资源台账：qram_queries =", dict(estimate.qram_queries),
        "qram_writes =", dict(estimate.qram_writes),
        "gate_total =", estimate.gate_total,
    )
    print(dumps(writer.program()))


if __name__ == "__main__":
    main()
