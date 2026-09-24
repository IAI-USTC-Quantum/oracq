"""生成具有嵌套模块、QRAM 参数和寄存器视图的演示程序。"""

from __future__ import annotations

from oracq import QRAM, Builder, Operation, Program, UInt, dumps, export_originir, simulate


def make_lookup(address_width: int = 2, data_width: int = 3) -> Operation:
    """构造单次 QRAM 查表模块：``address`` 进、``data`` 出。"""
    b = Builder(
        f"lookup_{address_width}_{data_width}",
        {"address": UInt(address_width), "data": UInt(data_width)},
        {"table": QRAM(address_width, data_width)},
    )
    b.qram("table", b["address"], b["data"])
    return b.finish()


def make_program() -> Program:
    """构造两层模块调用包装的查表演示程序。"""
    lookup = make_lookup()
    wrapper = Builder("wrapped_lookup", {"a": UInt(2), "d": UInt(3)}, {"memory": QRAM(2, 3)})
    wrapper.call(lookup, address=wrapper["a"], data=wrapper["d"], resources={"table": "memory"})
    wrapped = wrapper.finish()
    b = Builder("demo", {"address": UInt(2), "data": UInt(3)}, {"values": QRAM(2, 3)})
    b.h(b["address"])
    b.x(b["data"][0])
    b.call(wrapped, a=b["address"], d=b["data"], resources={"memory": "values"})
    return b.finish().program()


if __name__ == "__main__":
    p = make_program()
    print(export_originir(p).text)
    print(simulate(p, {"values": [1, 2, 4, 7]}).amplitudes)
    print(dumps(p))
