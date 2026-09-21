"""Deutsch–Jozsa、Bernstein–Vazirani 与 Simon 查询算法。"""

from __future__ import annotations

from collections.abc import Iterable

from pyqecclang.algorithms.input_model.contracts import positive_integer
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.algorithms.input_model.oracles import XorDatabase, annotate, invoke, resources_for
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, ValidationError


def deutsch_jozsa(function: XorDatabase) -> Operation:
    """生成 D-J 查询电路；调用者声明函数是常量或平衡函数。

    Args:
        function: 一位结果的 XOR database，须满足常量或平衡承诺。

    Returns:
        Operation: 寄存器 input 与 answer 的查询电路；读出 input 判定常量/平衡。
    """
    if function.data_width != 1:
        raise ValidationError("D-J oracle 必须只有一个结果位")
    b = Builder(
        _name("deutsch_jozsa", function.operation),
        {"input": Bits(function.address_width), "answer": Bits(1)},
        resources_for(("function", function.operation)),
        attributes={
            "algorithm": "deutsch_jozsa",
            "readout_register": "input",
            "validation_stage": "paradigm",
        },
    )
    b.x(b["answer"])
    b.h(b["answer"])
    b.h(b["input"])
    invoke(b, function.operation, "function", address=b["input"], data=b["answer"])
    b.h(b["input"])
    return b.finish()


def bernstein_vazirani(function: XorDatabase) -> Operation:
    """生成 Bernstein–Vazirani 秘密字符串读出电路。

    Args:
        function: 一位结果的 XOR database，调用者保证 f(x)=s·x XOR c。

    Returns:
        Operation: 公开 input 和 answer。满足输入前提时，读取 input 得到 s。

    Oracle 可以保留为开放声明，后续绑定 gate 或 QRAM 实现。"""
    from dataclasses import replace

    from pyqecclang.infrastructure.builder import Operation

    operation = deutsch_jozsa(function)
    attrs = dict(operation.module.attributes)
    attrs.update(algorithm="bernstein_vazirani", input_promise="f(x)=dot(s,x) XOR c over GF(2)")
    return Operation(
        replace(
            operation.module,
            name=_name("bernstein_vazirani", function.operation),
            attributes=tuple(sorted(attrs.items())),
        ),
        operation.dependencies,
    )


def affine_boolean_oracle(width: int, secret: int, *, bias: int = 0) -> XorDatabase:
    """构造 BV 的普通门 oracle，secret 的第 i 位对应地址第 i 位。

    Args:
        width: 地址寄存器位宽，取 1..64。
        secret: 秘密字符串的整数编码，取 0..2^width−1。
        bias: 常数偏置 c，取 0 或 1。

    Returns:
        XorDatabase: 实现 f(x)=s·x XOR c 的门级 XOR database 句柄。
    """
    positive_integer(width, "affine_boolean_oracle.width", maximum=64)
    positive_integer(secret, "affine_boolean_oracle.secret", minimum=0, maximum=(1 << width) - 1)
    positive_integer(bias, "affine_boolean_oracle.bias", minimum=0, maximum=1)
    b = Builder(
        _name("affine_boolean", width, secret, bias), {"address": Bits(width), "data": Bits(1)}
    )
    if bias:
        b.x(b["data"])
    for bit in range(width):
        if (secret >> bit) & 1:
            b.xor(b["address"][bit], b["data"])
    return XorDatabase(annotate(b.finish(), "database_xor"))


def simon_sample(function: XorDatabase) -> Operation:
    """生成一次 Simon 采样电路。

    Args:
        function: XOR database，调用者保证具有非零 XOR 周期的二对一映射。

    Returns:
        Operation: 公开 input 和 output。只读取 input 即可获得周期正交约束。

    输出寄存器无需在电路中提前测量；重复采样后的消元在经典侧进行。"""
    if not isinstance(function, XorDatabase):
        raise ValidationError("Simon 需要 XOR database")
    b = Builder(
        _name("simon_sample", function.operation),
        {"input": Bits(function.address_width), "output": Bits(function.data_width)},
        resources_for(("function", function.operation)),
        attributes={
            "algorithm": "simon_sample",
            "input_promise": "two-to-one XOR period",
            "readout_register": "input",
        },
    )
    b.h(b["input"])
    invoke(b, function.operation, "function", address=b["input"], data=b["output"])
    b.h(b["input"])
    return b.finish()


def simon_nullspace(samples: Iterable[int], width: int) -> tuple[int, ...]:
    """求 Simon 样本约束在 GF(2) 上的零空间。

    Args:
        samples: 按 little endian 编码的整数样本。
        width: 秘密字符串的位宽。

    Returns:
        tuple: 零空间基向量。只有一维零空间时才能直接确定非零周期；样本不足时保留多个基向量。"""
    positive_integer(width, "simon.width", maximum=64)
    rows = []
    for sample in samples:
        positive_integer(sample, "simon.sample", minimum=0, maximum=(1 << width) - 1)
        rows.append(sample)
    pivot_rows = {}
    cursor = 0
    for bit in range(width):
        candidate = next((i for i in range(cursor, len(rows)) if (rows[i] >> bit) & 1), None)
        if candidate is None:
            continue
        rows[cursor], rows[candidate] = rows[candidate], rows[cursor]
        for i in range(len(rows)):
            if i != cursor and (rows[i] >> bit) & 1:
                rows[i] ^= rows[cursor]
        pivot_rows[bit] = cursor
        cursor += 1
    basis = []
    for free in range(width):
        if free in pivot_rows:
            continue
        vector = 1 << free
        for pivot, row in pivot_rows.items():
            if (rows[row] >> free) & 1:
                vector |= 1 << pivot
        basis.append(vector)
    return tuple(basis)
