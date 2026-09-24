"""稀疏位置与元素访问，以及数值字到 block encoding 的显式适配。"""

from __future__ import annotations

import math
from collections.abc import Mapping
from functools import lru_cache

from oracq.algorithms.common.arithmetic import BooleanNetwork, FixedFormat
from oracq.algorithms.input_model.block_encoding import reflect_zero
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import (
    SparseAccess,
    XorDatabase,
    annotate,
    declare,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Ref, RegType, ValidationError, fuse


def reversible_lookup(
    inputs: Mapping[str, int],
    outputs: Mapping[str, int],
    database: XorDatabase,
    *,
    name: str | None = None,
) -> Operation:
    """把 ``XorDatabase`` 的 address/data 接口适配到命名的输入/输出寄存器组。

    输入寄存器按声明顺序融合为地址视图，输出寄存器融合为数据视图，查询语义
    仍是 XOR 数据库的 ``data ^= memory[address]``。

    Args:
        inputs: 名字到位宽的映射，总宽度须等于 ``database.address_width``。
        outputs: 名字到位宽的映射，总宽度须等于 ``database.data_width``。
        database: 被适配的 XOR 数据库。
        name: 可选模块名，省略时由内容生成。

    Returns:
        Operation: 寄存器为 inputs 与 outputs 的并集，范式为 ``reversible_function``。

    Raises:
        ValidationError: 输入/输出总宽度与数据库不符，或两边名字冲突。
    """
    if (
        sum(inputs.values()) != database.address_width
        or sum(outputs.values()) != database.data_width
    ):
        raise ValidationError("可逆查询的输入/输出总宽度不符")
    if set(inputs) & set(outputs):
        raise ValidationError("可逆查询的输入与输出名字冲突")
    b = Builder(
        name
        or _name(
            "reversible_lookup", database.operation, tuple(inputs.items()), tuple(outputs.items())
        ),
        {k: Bits(v) for k, v in {**inputs, **outputs}.items()},
        resources_for(("db", database.operation)),
    )
    invoke(
        b,
        database.operation,
        "db",
        address=fuse(*(b[k] for k in inputs)),
        data=fuse(*(b[k] for k in outputs)),
    )
    return annotate(b.finish(), "reversible_function")


def word_rotation(
    value_width: int, *, scale: float | None = None, name: str | None = None
) -> Operation:
    """把数值字的整数值线性转成 ``amplitude`` 比特上 Ry 角度的可逆转导。

    逐位受控叠加后总旋转角为 ``scale * value``（value 取无符号整数）；
    缺省 ``scale = 2*pi / 2**value_width``，使整个取值域恰好转满一周。

    Args:
        value_width: value 寄存器位宽。
        scale: 每单位整数值的旋转角；省略时取 ``2*pi / 2**value_width``。
        name: 可选模块名，省略时由内容生成。

    Returns:
        Operation: 寄存器为 value 与单比特 amplitude，范式为 ``reversible_function``，
        模块属性 ``angle_scale`` 记录所用 scale。
    """
    scale = 2 * math.pi / (1 << value_width) if scale is None else scale
    b = Builder(
        name or _name("word_rotation", value_width, scale),
        {"value": Bits(value_width), "amplitude": Bits(1)},
    )
    for bit in range(value_width):
        with b.control(b["value"][bit]):
            b.ry(b["amplitude"], scale * (1 << bit))
    return annotate(b.finish(), "reversible_function", angle_scale=scale)


def sparse_block_encoding(
    access: SparseAccess, transducer: Operation | None = None, *, alpha: float | None = None
) -> BlockEncoding:
    """历史候选，仅供旧目录描述；正式稀疏适配见 real_symmetric_sparse_encoding。

    Args:
        access: CKS 稀疏访问束。
        transducer: 值到幅度的转导操作；缺省为开放声明。
        alpha: 显式覆盖的归一化常数；缺省取稀疏度。

    Returns:
        BlockEncoding: 遗留转导构造、契约未指定的块编码。
    """
    n, v = access.width, access.value_width
    lw = next(r.type.width for r in access.location.module.registers if r.name == "work")
    transducer = transducer or declare(
        _name("sparse_amplitude", n, v),
        {"value": Bits(v), "amplitude": Bits(1)},
        paradigm="reversible_function",
    )
    uniform = gate_state_prep([1 if i < access.sparsity else 0 for i in range(1 << n)])
    width = n + lw + v + 1
    b = Builder(
        _name("sparse_prepare", access.location, access.entry, transducer, access.sparsity),
        {"target": Bits(n), "signal": Bits(width)},
        resources_for(
            ("position", access.location), ("entry", access.entry), ("amplitude", transducer)
        ),
        attributes={"algorithm": "sparse_isometry_extension", "validation_stage": "paradigm"},
    )
    neighbor, work = b["signal"][:n], b["signal"][n : n + lw]
    value, amplitude = b["signal"][n + lw : n + lw + v], b["signal"][width - 1 :]
    invoke(b, uniform.operation, target=neighbor, work=neighbor[:0])
    invoke(b, access.location, "position", column=b["target"], index=neighbor, work=work)
    invoke(b, access.entry, "entry", row=neighbor, column=b["target"], data=value)
    invoke(b, transducer, "amplitude", value=value, amplitude=amplitude)
    with b.adjoint():
        invoke(b, access.entry, "entry", row=neighbor, column=b["target"], data=value)
    prep = b.finish()
    out = Builder(
        _name("sparse_be", prep),
        {"target": Bits(n), "signal": Bits(width)},
        resources_for(("prep", prep)),
    )
    invoke(out, prep, "prep", target=out["target"], signal=out["signal"])
    out.swap(out["target"], out["signal"][:n])
    with out.adjoint():
        invoke(out, prep, "prep", target=out["target"], signal=out["signal"])
    return BlockEncoding(
        annotate(
            out.finish(),
            "block_encoding",
            be_alpha=float(access.sparsity if alpha is None else alpha),
            construction="Tdag_SWAP_T",
            validation_stage="paradigm",
            matrix_contract="legacy transduction unspecified; not a general sparse input adapter",
            legacy_input_model=True,
        )
    )


def batch_lookup(database: XorDatabase, count: int) -> Operation:
    """对同一 XOR 数据库做多路并发查询的批量线路。

    Args:
        database: XOR 数据库。
        count: 查询路数。

    Returns:
        Operation: 寄存器为 ``address{i}`` 与 ``data{i}``，i 从 0 到 count-1，
        位宽分别等于数据库的 address/data 宽度，逐路调用同一数据库操作。
    """
    registers: dict[str, RegType] = {}
    for i in range(count):
        registers[f"address{i}"] = Bits(database.address_width)
        registers[f"data{i}"] = Bits(database.data_width)
    b = Builder(
        _name("batch_lookup", database.operation, count),
        registers,
        resources_for(("db", database.operation)),
    )
    for i in range(count):
        invoke(b, database.operation, "db", address=b[f"address{i}"], data=b[f"data{i}"])
    return b.finish()


@lru_cache(maxsize=64)
def compare_words(width: int, kind: str = "eq") -> Operation:
    """两个字宽度无符号整数的相等或小于比较网络。

    Args:
        width: 每个输入字的位宽。
        kind: ``"eq"`` 生成相等判定，``"lt"`` 生成无符号小于判定。

    Returns:
        Operation: 输入寄存器 a、b，输出单比特 flag，比较成立时为 1。
        布尔网络经 compute/copy/uncompute 编译为可逆量子操作，私有 bank
        零进零出。结果按 ``(width, kind)`` 缓存复用。
    """
    net = BooleanNetwork()
    a, b = net.input("a", width), net.input("b", width)
    bit = (
        net.inv(net.any([net.xor(x, y) for x, y in zip(a, b, strict=True)]))
        if kind == "eq"
        else net.lt(a, b)
    )
    net.outputs = {"flag": [bit]}
    return net.operation()


@lru_cache(maxsize=64)
def value_transposition(width: int) -> Operation:
    """在 index 中交换 a/b 两个位模式；a、b 保留，适用于量子地址。

    Args:
        width: index 与 a、b 寄存器各自的位宽。

    Returns:
        Operation: 把 index 中等于 a 或 b 的基态互换、其余基态不变的操作。
    """
    net = BooleanNetwork()
    x, a, c = net.input("index", width), net.input("a", width), net.input("b", width)

    def eq(y: list[int]) -> int:
        """输出 1 当且仅当 ``y`` 与输入 ``x`` 逐位相等。"""
        return net.inv(net.any([net.xor(v, w) for v, w in zip(x, y, strict=True)]))

    net.outputs = {"flag": [net.or_(eq(a), eq(c))]}
    predicate = net.operation()
    b = Builder(
        "value_transposition_" + str(width),
        {"index": Bits(width), "a": Bits(width), "b": Bits(width)},
    )
    flag = b.local("membership", Bits(1))
    b.call(predicate, index=b["index"], a=b["a"], b=b["b"], flag=flag)
    with b.control(flag):
        b.xor(b["a"], b["index"])
        b.xor(b["b"], b["index"])
    b.call(predicate, index=b["index"], a=b["a"], b=b["b"], flag=flag)
    return b.finish()


def prefix_state(width: int, count: int) -> Operation:
    """在前 count 个基态上制备均匀叠加态。

    Args:
        width: target 寄存器位宽。
        count: 叠加覆盖的基态个数，范围为 ``1..2**width``。

    Returns:
        Operation: 单个 target 寄存器上的零输入制备，支撑集为 ``|0>`` 到
        ``|count-1>`` 且幅度相等。

    Raises:
        ValidationError: count 不在 ``1..2**width`` 范围内。
    """
    if not 1 <= count <= 1 << width:
        raise ValidationError("均匀前缀范围无效")
    b = Builder("uniform_prefix_" + str(width) + "_" + str(count), {"target": Bits(width)})

    def prepare(ref: Ref, size: int) -> None:
        """在 ``ref`` 的前 ``size`` 个基态上递归制备均匀叠加。"""
        if not ref.width:
            return
        if size == 1 << ref.width:
            b.h(ref)
            return
        half = 1 << (ref.width - 1)
        left, right = min(size, half), max(0, size - half)
        b.ry(ref[ref.width - 1], 2 * math.asin(math.sqrt(right / size)))
        with b.control(ref[ref.width - 1], 0):
            prepare(ref[: ref.width - 1], left)
        if right:
            with b.control(ref[ref.width - 1]):
                prepare(ref[: ref.width - 1], right)

    prepare(b["target"], count)
    return b.finish()


def magnitude_rotation(fmt: FixedFormat, amax: float) -> Operation:
    """数值字的小型普通实现；超过 12 位保留显式待绑定 transducer。

    Args:
        fmt: 元素值的定点格式。
        amax: 元素幅值上界，正有限实数。

    Returns:
        Operation: 幅度 ``sqrt(|value|/amax)`` 的受控 Ry 转导；格式超过
        12 位时返回待绑定的开放声明。
    """
    name = _name("sparse_sqrt_rotation", fmt, amax)
    registers = {"value": Bits(fmt.width), "amplitude": Bits(1)}
    attrs: dict[str, str | int | float | bool] = {
        "entry_bound": float(amax),
        "value_fraction": fmt.fraction,
        "amplitude_contract": "good amplitude sqrt(abs(value)/entry_bound)",
    }
    if fmt.width > 12:
        return declare(name, registers, paradigm="reversible_function", attributes=attrs)
    b = Builder(name, registers, attributes=attrs)
    for raw in range(1 << fmt.width):
        magnitude = abs(fmt.decode(raw)) / amax
        angle = 2 * math.acos(math.sqrt(min(1, magnitude)))
        if angle:
            with b.control(b["value"], raw):
                b.ry(b["amplitude"], angle)
    return b.finish()


def real_symmetric_sparse_encoding(
    access: SparseAccess,
    fmt: FixedFormat,
    amax: float,
    *,
    diagonal_nonnegative: bool = False,
    rotation: Operation | None = None,
) -> BlockEncoding:
    """CKS 型 T†ST：实 Hermitian、非负对角；交换两侧坐标及失败旗标。

    Args:
        access: CKS 稀疏访问束。
        fmt: 元素值的定点格式，位宽须与访问束的值宽一致。
        amax: 元素幅值上界，正有限实数。
        diagonal_nonnegative: 须为 True；当前仅支持非负对角的矩阵。
        rotation: 可选幅度转导操作；缺省由 magnitude_rotation 生成。

    Returns:
        BlockEncoding: T†ST 型实对称稀疏矩阵的块编码。
    """
    if not diagonal_nonnegative:
        raise ValidationError("当前对称稀疏适配要求非负对角；一般矩阵请显式 Hermitian dilation")
    if not math.isfinite(amax) or amax <= 0 or fmt.width != access.value_width:
        raise ValidationError("稀疏值格式/元素上界无效")
    n = access.width
    rotation = rotation or magnitude_rotation(fmt, amax)
    location_work = next(r.type.width for r in access.location.module.registers if r.name == "work")
    operands = (("location", access.location), ("entry", access.entry), ("rotation", rotation))
    b = Builder(
        _name("cks_isometry", access.location, access.entry, rotation, amax),
        {"target": Bits(n), "signal": Bits(n + 2)},
        resources_for(*operands),
        attributes={
            "algorithm": "cks_real_symmetric_isometry",
            "entry_bound": float(amax),
            "sparsity": access.sparsity,
            "correctness": "pending",
            "sign_convention": "negative offdiagonal phase pi only for target < neighbor",
        },
    )
    neighbor, _first_flag, second_flag = b["signal"][:n], b["signal"][n], b["signal"][n + 1]
    value = b.local("entry_value", Bits(fmt.width))
    work = b.local("location_work", Bits(location_work))
    order = b.local("ordered", Bits(1))
    b.call(prefix_state(n, access.sparsity), target=neighbor)
    invoke(b, access.location, "location", column=b["target"], index=neighbor, work=work)
    invoke(b, access.entry, "entry", row=neighbor, column=b["target"], data=value)
    invoke(b, rotation, "rotation", value=value, amplitude=second_flag)
    if fmt.signed:
        b.call(compare_words(n, "lt"), a=b["target"], b=neighbor, flag=order)
        with b.control(fuse(order, value[fmt.width - 1])):
            b.global_phase(math.pi)
        b.call(compare_words(n, "lt"), a=b["target"], b=neighbor, flag=order)
    invoke(b, access.entry, "entry", row=neighbor, column=b["target"], data=value)
    preparation = b.finish()
    out = Builder(
        _name("cks_sparse_be", preparation),
        {"target": Bits(n), "signal": Bits(n + 2)},
        resources_for(("prep", preparation)),
    )
    invoke(out, preparation, "prep", target=out["target"], signal=out["signal"])
    out.swap(out["target"], out["signal"][:n])
    out.swap(out["signal"][n], out["signal"][n + 1])
    with out.adjoint():
        invoke(out, preparation, "prep", target=out["target"], signal=out["signal"])
    return BlockEncoding(
        annotate(
            out.finish(),
            "block_encoding",
            be_alpha=access.sparsity * amax,
            construction="CKS_Tdag_S_T",
            self_adjoint_extension=True,
            correctness="pending",
        )
    )


def chebyshev_block(a: BlockEncoding, degree: int) -> BlockEncoding:
    """Chebyshev 行走幂：零信号块实现缩放矩阵的第 degree 阶 Chebyshev 多项式。

    阶数以 Repeat 保存，每步交替信号零态正反射与调用 ``a``；被编码矩阵按
    ``a.alpha`` 缩放后进入多项式。

    Args:
        a: 模块属性显式声明 ``self_adjoint_extension`` 的 ``BlockEncoding``。
        degree: Chebyshev 阶数，非负整数。

    Returns:
        BlockEncoding: ``be_alpha`` 为 1.0，``argument_scale`` 记录 ``a.alpha``，
        零信号块为 ``T_degree(A / a.alpha)``。

    Raises:
        ValidationError: 输入块编码缺少自伴酉扩张声明，或阶数为负。
    """
    if not dict(a.operation.module.attributes).get("self_adjoint_extension"):
        raise ValidationError("Chebyshev walk 需要显式的自伴酉扩张，普通 BE 不足以保证")
    if degree < 0:
        raise ValidationError("Chebyshev 阶数不能为负")
    b = Builder(
        _name("chebyshev_walk_power", a.operation, degree),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        resources_for(("a", a.operation)),
    )
    with b.repeat(degree):
        reflect_zero(b, b["signal"], positive=True)
        invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
    return BlockEncoding(
        annotate(
            b.finish(),
            "block_encoding",
            be_alpha=1.0,
            polynomial="Chebyshev",
            degree=degree,
            argument_scale=a.alpha,
            correctness="pending",
        )
    )
