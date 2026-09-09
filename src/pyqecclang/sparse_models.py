"""稀疏输入的显式适配；位置原地置换与元素 XOR 分开。"""

import math
from functools import lru_cache

from .arithmetic import BooleanNetwork
from .builder import Builder
from .combinators import reflect_zero
from .ir import Bits, ValidationError, fuse
from .library import BlockEncoding, _name
from .oracles import annotate, declare, invoke, resources_for


@lru_cache(maxsize=64)
def compare_words(width, kind="eq"):
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
def value_transposition(width):
    """在 index 中交换 a/b 两个位模式；a、b 保留，适用于量子地址。"""
    net = BooleanNetwork()
    x, a, c = net.input("index", width), net.input("a", width), net.input("b", width)

    def eq(y):
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


def prefix_state(width, count):
    if not 1 <= count <= 1 << width:
        raise ValidationError("均匀前缀范围无效")
    b = Builder("uniform_prefix_" + str(width) + "_" + str(count), {"target": Bits(width)})

    def prepare(ref, size):
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


def magnitude_rotation(fmt, amax):
    """数值字的小型普通实现；超过 12 位保留显式待绑定 transducer。"""
    name = _name("sparse_sqrt_rotation", fmt, amax)
    registers = {"value": Bits(fmt.width), "amplitude": Bits(1)}
    attrs = {
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


def real_symmetric_sparse_encoding(access, fmt, amax, *, diagonal_nonnegative=False, rotation=None):
    """CKS 型 T†ST：实 Hermitian、非负对角；交换两侧坐标及失败旗标。"""
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


def chebyshev_block(a, degree):
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
