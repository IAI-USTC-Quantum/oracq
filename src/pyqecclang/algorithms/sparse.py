"""稀疏位置与元素访问，以及数值字到 block encoding 的显式适配。"""

from __future__ import annotations

import math
from functools import lru_cache

from pyqecclang.algorithms.arithmetic import BooleanNetwork
from pyqecclang.algorithms.block_encoding import reflect_zero
from pyqecclang.algorithms.operators import BlockEncoding, _name
from pyqecclang.algorithms.oracles import (
    SparseAccess,
    XorDatabase,
    annotate,
    declare,
    gate_state_prep,
    invoke,
    resources_for,
)
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse


def reversible_lookup(inputs, outputs, database: XorDatabase, *, name=None):
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


def word_rotation(value_width, *, scale=None, name=None):
    scale = 2 * math.pi / (1 << value_width) if scale is None else scale
    b = Builder(
        name or _name("word_rotation", value_width, scale),
        {"value": Bits(value_width), "amplitude": Bits(1)},
    )
    for bit in range(value_width):
        with b.control(b["value"][bit]):
            b.ry(b["amplitude"], scale * (1 << bit))
    return annotate(b.finish(), "reversible_function", angle_scale=scale)


def sparse_block_encoding(access: SparseAccess, transducer=None, *, alpha=None):
    """历史候选，仅供旧目录描述；正式稀疏适配见 real_symmetric_sparse_encoding。"""
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


def batch_lookup(database, count):
    registers = {}
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
