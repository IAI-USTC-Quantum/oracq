"""数据库、稀疏访问和领域可逆计算的组装适配。"""

from __future__ import annotations

import math

from .builder import Builder
from .ir import Bits, ValidationError, fuse
from .library import BlockEncoding, _name
from .oracles import (
    SparseAccess,
    XorDatabase,
    annotate,
    declare,
    gate_state_prep,
    invoke,
    resources_for,
)


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
