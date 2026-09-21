"版本化 JSON 编码；模块定义和调用边在往返过程中原样保留。"

from __future__ import annotations

import json
from dataclasses import fields, is_dataclass, replace
from typing import cast

from pyqecclang.infrastructure import ir
from pyqecclang.infrastructure.validation import validate

TYPES = {
    cls.__name__: cls
    for cls in (
        ir.RegType,
        ir.Register,
        ir.Span,
        ir.Ref,
        ir.QRAM,
        ir.Resource,
        ir.Primitive,
        ir.Load,
        ir.Store,
        ir.Call,
        ir.Repeat,
        ir.Control,
        ir.Adjoint,
        ir.Module,
        ir.Program,
    )
}
"""JSON ``tag`` 到 ``ir`` 数据类的映射，``decode`` 据此还原节点类型。"""


def encode(value: object) -> dict[str, object] | list[object] | int | float | str | bool | None:
    """把 RIR 值编码为 JSON 兼容结构。

    数据类编码为携带 ``tag`` 字段的对象，元组编码为数组，标量与 ``None``
    原样传递。

    Args:
        value: 待编码的数据类节点、元组或标量。

    Returns:
        由 dict、list 和标量构成的 JSON 兼容结构。

    Raises:
        ValidationError: 值属于不可序列化的类型。
    """
    if is_dataclass(value):
        return {
            "tag": type(value).__name__,
            **{field.name: encode(getattr(value, field.name)) for field in fields(value)},
        }
    if isinstance(value, tuple):
        return [encode(item) for item in value]
    if value is None or type(value) in (int, float, str, bool):
        return cast("int | float | str | bool | None", value)
    raise ir.ValidationError(f"不可序列化的值：{type(value).__name__}")


def decode(value: object) -> object:
    """把 JSON 兼容结构重建为 RIR 数据类树。

    Args:
        value: ``encode`` 产生的结构；JSON 数组恢复为元组。

    Returns:
        重建的数据类节点、元组或标量。

    Raises:
        ValidationError: 遇到未知 ``tag``，或节点字段缺失、多余。
    """
    if isinstance(value, list):
        return tuple(decode(item) for item in value)
    if isinstance(value, dict):
        tag = value.get("tag")
        if tag not in TYPES:
            raise ir.ValidationError(f"未知 JSON 节点：{tag}")
        cls = TYPES[tag]
        expected = {field.name for field in fields(cls)}
        if set(value) != expected | {"tag"}:
            raise ir.ValidationError(f"{tag} 的字段缺失或多余")
        # 按 tag 泛型分发到各数据类构造器，mypy 无法静态验证 **payload 字段类型。
        return cls(**{key: decode(value[key]) for key in expected})  # type: ignore[arg-type]
    return value


def dumps(program: ir.Program) -> str:
    """把程序序列化为规范 JSON 文本。

    输出按键排序、两空格缩进、末尾一个换行；模块定义按模块名排序，签名
    参数与指令顺序保留。版本 ``0.1`` 与 ``0.2`` 的模块节点省略 ``locals``
    字段。开放主体保持为 null，binding_captures 等绑定来源属性原样保存；
    绑定与资源分析报告不混入可执行节点。

    Args:
        program: 待序列化的程序。

    Returns:
        str: 规范 JSON 文本。

    Raises:
        ValidationError: 程序未通过结构或语义验证，或含不可序列化的值。
    """
    validate(program)
    canonical = replace(program, modules=tuple(sorted(program.modules, key=lambda m: m.name)))
    data = encode(canonical)
    if program.version in {"0.1", "0.2"}:
        # encode(Program) 的规范形状已知：顶层为 dict，modules 为模块 dict 列表。
        for module in cast("list[dict[str, object]]", cast("dict[str, object]", data)["modules"]):
            module.pop("locals")
    return json.dumps(data, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2) + "\n"


def _unique(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """把 JSON 对象的键值对列表组装成字典，遇重复键即抛 ``ValidationError``。"""
    result = {}
    for key, value in pairs:
        if key in result:
            raise ir.ValidationError(f"JSON 键重复：{key}")
        result[key] = value
    return result


def loads(text: str) -> ir.Program:
    """解析 JSON 文本并重建通过验证的程序。

    版本 ``0.1`` 与 ``0.2`` 的模块节点自动补全空 ``locals``。解析阶段拒绝
    重复 JSON 键与非有限数值。

    Args:
        text: ``dumps`` 产生的 JSON 文本。

    Returns:
        Program: 重建并通过语义验证的程序。

    Raises:
        ValidationError: JSON 语法或节点结构非法、根节点不是 ``Program``，
            或验证未通过。
    """
    try:
        data = json.loads(
            text,
            object_pairs_hook=_unique,
            parse_constant=lambda s: (_ for _ in ()).throw(ir.ValidationError(f"非法数值：{s}")),
        )
        if isinstance(data, dict) and data.get("version") in {"0.1", "0.2"}:
            for module in data.get("modules", []):
                if isinstance(module, dict):
                    module.setdefault("locals", [])
        program = decode(data)
        if not isinstance(program, ir.Program):
            raise ir.ValidationError("JSON 根节点必须为 Program")
        return validate(program)
    except (TypeError, KeyError, AttributeError, RecursionError, json.JSONDecodeError) as exc:
        raise ir.ValidationError(f"非法 RIR JSON：{exc}") from exc
