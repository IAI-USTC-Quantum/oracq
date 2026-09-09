"版本化 JSON 编码；模块定义和调用边在往返过程中原样保留。"

from __future__ import annotations

import json
from dataclasses import fields, is_dataclass, replace

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
        ir.Call,
        ir.Repeat,
        ir.Control,
        ir.Adjoint,
        ir.Module,
        ir.Program,
    )
}


def encode(value):
    if is_dataclass(value):
        return {
            "tag": type(value).__name__,
            **{field.name: encode(getattr(value, field.name)) for field in fields(value)},
        }
    if isinstance(value, tuple):
        return [encode(item) for item in value]
    if value is None or type(value) in (int, float, str, bool):
        return value
    raise ir.ValidationError(f"不可序列化的值：{type(value).__name__}")


def decode(value):
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
        return cls(**{key: decode(value[key]) for key in expected})
    return value


def dumps(program: ir.Program) -> str:
    validate(program)
    canonical = replace(program, modules=tuple(sorted(program.modules, key=lambda m: m.name)))
    data = encode(canonical)
    if program.version in {"0.1", "0.2"}:
        for module in data["modules"]:
            module.pop("locals")
    return json.dumps(data, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2) + "\n"


def _unique(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ir.ValidationError(f"JSON 键重复：{key}")
        result[key] = value
    return result


def loads(text: str) -> ir.Program:
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
