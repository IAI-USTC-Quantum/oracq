"Versioned dual-format YAML/JSON encoding, YAML by default; module definitions and call edges survive the round trip unchanged."

from __future__ import annotations

import json
import math
import re
from collections.abc import Hashable
from dataclasses import fields, is_dataclass, replace
from typing import Literal, cast

import yaml
from yaml import MappingNode, Node
from yaml.constructor import ConstructorError

from oracq.infrastructure import ir
from oracq.infrastructure.validation import validate

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
"""Mapping from node ``tag`` values to ``ir`` dataclasses; ``decode`` restores node types from it."""


def encode(value: object) -> dict[str, object] | list[object] | int | float | str | bool | None:
    """Encode an RIR value into a structure compatible with both YAML and JSON.

    Dataclasses encode as objects carrying a ``tag`` field, tuples encode as
    arrays, and scalars and ``None`` pass through unchanged.

    Args:
        value: dataclass node, tuple or scalar to encode.

    Returns:
        text-format-compatible structure built from dict, list and scalars.

    Raises:
        ValidationError: the value belongs to a non-serializable type.
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
    raise ir.ValidationError(f"non-serializable value: {type(value).__name__}")


def decode(value: object) -> object:
    """Rebuild an encoded structure into an RIR dataclass tree.

    Args:
        value: structure produced by ``encode``; arrays are restored to tuples.

    Returns:
        rebuilt dataclass node, tuple or scalar.

    Raises:
        ValidationError: an unknown ``tag`` is encountered, or node fields are missing or extra.
    """
    if isinstance(value, list):
        return tuple(decode(item) for item in value)
    if isinstance(value, dict):
        tag = value.get("tag")
        if tag not in TYPES:
            raise ir.ValidationError(f"unknown node: {tag}")
        cls = TYPES[tag]
        expected = {field.name for field in fields(cls)}
        if set(value) != expected | {"tag"}:
            raise ir.ValidationError(f"missing or extra fields for {tag}")
        # Dispatch generically by tag to the dataclass constructors; mypy cannot statically verify **payload field types.
        return cls(**{key: decode(value[key]) for key in expected})  # type: ignore[arg-type]
    return value


_EXPONENT = re.compile(r"^[-+]?[0-9][0-9_]*[eE][-+]?[0-9]+$")
"""Scientific notation without a decimal point (e.g. ``1e+16``); not covered by PyYAML 1.1 implicit floats, so extra recognition is needed."""


class _StrictLoader(yaml.SafeLoader):
    """Strict YAML loader.

    Rejects duplicate keys and merge keys; beyond YAML 1.1 implicit floats it
    accepts scientific notation without a decimal point, ensuring that
    hand-written text of Python ``repr``-style floats can be restored to floats.
    """

    def construct_mapping(self, node: Node, deep: bool = False) -> dict[object, object]:
        if not isinstance(node, MappingNode):
            raise ConstructorError(
                None, None, f"expected a mapping node, got {type(node).__name__}", node.start_mark
            )
        mapping: dict[object, object] = {}
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            if not isinstance(key, Hashable):
                raise ConstructorError("while constructing a mapping", node.start_mark, "found unhashable key", key_node.start_mark)
            if key in mapping:
                raise ir.ValidationError(f"duplicate YAML key: {key}")
            mapping[key] = self.construct_object(value_node, deep=deep)
        return mapping


_StrictLoader.add_implicit_resolver("tag:yaml.org,2002:float", _EXPONENT, list("-+0123456789"))


class _CanonicalDumper(yaml.SafeDumper):
    """Canonical YAML dumper: disables anchors and aliases, indents nested sequences relative to their owning key."""

    def ignore_aliases(self, data: object) -> bool:
        return True

    def increase_indent(self, flow: bool = False, indentless: bool = True) -> None:
        super().increase_indent(flow, indentless=False)


def _check_values(value: object) -> None:
    """Recursively check that an encoded structure contains only RIR-allowed scalars and containers, with floats required to be finite.

    The YAML path uses this to reject implicit scalar types such as dates and
    ``.nan``/``.inf``; before YAML output it likewise stands in for the
    ``allow_nan=False`` semantics of ``json.dumps``.

    Raises:
        ValidationError: a non-finite float or a scalar type outside RIR appears.
    """
    if isinstance(value, dict):
        for item in value.values():
            _check_values(item)
    elif isinstance(value, list):
        for item in value:
            _check_values(item)
    elif isinstance(value, float):
        if not math.isfinite(value):
            raise ir.ValidationError(f"invalid number: {value}")
    elif not isinstance(value, (str, int, bool)) and value is not None:
        raise ir.ValidationError(f"invalid YAML scalar: {type(value).__name__}")


def dumps(program: ir.Program, *, format: Literal["yaml", "json"] = "yaml") -> str:
    """Serialize a program to canonical YAML (default) or canonical JSON text.

    Both formats use sorted keys, two-space indentation, UTF-8 and one trailing
    newline; module definitions are sorted by module name, with signature
    parameters and instruction order preserved. YAML uses block style, indents
    nested sequences relative to their owning key, and uses no anchors or
    aliases; floats must be finite. Open bodies stay null, and binding-origin
    attributes such as binding_captures are saved as-is; binding and resource
    analysis reports are not mixed into executable nodes.

    Args:
        program: program to serialize.
        format: output text format, ``yaml`` (default) or ``json``.

    Returns:
        str: canonical text.

    Raises:
        ValidationError: the program failed structural or semantic validation, or contains non-serializable values.
    """
    validate(program)
    canonical = replace(program, modules=tuple(sorted(program.modules, key=lambda m: m.name)))
    data = encode(canonical)
    if format == "json":
        return (
            json.dumps(data, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2) + "\n"
        )
    _check_values(data)
    return yaml.dump(
        data,
        Dumper=_CanonicalDumper,
        allow_unicode=True,
        sort_keys=True,
        default_flow_style=False,
        indent=2,
        width=1_000_000,
    )


def _unique(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Assemble a JSON object's key-value pair list into a dict, raising ``ValidationError`` on duplicate keys."""
    result = {}
    for key, value in pairs:
        if key in result:
            raise ir.ValidationError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def loads(text: str) -> ir.Program:
    """Parse YAML or JSON text and rebuild a validated program.

    The format is auto-detected from the text: text starting with ``{`` is
    first parsed as strict JSON, falling back to YAML on failure (YAML is a
    superset of JSON); other text is parsed as YAML. Both paths reject
    duplicate keys and non-finite numbers; the YAML path additionally rejects
    implicit scalar types such as dates.

    Args:
        text: YAML or JSON text produced by ``dumps``.

    Returns:
        Program: rebuilt program that passed semantic validation.

    Raises:
        ValidationError: the text syntax or node structure is invalid, the root node is not a ``Program``,
            or validation failed.
    """
    try:
        stripped = text.lstrip("\ufeff \t\r\n")
        if stripped[:1] == "{":
            try:
                data = json.loads(
                    stripped,
                    object_pairs_hook=_unique,
                    parse_constant=lambda s: (_ for _ in ()).throw(
                        ir.ValidationError(f"invalid number: {s}")
                    ),
                )
            except json.JSONDecodeError:
                data = yaml.load(stripped, Loader=_StrictLoader)
        else:
            data = yaml.load(stripped, Loader=_StrictLoader)
        _check_values(data)
        program = decode(data)
        if not isinstance(program, ir.Program):
            raise ir.ValidationError("the root node must be a Program")
        return validate(program)
    except (yaml.YAMLError, TypeError, KeyError, AttributeError, RecursionError) as exc:
        raise ir.ValidationError(f"invalid RIR YAML: {exc}") from exc
