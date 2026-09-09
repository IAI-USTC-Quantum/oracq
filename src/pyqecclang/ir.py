"""RIR v0.3：不可变、可序列化的寄存器级操作图。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

VERSION = "0.3"


class ValidationError(ValueError):
    """输入违反 RIR 的结构或语义规则。"""


@dataclass(frozen=True)
class RegType:
    kind: str
    width: int


def Bits(width: int) -> RegType:
    return RegType("bits", width)


def UInt(width: int) -> RegType:
    return RegType("uint", width)


def SInt(width: int) -> RegType:
    return RegType("sint", width)


def Rational(width: int) -> RegType:
    return RegType("rational", width)


@dataclass(frozen=True)
class Register:
    name: str
    type: RegType


@dataclass(frozen=True)
class Span:
    register: str
    start: int
    width: int


@dataclass(frozen=True)
class Ref:
    """按低位到高位拼接的寄存器视图，不含物理量子位编号。"""

    parts: tuple[Span, ...]
    type: RegType

    @property
    def width(self) -> int:
        return self.type.width

    def __getitem__(self, key: int | slice) -> Ref:
        if isinstance(key, int):
            start, stop = key, key + 1
        elif isinstance(key, slice) and key.step in (None, 1):
            start = 0 if key.start is None else key.start
            stop = self.width if key.stop is None else key.stop
        else:
            raise ValidationError("视图只允许整数下标和连续切片")
        if type(start) is not int or type(stop) is not int or not 0 <= start <= stop <= self.width:
            raise ValidationError("寄存器切片越界")
        result = []
        offset = 0
        for part in self.parts:
            lo, hi = max(start, offset), min(stop, offset + part.width)
            if lo < hi:
                result.append(Span(part.register, part.start + lo - offset, hi - lo))
            offset += part.width
        return Ref(tuple(result), Bits(stop - start))

    def reinterpret(self, kind: str) -> Ref:
        return Ref(self.parts, RegType(kind, self.width))


def fuse(*refs: Ref) -> Ref:
    return Ref(tuple(part for ref in refs for part in ref.parts), Bits(sum(r.width for r in refs)))


@dataclass(frozen=True)
class QRAM:
    address_width: int
    data_width: int


@dataclass(frozen=True)
class Resource:
    name: str
    type: QRAM


@dataclass(frozen=True)
class Primitive:
    op: str
    operands: tuple[Ref, ...]
    angle: float | None = None
    value: int | None = None


@dataclass(frozen=True)
class Load:
    resource: str
    address: Ref
    data: Ref


@dataclass(frozen=True)
class Call:
    module: str
    arguments: tuple[Ref, ...]
    resources: tuple[str, ...] = ()


@dataclass(frozen=True)
class Repeat:
    count: int
    body: tuple[Instruction, ...]


@dataclass(frozen=True)
class Control:
    register: Ref
    value: int
    body: tuple[Instruction, ...]


@dataclass(frozen=True)
class Adjoint:
    body: tuple[Instruction, ...]


Instruction: TypeAlias = Primitive | Load | Call | Repeat | Control | Adjoint


@dataclass(frozen=True)
class Module:
    name: str
    registers: tuple[Register, ...]
    resources: tuple[Resource, ...]
    body: tuple[Instruction, ...] | None
    attributes: tuple[tuple[str, str | int | float | bool], ...] = ()
    locals: tuple[Register, ...] = ()


@dataclass(frozen=True)
class Program:
    entry: str
    modules: tuple[Module, ...]
    version: str = VERSION

    @property
    def module_map(self) -> dict[str, Module]:
        return {module.name: module for module in self.modules}

    @property
    def main(self) -> Module:
        return self.module_map[self.entry]
