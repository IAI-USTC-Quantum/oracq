"RIR v0.3：不可变、可序列化的寄存器级操作图。"

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

VERSION = "0.3"
"""当前 RIR 规范版本；新建 ``Program`` 的默认 ``version`` 值。"""


class ValidationError(ValueError):
    """输入违反 RIR 的结构或语义规则。"""


@dataclass(frozen=True)
class RegType:
    """寄存器的存储类型，由位模式解释与位宽构成。

    Attributes:
        kind: 位模式解释，取 ``bits``、``uint``、``sint`` 或 ``rational``。
        width: 位宽，范围为 0..64；零宽度表示空接口。
    """

    kind: str
    width: int


def Bits(width: int) -> RegType:
    """返回 ``bits`` 解释的 ``RegType``，表示不指定数值意义的位串。"""
    return RegType("bits", width)


def UInt(width: int) -> RegType:
    """返回 ``uint`` 解释的 ``RegType``，表示 0 到 2^width-1 的无符号整数。"""
    return RegType("uint", width)


def SInt(width: int) -> RegType:
    """返回 ``sint`` 解释的 ``RegType``，表示二补码整数；位操作仍作用于原始字。"""
    return RegType("sint", width)


def Rational(width: int) -> RegType:
    """返回 ``rational`` 解释的 ``RegType``，表示无符号字除以 2^width。"""
    return RegType("rational", width)


@dataclass(frozen=True)
class Register:
    """模块的形式量子寄存器声明。

    Attributes:
        name: 寄存器名，在模块内唯一且不与资源名冲突。
        type: 寄存器的存储类型。
    """

    name: str
    type: RegType


@dataclass(frozen=True)
class Span:
    """根寄存器内一段连续的量子位范围。

    Attributes:
        register: 根寄存器名。
        start: 范围起点的下标；下标零表示最低位。
        width: 范围的位宽。
    """

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
        """视图的总位宽，即 ``type`` 的位宽。"""
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
        """在不改变位宽和量子状态的前提下改变视图的解释。

        Args:
            kind: 目标解释，取 ``bits``、``uint``、``sint`` 或 ``rational``。

        Returns:
            Ref: 复用原 ``parts``、按新解释查看的视图。
        """
        return Ref(self.parts, RegType(kind, self.width))


def fuse(*refs: Ref) -> Ref:
    """按低位到高位拼接多个视图，生成可跨根寄存器的复合视图。

    Args:
        *refs: 依次占据低位到高位的视图序列。

    Returns:
        Ref: ``parts`` 依次连接、解释为 ``bits`` 的新视图。
    """
    return Ref(tuple(part for ref in refs for part in ref.parts), Bits(sum(r.width for r in refs)))


@dataclass(frozen=True)
class QRAM:
    """QRAM 资源的类型声明。

    内存是由地址到无符号数据字的映射，未指定单元为零；内容不写入程序
    JSON，作为执行输入另行绑定。

    Attributes:
        address_width: 地址字的位宽，范围为 1..64。
        data_width: 数据字的位宽，范围为 1..64。
    """

    address_width: int
    data_width: int


@dataclass(frozen=True)
class Resource:
    """模块的 QRAM 形式参数声明。

    Attributes:
        name: 资源名，在模块内唯一且不与寄存器名冲突。
        type: 资源的 ``QRAM`` 类型。
    """

    name: str
    type: QRAM


@dataclass(frozen=True)
class Primitive:
    """寄存器级基元指令。

    广播门在 IR 中保持为单条寄存器级操作，不因内部有多个物理门而拆分。

    Attributes:
        op: 门名，如 ``h``、``rx``、``phase``、``xor``、``swap``、``add_const``。
        operands: 操作数视图元组。
        angle: 旋转或相位角（弧度）；未使用的门为 ``None``。
        value: ``add_const`` 的模加数；未使用的门为 ``None``。
    """

    op: str
    operands: tuple[Ref, ...]
    angle: float | None = None
    value: int | None = None


@dataclass(frozen=True)
class Load:
    """QRAM 查询指令：``data`` 异或所寻址单元的内容。

    该指令自逆，不要求目标初始为零，也不修改地址或经典内存。

    Attributes:
        resource: 被查询的 QRAM 形式资源名。
        address: 地址视图，宽度须与资源声明一致。
        data: 数据视图，宽度须与资源声明一致，且与地址不重叠。
    """

    resource: str
    address: Ref
    data: Ref


@dataclass(frozen=True)
class Store:
    """QRAM 随机写指令：在执行时刻把经典单元赋值为数据视图。

    地址与数据视图须处于确定基矢；指令不改变任何量子位，只能出现在模块
    体或 ``Repeat`` 体内，``Control`` 与 ``Adjoint`` 体内禁止出现。

    Attributes:
        resource: 被写入的 QRAM 形式资源名。
        address: 地址视图，宽度须与资源声明一致。
        data: 数据视图，宽度须与资源声明一致，且与地址不重叠。
    """

    resource: str
    address: Ref
    data: Ref


@dataclass(frozen=True)
class Call:
    """模块调用指令；被调定义在 ``Program`` 中共享，不在调用点内联。

    Attributes:
        module: 被调模块名。
        arguments: 按被调签名顺序绑定的量子实参视图。
        resources: 按被调签名顺序绑定的资源实参名。
    """

    module: str
    arguments: tuple[Ref, ...]
    resources: tuple[str, ...] = ()


@dataclass(frozen=True)
class Repeat:
    """静态重复指令：顺序施加体指令共 ``count`` 次。

    生成与序列化阶段不按 ``count`` 复制指令体。

    Attributes:
        count: 重复次数，范围为 0..2^63-1；零次为恒等操作。
        body: 被重复的指令元组，即使 ``count`` 为零也必须结构合法。
    """

    count: int
    body: tuple[Instruction, ...]


@dataclass(frozen=True)
class Control:
    """相干控制指令：控制视图等于给定值时施加体指令，否则为恒等。

    控制寄存器作为量子条件参与相干控制，不被测量；控制位在整个 ``body``
    内受到保护。

    Attributes:
        register: 非空控制视图。
        value: 控制视图的无符号比较值。
        body: 受控执行的指令元组。
    """

    register: Ref
    value: int
    body: tuple[Instruction, ...]


@dataclass(frozen=True)
class Adjoint:
    """伴随指令：将体指令顺序反转，并对每条操作取伴随。

    创建时不改写或展开其主体；双重伴随恢复原操作。

    Attributes:
        body: 被取伴随的指令元组。
    """

    body: tuple[Instruction, ...]


Instruction: TypeAlias = Primitive | Load | Store | Call | Repeat | Control | Adjoint
"""RIR 指令节点的联合类型别名。"""


@dataclass(frozen=True)
class Module:
    """模块定义，由调用签名、QRAM 资源与有序指令体构成。

    Attributes:
        name: 模块名，在 ``Program`` 内唯一。
        registers: 有序公开量子参数列表。
        resources: 有序 QRAM 形式参数列表。
        body: 有序指令体；``None`` 表示开放声明，即尚未绑定实现的 oracle 槽位。
        attributes: 有序属性二元组，值限于字符串、整数、有限浮点数或布尔值。
        locals: 模块私有工作寄存器，零输入零输出，不属于公开调用签名。
    """

    name: str
    registers: tuple[Register, ...]
    resources: tuple[Resource, ...]
    body: tuple[Instruction, ...] | None
    attributes: tuple[tuple[str, str | int | float | bool], ...] = ()
    locals: tuple[Register, ...] = ()


@dataclass(frozen=True)
class Program:
    """完整程序：模块表加上入口模块名。

    Attributes:
        entry: 入口模块名，必须指向 ``modules`` 中已定义的模块。
        modules: 全部模块定义；定义在调用点间共享，不复制主体。
        version: RIR 版本字符串；读取器另接受历史版本 ``0.1`` 与 ``0.2``。
    """

    entry: str
    modules: tuple[Module, ...]
    version: str = VERSION

    @property
    def module_map(self) -> dict[str, Module]:
        """以模块名为键的模块定义索引。"""
        return {module.name: module for module in self.modules}

    @property
    def main(self) -> Module:
        """入口模块的定义。"""
        return self.module_map[self.entry]
