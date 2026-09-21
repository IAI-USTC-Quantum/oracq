"QRAM 数据结构的指针式读写：C 风格基地址、偏移、多维视图与随机写。"

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import QRAM, Ref, UInt, ValidationError


@dataclass(frozen=True)
class _Addr:
    """地址表达式：constant + Σ(寄存器视图 × 正系数)，按 2^地址宽度取模。"""

    constant: int = 0
    terms: tuple[tuple[Ref, int], ...] = ()

    def plus_const(self, value: int) -> _Addr:
        """返回常量项平移 ``value`` 后的新地址表达式。

        Args:
            value: 加到常量项上的整数偏移，可正可负。

        Returns:
            _Addr: 常量项平移后的新地址表达式；原表达式不变。
        """
        return _Addr(self.constant + value, self.terms)

    def plus_ref(self, ref: Ref, coefficient: int = 1) -> _Addr:
        """返回追加一项 ``ref × coefficient`` 后的新地址表达式。

        Args:
            ref: 参与求地址和的寄存器视图。
            coefficient: 视图乘上的整数系数；缺省为 1。

        Returns:
            _Addr: 追加该项后的新地址表达式；原表达式不变。
        """
        return _Addr(self.constant, self.terms + ((ref, coefficient),))


def _fresh_local(builder: Builder, prefix: str, width: int) -> Ref:
    """分配首个空闲的 ``{prefix}_{i}`` ``UInt(width)`` 局部寄存器并返回其视图。"""
    if getattr(builder, "_closed", False):
        raise ValidationError("构造器已经结束")
    existing = {r.name for r in builder.locals} | {
        r.name for r in getattr(builder, "registers", ())
    }
    index = 1
    while f"{prefix}_{index}" in existing:
        index += 1
    return builder.local(f"{prefix}_{index}", UInt(width))


def _ripple_add(builder: Builder, target: Ref, source: Ref, carry: Ref | None) -> None:
    """target += source（模 2^w），source 保持不变；carry 为零初始工作区。

    逐位进位链：c_{i+1} = maj(target_i, source_i, c_i)，sum_i = target_i ⊕ source_i ⊕ c_i。
    进位位留在 carry 中，由整体 Adjoint 回放复净。
    """
    width = target.width
    for i in range(width):
        bit_t, bit_s = target[i], source[i]
        if carry is not None and i + 1 < width:
            with builder.control(bit_t, 1), builder.control(bit_s, 1):
                builder.x(carry[i])
            if i:
                with builder.control(bit_t, 1), builder.control(carry[i - 1], 1):
                    builder.x(carry[i])
                with builder.control(bit_s, 1), builder.control(carry[i - 1], 1):
                    builder.x(carry[i])
        builder.xor(bit_s, bit_t)
        if i:
            # i>0 蕴含 width>1，此时调用方保证 carry 已分配（span>1）。
            builder.xor(cast(Ref, carry)[i - 1], bit_t)


def _run_steps(
    builder: Builder, steps: list[tuple[str, Ref, Ref | None, Ref | int | None]]
) -> None:
    """按步骤序列合成地址：copy 为异或拷贝、add 为逐位加法、const 为常量加。"""
    for kind, target, source, extra in steps:
        if kind == "copy":
            builder.xor(cast(Ref, source), target)
        elif kind == "add":
            _ripple_add(builder, target, cast(Ref, source), cast("Ref | None", extra))
        else:
            builder.add_const(target, cast(int, extra))


class QMem:
    """QRAM 资源的多维数组视图：row-major 展平，下标可为整数或寄存器视图。"""

    def __init__(
        self,
        builder: Builder,
        resource: str,
        *,
        shape: tuple[int, ...] | None = None,
        _base: _Addr | None = None,
    ) -> None:
        """初始化 QRAM 数组视图。

        Args:
            builder: 已声明该 QRAM 资源的构造器。
            resource: QRAM 资源名。
            shape: row-major 形状元组；缺省为覆盖整个地址空间的一维数组。
            _base: 视图共享的基地址表达式；缺省从零地址开始。

        Raises:
            ValidationError: 资源未在构造器中声明、shape 非法或单元总数超过地址空间。
        """
        specs = {r.name: r.type for r in builder.resources}
        if resource not in specs:
            raise ValidationError(f"构造器没有声明 QRAM 资源：{resource}")
        self.builder: Builder = builder
        self.resource: str = resource
        self.spec: QRAM = specs[resource]
        if shape is None:
            shape = (1 << self.spec.address_width,)
        if not isinstance(shape, tuple) or not shape:
            raise ValidationError("shape 必须是非空正整数元组")
        cells = 1
        for dim in shape:
            if type(dim) is not int or dim <= 0:
                raise ValidationError("shape 必须是非空正整数元组")
            cells *= dim
        if cells > 1 << self.spec.address_width:
            raise ValidationError("shape 的单元总数超过 QRAM 地址空间")
        self.shape: tuple[int, ...] = shape
        self._base: _Addr = _base if _base is not None else _Addr()

    @property
    def strides(self) -> tuple[int, ...]:
        """row-major 各维跨步，即沿该维前进一个下标对应的地址增量。"""
        result = [1] * len(self.shape)
        for axis in range(len(self.shape) - 2, -1, -1):
            result[axis] = result[axis + 1] * self.shape[axis + 1]
        return tuple(result)

    @property
    def address_width(self) -> int:
        """底层 QRAM 资源的地址宽度。"""
        return self.spec.address_width

    @property
    def data_width(self) -> int:
        """底层 QRAM 资源的数据字宽度。"""
        return self.spec.data_width

    def ptr(self, base: int | Ref | None = None) -> QPtr:
        """返回指针；base 为整数偏移或持有地址的寄存器视图（量子指针）。

        Args:
            base: 整数偏移，或位宽不超过地址宽度的寄存器视图；省略时
                沿用当前视图的基地址。

        Returns:
            QPtr: 指向平移后基地址的量子内存指针。

        Raises:
            ValidationError: 基地址既不是整数、寄存器视图，或位宽超限。
        """
        if base is None:
            return QPtr(self, self._base)
        if isinstance(base, Ref):
            if base.width > self.spec.address_width:
                raise ValidationError("指针寄存器宽度超过地址宽度")
            return QPtr(self, self._base.plus_ref(base))
        if type(base) is int:
            return QPtr(self, self._base.plus_const(base))
        raise ValidationError("指针基地址必须是整数或寄存器视图")

    def __getitem__(
        self, key: int | slice | Ref | tuple[int | slice | Ref, ...]
    ) -> QMem | QPtr:
        """按多维下标访问视图。

        Args:
            key: 每维下标为整数、寄存器视图或切片；缺省维度按整维切片补全。

        Returns:
            QMem | QPtr: 含切片维时的子数组视图，或全部维度精确定位时的指针。

        Raises:
            ValidationError: 维数与 shape 不符、切片步长非 1、下标越界或类型非法。
        """
        if not isinstance(key, tuple):
            key = (key,)
        if len(key) < len(self.shape):
            key = key + (slice(None),) * (len(self.shape) - len(key))
        if len(key) != len(self.shape):
            raise ValidationError("下标维数与 shape 不符")
        addr, sliced, new_shape = self._base, False, list(self.shape)
        for axis, (index, size, stride) in enumerate(
            zip(key, self.shape, self.strides, strict=True)
        ):
            if isinstance(index, slice):
                if index.step not in (None, 1):
                    raise ValidationError("切片步长只支持 1")
                start = 0 if index.start is None else index.start
                stop = size if index.stop is None else index.stop
                if type(start) is not int or type(stop) is not int or not 0 <= start < stop <= size:
                    raise ValidationError("切片边界非法")
                addr = addr.plus_const(start * stride)
                new_shape[axis] = stop - start
                sliced = True
            elif isinstance(index, Ref):
                if (1 << index.width) > size:
                    raise ValidationError(f"第 {axis} 维量子下标宽度超过维度长度")
                addr = addr.plus_ref(index, stride)
            elif type(index) is int:
                if not 0 <= index < size:
                    raise ValidationError("下标越界")
                addr = addr.plus_const(index * stride)
            else:
                raise ValidationError("下标必须是整数、寄存器视图或切片")
        if sliced:
            return QMem(self.builder, self.resource, shape=tuple(new_shape), _base=addr)
        return QPtr(self, addr)

    def _emit(self, addr: _Addr, data: Ref, kind: str) -> None:
        """解引用：物化地址表达式后发射 Load（XOR 读）或 Store（随机写）。"""
        builder, spec = self.builder, self.spec
        if data.width != spec.data_width:
            raise ValidationError("数据寄存器宽度与 QRAM 数据宽度不符")
        address_width = spec.address_width
        constant = addr.constant % (1 << address_width)
        direct = (
            len(addr.terms) == 1
            and addr.terms[0][1] == 1
            and addr.terms[0][0].width == address_width
            and constant == 0
        )
        if direct:
            ref = addr.terms[0][0]
            if kind == "load":
                builder.qram(self.resource, ref, data)
            else:
                builder.store(self.resource, ref, data)
            return
        temp = _fresh_local(builder, "qmem_addr", address_width)
        steps: list[tuple[str, Ref, Ref | None, Ref | int | None]] = []
        occupied: list[tuple[int, int]] = []
        for ref, coefficient in addr.terms:
            if coefficient < 1:
                raise ValidationError("地址系数必须为正整数")
            for shift in range(coefficient.bit_length()):
                if not (coefficient >> shift) & 1:
                    continue
                span = min(ref.width, address_width - shift)
                if span <= 0:
                    break
                target, source = temp[shift : shift + span], ref[:span]
                # 位移片段落在尚未写过的区间时，temp 该段仍为零，XOR 即赋值，无需加法器；
                # 与已占用区间重叠（如指针 + 量子偏移的和）才合成逐位进位加法。
                if all(shift >= end or shift + span <= start for start, end in occupied):
                    steps.append(("copy", target, source, None))
                else:
                    carry = _fresh_local(builder, "qmem_carry", span - 1) if span > 1 else None
                    steps.append(("add", target, source, carry))
                occupied.append((shift, shift + span))
        if constant:
            steps.append(("const", temp, None, constant))
        _run_steps(builder, steps)
        if kind == "load":
            builder.qram(self.resource, temp, data)
        else:
            builder.store(self.resource, temp, data)
        if steps:
            with builder.adjoint():
                _run_steps(builder, steps)


class QPtr:
    """QRAM 指针：基地址 + 偏移，可解引用读（XOR-Load）或随机写（Store）。"""

    def __init__(self, mem: QMem, addr: _Addr) -> None:
        """绑定指针所属的内存视图与地址表达式。

        Args:
            mem: 指针指向的数组视图。
            addr: 指针持有的地址表达式。
        """
        self.mem: QMem = mem
        self.addr: _Addr = addr

    def __add__(self, other: int | Ref) -> QPtr:
        """返回沿地址前进 ``other`` 后的新指针。

        Args:
            other: 整数偏移或持有偏移的寄存器视图。

        Returns:
            QPtr: 偏移后的新指针，原指针不变。

        Raises:
            ValidationError: 偏移既不是整数也不是寄存器视图。
        """
        if type(other) is int:
            return QPtr(self.mem, self.addr.plus_const(other))
        if isinstance(other, Ref):
            return QPtr(self.mem, self.addr.plus_ref(other))
        raise ValidationError("指针偏移必须是整数或寄存器视图")

    def __radd__(self, other: int | Ref) -> QPtr:
        """支持 ``整数 + 指针`` 写法，转发到 ``__add__``。

        Args:
            other: 整数偏移或持有偏移的寄存器视图。

        Returns:
            QPtr: 偏移后的新指针。
        """
        return self.__add__(other)

    def __sub__(self, other: int) -> QPtr:
        """返回沿地址回退整数 ``other`` 后的新指针。

        Args:
            other: 整数偏移。

        Returns:
            QPtr: 基地址减去偏移后的新指针。

        Raises:
            ValidationError: 偏移不是整数。
        """
        if type(other) is not int:
            raise ValidationError("指针减法只支持整数偏移")
        return QPtr(self.mem, self.addr.plus_const(-other))

    def load(self, data: Ref) -> None:
        """XOR-Load：data ^= M[addr]；地址计算自动复净，叠加地址天然支持。

        Args:
            data: 与数据字等宽、异或进命中内容的视图。
        """
        self.mem._emit(self.addr, data, "load")

    def store(self, data: Ref) -> None:
        """随机写：M[addr] := d。执行时地址与数据须处于确定基矢；不计门成本。

        Args:
            data: 与数据字等宽、提供写入内容的视图。
        """
        self.mem._emit(self.addr, data, "store")
