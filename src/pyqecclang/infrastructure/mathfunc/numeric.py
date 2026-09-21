"固定点实数/复数数学分解与按路径传播的数值状态。"

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from typing import cast

from pyqecclang.algorithms.common.arithmetic import (
    BooleanNetwork,
    FixedFormat,
    fixed_arithmetic,
)
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Adjoint, Bits, Ref, ValidationError

BOUNDS = {
    "sin": (-math.pi, math.pi),
    "cos": (-math.pi, math.pi),
    "tan": (-1.0, 1.0),
    "exp": (-2.0, 2.0),
    "log": (0.25, 4.0),
    "log10": (0.25, 4.0),
    "asin": (-1.0, 1.0),
    "acos": (-1.0, 1.0),
    "atan": (-1.0, 1.0),
    "sinh": (-2.0, 2.0),
    "cosh": (-2.0, 2.0),
    "tanh": (-2.0, 2.0),
    "asinh": (-4.0, 4.0),
    "acosh": (1.0, 4.0),
    "atanh": (-0.75, 0.75),
}
"""各数学核的默认 Chebyshev 近似区间，函数名映射到 (下界, 上界)。"""


@dataclass(frozen=True)
class MathConfig:
    """数学核生成配置：多项式阶数与可替换的近似区间。

    Attributes:
        degree: Chebyshev 多项式阶数，范围为 1..32。
        intervals: 每项为 (函数名, 下界, 上界)，覆盖指定数学核的近似区间；
            未列出的数学核使用 BOUNDS 中的默认区间。

    Raises:
        ValidationError: 阶数不在 1..32，或区间引用未知函数、重复、含非有限值或下界不小于上界。
    """

    degree: int = 6
    intervals: tuple[tuple[str, float, float], ...] = ()

    def __post_init__(self) -> None:
        """校验多项式阶数范围与近似区间覆盖表。"""
        if type(self.degree) is not int or not 1 <= self.degree <= 32:
            raise ValidationError("数学核多项式阶数必须为 1..32")
        seen: set[str] = set()
        for kind, lo, hi in self.intervals:
            if (
                kind not in BOUNDS
                or kind in seen
                or not all(math.isfinite(x) for x in (lo, hi))
                or lo >= hi
            ):
                raise ValidationError("数学核近似区间无效或重复")
            seen.add(kind)

    def bounds(self, kind: str) -> tuple[float, float]:
        """查询指定数学核的近似区间。

        Args:
            kind: 数学函数名，如 sin 或 exp。

        Returns:
            tuple: ``(下界, 上界)``；优先取 ``intervals`` 中的覆盖项，否则取默认区间。

        Raises:
            KeyError: 名字既不在 ``intervals`` 也不在 BOUNDS 中。
        """
        return next(((lo, hi) for name, lo, hi in self.intervals if name == kind), BOUNDS[kind])


@dataclass(frozen=True)
class Numeric:
    """数学值在生成期的电路表示。

    Attributes:
        kind: 值类型 real、complex 或 bool。
        parts: 承载该值的寄存器引用 tuple；复数为 (实部, 虚部)，其余为一元。
        status: 2 位状态寄存器引用，第 0 位汇总定义域失效，第 1 位汇总值域/字长越界。
    """

    kind: str
    parts: tuple[Ref, ...]
    status: Ref


@lru_cache(maxsize=128)
def logic_operation(kind: str, width: int = 1) -> Operation:
    """构造逐位逻辑/选择运算的布尔网络。

    Args:
        kind: 运算名：not、flag、range_flag、and、or、xor 或 select。
        width: 输入寄存器位宽。

    Returns:
        Operation: 接口为 a、（按需的）b 与 select，输出 out；
        flag 与 range_flag 把条件位分别复制到输出的第 0 位与第 1 位。
    """
    net = BooleanNetwork()
    a = net.input("a", width)
    if kind == "not":
        out = [net.inv(x) for x in a]
    elif kind == "flag":
        out = [a[0], 0]
    elif kind == "range_flag":
        out = [0, a[0]]
    else:
        b = net.input("b", width)
        if kind == "select":
            select = net.input("select", 1)[0]
            out = net.mux(select, a, b)
        else:
            fn = {"and": net.and_, "or": net.or_, "xor": net.xor}[kind]
            out = [fn(x, y) for x, y in zip(a, b, strict=True)]
    net.outputs = {"out": out}
    return net.operation(attributes={"arithmetic_kind": "logic_" + kind})


class NumericEmitter:
    """在 Builder 上把 MIR 数值操作展开为定点可逆电路。

    以 Numeric 值为中介执行定点算术、逐位逻辑与数学核调用，逐值追踪
    status 位，并在收尾时复制输出、反算全部临时寄存器。

    Args:
        builder: 目标电路 Builder。
        fmt: 定点格式 FixedFormat，当前要求有符号。
        config: 数学核生成配置 MathConfig。

    Raises:
        ValidationError: fmt 不是有符号定点格式。
    """

    def __init__(self, builder: Builder, fmt: FixedFormat, config: MathConfig) -> None:
        """绑定构造器、定点格式与数学核配置，初始化常量缓存。"""
        if not fmt.signed:
            raise ValidationError("数学函数编译目前使用有符号定点格式")
        self.b: Builder = builder
        self.fmt: FixedFormat = fmt
        self.config: MathConfig = config
        self.counter: int = 0
        self.constants: dict[tuple[str, bool | int | float | complex], Numeric] = {}
        self.zero_status: Ref = self.local(2)

    def local(self, width: int | None = None) -> Ref:
        """分配一个数学临时寄存器。

        Args:
            width: 位宽；省略时使用定点格式全宽。

        Returns:
            Ref: 新分配的寄存器引用，名字按 ``math_tmp_`` 前缀递增编号。
        """
        self.counter += 1
        return self.b.local(
            "math_tmp_" + str(self.counter), Bits(self.fmt.width if width is None else width)
        )

    def invoke(self, op: Operation, arguments: dict[str, Ref]) -> None:
        """复制重叠的只读参数，满足 RIR 无别名的调用 ABI。

        Args:
            op: 待调用的已实现操作。
            arguments: 形参寄存器名到实参视图的映射；互相重叠的只读视图
                会被临时复制并在调用后复净。
        """
        seen: set[tuple[str, int]] = set()
        clones: list[tuple[Ref, Ref]] = []
        from pyqecclang.infrastructure.validation import locations

        actual: dict[str, Ref] = {}
        for name, ref in arguments.items():
            loc = set(locations(ref))
            if loc & seen:
                clone = self.local(ref.width)
                self.b.xor(ref, clone)
                clones.append((ref, clone))
                actual[name] = clone
            else:
                actual[name] = ref
            seen |= loc
        self.b.call(op, **actual)  # type: ignore[arg-type]
        for ref, clone in reversed(clones):
            self.b.xor(ref, clone)

    def logic(self, kind: str, *refs: Ref, select: Ref | None = None) -> Ref:
        """调用逐位逻辑/选择网络并返回结果寄存器。

        Args:
            kind: 运算名，见 logic_operation。
            refs: 参与运算的寄存器引用；二元运算取前两个。
            select: select 运算的条件位引用。

        Returns:
            Ref: 写入运算结果的新寄存器引用。
        """
        operation = logic_operation(kind, refs[0].width)
        width = next(r.type.width for r in operation.module.registers if r.name == "out")
        out = self.local(width)
        arguments = dict(zip(("a", "b"), refs, strict=False))
        if select is not None:
            arguments["select"] = select
        self.invoke(operation, {**arguments, "out": out})
        return out

    def flags(self, *refs: Ref) -> Ref:
        """合并若干 2 位状态寄存器。

        全零的 zero_status 不参与合并；没有其他输入时直接返回 zero_status。

        Args:
            refs: 状态寄存器引用。

        Returns:
            Ref: 各输入按位或得到的状态寄存器引用。
        """
        unique: list[Ref] = []
        for ref in refs:
            if ref != self.zero_status and ref not in unique:
                unique.append(ref)
        if not unique:
            return self.zero_status
        result = unique[0]
        for ref in unique[1:]:
            result = self.logic("or", result, ref)
        return result

    def constant(self, value: bool | int | float | complex, kind: str | None = None) -> Numeric:
        """把生成期数值固化为电路常量并缓存。

        复数拆为两个实常量并合并状态；实数超出格式的可表示幅值时置
        值域越界位。相同 (kind, value) 只编码一次。

        Args:
            value: bool、int、float 或 complex 常量。
            kind: 指定值类型；省略时按 Python 类型推导。

        Returns:
            Numeric: 常量的电路表示。
        """
        kind = kind or (
            "bool" if type(value) is bool else "complex" if type(value) is complex else "real"
        )
        key = (kind, value)
        if key in self.constants:
            return self.constants[key]
        if kind == "complex":
            re, im = self.constant(value.real), self.constant(value.imag)
            result = Numeric(
                "complex", (re.parts[0], im.parts[0]), self.flags(re.status, im.status)
            )
        else:
            width = 1 if kind == "bool" else self.fmt.width
            raw = (
                int(cast("bool", value))
                if kind == "bool"
                else self.fmt.encode(cast("int | float", value))
            )
            ref = self.local(width)
            for bit in range(width):
                if (raw >> bit) & 1:
                    self.b.x(ref[bit])
            status = self.zero_status
            limit = 2 ** (self.fmt.width - self.fmt.fraction - 1)
            if kind == "real" and not -limit <= value < limit:
                status = self.local(2)
                self.b.x(status[1])
            result = Numeric(kind, (ref,), status)
        self.constants[key] = result
        return result

    def as_complex(self, value: Numeric) -> Numeric:
        """把实数值提升为复数表示。

        Args:
            value: real 类型的 Numeric。

        Returns:
            Numeric: 虚部绑定常量 0 的复数值，状态沿用原值。

        Raises:
            ValidationError: 输入为 bool，不允许隐式转换为复数。
        """
        if value.kind == "complex":
            return value
        if value.kind != "real":
            raise ValidationError("布尔值不能隐式转换为复数")
        return Numeric("complex", (value.parts[0], self.constant(0).parts[0]), value.status)

    def component(self, value: Numeric, index: int) -> Numeric:
        """取出复数的实部或虚部分量。

        Args:
            value: real 或 complex 的 Numeric。
            index: 0 取实部，1 取虚部。

        Returns:
            Numeric: 实数分量；real 输入取虚部时返回常量 0。
        """
        return (
            Numeric("real", (value.parts[index],), value.status)
            if value.kind == "complex"
            else (value if index == 0 else Numeric("real", self.constant(0).parts, value.status))
        )

    def arithmetic(self, kind: str, *values: Numeric) -> Numeric:
        """调用定点算术网络执行一元或二元运算。

        Args:
            kind: 算术名，如 add、sub、mul、div、sqrt、lt、eq。
            values: 参与运算的 Numeric 操作数。

        Returns:
            Numeric: 比较运算为 bool，其余为 real；状态合并操作数状态与新产生的标志位。
        """
        operation = fixed_arithmetic(kind, self.fmt)
        out_width = next(r.type.width for r in operation.module.registers if r.name == "out")
        out, flag = self.local(out_width), self.local(2)
        args = {name: value.parts[0] for name, value in zip(("a", "b"), values, strict=False)}
        self.invoke(operation, {**args, "out": out, "status": flag})
        return Numeric(
            "bool" if kind in {"lt", "eq"} else "real",
            (out,),
            self.flags(*(v.status for v in values), flag),
        )

    def unary(self, kind: str, value: Numeric) -> Numeric:
        """展开一元数值操作。

        real/imag 取复数分量，not 按位取反，conj 将虚部取负，复数 abs
        分解为平方和开方；复数输入的其余算术逐分量分解后重装。

        Args:
            kind: 操作名：real、imag、not、conj、abs、neg 或定点一元算术名。
            value: 操作数 Numeric。

        Returns:
            Numeric: 操作结果。
        """
        if kind in {"real", "imag"}:
            return self.component(value, int(kind == "imag"))
        if kind == "not":
            return Numeric("bool", (self.logic("not", value.parts[0]),), value.status)
        if kind == "conj":
            if value.kind != "complex":
                return value
            return self.complex(
                self.component(value, 0), self.unary("neg", self.component(value, 1))
            )
        if kind == "abs" and value.kind == "complex":
            re, im = self.component(value, 0), self.component(value, 1)
            return self.arithmetic(
                "sqrt", self.binary("add", self.binary("mul", re, re), self.binary("mul", im, im))
            )
        if value.kind == "complex":
            return self.complex(*(self.unary(kind, self.component(value, j)) for j in range(2)))
        return self.arithmetic(kind, value)

    def complex(self, re: Numeric, im: Numeric) -> Numeric:
        """由两个实数值组装复数。

        Args:
            re: 实部分量。
            im: 虚部分量。

        Returns:
            Numeric: 状态为两分量状态合并的复数值。

        Raises:
            ValidationError: 任一分量不是 real。
        """
        if re.kind != "real" or im.kind != "real":
            raise ValidationError("complex(real, imag) 要求两个实数")
        return Numeric("complex", (re.parts[0], im.parts[0]), self.flags(re.status, im.status))

    def binary(self, kind: str, a: Numeric, b: Numeric) -> Numeric:
        """展开二元运算的定点与复数分解。

        布尔逻辑走逐位网络，实数走定点算术；复数按代数关系分解为实数
        分量运算，除法以分母模方归一。

        Args:
            kind: 运算名：and、or、eq、add、sub、mul 或 div。
            a: 左操作数 Numeric。
            b: 右操作数 Numeric。

        Returns:
            Numeric: 运算结果。

        Raises:
            ValidationError: 复数分解遇到未知运算名。
        """
        if kind in {"and", "or"}:
            return Numeric(
                "bool", (self.logic(kind, a.parts[0], b.parts[0]),), self.flags(a.status, b.status)
            )
        if "complex" not in (a.kind, b.kind):
            if kind == "eq" and a.kind == b.kind == "bool":
                ref = self.logic("not", self.logic("xor", a.parts[0], b.parts[0]))
                return Numeric("bool", (ref,), self.flags(a.status, b.status))
            return self.arithmetic(kind, a, b)
        a, b = self.as_complex(a), self.as_complex(b)
        ar, ai, br, bi = (
            self.component(a, 0),
            self.component(a, 1),
            self.component(b, 0),
            self.component(b, 1),
        )
        if kind == "eq":
            return self.binary("and", self.binary("eq", ar, br), self.binary("eq", ai, bi))
        if kind in {"add", "sub"}:
            return self.complex(self.binary(kind, ar, br), self.binary(kind, ai, bi))

        def mul(x: Numeric, y: Numeric) -> Numeric:
            """两值相乘的简写。"""
            return self.binary("mul", x, y)

        if kind == "mul":
            return self.complex(
                self.binary("sub", mul(ar, br), mul(ai, bi)),
                self.binary("add", mul(ar, bi), mul(ai, br)),
            )
        if kind == "div":
            den = self.binary("add", mul(br, br), mul(bi, bi))
            return self.complex(
                self.binary("div", self.binary("add", mul(ar, br), mul(ai, bi)), den),
                self.binary("div", self.binary("sub", mul(ai, br), mul(ar, bi)), den),
            )
        raise ValidationError("未知复数二元分解：" + kind)

    def choose(self, test: Numeric, yes: Numeric, no: Numeric) -> Numeric:
        """按布尔条件逐位选择两分支的电路值。

        real 与 complex 混合的分支先统一提升为复数；选择同时作用于状态位。

        Args:
            test: bool 类型的条件值。
            yes: 条件为真的分支值。
            no: 条件为假的分支值。

        Returns:
            Numeric: 条件状态与被选分支状态合并后的结果。
        """
        if yes.kind != no.kind and {yes.kind, no.kind} <= {"real", "complex"}:
            yes, no = self.as_complex(yes), self.as_complex(no)
        values = tuple(
            self.logic("select", a, b, select=test.parts[0])
            for a, b in zip(yes.parts, no.parts, strict=True)
        )
        selected_status = self.logic("select", yes.status, no.status, select=test.parts[0])
        return Numeric(yes.kind, values, self.flags(test.status, selected_status))

    def add_flag(self, value: Numeric, test: Numeric, *, domain: bool = False) -> Numeric:
        """把一个条件位并入值的状态。

        Args:
            value: 待标记的值；部件不变，仅更新状态。
            test: bool 条件值。
            domain: 为 True 时条件进入 status 第 0 位（定义域失效），
                否则进入第 1 位（值域越界）。

        Returns:
            Numeric: 状态并位后的新值。
        """
        flag = self.logic("flag" if domain else "range_flag", test.parts[0])
        return Numeric(value.kind, value.parts, self.flags(value.status, test.status, flag))

    def kernel(self, name: str, value: Numeric) -> Numeric:
        """调用初等数学核并合并其状态标志。

        Args:
            name: 数学函数名，如 sin、exp 或 log。
            value: real 输入值。

        Returns:
            Numeric: 核输出的实数值，状态合并输入状态与核标志。
        """
        op = elementary_kernel(name, self.fmt, self.config)
        out, flag = self.local(), self.local(2)
        self.invoke(op, {"a": value.parts[0], "out": out, "status": flag})
        return Numeric("real", (out,), self.flags(value.status, flag))

    def atan2(self, y: Numeric, x: Numeric) -> Numeric:
        """以幅值比与象限选择分解二元反正切。

        用 ``|x|``、``|y|`` 中较小者作被除数调用 atan 核，再按大小交换、x 符号
        与 y 符号依次修正象限；两输入均为零时返回零。

        Args:
            y: 纵坐标实数值。
            x: 横坐标实数值。

        Returns:
            Numeric: 象限修正后的角度值。
        """
        zero = self.constant(0)
        ax, ay = self.unary("abs", x), self.unary("abs", y)
        swap = self.binary("lt", ax, ay)
        small, big = self.choose(swap, ax, ay), self.choose(swap, ay, ax)
        iszero = self.binary("eq", big, zero)
        ratio = self.binary("div", small, self.choose(iszero, self.constant(1), big))
        angle = self.kernel("atan", ratio)
        angle = self.choose(swap, self.binary("sub", self.constant(math.pi / 2), angle), angle)
        angle = self.choose(
            self.binary("lt", x, zero), self.binary("sub", self.constant(math.pi), angle), angle
        )
        angle = self.choose(self.binary("lt", y, zero), self.unary("neg", angle), angle)
        return self.choose(iszero, zero, angle)

    def integer_power(self, value: Numeric, exponent: int) -> Numeric:
        """用平方-乘算法计算整数常量幂。

        Args:
            value: 底数值。
            exponent: 整数指数，绝对值不超过 128。

        Returns:
            Numeric: 幂结果；负指数按倒数计算。

        Raises:
            ValidationError: 指数绝对值超过 128。
        """
        if abs(exponent) > 128:
            raise ValidationError("整数幂超过生成上限 128")
        result = self.constant(1, "complex" if value.kind == "complex" else "real")
        factor = value
        n = abs(exponent)
        while n:
            if n & 1:
                result = self.binary("mul", result, factor)
            n //= 2
            if n:
                factor = self.binary("mul", factor, factor)
        return self.binary("div", self.constant(1), result) if exponent < 0 else result

    def intrinsic(self, name: str, args: list[Numeric], kind: str) -> Numeric:
        """展开数学 intrinsic 的实数与复数分解。

        rect、phase、atan2、hypot 与双参 log 先行处理；实数路径走
        sqrt 算术或数学核，复数路径以实数核和代数恒等式组合，log 类
        核额外标记 log(0) 定义域失效。

        Args:
            name: intrinsic 函数名。
            args: 实参 Numeric 列表。
            kind: 期望的结果类型 real 或 complex。

        Returns:
            Numeric: 分解后的结果值。

        Raises:
            ValidationError: 缺少该函数的复数分解实现。
        """
        if name == "rect":
            r, angle = args
            return self.complex(
                self.binary("mul", r, self.kernel("cos", angle)),
                self.binary("mul", r, self.kernel("sin", angle)),
            )
        if name == "phase":
            z = self.as_complex(args[0])
            return self.atan2(self.component(z, 1), self.component(z, 0))
        if name == "atan2":
            return self.atan2(*args)
        if name == "hypot":
            return self.arithmetic(
                "sqrt", self.binary("add", *(self.binary("mul", v, v) for v in args))
            )
        if name == "log" and len(args) == 2:
            return self.binary(
                "div", self.intrinsic("log", args[:1], kind), self.intrinsic("log", args[1:], kind)
            )
        value = args[0]
        if kind == "real":
            return self.arithmetic("sqrt", value) if name == "sqrt" else self.kernel(name, value)
        z = self.as_complex(value)
        re, im = self.component(z, 0), self.component(z, 1)
        c = self.constant

        def add(a: Numeric, b: Numeric) -> Numeric:
            """两值相加的简写。"""
            return self.binary("add", a, b)

        def sub(a: Numeric, b: Numeric) -> Numeric:
            """两值相减的简写。"""
            return self.binary("sub", a, b)

        def mul(a: Numeric, b: Numeric) -> Numeric:
            """两值相乘的简写。"""
            return self.binary("mul", a, b)

        def div(a: Numeric, b: Numeric) -> Numeric:
            """两值相除的简写。"""
            return self.binary("div", a, b)

        def fn(k: str, v: Numeric) -> Numeric:
            """以复数路径递归调用 intrinsic 的简写。"""
            return self.intrinsic(k, [v], "complex")

        if name == "sqrt":
            radius = self.unary("abs", z)
            positive = self.arithmetic("sqrt", mul(c(0.5), add(radius, re)))
            negative = self.arithmetic("sqrt", mul(c(0.5), sub(radius, re)))
            negative = self.choose(
                self.binary("lt", im, c(0)), self.unary("neg", negative), negative
            )
            return self.complex(positive, negative)
        if name in {"log", "log10"}:
            logarithm = self.complex(self.kernel("log", self.unary("abs", z)), self.atan2(im, re))
            # log(0) 明确标为定义域失效，即便多项式总化后有输出。
            logarithm = self.add_flag(logarithm, self.binary("eq", z, c(0j)), domain=True)
            return div(logarithm, c(math.log(10))) if name == "log10" else logarithm
        if name == "exp":
            magnitude = self.kernel("exp", re)
            return self.complex(
                mul(magnitude, self.kernel("cos", im)), mul(magnitude, self.kernel("sin", im))
            )
        if name in {"sin", "cos"}:
            sr, cr = self.kernel("sin", re), self.kernel("cos", re)
            shi, chi = self.kernel("sinh", im), self.kernel("cosh", im)
            return (
                self.complex(mul(sr, chi), mul(cr, shi))
                if name == "sin"
                else self.complex(mul(cr, chi), self.unary("neg", mul(sr, shi)))
            )
        if name in {"sinh", "cosh"}:
            shr, chr_ = self.kernel("sinh", re), self.kernel("cosh", re)
            si, ci = self.kernel("sin", im), self.kernel("cos", im)
            return (
                self.complex(mul(shr, ci), mul(chr_, si))
                if name == "sinh"
                else self.complex(mul(chr_, ci), mul(shr, si))
            )
        if name == "tan":
            return div(fn("sin", z), fn("cos", z))
        if name == "tanh":
            return div(fn("sinh", z), fn("cosh", z))
        if name == "asin":
            return mul(c(-1j), fn("log", add(mul(c(1j), z), fn("sqrt", sub(c(1), mul(z, z))))))
        if name == "acos":
            return sub(c(math.pi / 2), fn("asin", z))
        if name == "atan":
            return mul(
                c(0.5j),
                sub(fn("log", sub(c(1), mul(c(1j), z))), fn("log", add(c(1), mul(c(1j), z)))),
            )
        if name == "asinh":
            return fn("log", add(z, fn("sqrt", add(mul(z, z), c(1)))))
        if name == "acosh":
            return fn("log", add(z, mul(fn("sqrt", add(z, c(1))), fn("sqrt", sub(z, c(1))))))
        if name == "atanh":
            return mul(c(0.5), sub(fn("log", add(c(1), z)), fn("log", sub(c(1), z))))
        raise ValidationError("缺少复数数学分解：" + name)

    def finish(self, outputs: list[tuple[tuple[Ref, ...], Numeric]], status: Ref) -> Operation:
        """复制输出并反算全部临时寄存器后结束构建。

        Args:
            outputs: 每项为 (目标引用 tuple, 输出值 Numeric) 的序列。
            status: 汇总状态写入的 2 位目标寄存器。

        Returns:
            Operation: 构建完成的完整可逆模块。
        """
        flag = self.flags(*(value.status for _, value in outputs))
        forward = tuple(self.b._frames[0])
        for refs, value in outputs:
            for target, source in zip(refs, value.parts, strict=True):
                self.b.xor(source, target)
        self.b.xor(flag, status)
        self.b.emit(Adjoint(forward))
        return self.b.finish()


@lru_cache(maxsize=128)
def elementary_kernel(name: str, fmt: FixedFormat, config: MathConfig) -> Operation:
    """生成单个初等实函数的 Chebyshev 近似核模块。

    在配置区间上取 degree+1 个采样点求 Chebyshev 系数，以 Clenshaw
    递推组合成多项式电路，并对区间越界与定义域错误置相应状态位。

    Args:
        name: 数学函数名，须能取到 math 同名函数并配置近似区间。
        fmt: 定点格式。
        config: MathConfig 生成配置。

    Returns:
        Operation: 接口为 a、out、status 的可逆模块；采样系数配方
        记录在属性 math_approximation 中。

    Raises:
        ValidationError: 近似区间越过函数定义域。
    """
    lo, hi = config.bounds(name)
    degree = config.degree
    count = degree + 1
    center, half = (lo + hi) / 2, (hi - lo) / 2
    theta = [math.pi * (j + 0.5) / count for j in range(count)]
    try:
        samples = [getattr(math, name)(center + half * math.cos(t)) for t in theta]
    except (ValueError, OverflowError) as exc:
        raise ValidationError("数学核区间越过函数定义域：" + name) from exc
    coefficients = [
        2 / count * sum(y * math.cos(k * t) for y, t in zip(samples, theta, strict=True))
        for k in range(count)
    ]
    coefficients[0] *= 0.5
    recipe = json.dumps(
        {
            "method": "chebyshev_clenshaw",
            "function": name,
            "interval": [lo, hi],
            "degree": degree,
            "coefficients": coefficients,
        },
        sort_keys=True,
    )
    b = Builder(
        _name("math_kernel", name, fmt, recipe),
        {"a": Bits(fmt.width), "out": Bits(fmt.width), "status": Bits(2)},
        attributes={"math_kernel": name, "math_approximation": recipe, "correctness": "pending"},
    )
    e = NumericEmitter(b, fmt, config)
    a = Numeric("real", (b["a"],), e.zero_status)
    t = e.binary("mul", e.binary("sub", a, e.constant(center)), e.constant(1 / half))
    next_, next2 = e.constant(0), e.constant(0)
    for coefficient in reversed(coefficients[1:]):
        current = e.binary(
            "add",
            e.binary("sub", e.binary("mul", e.constant(2), e.binary("mul", t, next_)), next2),
            e.constant(coefficient),
        )
        next2, next_ = next_, current
    value = e.binary(
        "add", e.binary("sub", e.binary("mul", t, next_), next2), e.constant(coefficients[0])
    )
    outside = e.binary("or", e.binary("lt", a, e.constant(lo)), e.binary("lt", e.constant(hi), a))
    value = e.add_flag(value, outside)
    if name in {"log", "log10"}:
        bad = e.unary("not", e.binary("lt", e.constant(0), a))
        value = e.add_flag(value, bad, domain=True)
    elif name in {"asin", "acos", "atanh"}:
        magnitude = e.unary("abs", a)
        bad = e.binary("lt", e.constant(1), magnitude)
        if name == "atanh":
            bad = e.unary("not", e.binary("lt", magnitude, e.constant(1)))
        value = e.add_flag(value, bad, domain=True)
    elif name == "acosh":
        value = e.add_flag(value, e.binary("lt", a, e.constant(1)), domain=True)
    return e.finish([((b["out"],), value)], b["status"])
