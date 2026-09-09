"""固定点实数/复数数学分解与按路径传播的数值状态。"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache

from ..arithmetic import BooleanNetwork, fixed_arithmetic
from ..builder import Builder
from ..ir import Adjoint, Bits, ValidationError
from ..library import _name

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


@dataclass(frozen=True)
class MathConfig:
    degree: int = 6
    intervals: tuple[tuple[str, float, float], ...] = ()

    def __post_init__(self):
        if type(self.degree) is not int or not 1 <= self.degree <= 32:
            raise ValidationError("数学核多项式阶数必须为 1..32")
        seen = set()
        for kind, lo, hi in self.intervals:
            if (
                kind not in BOUNDS
                or kind in seen
                or not all(math.isfinite(x) for x in (lo, hi))
                or lo >= hi
            ):
                raise ValidationError("数学核近似区间无效或重复")
            seen.add(kind)

    def bounds(self, kind):
        return next(((lo, hi) for name, lo, hi in self.intervals if name == kind), BOUNDS[kind])


@dataclass(frozen=True)
class Numeric:
    kind: str
    parts: tuple
    status: object


@lru_cache(maxsize=128)
def logic_operation(kind, width=1):
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
    def __init__(self, builder, fmt, config):
        if not fmt.signed:
            raise ValidationError("数学函数编译目前使用有符号定点格式")
        self.b, self.fmt, self.config = builder, fmt, config
        self.counter = 0
        self.constants = {}
        self.zero_status = self.local(2)

    def local(self, width=None):
        self.counter += 1
        return self.b.local(
            "math_tmp_" + str(self.counter), Bits(self.fmt.width if width is None else width)
        )

    def invoke(self, op, arguments):
        """复制重叠的只读参数，满足 RIR 无别名的调用 ABI。"""
        seen, clones = set(), []
        from ..validation import locations

        actual = {}
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
        self.b.call(op, **actual)
        for ref, clone in reversed(clones):
            self.b.xor(ref, clone)

    def logic(self, kind, *refs, select=None):
        operation = logic_operation(kind, refs[0].width)
        width = next(r.type.width for r in operation.module.registers if r.name == "out")
        out = self.local(width)
        arguments = dict(zip(("a", "b"), refs, strict=False))
        if select is not None:
            arguments["select"] = select
        self.invoke(operation, {**arguments, "out": out})
        return out

    def flags(self, *refs):
        unique = []
        for ref in refs:
            if ref != self.zero_status and ref not in unique:
                unique.append(ref)
        if not unique:
            return self.zero_status
        result = unique[0]
        for ref in unique[1:]:
            result = self.logic("or", result, ref)
        return result

    def constant(self, value, kind=None):
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
            raw = int(value) if kind == "bool" else self.fmt.encode(value)
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

    def as_complex(self, value):
        if value.kind == "complex":
            return value
        if value.kind != "real":
            raise ValidationError("布尔值不能隐式转换为复数")
        return Numeric("complex", (value.parts[0], self.constant(0).parts[0]), value.status)

    def component(self, value, index):
        return (
            Numeric("real", (value.parts[index],), value.status)
            if value.kind == "complex"
            else (value if index == 0 else Numeric("real", self.constant(0).parts, value.status))
        )

    def arithmetic(self, kind, *values):
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

    def unary(self, kind, value):
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

    def complex(self, re, im):
        if re.kind != "real" or im.kind != "real":
            raise ValidationError("complex(real, imag) 要求两个实数")
        return Numeric("complex", (re.parts[0], im.parts[0]), self.flags(re.status, im.status))

    def binary(self, kind, a, b):
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

        def mul(x, y):
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

    def choose(self, test, yes, no):
        if yes.kind != no.kind and {yes.kind, no.kind} <= {"real", "complex"}:
            yes, no = self.as_complex(yes), self.as_complex(no)
        values = tuple(
            self.logic("select", a, b, select=test.parts[0])
            for a, b in zip(yes.parts, no.parts, strict=True)
        )
        selected_status = self.logic("select", yes.status, no.status, select=test.parts[0])
        return Numeric(yes.kind, values, self.flags(test.status, selected_status))

    def add_flag(self, value, test, *, domain=False):
        flag = self.logic("flag" if domain else "range_flag", test.parts[0])
        return Numeric(value.kind, value.parts, self.flags(value.status, test.status, flag))

    def kernel(self, name, value):
        op = elementary_kernel(name, self.fmt, self.config)
        out, flag = self.local(), self.local(2)
        self.invoke(op, {"a": value.parts[0], "out": out, "status": flag})
        return Numeric("real", (out,), self.flags(value.status, flag))

    def atan2(self, y, x):
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

    def integer_power(self, value, exponent):
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

    def intrinsic(self, name, args, kind):
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

        def add(a, b):
            return self.binary("add", a, b)

        def sub(a, b):
            return self.binary("sub", a, b)

        def mul(a, b):
            return self.binary("mul", a, b)

        def div(a, b):
            return self.binary("div", a, b)

        def fn(k, v):
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

    def finish(self, outputs, status):
        flag = self.flags(*(value.status for _, value in outputs))
        forward = tuple(self.b._frames[0])
        for refs, value in outputs:
            for target, source in zip(refs, value.parts, strict=True):
                self.b.xor(source, target)
        self.b.xor(flag, status)
        self.b.emit(Adjoint(forward))
        return self.b.finish()


@lru_cache(maxsize=128)
def elementary_kernel(name, fmt, config):
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
