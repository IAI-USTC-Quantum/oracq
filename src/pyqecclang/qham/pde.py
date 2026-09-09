"""有限多项式演化 PDE 的不可变表示与 Python 构造面。"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, replace
from numbers import Number

from ..ir import ValidationError
from ..validation import name as check_name


def derivative_add(derivative, axis, order):
    check_name(axis)
    if type(order) is not int or order < 0:
        raise ValidationError("空间导数阶数必须为非负整数")
    values = dict(derivative)
    values[axis] = values.get(axis, 0) + order
    return tuple(sorted((a, n) for a, n in values.items() if n))


@dataclass(frozen=True)
class Atom:
    name: str
    derivative: tuple[tuple[str, int], ...] = ()


@dataclass(frozen=True)
class Monomial:
    coefficient: complex = 1.0
    fields: tuple[Atom, ...] = ()
    known: tuple[Atom, ...] = ()
    outer_derivative: tuple[tuple[str, int], ...] = ()


@dataclass(frozen=True)
class Expr:
    terms: tuple[Monomial, ...]

    def __add__(self, other):
        return Expr(self.terms + expression(other).terms).normalized()

    __radd__ = __add__

    def __neg__(self):
        return Expr(tuple(replace(t, coefficient=-t.coefficient) for t in self.terms))

    def __sub__(self, other):
        return self + -expression(other)

    def __rsub__(self, other):
        return expression(other) + -self

    def __mul__(self, other):
        other = expression(other)
        terms = []
        for a in self.terms:
            for b in other.terms:
                if a.outer_derivative or b.outer_derivative:
                    # 数值缩放可交换到导数外；未知场乘以已收缩的导数项不在本规则内。
                    if not b.fields and not b.known and not b.outer_derivative:
                        terms.append(replace(a, coefficient=a.coefficient * b.coefficient))
                        continue
                    if not a.fields and not a.known and not a.outer_derivative:
                        terms.append(replace(b, coefficient=a.coefficient * b.coefficient))
                        continue
                    raise ValidationError(
                        "带外导数的多项式项不能直接成为乘积因子；请写成规则内的字段导数乘积"
                    )
                terms.append(
                    Monomial(a.coefficient * b.coefficient, a.fields + b.fields, a.known + b.known)
                )
        return Expr(tuple(terms)).normalized()

    __rmul__ = __mul__

    def __truediv__(self, value):
        if not isinstance(value, Number) or value == 0:
            raise ValidationError("PDE 只允许除以非零数值常量，未知场分母不在多项式规则内")
        return self * (1 / value)

    def __rtruediv__(self, value):
        raise ValidationError("未知场分母不在多项式 PDE 规则内")

    def __pow__(self, power):
        if type(power) is not int or power < 0:
            raise ValidationError("PDE 多项式幂必须是非负整数")
        result = expression(1)
        for _ in range(power):
            result = result * self
        return result

    def d(self, axis, order=1):
        result = []
        for term in self.terms:
            if len(term.fields) == 1 and not term.known and not term.outer_derivative:
                atom = term.fields[0]
                result.append(
                    replace(
                        term,
                        fields=(
                            replace(atom, derivative=derivative_add(atom.derivative, axis, order)),
                        ),
                    )
                )
            elif len(term.known) == 1 and not term.fields and not term.outer_derivative:
                atom = term.known[0]
                result.append(
                    replace(
                        term,
                        known=(
                            replace(atom, derivative=derivative_add(atom.derivative, axis, order)),
                        ),
                    )
                )
            else:
                result.append(
                    replace(
                        term, outer_derivative=derivative_add(term.outer_derivative, axis, order)
                    )
                )
        return Expr(tuple(result))

    def normalized(self):
        terms = {}
        for term in self.terms:
            key = (term.fields, term.known, term.outer_derivative)
            terms[key] = terms.get(key, 0) + term.coefficient
        return Expr(
            tuple(
                Monomial(c, *key)
                for key, c in sorted(terms.items(), key=lambda item: repr(item[0]))
                if c != 0
            )
        )


def expression(value):
    if isinstance(value, Expr):
        return value
    if isinstance(value, Number):
        value = complex(value)
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            raise ValidationError("PDE 系数必须有限")
        return Expr(()) if value == 0 else Expr((Monomial(value),))
    raise ValidationError("需要 PDE 表达式或数值常量")


def Field(name):
    check_name(name)
    return Expr((Monomial(fields=(Atom(name),)),))


def Known(name):
    check_name(name)
    return Expr((Monomial(known=(Atom(name),)),))


@dataclass(frozen=True)
class EquationTerm:
    output: str
    monomial: Monomial


@dataclass(frozen=True)
class OperatorPort:
    name: str
    arity: int
    terms: tuple[EquationTerm, ...]


@dataclass(frozen=True)
class PolynomialPDE:
    fields: tuple[str, ...]
    axes: tuple[str, ...]
    terms: tuple[EquationTerm, ...]
    label: str = "polynomial_pde"
    version: str = "0.1"

    @classmethod
    def from_equations(cls, equations, *, axes=("x",), label="polynomial_pde"):
        result = cls(
            tuple(equations),
            tuple(axes),
            tuple(
                EquationTerm(output, term)
                for output, value in equations.items()
                for term in expression(value).normalized().terms
            ),
            label,
        )
        return result.validate()

    def validate(self):
        if not isinstance(self.label, str):
            raise ValidationError("PDE label 必须为字符串")
        if (
            self.version != "0.1"
            or not self.fields
            or len(set(self.fields)) != len(self.fields)
            or len(set(self.axes)) != len(self.axes)
        ):
            raise ValidationError("PDE 字段/空间轴/版本无效")
        if any(type(x) is not tuple for x in (self.fields, self.axes, self.terms)):
            raise ValidationError("PDE 表示必须不可变")
        for key in (*self.fields, *self.axes):
            check_name(key)
        for term in self.terms:
            if term.output not in self.fields:
                raise ValidationError("PDE 输出分量未声明")
            m = term.monomial
            if any(type(value) is not tuple for value in (m.fields, m.known, m.outer_derivative)):
                raise ValidationError("PDE 单项式必须不可变")
            if not math.isfinite(m.coefficient.real) or not math.isfinite(m.coefficient.imag):
                raise ValidationError("PDE 系数必须有限")
            for atom in m.fields:
                if atom.name not in self.fields:
                    raise ValidationError("PDE 未知字段：" + atom.name)
            for atom in (*m.fields, *m.known):
                if type(atom.derivative) is not tuple:
                    raise ValidationError("PDE 导数说明必须不可变")
            for atom in m.known:
                check_name(atom.name)
            for derivative in (m.outer_derivative, *(a.derivative for a in (*m.fields, *m.known))):
                if len(dict(derivative)) != len(derivative) or any(
                    axis not in self.axes or type(order) is not int or order < 1
                    for axis, order in derivative
                ):
                    raise ValidationError("PDE 空间导数无效")
        return self

    @property
    def degree(self):
        return max((len(t.monomial.fields) for t in self.terms), default=0)

    @property
    def ports(self):
        linear = tuple(t for t in self.terms if len(t.monomial.fields) == 1)
        forcing = tuple(t for t in self.terms if not t.monomial.fields)
        result = []
        if linear:
            result.append(OperatorPort("L", 1, linear))
        if forcing:
            result.append(OperatorPort("F", 0, forcing))
        result += [
            OperatorPort("B_" + str(i), len(t.monomial.fields), (t,))
            for i, t in enumerate(self.terms)
            if len(t.monomial.fields) >= 2
        ]
        return tuple(result)

    def dumps(self):
        self.validate()
        data = {
            "version": self.version,
            "label": self.label,
            "fields": self.fields,
            "axes": self.axes,
            "terms": [],
        }
        for term in self.terms:
            item = asdict(term)
            c = term.monomial.coefficient
            item["monomial"]["coefficient"] = [c.real, c.imag]
            data["terms"].append(item)
        return (
            json.dumps(data, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n"
        )

    @classmethod
    def loads(cls, text):
        try:
            data = json.loads(text)
            if not isinstance(data, dict) or any(
                type(data.get(key)) is not list for key in ("fields", "axes", "terms")
            ):
                raise ValidationError("PDE 字段/空间轴/项必须为 JSON 数组")
            if set(data) != {"version", "label", "fields", "axes", "terms"}:
                raise ValidationError("未知 PDE 字段")

            def atom(raw):
                if set(raw) != {"name", "derivative"} or type(raw["derivative"]) is not list:
                    raise ValidationError("未知或无效 PDE atom 字段")
                return Atom(raw["name"], tuple(tuple(d) for d in raw["derivative"]))

            terms = []
            for item in data["terms"]:
                if set(item) != {"output", "monomial"}:
                    raise ValidationError("未知 PDE 项字段")
                raw = item["monomial"]
                if set(raw) != {"coefficient", "fields", "known", "outer_derivative"} or any(
                    type(raw[key]) is not list for key in raw
                ):
                    raise ValidationError("未知或无效 PDE 单项式字段")
                terms.append(
                    EquationTerm(
                        item["output"],
                        Monomial(
                            complex(*raw["coefficient"]),
                            tuple(atom(x) for x in raw["fields"]),
                            tuple(atom(x) for x in raw["known"]),
                            tuple(tuple(d) for d in raw["outer_derivative"]),
                        ),
                    )
                )
            return cls(
                tuple(data["fields"]),
                tuple(data["axes"]),
                tuple(terms),
                data["label"],
                data["version"],
            ).validate()
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError("无效 PDE JSON：" + str(exc)) from exc
