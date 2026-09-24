"有限多项式演化 PDE 的不可变表示与 Python 构造面。"

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from numbers import Number
from typing import Never, cast

from oracq.infrastructure.ir import ValidationError
from oracq.infrastructure.validation import name as check_name


def derivative_add(
    derivative: Sequence[tuple[str, int]], axis: str, order: int
) -> tuple[tuple[str, int], ...]:
    """返回在指定轴上累加导数阶数后的新导数说明。

    同轴阶数相加，零阶条目被去除，结果按轴名排序。

    Args:
        derivative: 现有的 ``(axis, order)`` 二元组序列。
        axis: 空间轴名，须是合法标识符。
        order: 追加的导数阶数，必须为非负整数。

    Returns:
        tuple: 按轴名排序、不含零阶条目的二元组元组。

    Raises:
        ValidationError: 轴名非法或阶数不是非负整数。
    """
    check_name(axis)
    if type(order) is not int or order < 0:
        raise ValidationError("空间导数阶数必须为非负整数")
    values = dict(derivative)
    values[axis] = values.get(axis, 0) + order
    return tuple(sorted((a, n) for a, n in values.items() if n))


@dataclass(frozen=True)
class Atom:
    """按名字引用一个未知场或已知系数的多项式原子，附带已作用的空间导数说明。

    Attributes:
        name: 场或已知系数的名字，须是合法标识符。
        derivative: ``(axis, order)`` 二元组序列，按轴累计的空间导数阶。
    """
    name: str
    derivative: tuple[tuple[str, int], ...] = ()


@dataclass(frozen=True)
class Monomial:
    """多项式单项式：系数乘以未知场原子与已知系数原子的乘积，可再施加外层导数。

    Attributes:
        coefficient: 复数系数。
        fields: 参与乘积的未知场原子。
        known: 参与乘积的已知系数原子，如强迫数据。
        outer_derivative: 作用在完整单项式上的外层空间导数说明。
    """
    coefficient: complex = 1.0
    fields: tuple[Atom, ...] = ()
    known: tuple[Atom, ...] = ()
    outer_derivative: tuple[tuple[str, int], ...] = ()


@dataclass(frozen=True)
class Expr:
    """未知场的有限多项式表达式，各项按加法意义求和。

    支持与 ``Expr`` 或数值常量做加减乘、除以非零数值常量、非负整数幂，
    以及用 ``d`` 施加空间导数。

    Attributes:
        terms: 单项式元组。
    """
    terms: tuple[Monomial, ...]

    def __add__(self, other: Expr | float | complex) -> Expr:
        """返回与 ``other`` 之和的规范化 ``Expr``。"""
        return Expr(self.terms + expression(other).terms).normalized()

    __radd__ = __add__

    def __neg__(self) -> Expr:
        """返回各项系数取反后的 ``Expr``。"""
        return Expr(tuple(replace(t, coefficient=-t.coefficient) for t in self.terms))

    def __sub__(self, other: Expr | float | complex) -> Expr:
        """返回减去 ``other`` 后的规范化 ``Expr``。"""
        return self + -expression(other)

    def __rsub__(self, other: Expr | float | complex) -> Expr:
        """返回 ``other`` 减去 ``self`` 后的规范化 ``Expr``。"""
        return expression(other) + -self

    def __mul__(self, other: Expr | float | complex) -> Expr:
        """返回两个表达式的多项式乘积并规范化。"""
        other = expression(other)
        terms: list[Monomial] = []
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

    def __truediv__(self, value: float | complex) -> Expr:
        """返回除以非零数值常量后的 ``Expr``。"""
        if not isinstance(value, Number) or value == 0:
            raise ValidationError("PDE 只允许除以非零数值常量，未知场分母不在多项式规则内")
        return self * (1 / value)

    def __rtruediv__(self, value: float | complex) -> Never:
        """未知场作分母不在多项式 PDE 规则内，总是抛出 ``ValidationError``。"""
        raise ValidationError("未知场分母不在多项式 PDE 规则内")

    def __pow__(self, power: int) -> Expr:
        """返回表达式的非负整数次幂。"""
        if type(power) is not int or power < 0:
            raise ValidationError("PDE 多项式幂必须是非负整数")
        result = expression(1)
        for _ in range(power):
            result = result * self
        return result

    def d(self, axis: str, order: int = 1) -> Expr:
        """对表达式的每个单项式施加空间导数。

        只含单个未知场原子或只含单个已知系数原子且无外导数的单项式，导数直接
        累加到该原子上；其余单项式把导数记入外层导数，保留算子作用语义而不
        改写成字段导数乘积。

        Args:
            axis: 空间轴名，须是合法标识符。
            order: 导数阶数，必须为非负整数，默认为 1。

        Returns:
            Expr: 施加导数后的新表达式；自身不变。

        Raises:
            ValidationError: 轴名非法或阶数不是非负整数。
        """
        result: list[Monomial] = []
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

    def normalized(self) -> Expr:
        """合并同结构单项式并去除零系数项，返回规范化后的新表达式。

        以 ``(fields, known, outer_derivative)`` 为键累加系数，按键的 ``repr``
        字符串排序，丢弃系数为零的项。

        Returns:
            Expr: 规范化后的新表达式；自身不变。
        """
        terms: dict[
            tuple[tuple[Atom, ...], tuple[Atom, ...], tuple[tuple[str, int], ...]], complex
        ] = {}
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


def expression(value: Expr | float | complex) -> Expr:
    """把 ``value`` 归一为 ``Expr``：数值常量包装为常量表达式，``Expr`` 原样返回。

    零数值对应不含任何单项式的空表达式。

    Args:
        value: ``Expr`` 实例或数值常量。

    Returns:
        Expr: 与 ``value`` 等价的常量 ``Expr``，或 ``value`` 本身。

    Raises:
        ValidationError: ``value`` 既不是 ``Expr`` 也不是数值，或数值非有限。
    """
    if isinstance(value, Expr):
        return value
    if isinstance(value, Number):
        value = complex(value)
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            raise ValidationError("PDE 系数必须有限")
        return Expr(()) if value == 0 else Expr((Monomial(value),))
    raise ValidationError("需要 PDE 表达式或数值常量")


def Field(name: str) -> Expr:
    """构造引用指定未知场的单位 ``Expr``。

    Args:
        name: 未知场名，须是合法标识符。

    Returns:
        Expr: 只含一个无导数未知场原子的表达式。

    Raises:
        ValidationError: 名字非法。
    """
    check_name(name)
    return Expr((Monomial(fields=(Atom(name),)),))


def Known(name: str) -> Expr:
    """构造引用指定已知系数（如强迫数据）的单位 ``Expr``。

    Args:
        name: 已知系数名，须是合法标识符。

    Returns:
        Expr: 只含一个无导数已知系数原子的表达式。

    Raises:
        ValidationError: 名字非法。
    """
    check_name(name)
    return Expr((Monomial(known=(Atom(name),)),))


@dataclass(frozen=True)
class EquationTerm:
    """PDE 方程组中的一项：一个单项式对指定输出分量的贡献。

    Attributes:
        output: 该项贡献到的未知场名字。
        monomial: 单项式内容。
    """
    output: str
    monomial: Monomial


@dataclass(frozen=True)
class OperatorPort:
    """线性化端口：把同类方程项归组为一个待绑定的算子。

    Attributes:
        name: 端口名；线性端口为 ``L``，强迫端口为 ``F``，非线性端口形如 ``B_0``。
        arity: 端口输入重数；线性为 1，强迫为 0，非线性为对应单项式的未知场原子个数。
        terms: 归入该端口的方程项。
    """
    name: str
    arity: int
    terms: tuple[EquationTerm, ...]


@dataclass(frozen=True)
class PolynomialPDE:
    """一阶自治时间演化多项式 PDE 系统的不可变表示（PDE 0.1）。

    右端拆分为已知线性映射，即 ``u' = f + L u + sum_tau B_tau(u, ..., u)``。
    一般经由 ``from_equations`` 构造，或用 ``loads`` 从 JSON 重建。

    Attributes:
        fields: 全部未知场名，非空且不得重复。
        axes: 全部空间轴名，不得重复。
        terms: 全部方程项；每项的 ``output`` 必须在 ``fields`` 中声明。
        label: 标识该 PDE 的标签字符串。
        version: 表示版本，当前为 ``"0.1"``。
    """
    fields: tuple[str, ...]
    axes: tuple[str, ...]
    terms: tuple[EquationTerm, ...]
    label: str = "polynomial_pde"
    version: str = "0.1"

    @classmethod
    def from_equations(
        cls,
        equations: Mapping[str, Expr | float | complex],
        *,
        axes: Sequence[str] = ("x",),
        label: str = "polynomial_pde",
    ) -> PolynomialPDE:
        """从方程右端字典构造并校验 PDE。

        每个右端经 ``expression`` 转为 ``Expr``，用 ``normalized`` 合并同类项后
        冻结为方程项集合。

        Args:
            equations: 从输出分量名到右端表达式的映射；右端可为 ``Expr`` 或数值常量。
            axes: 空间轴名序列，默认为 ``("x",)``。
            label: PDE 标签，默认为 ``"polynomial_pde"``。

        Returns:
            PolynomialPDE: 已通过 ``validate`` 的不可变实例。

        Raises:
            ValidationError: 名字、表达式或拆分后的结构不符合 PDE 规则。
        """
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

    def validate(self) -> PolynomialPDE:
        """校验 PDE 表示的结构、语义与不可变性。

        检查标签与版本、字段非空且不重复、各方程项输出分量已声明、系数有限、
        未知场原子均已声明，以及全部导数说明不可变、轴已声明且阶数为正整数。

        Returns:
            PolynomialPDE: 返回自身，便于链式构造。

        Raises:
            ValidationError: 任一检查不通过。
        """
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
    def degree(self) -> int:
        """所有方程项中未知场原子个数的最大值，即推导使用的非线性次数。

        强迫项计 0，线性项计 1；没有任何项时返回 0。
        """
        return max((len(t.monomial.fields) for t in self.terms), default=0)

    @property
    def ports(self) -> tuple[OperatorPort, ...]:
        """把方程项归组为线性化端口。

        恰含一个未知场原子的项归入线性端口 ``L``；不含未知场原子的项归入强迫
        端口 ``F``；含两个及以上未知场原子的项各成一个多线性端口 ``B_i``，
        ``i`` 为该项在 ``terms`` 中的序号。

        Returns:
            tuple[OperatorPort, ...]: 依次为 ``L``、``F`` 和各 ``B_i``，仅包含存在对应项的端口。
        """
        linear = tuple(t for t in self.terms if len(t.monomial.fields) == 1)
        forcing = tuple(t for t in self.terms if not t.monomial.fields)
        result: list[OperatorPort] = []
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

    def dumps(self) -> str:
        """先执行 ``validate``，再把 PDE 序列化为符合 PDE 0.1 格式的 JSON 文本。

        复系数以 ``[实部, 虚部]`` 数组表示；键排序、两格缩进，文本以换行结尾。

        Returns:
            str: JSON 文本。

        Raises:
            ValidationError: 表示未通过 ``validate``。
        """
        self.validate()
        data: dict[str, object] = {
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
            cast(list[object], data["terms"]).append(item)
        return (
            json.dumps(data, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n"
        )

    @classmethod
    def loads(cls, text: str) -> PolynomialPDE:
        """从 JSON 文本重建 PDE。

        严格检查顶层与各项的键集合，复系数按 ``[实部, 虚部]`` 读回，重建结果
        返回前先经 ``validate``。

        Args:
            text: ``dumps`` 产生的 JSON 文本。

        Returns:
            PolynomialPDE: 重建并通过校验的实例。

        Raises:
            ValidationError: JSON 无效、包含未知字段或重建结果未通过校验。
        """
        try:
            data = json.loads(text)
            if not isinstance(data, dict) or any(
                type(data.get(key)) is not list for key in ("fields", "axes", "terms")
            ):
                raise ValidationError("PDE 字段/空间轴/项必须为 JSON 数组")
            if set(data) != {"version", "label", "fields", "axes", "terms"}:
                raise ValidationError("未知 PDE 字段")

            def atom(raw: Mapping[str, object]) -> Atom:
                """从 JSON 对象重建单个 ``Atom``，字段非法时抛出 ``ValidationError``。"""
                if set(raw) != {"name", "derivative"} or type(raw["derivative"]) is not list:
                    raise ValidationError("未知或无效 PDE atom 字段")
                return Atom(cast(str, raw["name"]), tuple(tuple(d) for d in raw["derivative"]))

            terms: list[EquationTerm] = []
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
