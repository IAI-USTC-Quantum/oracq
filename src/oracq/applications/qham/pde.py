"Immutable representation and Python construction surface of finite polynomial evolution PDEs."

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
    """Return a new derivative specification with the order accumulated on the given axis.

    Orders on the same axis are added, zero-order entries are removed, and the
    result is sorted by axis name.

    Args:
        derivative: Existing sequence of ``(axis, order)`` pairs.
        axis: Spatial axis name, must be a valid identifier.
        order: Derivative order to append, must be a non-negative integer.

    Returns:
        tuple: Tuple of pairs sorted by axis name without zero-order entries.

    Raises:
        ValidationError: The axis name is invalid or the order is not a
            non-negative integer.
    """
    check_name(axis)
    if type(order) is not int or order < 0:
        raise ValidationError("spatial derivative order must be a non-negative integer")
    values = dict(derivative)
    values[axis] = values.get(axis, 0) + order
    return tuple(sorted((a, n) for a, n in values.items() if n))


@dataclass(frozen=True)
class Atom:
    """Polynomial atom referencing an unknown field or known coefficient by name, carrying the applied spatial derivative specification.

    Attributes:
        name: Name of the field or known coefficient, must be a valid
            identifier.
        derivative: Sequence of ``(axis, order)`` pairs, the accumulated
            spatial derivative orders per axis.
    """
    name: str
    derivative: tuple[tuple[str, int], ...] = ()


@dataclass(frozen=True)
class Monomial:
    """Polynomial monomial: a coefficient times the product of unknown-field atoms and known-coefficient atoms, optionally under an outer derivative.

    Attributes:
        coefficient: Complex coefficient.
        fields: Unknown-field atoms participating in the product.
        known: Known-coefficient atoms participating in the product, such as
            forcing data.
        outer_derivative: Outer spatial derivative specification applied to the
            complete monomial.
    """
    coefficient: complex = 1.0
    fields: tuple[Atom, ...] = ()
    known: tuple[Atom, ...] = ()
    outer_derivative: tuple[tuple[str, int], ...] = ()


@dataclass(frozen=True)
class Expr:
    """Finite polynomial expression of unknown fields, with terms summed additively.

    Supports addition, subtraction, and multiplication with ``Expr`` or numeric
    constants, division by a non-zero numeric constant, non-negative integer
    powers, and applying spatial derivatives with ``d``.

    Attributes:
        terms: Tuple of monomials.
    """
    terms: tuple[Monomial, ...]

    def __add__(self, other: Expr | float | complex) -> Expr:
        """Return the normalized ``Expr`` of the sum with ``other``."""
        return Expr(self.terms + expression(other).terms).normalized()

    __radd__ = __add__

    def __neg__(self) -> Expr:
        """Return the ``Expr`` with every coefficient negated."""
        return Expr(tuple(replace(t, coefficient=-t.coefficient) for t in self.terms))

    def __sub__(self, other: Expr | float | complex) -> Expr:
        """Return the normalized ``Expr`` after subtracting ``other``."""
        return self + -expression(other)

    def __rsub__(self, other: Expr | float | complex) -> Expr:
        """Return the normalized ``Expr`` of ``other`` minus ``self``."""
        return expression(other) + -self

    def __mul__(self, other: Expr | float | complex) -> Expr:
        """Return the normalized polynomial product of two expressions."""
        other = expression(other)
        terms: list[Monomial] = []
        for a in self.terms:
            for b in other.terms:
                if a.outer_derivative or b.outer_derivative:
                    # Numeric scaling commutes out of the derivative; an unknown field
                    # multiplied by an already-contracted derivative term is outside this rule.
                    if not b.fields and not b.known and not b.outer_derivative:
                        terms.append(replace(a, coefficient=a.coefficient * b.coefficient))
                        continue
                    if not a.fields and not a.known and not a.outer_derivative:
                        terms.append(replace(b, coefficient=a.coefficient * b.coefficient))
                        continue
                    raise ValidationError(
                        "a polynomial term with an outer derivative cannot directly become a product factor; write it as a product of field derivatives within the rules"
                    )
                terms.append(
                    Monomial(a.coefficient * b.coefficient, a.fields + b.fields, a.known + b.known)
                )
        return Expr(tuple(terms)).normalized()

    __rmul__ = __mul__

    def __truediv__(self, value: float | complex) -> Expr:
        """Return the ``Expr`` after division by a non-zero numeric constant."""
        if not isinstance(value, Number) or value == 0:
            raise ValidationError("a PDE only allows division by a non-zero numeric constant; unknown-field denominators are outside the polynomial rules")
        return self * (1 / value)

    def __rtruediv__(self, value: float | complex) -> Never:
        """An unknown field in the denominator is outside the polynomial PDE rules; always raises ``ValidationError``."""
        raise ValidationError("unknown-field denominators are outside the polynomial PDE rules")

    def __pow__(self, power: int) -> Expr:
        """Return the expression raised to a non-negative integer power."""
        if type(power) is not int or power < 0:
            raise ValidationError("PDE polynomial powers must be non-negative integers")
        result = expression(1)
        for _ in range(power):
            result = result * self
        return result

    def d(self, axis: str, order: int = 1) -> Expr:
        """Apply a spatial derivative to every monomial of the expression.

        For monomials containing exactly one unknown-field atom or exactly one
        known-coefficient atom and no outer derivative, the derivative
        accumulates directly on that atom; all other monomials record the
        derivative in the outer derivative, preserving operator-action
        semantics instead of rewriting into a product of field derivatives.

        Args:
            axis: Spatial axis name, must be a valid identifier.
            order: Derivative order, must be a non-negative integer; defaults
                to 1.

        Returns:
            Expr: New expression with the derivative applied; self is
            unchanged.

        Raises:
            ValidationError: The axis name is invalid or the order is not a
                non-negative integer.
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
        """Merge structurally identical monomials and drop zero-coefficient terms, returning the normalized new expression.

        Coefficients are accumulated keyed by ``(fields, known,
        outer_derivative)``, sorted by the ``repr`` string of the key, and
        zero-coefficient terms are discarded.

        Returns:
            Expr: The normalized new expression; self is unchanged.
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
    """Normalize ``value`` to an ``Expr``: numeric constants are wrapped as constant expressions and an ``Expr`` is returned unchanged.

    The numeric zero maps to the empty expression containing no monomials.

    Args:
        value: An ``Expr`` instance or a numeric constant.

    Returns:
        Expr: A constant ``Expr`` equivalent to ``value``, or ``value``
        itself.

    Raises:
        ValidationError: ``value`` is neither an ``Expr`` nor a number, or the
            number is not finite.
    """
    if isinstance(value, Expr):
        return value
    if isinstance(value, Number):
        value = complex(value)
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            raise ValidationError("PDE coefficients must be finite")
        return Expr(()) if value == 0 else Expr((Monomial(value),))
    raise ValidationError("a PDE expression or numeric constant is required")


def Field(name: str) -> Expr:
    """Build the unit ``Expr`` referencing the given unknown field.

    Args:
        name: Unknown field name, must be a valid identifier.

    Returns:
        Expr: Expression containing only one derivative-free unknown-field
        atom.

    Raises:
        ValidationError: The name is invalid.
    """
    check_name(name)
    return Expr((Monomial(fields=(Atom(name),)),))


def Known(name: str) -> Expr:
    """Build the unit ``Expr`` referencing the given known coefficient, such as forcing data.

    Args:
        name: Known coefficient name, must be a valid identifier.

    Returns:
        Expr: Expression containing only one derivative-free known-coefficient
        atom.

    Raises:
        ValidationError: The name is invalid.
    """
    check_name(name)
    return Expr((Monomial(known=(Atom(name),)),))


@dataclass(frozen=True)
class EquationTerm:
    """One term of the PDE system: the contribution of a monomial to a given output component.

    Attributes:
        output: Name of the unknown field this term contributes to.
        monomial: The monomial content.
    """
    output: str
    monomial: Monomial


@dataclass(frozen=True)
class OperatorPort:
    """Linearization port: groups same-class equation terms into one operator awaiting binding.

    Attributes:
        name: Port name; ``L`` for the linear port, ``F`` for the forcing
            port, and ``B_0``-style names for nonlinear ports.
        arity: Number of port inputs; 1 for linear, 0 for forcing, and for
            nonlinear the count of unknown-field atoms in the corresponding
            monomial.
        terms: Equation terms assigned to this port.
    """
    name: str
    arity: int
    terms: tuple[EquationTerm, ...]


@dataclass(frozen=True)
class PolynomialPDE:
    """Immutable representation of a first-order autonomous time-evolution polynomial PDE system (PDE 0.1).

    The right-hand side is split into known linear maps, i.e.
    ``u' = f + L u + sum_tau B_tau(u, ..., u)``. Usually constructed via
    ``from_equations``, or rebuilt from JSON with ``loads``.

    Attributes:
        fields: All unknown field names, non-empty and without duplicates.
        axes: All spatial axis names, without duplicates.
        terms: All equation terms; each term's ``output`` must be declared in
            ``fields``.
        label: Label string identifying this PDE.
        version: Representation version, currently ``"0.1"``.
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
        """Construct and validate a PDE from a dictionary of equation right-hand sides.

        Each right-hand side is converted to an ``Expr`` by ``expression`` and
        frozen into the equation-term set after merging like terms with
        ``normalized``.

        Args:
            equations: Mapping from output component name to right-hand-side
                expression; a right-hand side may be an ``Expr`` or a numeric
                constant.
            axes: Sequence of spatial axis names, defaulting to ``("x",)``.
            label: PDE label, defaulting to ``"polynomial_pde"``.

        Returns:
            PolynomialPDE: Immutable instance that has passed ``validate``.

        Raises:
            ValidationError: A name, an expression, or the structure after
                splitting violates the PDE rules.
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
        """Validate the structure, semantics, and immutability of the PDE representation.

        Checks the label and version, that the fields are non-empty and
        duplicate-free, that each equation term's output component is
        declared, that coefficients are finite, that all unknown-field atoms
        are declared, and that every derivative specification is immutable with
        declared axes and positive integer orders.

        Returns:
            PolynomialPDE: Returns self for chained construction.

        Raises:
            ValidationError: Any check fails.
        """
        if not isinstance(self.label, str):
            raise ValidationError("the PDE label must be a string")
        if (
            self.version != "0.1"
            or not self.fields
            or len(set(self.fields)) != len(self.fields)
            or len(set(self.axes)) != len(self.axes)
        ):
            raise ValidationError("invalid PDE fields, spatial axes, or version")
        if any(type(x) is not tuple for x in (self.fields, self.axes, self.terms)):
            raise ValidationError("the PDE representation must be immutable")
        for key in (*self.fields, *self.axes):
            check_name(key)
        for term in self.terms:
            if term.output not in self.fields:
                raise ValidationError("undeclared PDE output component")
            m = term.monomial
            if any(type(value) is not tuple for value in (m.fields, m.known, m.outer_derivative)):
                raise ValidationError("PDE monomials must be immutable")
            if not math.isfinite(m.coefficient.real) or not math.isfinite(m.coefficient.imag):
                raise ValidationError("PDE coefficients must be finite")
            for atom in m.fields:
                if atom.name not in self.fields:
                    raise ValidationError("undeclared PDE unknown field: " + atom.name)
            for atom in (*m.fields, *m.known):
                if type(atom.derivative) is not tuple:
                    raise ValidationError("PDE derivative specifications must be immutable")
            for atom in m.known:
                check_name(atom.name)
            for derivative in (m.outer_derivative, *(a.derivative for a in (*m.fields, *m.known))):
                if len(dict(derivative)) != len(derivative) or any(
                    axis not in self.axes or type(order) is not int or order < 1
                    for axis, order in derivative
                ):
                    raise ValidationError("invalid PDE spatial derivative")
        return self

    @property
    def degree(self) -> int:
        """Maximum count of unknown-field atoms over all equation terms, i.e. the nonlinear degree used by the derivation.

        Forcing terms count 0 and linear terms count 1; returns 0 when there
        are no terms at all.
        """
        return max((len(t.monomial.fields) for t in self.terms), default=0)

    @property
    def ports(self) -> tuple[OperatorPort, ...]:
        """Group equation terms into linearization ports.

        Terms with exactly one unknown-field atom go to the linear port ``L``;
        terms with none go to the forcing port ``F``; each term with two or
        more unknown-field atoms forms its own multilinear port ``B_i``, with
        ``i`` the term's index in ``terms``.

        Returns:
            tuple[OperatorPort, ...]: ``L``, ``F``, then each ``B_i``, in
            order, including only ports with corresponding terms.
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
        """Run ``validate`` first, then serialize the PDE to JSON text following the PDE 0.1 format.

        Complex coefficients are written as ``[real part, imaginary part]``
        arrays; keys are sorted, indentation is two spaces, and the text ends
        with a newline.

        Returns:
            str: JSON text.

        Raises:
            ValidationError: The representation fails ``validate``.
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
        """Rebuild a PDE from JSON text.

        Key sets at the top level and of each item are strictly checked,
        complex coefficients are read back as ``[real part, imaginary part]``,
        and the rebuilt result passes ``validate`` before being returned.

        Args:
            text: JSON text produced by ``dumps``.

        Returns:
            PolynomialPDE: Rebuilt and validated instance.

        Raises:
            ValidationError: The JSON is invalid, contains unknown fields, or
                the rebuilt result fails validation.
        """
        try:
            data = json.loads(text)
            if not isinstance(data, dict) or any(
                type(data.get(key)) is not list for key in ("fields", "axes", "terms")
            ):
                raise ValidationError("PDE fields, spatial axes, and terms must be JSON arrays")
            if set(data) != {"version", "label", "fields", "axes", "terms"}:
                raise ValidationError("unknown PDE field")

            def atom(raw: Mapping[str, object]) -> Atom:
                """Rebuild a single ``Atom`` from a JSON object, raising ``ValidationError`` on invalid fields."""
                if set(raw) != {"name", "derivative"} or type(raw["derivative"]) is not list:
                    raise ValidationError("unknown or invalid PDE atom field")
                return Atom(cast(str, raw["name"]), tuple(tuple(d) for d in raw["derivative"]))

            terms: list[EquationTerm] = []
            for item in data["terms"]:
                if set(item) != {"output", "monomial"}:
                    raise ValidationError("unknown PDE term field")
                raw = item["monomial"]
                if set(raw) != {"coefficient", "fields", "known", "outer_derivative"} or any(
                    type(raw[key]) is not list for key in raw
                ):
                    raise ValidationError("unknown or invalid PDE monomial field")
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
            raise ValidationError("invalid PDE JSON: " + str(exc)) from exc
