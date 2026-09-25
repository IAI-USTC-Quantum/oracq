"Fixed-point real/complex math decomposition and per-path propagated numeric status."

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from typing import cast

from oracq.algorithms.common.arithmetic import (
    BooleanNetwork,
    FixedFormat,
    fixed_arithmetic,
)
from oracq.algorithms.input_model.operators import _name
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Adjoint, Bits, Ref, ValidationError

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
"""Default Chebyshev approximation intervals for each math kernel, mapping function names to (lower bound, upper bound)."""


@dataclass(frozen=True)
class MathConfig:
    """Math kernel generation configuration: polynomial degree and overridable approximation intervals.

    Attributes:
        degree: Chebyshev polynomial degree, in the range 1..32.
        intervals: Each entry is (function name, lower bound, upper bound), overriding the approximation interval of the named math kernel;
            kernels not listed use the default interval from BOUNDS.

    Raises:
        ValidationError: The degree is outside 1..32, or an interval references an unknown function, is duplicated, contains non-finite values, or has a lower bound not below the upper bound.
    """

    degree: int = 6
    intervals: tuple[tuple[str, float, float], ...] = ()

    def __post_init__(self) -> None:
        """Validate the polynomial degree range and the approximation interval override table."""
        if type(self.degree) is not int or not 1 <= self.degree <= 32:
            raise ValidationError("math kernel polynomial degree must be between 1 and 32")
        seen: set[str] = set()
        for kind, lo, hi in self.intervals:
            if (
                kind not in BOUNDS
                or kind in seen
                or not all(math.isfinite(x) for x in (lo, hi))
                or lo >= hi
            ):
                raise ValidationError("invalid or duplicate math kernel approximation interval")
            seen.add(kind)

    def bounds(self, kind: str) -> tuple[float, float]:
        """Query the approximation interval of a named math kernel.

        Args:
            kind: Math function name, such as sin or exp.

        Returns:
            tuple: ``(lower bound, upper bound)``; the override entry in ``intervals`` wins, otherwise the default interval is used.

        Raises:
            KeyError: The name is neither in ``intervals`` nor in BOUNDS.
        """
        return next(((lo, hi) for name, lo, hi in self.intervals if name == kind), BOUNDS[kind])


@dataclass(frozen=True)
class Numeric:
    """Generation-time circuit representation of a math value.

    Attributes:
        kind: Value kind: real, complex or bool.
        parts: Tuple of register references carrying the value; complex uses (real part, imaginary part), others are unary.
        status: Reference to the 2-bit status register; bit 0 aggregates domain failures, bit 1 aggregates range/word-length overflows.
    """

    kind: str
    parts: tuple[Ref, ...]
    status: Ref


@lru_cache(maxsize=128)
def logic_operation(kind: str, width: int = 1) -> Operation:
    """Construct the boolean network of a bitwise logic/select operation.

    Args:
        kind: Operation name: not, flag, range_flag, and, or, xor or select.
        width: Input register bit width.

    Returns:
        Operation: The interface is a, b as needed and select, with output out;
            flag and range_flag copy the condition bit into output bit 0 and bit 1 respectively.
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
    """Expands MIR numeric operations into fixed-point reversible circuits on a Builder.

    Performs fixed-point arithmetic, bitwise logic and math kernel calls
    through Numeric values, tracks status bits per value, and at the end
    copies outputs and uncomputes all temporary registers.

    Args:
        builder: The target circuit Builder.
        fmt: The fixed-point format FixedFormat; currently required to be signed.
        config: The math kernel generation configuration MathConfig.

    Raises:
        ValidationError: fmt is not a signed fixed-point format.
    """

    def __init__(self, builder: Builder, fmt: FixedFormat, config: MathConfig) -> None:
        """Bind the builder, fixed-point format and math kernel configuration, and initialize the constant cache."""
        if not fmt.signed:
            raise ValidationError("math function compilation currently uses a signed fixed-point format")
        self.b: Builder = builder
        self.fmt: FixedFormat = fmt
        self.config: MathConfig = config
        self.counter: int = 0
        self.constants: dict[tuple[str, bool | int | float | complex], Numeric] = {}
        self.zero_status: Ref = self.local(2)

    def local(self, width: int | None = None) -> Ref:
        """Allocate one math temporary register.

        Args:
            width: Bit width; the full width of the fixed-point format is used when omitted.

        Returns:
            Ref: The newly allocated register reference, named with an increasing number after the ``math_tmp_`` prefix.
        """
        self.counter += 1
        return self.b.local(
            "math_tmp_" + str(self.counter), Bits(self.fmt.width if width is None else width)
        )

    def invoke(self, op: Operation, arguments: dict[str, Ref]) -> None:
        """Copy overlapping read-only arguments to satisfy the alias-free RIR call ABI.

        Args:
            op: The implemented operation to call.
            arguments: Mapping from formal parameter register names to actual argument views; read-only views that overlap
                each other are temporarily copied and uncomputed after the call.
        """
        seen: set[tuple[str, int]] = set()
        clones: list[tuple[Ref, Ref]] = []
        from oracq.infrastructure.validation import locations

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
        """Invoke a bitwise logic/select network and return the result register.

        Args:
            kind: Operation name, see logic_operation.
            refs: Register references participating in the operation; binary operations take the first two.
            select: The condition bit reference for the select operation.

        Returns:
            Ref: The new register reference holding the operation result.
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
        """Merge several 2-bit status registers.

        The all-zero zero_status does not participate in the merge; with no
        other input, zero_status itself is returned.

        Args:
            refs: Status register references.

        Returns:
            Ref: The status register reference obtained by bitwise-OR of the inputs.
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
        """Freeze a generation-time value into a circuit constant and cache it.

        Complex values are split into two real constants with merged status; a
        real value beyond the representable magnitude of the format sets the
        range overflow bit. Identical (kind, value) pairs are encoded once.

        Args:
            value: A bool, int, float or complex constant.
            kind: Explicit value kind; inferred from the Python type when omitted.

        Returns:
            Numeric: The circuit representation of the constant.
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
        """Promote a real value to the complex representation.

        Args:
            value: A Numeric of real kind.

        Returns:
            Numeric: The complex value with the imaginary part bound to constant 0, reusing the original status.

        Raises:
            ValidationError: The input is bool; implicit conversion to complex is not allowed.
        """
        if value.kind == "complex":
            return value
        if value.kind != "real":
            raise ValidationError("boolean values cannot be implicitly converted to complex")
        return Numeric("complex", (value.parts[0], self.constant(0).parts[0]), value.status)

    def component(self, value: Numeric, index: int) -> Numeric:
        """Extract the real or imaginary component of a complex value.

        Args:
            value: A Numeric of real or complex kind.
            index: 0 selects the real part, 1 the imaginary part.

        Returns:
            Numeric: The real component; taking the imaginary part of a real input returns constant 0.
        """
        return (
            Numeric("real", (value.parts[index],), value.status)
            if value.kind == "complex"
            else (value if index == 0 else Numeric("real", self.constant(0).parts, value.status))
        )

    def arithmetic(self, kind: str, *values: Numeric) -> Numeric:
        """Invoke a fixed-point arithmetic network for a unary or binary operation.

        Args:
            kind: Arithmetic name, such as add, sub, mul, div, sqrt, lt, eq.
            values: The Numeric operands of the operation.

        Returns:
            Numeric: Comparisons yield bool, the rest real; the status merges the operand statuses with the newly produced flag bits.
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
        """Expand a unary numeric operation.

        real/imag extract complex components, not inverts bitwise, conj
        negates the imaginary part, and complex abs decomposes into a
        square-sum square root; the remaining arithmetic on complex inputs
        is decomposed per component and reassembled.

        Args:
            kind: Operation name: real, imag, not, conj, abs, neg, or a fixed-point unary arithmetic name.
            value: The operand Numeric.

        Returns:
            Numeric: The operation result.
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
        """Assemble a complex value from two real values.

        Args:
            re: The real component.
            im: The imaginary component.

        Returns:
            Numeric: The complex value whose status merges the two component statuses.

        Raises:
            ValidationError: Either component is not real.
        """
        if re.kind != "real" or im.kind != "real":
            raise ValidationError("complex requires two real values")
        return Numeric("complex", (re.parts[0], im.parts[0]), self.flags(re.status, im.status))

    def binary(self, kind: str, a: Numeric, b: Numeric) -> Numeric:
        """Expand the fixed-point and complex decompositions of a binary operation.

        Boolean logic goes through the bitwise network, real values through
        fixed-point arithmetic; complex values decompose into real-component
        operations by algebraic identities, with division normalized by the
        squared modulus of the denominator.

        Args:
            kind: Operation name: and, or, eq, add, sub, mul or div.
            a: Left operand Numeric.
            b: Right operand Numeric.

        Returns:
            Numeric: The operation result.

        Raises:
            ValidationError: The complex decomposition hits an unknown operation name.
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
            """Shorthand for multiplying two values."""
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
        raise ValidationError("unknown complex binary decomposition: " + kind)

    def choose(self, test: Numeric, yes: Numeric, no: Numeric) -> Numeric:
        """Select between the circuit values of two branches bitwise by a boolean condition.

        Branches mixing real and complex are first uniformly promoted to
        complex; the selection also applies to the status bits.

        Args:
            test: The condition value of bool kind.
            yes: The branch value chosen when the condition is true.
            no: The branch value chosen when the condition is false.

        Returns:
            Numeric: The result with the condition status merged with the selected branch status.
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
        """Merge one condition bit into the status of a value.

        Args:
            value: The value to flag; its parts stay unchanged, only the status is updated.
            test: The bool condition value.
            domain: When True the condition enters status bit 0 for domain failure,
                otherwise bit 1 for range overflow.

        Returns:
            Numeric: The new value with the bit merged into its status.
        """
        flag = self.logic("flag" if domain else "range_flag", test.parts[0])
        return Numeric(value.kind, value.parts, self.flags(value.status, test.status, flag))

    def kernel(self, name: str, value: Numeric) -> Numeric:
        """Invoke an elementary math kernel and merge its status flag.

        Args:
            name: Math function name, such as sin, exp or log.
            value: The real input value.

        Returns:
            Numeric: The real output value of the kernel, with status merging the input status and the kernel flag.
        """
        op = elementary_kernel(name, self.fmt, self.config)
        out, flag = self.local(), self.local(2)
        self.invoke(op, {"a": value.parts[0], "out": out, "status": flag})
        return Numeric("real", (out,), self.flags(value.status, flag))

    def atan2(self, y: Numeric, x: Numeric) -> Numeric:
        """Decompose the two-argument arctangent via magnitude ratio and quadrant selection.

        Calls the atan kernel with the smaller of ``|x|`` and ``|y|`` as the
        dividend, then corrects the quadrant in turn by magnitude swap, the
        sign of x and the sign of y; returns zero when both inputs are zero.

        Args:
            y: The ordinate real value.
            x: The abscissa real value.

        Returns:
            Numeric: The angle value after quadrant correction.
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
        """Compute a constant integer power by the square-and-multiply algorithm.

        Args:
            value: The base value.
            exponent: The integer exponent, with absolute value at most 128.

        Returns:
            Numeric: The power result; negative exponents are computed as reciprocals.

        Raises:
            ValidationError: The absolute value of the exponent exceeds 128.
        """
        if abs(exponent) > 128:
            raise ValidationError("integer power exceeds the generation limit of 128")
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
        """Expand the real and complex decompositions of a math intrinsic.

        rect, phase, atan2, hypot and two-argument log are handled first; the
        real path goes through sqrt arithmetic or a math kernel, the complex
        path combines real kernels with algebraic identities, and log-type
        kernels additionally flag log(0) as a domain failure.

        Args:
            name: The intrinsic function name.
            args: List of actual Numeric arguments.
            kind: The expected result kind, real or complex.

        Returns:
            Numeric: The decomposed result value.

        Raises:
            ValidationError: A complex decomposition for the function is missing.
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
            """Shorthand for adding two values."""
            return self.binary("add", a, b)

        def sub(a: Numeric, b: Numeric) -> Numeric:
            """Shorthand for subtracting two values."""
            return self.binary("sub", a, b)

        def mul(a: Numeric, b: Numeric) -> Numeric:
            """Shorthand for multiplying two values."""
            return self.binary("mul", a, b)

        def div(a: Numeric, b: Numeric) -> Numeric:
            """Shorthand for dividing two values."""
            return self.binary("div", a, b)

        def fn(k: str, v: Numeric) -> Numeric:
            """Shorthand for recursively invoking an intrinsic on the complex path."""
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
            # log(0) is explicitly flagged as a domain failure even though the totalized polynomial still has an output.
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
        raise ValidationError("missing complex math decomposition: " + name)

    def finish(self, outputs: list[tuple[tuple[Ref, ...], Numeric]], status: Ref) -> Operation:
        """Finish the build after copying outputs and uncomputing all temporary registers.

        Args:
            outputs: A sequence of (target reference tuple, output value Numeric) pairs.
            status: The 2-bit target register the aggregated status is written to.

        Returns:
            Operation: The fully built reversible module.
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
    """Generate the Chebyshev approximation kernel module of one elementary real function.

    Samples degree+1 points on the configured interval to obtain the
    Chebyshev coefficients, combines them into a polynomial circuit via the
    Clenshaw recurrence, and sets the corresponding status bits for interval
    overflow and domain errors.

    Args:
        name: Math function name; the same-named math function must exist and have a configured approximation interval.
        fmt: The fixed-point format.
        config: The MathConfig generation configuration.

    Returns:
        Operation: The reversible module with interface a, out, status; the sampling coefficient recipe
            is recorded in the math_approximation attribute.

    Raises:
        ValidationError: The approximation interval crosses the function domain.
    """
    lo, hi = config.bounds(name)
    degree = config.degree
    count = degree + 1
    center, half = (lo + hi) / 2, (hi - lo) / 2
    theta = [math.pi * (j + 0.5) / count for j in range(count)]
    try:
        samples = [getattr(math, name)(center + half * math.cos(t)) for t in theta]
    except (ValueError, OverflowError) as exc:
        raise ValidationError("math kernel interval crosses the function domain: " + name) from exc
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
