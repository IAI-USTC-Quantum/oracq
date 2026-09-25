"Reversible arithmetic generation of one-dimensional Euler frozen-Roe matrix elements; physical and fixed-point correctness pending verification."

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.common.arithmetic import BooleanNetwork, FixedFormat, fixed_arithmetic
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Adjoint, Bits, Ref


@dataclass(frozen=True)
class Word:
    """A word in fixed-point arithmetic expressions; overloaded operators delegate to the owning builder.

    Attributes:
        algebra: The ``ArithmeticBuilder`` that produced this word and emits the actual arithmetic nodes.
        ref: The register reference carrying the value.
    """

    algebra: ArithmeticBuilder
    ref: Ref

    def __add__(self, other: float | Word) -> Word:
        """Emit an arithmetic node adding two words."""
        return self.algebra.calc("add", self, other)

    def __radd__(self, other: float | Word) -> Word:
        """Reflected addition, delegated to ``self + other``."""
        return self + other

    def __sub__(self, other: float | Word) -> Word:
        """Emit an arithmetic node subtracting two words."""
        return self.algebra.calc("sub", self, other)

    def __rsub__(self, other: float | Word) -> Word:
        """Reflected subtraction, returning the word for left operand minus ``self``."""
        return self.algebra.word(other) - self

    def __mul__(self, other: float | Word) -> Word:
        """Emit an arithmetic node multiplying two words."""
        return self.algebra.calc("mul", self, other)

    def __rmul__(self, other: float | Word) -> Word:
        """Reflected multiplication, delegated to ``self * other``."""
        return self * other

    def __truediv__(self, other: float | Word) -> Word:
        """Emit an arithmetic node dividing two words."""
        return self.algebra.calc("div", self, other)

    def __neg__(self) -> Word:
        """Emit an arithmetic node negating this word."""
        return self.algebra.calc("neg", self)

    def sqrt(self) -> Word:
        """Emit an arithmetic node taking the square root of this word and return the word carrying the result.

        Returns:
            Word: The new arithmetic word carrying the square root result.
        """
        return self.algebra.calc("sqrt", self)

    def abs(self) -> Word:
        """Emit an arithmetic node taking the absolute value of this word and return the word carrying the result.

        Returns:
            Word: The new arithmetic word carrying the absolute value result.
        """
        return self.algebra.calc("abs", self)


class ArithmeticBuilder:
    """Host-side expression builder; on finish only module calls or gates and private registers remain."""

    def __init__(self, builder: Builder, fmt: FixedFormat) -> None:
        """Initialize the builder and record the start of the host instruction stream.

        Args:
            builder: The host instruction builder that arithmetic nodes are emitted into.
            fmt: The fixed-point format used by the fixed-point arithmetic.
        """
        self.b: Builder = builder
        self.fmt: FixedFormat = fmt
        self.count: int = 0
        self.constants: dict[int, Word] = {}
        self.flags: list[Ref] = []
        self.start: int = len(builder._frames[0])

    def local(self, width: int | None = None) -> Ref:
        """Allocate a new private register for intermediate values or constants.

        Args:
            width: Bit width; defaults to the word width of the current fixed-point format.

        Returns:
            The newly allocated register reference, named ``v_`` plus an incrementing counter.
        """
        self.count += 1
        return self.b.local(
            "v_" + str(self.count), Bits(self.fmt.width if width is None else width)
        )

    def word(self, value: Word | Ref | float) -> Word:
        """Wrap various inputs into a ``Word`` uniformly.

        A ``Word`` is returned as-is; a register reference with a ``parts``
        attribute is wrapped directly; numeric values are encoded in the current
        fixed-point format, and constants with the same encoding are materialized
        once and share one register.

        Args:
            value: A ``Word``, a register reference, or a numeric value encodable in the fixed-point format.

        Returns:
            Word: The arithmetic word equivalent to the input.
        """
        if isinstance(value, Word):
            return value
        if hasattr(value, "parts"):
            return Word(self, cast(Ref, value))
        raw = self.fmt.encode(value)
        if raw not in self.constants:
            ref = self.local()
            for i in range(ref.width):
                if (raw >> i) & 1:
                    self.b.x(ref[i])
            self.constants[raw] = Word(self, ref)
        return self.constants[raw]

    def calc(self, kind: str, *args: Word | Ref | float, select: Ref | None = None) -> Word:
        """Emit one fixed-point arithmetic node and return its output word.

        When the same register is reused as an operand, it is first copied with
        XOR into a fresh register and restored after the call, supporting cases
        such as squaring where a value participates in its own operation; the
        node status register is recorded in ``flags`` and aggregated by
        ``finish``.

        Args:
            kind: Operation type name, e.g. ``add``, ``mul``, ``sqrt``.
            *args: Operands; each may be a ``Word``, a register, or a numeric value.
            select: The control bit for the ``select`` operation; unused by other operations.

        Returns:
            Word: The word carrying the operation result.
        """
        op = fixed_arithmetic(kind, self.fmt)
        refs: list[Ref] = []
        copied: list[tuple[Ref, Ref]] = []
        for arg in args:
            ref = self.word(arg).ref
            if ref in refs:
                clone = self.local()
                self.b.xor(ref, clone)
                copied.append((ref, clone))
                ref = clone
            refs.append(ref)
        out_width = next(r.type.width for r in op.module.registers if r.name == "out")
        out, status = self.local(out_width), self.local(2)
        kwargs = dict(zip(("a", "b"), refs, strict=False))
        kwargs.update(out=out, status=status)
        if select is not None:
            kwargs["select"] = select
        self.b.call(op, **kwargs)  # type: ignore[arg-type]  # dynamic keyword dispatch: mypy cannot rule out the resources parameter
        for source, clone in reversed(copied):
            self.b.xor(source, clone)
        self.flags.append(status)
        return Word(self, out)

    def choose(self, bit: Ref, yes: Word | Ref | float, no: Word | Ref | float) -> Word:
        """Select between two words by a single control bit: value 1 picks ``yes`` and 0 picks ``no``.

        Args:
            bit: A one-bit control register view.
            yes: The word, register or numeric value selected when the control bit is 1.
            no: The word, register or numeric value selected when the control bit is 0.

        Returns:
            Word: The arithmetic word carrying the selected value.
        """
        return self.calc("select", yes, no, select=bit)

    def choose3(self, index: Ref, values: Sequence[Word | Ref | float]) -> Word:
        """Select among three candidate words with a two-bit index.

        Args:
            index: The two-bit control sequence, least significant bit first.
            values: The three candidate words, corresponding to indices 0..2 in order.

        Returns:
            Word: The selected word; the result is zero when the index is 3.
        """
        low = self.choose(index[0], values[1], values[0])
        high = self.choose(index[0], 0, values[2])
        return self.choose(index[1], high, low)

    def finish(
        self, outputs: Iterable[tuple[Ref, Word | Ref | float]], status: Ref
    ) -> Operation:
        """Finalize and return the resulting operation.

        Each output word is XORed into its target register; the status bits of
        the arithmetic nodes are bitwise-ORed into ``status``; the forward
        instructions emitted during construction are replayed in reverse via
        ``Adjoint``, erasing all intermediate values.

        Args:
            outputs: Pairs of target register and result word.
            status: The two-bit register receiving the aggregated status.

        Returns:
            Operation: The complete generated operation.
        """
        outputs = [(target, self.word(word)) for target, word in outputs]
        forward = tuple(self.b._frames[0][self.start :])
        for target, word in cast(list[tuple[Ref, Word]], outputs):
            self.b.xor(word.ref, target)
        # The output status is the OR of the arithmetic node statuses; XOR would not cancel duplicated errors.
        if self.flags:
            net = BooleanNetwork()
            all_flags = [net.input("f_" + str(i), 2) for i in range(len(self.flags))]
            net.outputs = {"status": [net.any([f[j] for f in all_flags]) for j in range(2)]}
            op = net.operation()
            self.b.call(
                op,
                status=status,
                # dynamic keyword dispatch: mypy cannot rule out the resources parameter
                **{"f_" + str(i): flag for i, flag in enumerate(self.flags)},  # type: ignore[arg-type]
            )
        self.b.emit(Adjoint(forward))
        return self.b.finish()


def roe_face(
    *,
    fmt: FixedFormat | None = None,
    gamma: float = 1.4,
    entropy_delta: float = 0.125,
) -> Operation:
    """Automatic compilation of a plain Python Roe function, keeping the existing six-field plus row and column input ABI.

    Args:
        fmt: The fixed-point format for field values and fluxes; defaults to ``FixedFormat(10, 5)``.
        gamma: Ratio of specific heats.
        entropy_delta: Harten entropy fix threshold.

    Returns:
        Operation: The reversible face-flux computation module auto-compiled from
        ``frozen_roe_face`` with fixed-point arithmetic.
    """
    from dataclasses import replace

    from oracq.applications.roe_formulas import frozen_roe_face
    from oracq.infrastructure.mathfunc import Index, compile_function

    fmt = fmt or FixedFormat(10, 5)
    compiled = compile_function(
        frozen_roe_face,
        fmt=fmt,
        inputs={
            **{key: "real" for key in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r")},
            "row": Index(2),
            "col": Index(2),
        },
        constants={"gamma": gamma, "entropy_delta": entropy_delta},
        output_names=("left", "right"),
    )
    # CompiledFunction.operation is annotated as object on the infrastructure side; narrowed here by its runtime structure.
    attrs = dict(cast(Operation, compiled.operation).module.attributes)
    attrs.update(
        algorithm="frozen_roe_1d_euler",
        matrix_entry_source="compiled_pure_function_on_raw_conserved_variables",
    )
    return replace(
        cast(Operation, compiled.operation),
        module=replace(
            cast(Operation, compiled.operation).module,
            attributes=tuple(sorted(attrs.items())),
        ),
    )
