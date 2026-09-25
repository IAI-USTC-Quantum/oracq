"""Hamiltonian decomposition, Pauli evolution, Trotter and truncated Taylor encoding."""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import Protocol, cast, runtime_checkable

from oracq.algorithms.input_model.block_encoding import lcu, pauli_word
from oracq.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
    requires,
)
from oracq.algorithms.input_model.interfaces import BlockEncodingProtocol, as_block_encoding
from oracq.algorithms.input_model.operators import (
    BlockEncoding,
    _name,
    block_encoding,
    identity,
    product,
)
from oracq.algorithms.input_model.oracles import (
    annotate,
    invoke,
    resources_for,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError


def trotter_hamsim(
    terms: Iterable[tuple[float, str]], final_time: float, *, steps: int = 2
) -> Operation:
    """First-order Trotter product-formula simulation of the evolution under a
    Pauli-decomposed Hamiltonian.

    For H = sum_j c_j P_j it synthesizes exp(-i*H*t) as the steps-th power of
    prod_j exp(-i*c_j*P_j*t/steps): each Pauli word is folded, through basis
    changes and a CNOT chain, into a phase rotation on the last active qubit
    and then restored in reverse order. The circuit body is a single
    Repeat(steps) block and is not expanded per step at generation or
    serialization time; an all-I word consumes no qubits and degenerates to a
    global phase.

    Args:
        terms: Sequence of ``(coefficient, Pauli word)`` pairs; the coefficients are real
            numbers, converted through ``float``, and the Pauli words are equal-width
            strings of I/X/Y/Z.
        final_time: Total evolution time t.
        steps: Number of repeated steps, at least 1; more steps give a smaller
            product-formula error.

    Returns:
        Operation: A bare unitary operation rather than a block encoding; the registers are
        target, of word width, and a zero-width signal, and the whole approximates
        exp(-i*H*t).

    Raises:
        ValidationError: The term sequence is empty, steps is below 1, or the Pauli word
            widths differ.
    """
    terms = tuple((float(c), word) for c, word in terms)
    if not terms or steps < 1:
        raise ValidationError("Trotter requires a nonempty term list and a positive step count")
    width = len(terms[0][1])
    if any(len(word) != width for _, word in terms):
        raise ValidationError("Pauli terms have different widths")
    b = Builder(
        _name("trotter", terms, final_time, steps),
        {"target": Bits(width), "signal": Bits(0)},
        attributes={"algorithm": "trotter_hamsim", "validation_stage": "paradigm"},
    )
    with b.repeat(steps):
        for coefficient, word in terms:
            active = [i for i, char in enumerate(word) if char != "I"]
            if not active:
                b.global_phase(-coefficient * final_time / steps)
                continue
            for bit in active:
                if word[bit] == "Y":
                    b.gate("phase", b["target"][bit], -math.pi / 2)
                if word[bit] in "XY":
                    b.h(b["target"][bit])
            for bit in active[:-1]:
                b.xor(b["target"][bit], b["target"][active[-1]])
            b.rz(b["target"][active[-1]], 2 * coefficient * final_time / steps)
            for bit in reversed(active[:-1]):
                b.xor(b["target"][bit], b["target"][active[-1]])
            for bit in reversed(active):
                if word[bit] in "XY":
                    b.h(b["target"][bit])
                if word[bit] == "Y":
                    b.gate("phase", b["target"][bit], math.pi / 2)
    return b.finish()


def taylor_hamiltonian(
    hamiltonian: BlockEncoding, time: float, *, degree: int = 2
) -> BlockEncoding:
    """A closable plain Hamiltonian-function block encoding; replaceable by a
    QSP/HamSim protocol.

    Args:
        hamiltonian: Block encoding of the operator.
        time: Evolution time, a finite real number.
        degree: Taylor truncation order, an integer not less than 0; order 0 keeps only
            the identity term.

    Returns:
        BlockEncoding: LCU block encoding of the truncated Taylor series, with correctness
        marked pending.
    """
    require_instance(hamiltonian, BlockEncoding, "taylor_hamiltonian.H")
    positive_integer(degree, "taylor_hamiltonian.degree", minimum=0)
    finite_real(time, "taylor_hamiltonian.time")
    if degree < 0 or not math.isfinite(time):
        raise ValidationError("Invalid Taylor degree or evolution time")
    powers: list[tuple[complex, BlockEncoding]]
    current: BlockEncoding
    powers, current = [(1, identity(hamiltonian.width))], identity(hamiltonian.width)
    for k in range(1, degree + 1):
        current = product(hamiltonian, current)
        powers.append(((-1j * time) ** k / math.factorial(k), current))
    out = lcu(powers)
    return BlockEncoding(
        annotate(
            out.operation,
            "block_encoding",
            be_alpha=out.alpha,
            algorithm="truncated_taylor_hamiltonian_function",
            degree=degree,
            time=float(time),
            correctness="pending",
            success_condition="signal == 0",
        )
    )


@runtime_checkable
class HermitianProtocol(Protocol):
    """Access protocol by which the host declares the Hermiticity of an operator."""

    @property
    def hermitian(self) -> bool:
        """Whether the operator is declared Hermitian; hamiltonian_simulation accepts
        only True."""
        ...


@runtime_checkable
class EvolvableProtocol(Protocol):
    """Access protocol by which the host declares that an operator can produce its own
    unitary evolution."""

    def evolution(self, time: float) -> Operation:
        """Return the unitary Operation of this operator for the given evolution time.

        The Trotter path requires the returned evolution to need no
        postselection: apart from target, all public registers are zero width.

        Args:
            time: Evolution time, a finite real number.

        Returns:
            Operation: The unitary evolution operation of this operator for the given
            duration.
        """
        ...


@runtime_checkable
class TrotterizableProtocol(Protocol):
    """Access protocol by which the host declares that an operator decomposes into a
    Trotter term list."""

    def trotter_list(self) -> Sequence[TrotterTerm]:
        """Return the TrotterTerm sequence constituting the Hamiltonian; hamiltonian_simulation
        requires it to be nonempty.

        Returns:
            Sequence[TrotterTerm]: The product-formula term sequence constituting the
            Hamiltonian, required to be nonempty.
        """
        ...


@dataclass(frozen=True)
class TrotterTerm:
    """A single Hamiltonian term of a product formula: a real coefficient combined with
    an evolvable operator.

    Attributes:
        coefficient: Real coefficient of the term; the Trotter path calls
            ``operator.evolution`` with duration ``coefficient * time / steps``.
        operator: An operator satisfying EvolvableProtocol.

    Raises:
        ValidationError: coefficient is not a finite real number, or operator does not
            satisfy EvolvableProtocol.
    """

    coefficient: float
    operator: object

    def __post_init__(self) -> None:
        """Validate that the coefficient is a finite real number and the operator
        satisfies EvolvableProtocol."""
        finite_real(self.coefficient, "TrotterTerm.coefficient")
        requires(self.operator, EvolvableProtocol, path="TrotterTerm.operator")


@dataclass(frozen=True)
class PauliOperator:
    """A single Pauli word operator; satisfies HermitianProtocol and EvolvableProtocol.

    Attributes:
        word: A nonempty string of I/X/Y/Z with width at most 64.
        hermitian: Always True; this is the HermitianProtocol declaration.

    Raises:
        ValidationError: word is not a nonempty I/X/Y/Z string, or its width exceeds 64.
    """

    word: str
    hermitian = True

    def __post_init__(self) -> None:
        """Validate that word is an I/X/Y/Z string of width at most 64."""
        if (
            not isinstance(self.word, str)
            or not self.word
            or any(c not in "IXYZ" for c in self.word)
        ):
            raise ValidationError("A Pauli word must be a nonempty I/X/Y/Z string")
        positive_integer(len(self.word), "PauliOperator.width", maximum=64)

    def block_encoding(self) -> BlockEncoding:
        """Return the alpha=1.0 block encoding of this Pauli word: per-qubit single-qubit
        gates plus a zero-width signal.

        Returns:
            BlockEncoding: The Pauli word block encoding with scale 1.0.
        """
        return pauli_word(self.word)

    def evolution(self, time: float) -> Operation:
        """Return the exact unitary evolution of ``exp(-1j*word*time)`` (a single term in
        a single step, with no product-formula error).

        Args:
            time: Evolution time, a finite real number.

        Returns:
            Operation: The exact unitary evolution operation synthesized in a single term
            and single step.
        """
        finite_real(time, "PauliOperator.time")
        return trotter_hamsim(((1.0, self.word),), time, steps=1)


@dataclass(frozen=True)
class PauliHamiltonian:
    """A real-coefficient linear combination of equal-width Pauli words; supports both
    block-encoding and Trotter decomposition access.

    Attributes:
        terms: Tuple of ``(coefficient, Pauli word)`` pairs; the coefficients are finite
            real numbers and all words have equal width.
        hermitian: Always True; this is the HermitianProtocol declaration.

    Raises:
        ValidationError: The term list is empty, a coefficient is not a finite real number,
            a Pauli word is invalid, or the widths differ.
    """

    terms: tuple[tuple[float, str], ...]
    hermitian = True

    def __post_init__(self) -> None:
        """Normalize the term list and validate nonemptiness, finite coefficients and
        equal Pauli word widths."""
        object.__setattr__(self, "terms", tuple(tuple(t) for t in self.terms))
        if not self.terms:
            raise ValidationError("PauliHamiltonian requires a nonempty term list")
        for coefficient, word in self.terms:
            finite_real(coefficient, "PauliHamiltonian.coefficient")
            PauliOperator(word)
        if len({len(word) for _, word in self.terms}) != 1:
            raise ValidationError("Inconsistent Pauli term widths")

    def block_encoding(self) -> BlockEncoding:
        """Return the LCU block encoding of the nonzero-coefficient terms; with no
        nonzero term it degenerates to the zero-operator block encoding.

        Returns:
            BlockEncoding: LCU block encoding with scale equal to the sum of absolute
            coefficient values, or the zero operator with no nonzero term.
        """
        from oracq.algorithms.input_model.operators import zero

        terms = [(c, pauli_word(w)) for c, w in self.terms if c]
        return lcu(terms) if terms else zero(len(self.terms[0][1]))

    def trotter_list(self) -> tuple[TrotterTerm, ...]:
        """Wrap each ``(coefficient, word)`` into a TrotterTerm and return the tuple.

        Returns:
            tuple[TrotterTerm, ...]: Tuple of Trotter terms wrapping each coefficient and
            Pauli word one by one.
        """
        return tuple(TrotterTerm(c, PauliOperator(w)) for c, w in self.terms)


@dataclass(frozen=True)
class EncodedOperator:
    """Host declaration of matrix properties; may represent non-Hermitian operators and
    does not claim to itself be a unitary."""

    encoding: BlockEncoding
    hermitian: bool

    def __post_init__(self) -> None:
        """Adapt the encoding to ``BlockEncoding`` and validate that hermitian is a bool."""
        object.__setattr__(self, "encoding", as_block_encoding(self.encoding))
        if type(self.hermitian) is not bool:
            raise ValidationError("hermitian must be declared as a bool")

    def block_encoding(self) -> BlockEncoding:
        """Return the block encoding carried at construction.

        Returns:
            BlockEncoding: The block encoding passed at construction, adapted and
            validated.
        """
        return self.encoding


def hamiltonian_simulation(
    operator: HermitianProtocol,
    time: float,
    *,
    method: str = "auto",
    steps: int = 2,
    qsp: Callable[[BlockEncoding, float], BlockEncoding] | None = None,
) -> BlockEncoding:
    """Currently prefers the decomposable Trotter path; QSP requires the caller to
    supply an actual implementation.

    Args:
        operator: An operator declaring Hermitian, which must additionally satisfy the
            Trotter-decomposable or block-encodable protocol.
        time: Evolution time, a finite real number.
        method: Path selection, one of ``auto``, ``trotter`` or ``qsp``; ``auto`` picks
            by the operator protocol.
        steps: Number of repeated segments of the Trotter product formula, a positive
            integer.
        qsp: An actual implementation of the form qsp(BE, time) returning a
            BlockEncoding; needed only by the ``qsp`` path.

    Returns:
        BlockEncoding: Block encoding of the evolution exp(-iHt); the Trotter path has
        scale 1.0.
    """
    requires(operator, HermitianProtocol, path="HamSim.operator")
    if operator.hermitian is not True:
        raise ValidationError(
            "Hamiltonian simulation requires a Hermitian operator; use QODE for"
            " non-Hermitian dynamics"
        )
    finite_real(time, "HamSim.time")
    if method not in {"auto", "trotter", "qsp"}:
        raise ValidationError("Unknown Hamiltonian simulation method")
    if method == "auto":
        method = "trotter" if isinstance(operator, TrotterizableProtocol) else "qsp"
    if method == "qsp":
        requires(operator, BlockEncodingProtocol, path="QSP.operator")
        if not callable(qsp):
            raise ValidationError(
                "The QSP path requires injecting an actual qsp implementation of a block"
                " encoding and a time; this library has no generic QSP-HamSim kernel"
            )
        encoded = as_block_encoding(cast("BlockEncodingProtocol", operator))
        result = qsp(encoded, time)
        require_instance(result, BlockEncoding, "QSP.output")
        if result.width != encoded.width:
            raise ValidationError("QSP returned an inconsistent target width")
        return result
    requires(operator, TrotterizableProtocol, path="Trotter.operator")
    positive_integer(steps, "Trotter.steps")
    terms = tuple(cast("TrotterizableProtocol", operator).trotter_list())
    if not terms:
        raise ValidationError("trotter_list must return a nonempty term list")
    evolutions: list[Operation] = []
    for term in terms:
        require_instance(term, TrotterTerm, "Trotter.term")
        # Each term must provide a postselection-free unitary evolution; a generic
        # polynomial block encoding cannot substitute for it.
        evolution = cast("EvolvableProtocol", term.operator).evolution(
            term.coefficient * time / steps
        )
        from oracq.infrastructure.builder import Operation

        require_instance(evolution, Operation, "Trotter.term.evolution")
        if any(r.type.width for r in evolution.module.registers if r.name != "target"):
            raise ValidationError(
                "The current Trotter implementation requires each term evolution to expose"
                " only target; other public registers must have zero width"
            )
        if not any(r.name == "target" for r in evolution.module.registers):
            raise ValidationError("Each Trotter term evolution needs a target register")
        evolutions.append(evolution)
    widths = {
        next(r.type.width for r in op.module.registers if r.name == "target") for op in evolutions
    }
    if len(widths) != 1:
        raise ValidationError("Inconsistent target widths among Trotter term evolutions")
    width = widths.pop()
    b = Builder(
        _name("trotter_protocol", *evolutions, steps),
        {"target": Bits(width), "signal": Bits(0)},
        resources_for(*[(f"term{i}", op) for i, op in enumerate(evolutions)]),
        attributes={
            "algorithm": "trotter_hamiltonian_protocol",
            "steps": steps,
            "correctness": "product formula approximation pending",
        },
    )
    with b.repeat(steps):
        for i, op in enumerate(evolutions):
            arguments = {
                r.name: b["target"] if r.name == "target" else b["signal"]
                for r in op.module.registers
            }
            invoke(b, op, f"term{i}", **arguments)
    return block_encoding(b.finish(), alpha=1.0)
