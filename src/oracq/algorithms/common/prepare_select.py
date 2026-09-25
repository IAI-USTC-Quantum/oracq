"Standard PREPARE–SELECT decomposition of an LCU (Low & Chuang 2019; Babbush et al. 2018 alias sampling)."

from __future__ import annotations

import cmath
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from oracq.algorithms.common.arithmetic import FixedFormat, fixed_arithmetic
from oracq.algorithms.common.hamiltonian import PauliHamiltonian
from oracq.algorithms.input_model.block_encoding import pauli_word
from oracq.algorithms.input_model.contracts import positive_integer, require_instance
from oracq.algorithms.input_model.interfaces import (
    StatePreparationProtocol,
    as_state_preparation,
)
from oracq.algorithms.input_model.operators import BlockEncoding, _name, scale
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    XorDatabase,
    annotate,
    declare,
    gate_state_prep,
    invoke,
    qram_database,
    qram_state_angles,
    qram_state_prep,
    resources_for,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError


def _normalized(
    coefficients: Iterable[complex],
) -> tuple[tuple[complex, ...], float, int, list[float]]:
    """Validate the coefficients and return, coefficients, alpha equal to the l1 norm,
    the selector bit width, and the ``√|c|/√α`` amplitudes."""
    values = tuple(complex(c) for c in coefficients)
    if len(values) < 2:
        raise ValidationError("PREPARE requires at least two coefficients; use scale directly"
                              " for a single term")
    if not all(math.isfinite(c.real) and math.isfinite(c.imag) for c in values):
        raise ValidationError("LCU coefficients must be finite")
    alpha = sum(abs(c) for c in values)
    if not alpha:
        raise ValidationError("LCU coefficients must not all be zero")
    width = (len(values) - 1).bit_length()
    amplitudes = [math.sqrt(abs(c) / alpha) for c in values]
    amplitudes += [0.0] * ((1 << width) - len(amplitudes))
    return values, alpha, width, amplitudes


def _pauli_terms(
    terms: Iterable[tuple[complex, str]] | PauliHamiltonian,
) -> tuple[tuple[complex, str], ...]:
    """Accept either a sequence of, coefficient and Pauli word, pairs or a
    PauliHamiltonian; zero coefficients are kept to align indices."""
    from oracq.algorithms.common.hamiltonian import PauliHamiltonian

    if isinstance(terms, PauliHamiltonian):
        terms = terms.terms
    terms = tuple((complex(c), str(w)) for c, w in terms)
    if not terms:
        raise ValidationError("The LCU needs at least one Pauli term")
    for coefficient, word in terms:
        if not (math.isfinite(coefficient.real) and math.isfinite(coefficient.imag)):
            raise ValidationError("Pauli coefficients must be finite")
        if not word or any(letter not in "IXYZ" for letter in word):
            raise ValidationError("A Pauli word may only contain I/X/Y/Z")
    if len({len(word) for _, word in terms}) != 1:
        raise ValidationError("Inconsistent Pauli term widths")
    return terms


def _prepare_attributes(
    values: Sequence[complex], alpha: float, width: int
) -> dict[str, int | float]:
    """Collect the term count, alpha and selector bit width that PREPARE writes into
    module attributes."""
    return {
        "prepare_terms": len(values),
        "prepare_alpha": alpha,
        "selector_width": width,
    }


def abstract_prepare(
    coefficients: Iterable[complex], *, work_width: int = 0, name: str | None = None
) -> StatePreparation:
    """Open declaration of PREPARE: the coefficients go into declaration attributes, the
    body is empty, and bind defers the implementation.

    work_width must match the work width of the implementation bound later,
    0 for the gate variant and the selector bit width plus the angle bit width
    for the QRAM variant, consistently with the usage of
    abstract_state_prep in the catalog.

    Args:
        coefficients: LCU complex coefficient sequence, with at least two finite values,
            not all zero.
        work_width: Bit width of the work register, in 0..64; must match the work width
            of the implementation bound later.
        name: Declaration module name; generated automatically from the coefficients and
            width when omitted.

    Returns:
        StatePreparation: State preparation handle wrapping the open PREPARE declaration
        module.
    """
    values, alpha, width, _ = _normalized(coefficients)
    positive_integer(work_width, "abstract_prepare.work_width", minimum=0)
    return StatePreparation(
        declare(
            name or _name("prepare_abstract", values, work_width),
            {"target": Bits(width), "work": Bits(work_width)},
            paradigm="state_prep_isometry",
            attributes={
                "zero_input": True,
                "clean_work": True,
                **_prepare_attributes(values, alpha, width),
            },
        )
    )


def gate_prepare(
    coefficients: Iterable[complex], *, name: str | None = None
) -> StatePreparation:
    """Gate-level PREPARE: a multiplexed Ry with amplitudes proportional to ``√|c_i|``;
    coefficient phases enter SELECT following the repository convention.

    Args:
        coefficients: LCU complex coefficient sequence, with at least two finite values,
            not all zero.
        name: Name of the generated state preparation module; generated automatically
            from the coefficients when omitted.

    Returns:
        StatePreparation: Handle of the gate-level state preparation operation.
    """
    values, alpha, width, amplitudes = _normalized(coefficients)
    prep = gate_state_prep(amplitudes, name=name)
    return StatePreparation(
        annotate(
            prep.operation,
            "state_prep_isometry",
            zero_input=True,
            clean_work=True,
            **_prepare_attributes(values, alpha, width),
        )
    )


@dataclass(frozen=True)
class QramPreparation:
    """QRAM-variant PREPARE handle; the angle table is provided in memory and never
    enters the IR."""

    preparation: StatePreparation
    memory: dict[str, dict[int, int]]

    def state_preparation(self) -> StatePreparation:
        """Return the ``StatePreparation`` operation wrapped inside the handle.

        Returns:
            StatePreparation: The state preparation operation wrapped inside the handle.
        """
        return self.preparation


def qram_prepare(coefficients: Iterable[complex], *, angle_width: int = 8) -> QramPreparation:
    """QRAM-resource PREPARE: the rotation angle table is bound by the caller as QRAM
    data.

    Args:
        coefficients: LCU complex coefficient sequence, with at least two finite values,
            not all zero.
        angle_width: Quantization bit width of the rotation angles, in 2..32.

    Returns:
        QramPreparation: Handle holding the PREPARE operation and the default angle
        table QRAM data.
    """
    values, alpha, width, amplitudes = _normalized(coefficients)
    positive_integer(angle_width, "qram_prepare.angle_width", minimum=2, maximum=32)
    prep = qram_state_prep(width, angle_width)
    memory = {"angles": qram_state_angles(amplitudes, angle_width)}
    operation = annotate(
        prep.operation,
        "state_prep_isometry",
        zero_input=True,
        clean_work=True,
        **_prepare_attributes(values, alpha, width),
    )
    return QramPreparation(StatePreparation(operation), memory)


@dataclass(frozen=True)
class AliasTable:
    """Classical preprocessing result of Babbush et al. alias sampling; words are packed
    as ``keep | (alt << precision)``."""

    probabilities: tuple[float, ...]
    keep: tuple[float, ...]
    alt: tuple[int, ...]
    quantized: tuple[int, ...]
    precision: int
    table: dict[int, int]

    def distribution(self) -> tuple[float, ...]:
        """Classical sampling distribution under quantized keep and uniform drawing.

        Returns:
            tuple[float, ...]: Sampling probability of each slot, with length equal to
            the slot count and total sum one.
        """
        scale = 1 << self.precision
        size = len(self.keep)
        result = [0] * size
        for i in range(size):
            result[i] += self.quantized[i]
            result[self.alt[i]] += scale - self.quantized[i]
        return tuple(v / (size * scale) for v in result)


def alias_table(coefficients: Iterable[complex], *, precision: int = 8) -> AliasTable:
    """Vose alias preprocessing: padded to 2^selector, then split into keep and alt by
    the mean; zero-probability slots also work.

    Args:
        coefficients: LCU complex coefficient sequence, with at least two finite values,
            not all zero.
        precision: Quantization bit width of the keep probabilities, in 2..32; its sum
            with the selector bit width must not exceed 64.

    Returns:
        AliasTable: Preprocessing result holding the probabilities, the keep and alt
        tables, the quantized words and the packed data table.
    """
    values, alpha, width, _ = _normalized(coefficients)
    positive_integer(precision, "alias_table.precision", minimum=2, maximum=32)
    size = 1 << width
    if width + precision > 64:
        raise ValidationError("The alias data word exceeds the QRAM bit width limit")
    probabilities = [abs(c) / alpha for c in values] + [0.0] * (size - len(values))
    scaled = [p * size for p in probabilities]
    keep, alt = [0.0] * size, list(range(size))
    small = [i for i, v in enumerate(scaled) if v < 1]
    large = [i for i, v in enumerate(scaled) if v >= 1]
    while small and large:
        lo, hi = small.pop(), large.pop()
        keep[lo] = scaled[lo]
        alt[lo] = hi
        scaled[hi] += scaled[lo] - 1
        (small if scaled[hi] < 1 else large).append(hi)
    for i in small + large:
        keep[i] = 1.0
    quantized = [min((1 << precision) - 1, math.floor(k * (1 << precision))) for k in keep]
    table = {i: quantized[i] | (alt[i] << precision) for i in range(size)}
    return AliasTable(
        tuple(probabilities), tuple(keep), tuple(alt), tuple(quantized), precision, table
    )


@dataclass(frozen=True)
class AliasPreparation:
    """Alias-sampling PREPARE handle; table is the classical table and memory is the
    data of the default QRAM binding."""

    preparation: StatePreparation
    table: AliasTable
    memory: dict[str, dict[int, int]]

    def state_preparation(self) -> StatePreparation:
        """Return the ``StatePreparation`` operation wrapped inside the handle.

        Returns:
            StatePreparation: The state preparation operation wrapped inside the handle.
        """
        return self.preparation


def alias_prepare(
    coefficients: Iterable[complex],
    *,
    precision: int = 8,
    database: XorDatabase | None = None,
) -> AliasPreparation:
    """Babbush et al. alias-sampling PREPARE: uniform state plus keep and alt loading,
    a comparator and a controlled swap.

    The data table is injected through the XorDatabase interface, by default
    the QRAM resource but a gate_database also works, and is not embedded in
    the IR. The data and comparison bits in work entangle with the selector,
    forming the purification junk of ``Σ√p_i|i>|junk_i>``; the (0,0) block of
    the block encoding is unaffected by inner products over the junk, and the
    junk is cleaned by the adjoint of PREPARE. keep is quantized by rounding
    down to precision bits, with total variation from the exact distribution
    at most 2^selector·2^-precision.

    Args:
        coefficients: LCU complex coefficient sequence, with at least two finite values,
            not all zero.
        precision: Quantization bit width of the keep probabilities, in 2..32.
        database: XorDatabase handle carrying the keep alt packed data table; defaults
            to the QRAM resource, whose address bit width must equal the selector bit
            width and whose data bit width equals precision plus the selector bit width.

    Returns:
        AliasPreparation: Handle holding the PREPARE operation, the classical alias
        table and the default QRAM data.
    """
    values, alpha, width, _ = _normalized(coefficients)
    positive_integer(precision, "alias_prepare.precision", minimum=2, maximum=32)
    table = alias_table(coefficients, precision=precision)
    data_width = precision + width
    if database is None:
        database = qram_database(width, data_width)
        memory = {
            f"db__{r.name}": dict(table.table) for r in database.operation.module.resources
        }
    else:
        require_instance(database, XorDatabase, "alias_prepare.database")
        memory = {}
    if database.address_width != width or database.data_width != data_width:
        raise ValidationError("The alias data table requires the address width to equal the"
                              " selector bit width and the data width to equal the packed"
                              " keep alt word width")
    compare = fixed_arithmetic("lt", FixedFormat(precision, 0, signed=False))
    b = Builder(
        _name("alias_prepare", values, precision, database.operation),
        {"target": Bits(width), "work": Bits(data_width + precision)},
        resources_for(("db", database.operation), ("lt", compare)),
    )
    data, coin = b["work"][:data_width], b["work"][data_width:]
    keep, alt = data[:precision], data[precision:]
    flag, status = b.local("flag", Bits(1)), b.local("status", Bits(2))
    b.h(b["target"])
    b.h(coin)
    invoke(b, database.operation, "db", address=b["target"], data=data)
    invoke(b, compare, "lt", a=coin, b=keep, out=flag, status=status)
    with b.control(flag, 0):
        b.swap(b["target"], alt)
    # The compare inputs are unchanged by the swap, so calling it again cleans the flag;
    # data and coin stay in work as purification junk.
    invoke(b, compare, "lt", a=coin, b=keep, out=flag, status=status)
    operation = annotate(
        b.finish(),
        "state_prep_isometry",
        zero_input=True,
        clean_work=False,
        implementation="alias_sampling",
        alias_precision=precision,
        **_prepare_attributes(values, alpha, width),
    )
    return AliasPreparation(StatePreparation(operation), table, memory)


def select_pauli(terms: Iterable[tuple[complex, str]] | PauliHamiltonian) -> Operation:
    """SELECT: when selector equals i, apply the i-th Pauli word to target, and also
    apply the coefficient phase.

    The control condition is expressed with the RIR Control primitive on the
    binary value of the selector; when unary iteration is needed, the backend
    lowers the multi-bit control, and the generation stage does not expand it.

    Args:
        terms: Sequence of, coefficient and Pauli word, pairs or a PauliHamiltonian; at
            least two terms, whose Pauli words contain only I/X/Y/Z and have equal
            widths.

    Returns:
        Operation: SELECT operation applying the Pauli word and coefficient phase
        corresponding to the selector value.
    """
    terms = _pauli_terms(terms)
    if len(terms) < 2:
        raise ValidationError("SELECT requires at least two Pauli terms")
    width = len(terms[0][1])
    selector_width = (len(terms) - 1).bit_length()
    b = Builder(
        _name("select_pauli", terms),
        {"selector": Bits(selector_width), "target": Bits(width)},
    )
    for index, (coefficient, word) in enumerate(terms):
        with b.control(b["selector"], index):
            if cmath.phase(coefficient):
                b.global_phase(cmath.phase(coefficient))
            for bit, letter in enumerate(word):
                if letter != "I":
                    b.gate(letter.lower(), b["target"][bit])
    return annotate(
        b.finish(),
        "unitary",
        algorithm="select_pauli",
        select_terms=len(terms),
        selector_width=selector_width,
    )


def lcu_prepare_select(
    terms: Iterable[tuple[complex, str]] | PauliHamiltonian,
    *,
    prepare: StatePreparationProtocol | None = None,
) -> BlockEncoding:
    """Standard PREPARE–SELECT block encoding: (PREPARE†⊗I)·SELECT·(PREPARE⊗I).

    The low selector_width bits of signal are the selector and the high bits
    are the work of PREPARE, the purification junk for the alias variant;
    prepare defaults to gate_prepare, and an abstract, qram or alias handle
    also works. A single term degenerates to scale. The result can be fed
    directly to transforms.qubitization_walk.

    Args:
        terms: Sequence of, coefficient and Pauli word, pairs or a PauliHamiltonian;
            the Pauli words contain only I/X/Y/Z and have equal widths.
        prepare: PREPARE handle, defaulting to gate_prepare; it must satisfy the
            zero_input contract and its target width must equal the selector bit width.

    Returns:
        BlockEncoding: (PREPARE†⊗I)·SELECT·(PREPARE⊗I) block encoding with scale equal
        to the l1 norm of the coefficients.
    """
    terms = _pauli_terms(terms)
    if len(terms) == 1:
        return scale(terms[0][0], pauli_word(terms[0][1]))
    coefficients = tuple(c for c, _ in terms)
    _, alpha, selector_width, _ = _normalized(coefficients)
    prep = gate_prepare(coefficients) if prepare is None else as_state_preparation(prepare)
    if prep.width != selector_width:
        raise ValidationError("The PREPARE target width does not match the selector bit width")
    if dict(prep.operation.module.attributes).get("zero_input") is not True:
        raise ValidationError("PREPARE must satisfy the zero_input contract")
    select = select_pauli(terms)
    width = len(terms[0][1])
    junk = prep.work_width
    b = Builder(
        _name("lcu_prepare_select", select, prep.operation),
        {"target": Bits(width), "signal": Bits(selector_width + junk)},
        resources_for(("prep", prep.operation), ("select", select)),
    )
    selector, work = b["signal"][:selector_width], b["signal"][selector_width:]
    invoke(b, prep.operation, "prep", target=selector, work=work)
    invoke(b, select, "select", selector=selector, target=b["target"])
    with b.adjoint():
        invoke(b, prep.operation, "prep", target=selector, work=work)
    return BlockEncoding(
        annotate(
            b.finish(),
            "block_encoding",
            be_alpha=alpha,
            be_form="prepare_select",
            lcu_terms=len(terms),
            selector_width=selector_width,
            prepare_work=junk,
        )
    )
