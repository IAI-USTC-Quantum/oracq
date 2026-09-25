"""Phase estimation, amplitude estimation and overlap measurement circuits; readout
and statistics are done on the host side."""

from __future__ import annotations

import math
from collections.abc import Iterable

from oracq.algorithms.common.fourier import qft
from oracq.algorithms.input_model.contracts import positive_integer
from oracq.algorithms.input_model.interfaces import (
    BlockEncodingProtocol,
    StatePreparationProtocol,
    as_block_encoding,
    checked_state_preparation,
)
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import basis_state, invoke, resources_for
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError


def phase_estimation(operation: Operation, *, precision: int = 2) -> Operation:
    """Standard quantum phase estimation, QPE.

    Args:
        operation: A complete Operation supporting controlled invocation. The caller
            must prepare its input state.
        precision: Number of phase register bits, in the range 1..63.

    Returns:
        Operation: Keeps the original public registers and appends phase. For an
        eigenphase φ, the readout approximates 2**precision * φ.

    Raises:
        ValidationError: The precision is invalid, the input occupies the phase name,
            or the required controlled invocation is unsupported.

    The powers are kept in a Repeat; the generator does not expand the gate
    sequence per power."""
    positive_integer(precision, "phase_estimation.precision", maximum=63)
    registers = {r.name: r.type for r in operation.module.registers}
    if "phase" in registers:
        raise ValidationError("The invoked interface already occupies the parameter name phase")
    b = Builder(
        _name("qpe", operation, precision),
        {**registers, "phase": Bits(precision)},
        resources_for(("u", operation)),
        attributes={"algorithm": "qpe"},
    )
    b.h(b["phase"])
    for bit in range(precision):
        with b.control(b["phase"][bit]):
            with b.repeat(1 << bit):
                invoke(b, operation, "u", **{name: b[name] for name in registers})
    with b.adjoint():
        invoke(b, qft(precision), target=b["phase"])
    return b.finish()


def hadamard_test(
    unitary: BlockEncodingProtocol,
    preparation: StatePreparationProtocol | None = None,
    *,
    component: str = "real",
) -> Operation:
    """Generate a Hadamard test circuit for a complex expectation value.

    Args:
        unitary: A complete unitary operation, or a block encoding with no signal and
            alpha=1.
        preparation: Initial state preparation; when omitted the zero state of the
            target space is used. work must be cleaned.
        component: ``real`` or ``imag``, selecting the real or imaginary part of the
            expectation value to read out.

    Returns:
        Operation: Registers target, work and probe. The Z expectation of probe is the
        corresponding complex expectation component.

    This operation performs no measurement. Actual use requires sampling probe and
    computing the probability difference on the classical side."""
    encoded = as_block_encoding(unitary)
    if encoded.signal_qubits or encoded.alpha != 1:
        raise ValidationError("The Hadamard test currently requires a full unitary input"
                              " with no signal and alpha=1")
    prep = (
        basis_state(encoded.width)
        if preparation is None
        else checked_state_preparation(preparation)
    )
    if prep.width != encoded.width or component not in {"real", "imag"}:
        raise ValidationError("Invalid initial state width or component option for the"
                              " Hadamard test")
    b = Builder(
        _name("hadamard_test", encoded.operation, prep.operation, component),
        {"target": Bits(prep.width), "work": Bits(prep.work_width), "probe": Bits(1)},
        resources_for(("u", encoded.operation), ("prep", prep.operation)),
        attributes={
            "algorithm": "hadamard_test",
            "readout_register": "probe",
            "component": component,
        },
    )
    invoke(b, prep.operation, "prep", target=b["target"], work=b["work"])
    b.h(b["probe"])
    with b.control(b["probe"]):
        invoke(b, encoded.operation, "u", target=b["target"], signal=b["work"][:0])
    if component == "imag":
        b.gate("phase", b["probe"], -math.pi / 2)
    b.h(b["probe"])
    return b.finish()


def swap_test(first: StatePreparationProtocol, second: StatePreparationProtocol) -> Operation:
    """Generate an overlap measurement circuit for two pure states.

    Args:
        first: Preparation of the first state; requires zero input and a clean
            workspace.
        second: Preparation of the second state of the same width; requires zero input
            and a clean workspace.

    Returns:
        Operation: Keeps left, right, both work registers and probe. The probability
        that probe is zero equals the squared overlap modulus of the two states plus
        one, divided by two.

    The two preparations are only invoked forward; the conditional swap is generated
    by the algorithm."""
    a, c = checked_state_preparation(first), checked_state_preparation(second)
    if a.width != c.width:
        raise ValidationError("The two states in a swap test must have the same width")
    b = Builder(
        _name("swap_test", a.operation, c.operation),
        {
            "left": Bits(a.width),
            "right": Bits(c.width),
            "left_work": Bits(a.work_width),
            "right_work": Bits(c.work_width),
            "probe": Bits(1),
        },
        resources_for(("a", a.operation), ("b", c.operation)),
        attributes={"algorithm": "swap_test", "readout_register": "probe"},
    )
    invoke(b, a.operation, "a", target=b["left"], work=b["left_work"])
    invoke(b, c.operation, "b", target=b["right"], work=b["right_work"])
    b.h(b["probe"])
    with b.control(b["probe"]):
        b.swap(b["left"], b["right"])
    b.h(b["probe"])
    return b.finish()


def amplitude_estimation(
    preparation: StatePreparationProtocol, marked: Iterable[int], *, precision: int = 3
) -> Operation:
    """Generate a QPE circuit estimating the good-state probability.

    Args:
        preparation: A zero-input state preparation supporting adjoint and controlled
            invocation.
        marked: Set of integer indices of the good states in the target space.
        precision: Number of phase register bits, in the range 1..63.

    Returns:
        Operation: Registers target, work and phase. After reading phase, call
        amplitude_from_phase to decode.

    A finite-precision readout may correspond to several approximate probabilities.
    Sampling and statistics are the host's responsibility."""
    from oracq.algorithms.common.search import grover_iterate

    prep = checked_state_preparation(preparation, adjoint=True, controlled=True)
    iterate = grover_iterate(prep, marked)
    qpe = phase_estimation(iterate, precision=precision)
    b = Builder(
        _name("amplitude_estimation", prep.operation, iterate, precision),
        {"target": Bits(prep.width), "work": Bits(prep.work_width), "phase": Bits(precision)},
        resources_for(("prep", prep.operation), ("qpe", qpe)),
        attributes={
            "algorithm": "amplitude_estimation",
            "readout_register": "phase",
            "decoder": "sin(pi*phase/2**precision)**2",
        },
    )
    invoke(b, prep.operation, "prep", target=b["target"], work=b["work"])
    invoke(b, qpe, "qpe", target=b["target"], work=b["work"], phase=b["phase"])
    return b.finish()


def amplitude_from_phase(value: int, precision: int) -> float:
    """Convert one phase sample into an estimate of the good-state probability.

    Args:
        value: Unsigned integer readout of phase.
        precision: Bit width of phase.

    Returns:
        float: ``sin(pi*value/2**precision)**2``, between zero and one."""
    positive_integer(precision, "amplitude.precision", maximum=63)
    positive_integer(value, "amplitude.phase", minimum=0, maximum=(1 << precision) - 1)
    return math.sin(math.pi * value / (1 << precision)) ** 2
