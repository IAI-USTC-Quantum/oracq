"""Composition utilities for initial state extension, applying encodings and physical
subspace selection."""

from __future__ import annotations

from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    annotate,
    invoke,
    resources_for,
)
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, ValidationError, fuse


def extend_initial(prep: StatePreparation, extra_width: int) -> StatePreparation:
    """Extend the target space of a state preparation ``prep`` by ``extra_width`` high
    bits.

    The original preparation acts on the low bits of the extended target and
    the new high bits stay zero; the zero-input promise is preserved.

    Args:
        prep: The original ``StatePreparation``.
        extra_width: Number of high qubits to append.

    Returns:
        StatePreparation: Preparation whose target width is ``prep.width + extra_width``.
    """
    b = Builder(
        _name("extend_initial", prep.operation, extra_width),
        {"target": Bits(prep.width + extra_width), "work": Bits(prep.work_width)},
        resources_for(("initial", prep.operation)),
    )
    invoke(b, prep.operation, "initial", target=b["target"][: prep.width], work=b["work"])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def select_subspace(
    state: StateOracle,
    output_width: int,
    high_value: int = 0,
    *,
    label: str = "selection",
) -> StateOracle:
    """Select the physical subspace out of the target space of a state oracle.

    The low ``output_width`` bits serve as the output, while the remaining
    high bits are folded into signal and required to equal ``high_value``;
    the success condition is that the whole signal register returns to zero.

    Args:
        state: The input ``StateOracle``, with width not less than ``output_width``.
        output_width: The selected physical output width.
        high_value: The value that the discarded high bits must equal, zero by default.
        label: Label used for the generated module name and algorithm attribute.

    Returns:
        StateOracle: State oracle with output width ``output_width`` and a signal
        holding the original signal bits plus the appended high bits.

    Raises:
        ValidationError: ``output_width`` exceeds ``state.width``, or ``high_value``
            is out of range.
    """
    extra = state.width - output_width
    if extra < 0 or not 0 <= high_value < 1 << extra:
        raise ValidationError("Invalid output subspace layout")
    b = Builder(
        _name(label, state.operation, output_width, high_value),
        {"target": Bits(output_width), "signal": Bits(state.signal_qubits + extra)},
        resources_for(("state", state.operation)),
        attributes={
            "algorithm": label,
            "selected_high_value": high_value,
            "validation_stage": "paradigm",
            "success_condition": "signal == 0",
        },
    )
    old_signal, high = b["signal"][: state.signal_qubits], b["signal"][state.signal_qubits :]
    invoke(b, state.operation, "state", target=fuse(b["target"], high), signal=old_signal)
    for bit in range(extra):
        if (high_value >> bit) & 1:
            b.x(high[bit])
    return StateOracle(b.finish())


def apply_be_to_state(a: BlockEncoding, prep: StatePreparation) -> StateOracle:
    """Apply a block encoding ``a`` to an already prepared initial state, yielding a
    state oracle.

    The signal bits of ``a`` occupy the low bits of signal and the workspace
    of ``prep`` is folded into its high bits; the success condition is that
    the whole signal register returns to zero.

    Args:
        a: The ``BlockEncoding`` of the operator to apply.
        prep: A ``StatePreparation`` with the same target width as ``a``.

    Returns:
        StateOracle: State oracle with unchanged target width and an attached signal
        register.

    Raises:
        ValidationError: The two target widths differ.
    """
    if a.width != prep.width:
        raise ValidationError("The operator and preparation target widths do not match")
    b = Builder(
        _name("apply_be_state", a.operation, prep.operation),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits + prep.work_width)},
        resources_for(("a", a.operation), ("prep", prep.operation)),
    )
    invoke(b, prep.operation, "prep", target=b["target"], work=b["signal"][a.signal_qubits :])
    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"][: a.signal_qubits])
    return StateOracle(b.finish())
