"""Assembly of qubitization, explicit QSVT phase sequences and oblivious amplification."""

from __future__ import annotations

from collections.abc import Iterable

from oracq.algorithms.input_model.block_encoding import reflect_zero
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import invoke, resources_for
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits


def qubitization_walk(a: BlockEncoding) -> Operation:
    """Assemble the qubitization walk operation of a block encoding ``a``.

    It first invokes ``a`` forward, then applies the positive reflection of
    the signal about the zero subspace, yielding a walk operator of the form
    ``(2 P_0 - I) U``; the registers follow the target and signal signature
    of ``a``.

    Args:
        a: The input ``BlockEncoding``.

    Returns:
        Operation: Walk operation that keeps the module call to ``a``.
    """
    b = Builder(
        _name("qubitization_walk", a.operation),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        resources_for(("a", a.operation)),
        attributes={"algorithm": "qubitization_walk"},
    )
    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
    reflect_zero(b, b["signal"], positive=True)
    return b.finish()


def qsvt_sequence(a: BlockEncoding, phases: Iterable[float]) -> Operation:
    """Assemble a QSVT circuit for a block encoding ``a`` from an explicit phase
    sequence.

    The phases are ordered in time, ``phases[0]`` acting first; ``len(phases)``
    phases accompany ``len(phases) - 1`` calls to ``a``, alternating forward
    and inverse; each phase is composed from a global phase and a controlled
    phase on the zero subspace of signal. The phase convention matches the
    phase synthesis of the ``qsvt`` module.

    Args:
        a: The input ``BlockEncoding``.
        phases: Phase sequence in radians; elements may be any value convertible to
            ``float``.

    Returns:
        Operation: QSVT operation whose calls to ``a`` are kept as module calls.
    """
    phases = tuple(float(p) for p in phases)
    b = Builder(
        _name("qsvt", a.operation, phases),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        resources_for(("a", a.operation)),
        attributes={"algorithm": "qsvt_sequence", "validation_stage": "paradigm"},
    )
    for i, phase in enumerate(phases):
        b.global_phase(-phase)
        if b["signal"].width:
            with b.control(b["signal"], 0):
                b.global_phase(2 * phase)
        else:
            b.global_phase(2 * phase)
        if i + 1 < len(phases):
            if i % 2:
                with b.adjoint():
                    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
            else:
                invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
    return b.finish()


def oblivious_amplification(a: BlockEncoding, iterations: int = 1) -> Operation:
    """Assemble oblivious amplitude amplification for a block encoding ``a``.

    It first invokes ``a`` forward; each round then performs the zero
    reflection of signal, the inverse call to ``a``, another zero reflection
    and the forward call, amplifying the zero-signal projected component.

    Args:
        a: The input ``BlockEncoding``.
        iterations: Number of amplification rounds, one by default.

    Returns:
        Operation: Operation whose iterations are kept in an RIR ``Repeat``.
    """
    b = Builder(
        _name("oaa", a.operation, iterations),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        resources_for(("a", a.operation)),
        attributes={"algorithm": "oaa"},
    )
    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
    with b.repeat(iterations):
        reflect_zero(b, b["signal"])
        with b.adjoint():
            invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
        reflect_zero(b, b["signal"])
        invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
    return b.finish()
