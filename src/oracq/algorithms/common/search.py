"""Grover search, iteration operators and coherent amplitude amplification of the
success subspace."""

from __future__ import annotations

from collections.abc import Iterable

from oracq.algorithms.input_model.block_encoding import reflect_zero
from oracq.algorithms.input_model.contracts import positive_integer, require_instance
from oracq.algorithms.input_model.interfaces import (
    StatePreparationProtocol,
    checked_state_preparation,
)
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    XorDatabase,
    annotate,
    invoke,
    resources_for,
    uniform_state,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError, fuse


def grover(
    phase_oracle: Operation,
    width: int,
    *,
    iterations: int = 1,
    preparation: StatePreparationProtocol | None = None,
) -> StateOracle:
    """Generate a Grover search circuit with a replaceable initial state.

    Args:
        phase_oracle: Phase Operation marking the target basis states; its interface is
            target, possibly with an additional work.
        width: Target bit width of the search space.
        iterations: Nonnegative number of Grover iterations.
        preparation: Initial state preparation; when omitted the uniform state is used.
            Any provided work must be cleaned.

    Returns:
        StateOracle: target holds the search result, and signal keeps the phase query
        and the preparation workspace.

    The reflection about the initial state covers both target and the preparation
    workspace."""
    positive_integer(iterations, "grover.iterations", minimum=0)
    prep = (
        checked_state_preparation(preparation, adjoint=True)
        if preparation is not None
        else uniform_state(width)
    )
    if prep.width != width:
        raise ValidationError("The Grover initial state width does not match the search"
                              " space")
    phase_work = next((r.type.width for r in phase_oracle.module.registers if r.name == "work"), 0)
    b = Builder(
        _name("grover", phase_oracle, prep.operation, iterations),
        {"target": Bits(width), "signal": Bits(phase_work + prep.work_width)},
        resources_for(("phase", phase_oracle), ("prep", prep.operation)),
        attributes={
            "algorithm": "grover",
            "iterations": iterations,
            "validation_stage": "paradigm",
        },
    )
    pw, bw = b["signal"][:phase_work], b["signal"][phase_work:]
    invoke(b, prep.operation, "prep", target=b["target"], work=bw)
    with b.repeat(iterations):
        arguments = {"target": b["target"]}
        if any(r.name == "work" for r in phase_oracle.module.registers):
            arguments["work"] = pw
        invoke(b, phase_oracle, "phase", **arguments)
        with b.adjoint():
            invoke(b, prep.operation, "prep", target=b["target"], work=bw)
        reflect_zero(b, fuse(b["target"], bw), positive=True)
        invoke(b, prep.operation, "prep", target=b["target"], work=bw)
    return StateOracle(b.finish())


def phase_from_database(database: XorDatabase) -> Operation:
    """Convert a one-bit-output XOR database into a phase oracle.

    Args:
        database: A predicate database with a single result bit.

    Returns:
        Operation: Interface is target and a one-bit work; basis states with f(x)=1
        acquire a phase of -1, and work is cleaned.

    Raises:
        ValidationError: The database result bit width is not one.

    The input work must be zero; the database call and its inverse cancel, leaving only
    the phase flip."""
    if database.data_width != 1:
        raise ValidationError("The predicate database must have exactly one output bit")
    b = Builder(
        _name("phase_from_database", database.operation),
        {"target": Bits(database.address_width), "work": Bits(1)},
        resources_for(("db", database.operation)),
    )
    invoke(b, database.operation, "db", address=b["target"], data=b["work"])
    b.z(b["work"])
    with b.adjoint():
        invoke(b, database.operation, "db", address=b["target"], data=b["work"])
    return annotate(b.finish(), "phase_oracle")


def grover_iterate(preparation: StatePreparationProtocol, marked: Iterable[int]) -> Operation:
    """Return Q=A(2|0><0|-I)A†S_good; marked is the set of target basis state indices.

    Args:
        preparation: State preparation handle of the search space; must support adjoint
            invocation.
        marked: Set of target basis state indices, integers within 0..2^width−1.

    Returns:
        Operation: One round of the Grover amplification iterate, containing the target
        and work registers.
    """
    from oracq.algorithms.input_model.oracles import phase_marks

    prep = checked_state_preparation(preparation, adjoint=True)
    marker = phase_marks(prep.width, tuple(marked))
    b = Builder(
        _name("grover_iterate", prep.operation, marker),
        {"target": Bits(prep.width), "work": Bits(prep.work_width)},
        resources_for(("prep", prep.operation), ("marker", marker)),
    )
    invoke(b, marker, "marker", target=b["target"])
    with b.adjoint():
        invoke(b, prep.operation, "prep", target=b["target"], work=b["work"])
    reflect_zero(b, fuse(b["target"], b["work"]), positive=True)
    invoke(b, prep.operation, "prep", target=b["target"], work=b["work"])
    return b.finish()


def amplify_success(state: StateOracle, *, iterations: int = 1) -> StateOracle:
    """Coherently amplify the zero-signal success subspace of a state oracle.

    Args:
        state: A StateOracle supporting adjoint invocation.
        iterations: Nonnegative number of amplification rounds.

    Returns:
        StateOracle: Public interface identical to the input; the success condition
        remains signal==0.

    The round count must be chosen to match the input success probability; too many
    iterations can lower the success probability."""
    require_instance(state, StateOracle, "amplify_success.state")
    positive_integer(iterations, "amplify_success.iterations", minimum=0)
    b = Builder(
        _name("amplify_success", state.operation, iterations),
        {"target": Bits(state.width), "signal": Bits(state.signal_qubits)},
        resources_for(("state", state.operation)),
        attributes={"algorithm": "amplitude_amplification", "success_condition": "signal == 0"},
    )
    invoke(b, state.operation, "state", target=b["target"], signal=b["signal"])
    with b.repeat(iterations):
        reflect_zero(b, b["signal"])
        with b.adjoint():
            invoke(b, state.operation, "state", target=b["target"], signal=b["signal"])
        reflect_zero(b, fuse(b["target"], b["signal"]), positive=True)
        invoke(b, state.operation, "state", target=b["target"], signal=b["signal"])
    return StateOracle(b.finish())
