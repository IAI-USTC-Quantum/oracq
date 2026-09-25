"""Shared LCU branch assembly utilities for non-unitary evolution algorithms."""

from __future__ import annotations

from collections.abc import Callable, Sequence

from oracq.algorithms.common.state_preparation import apply_be_to_state
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.contracts import finite_real, require_instance
from oracq.algorithms.input_model.interfaces import (
    operator_state_contract,
)
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    annotate,
)
from oracq.algorithms.qode.ode_models import LinearODE
from oracq.infrastructure.builder import Operation
from oracq.infrastructure.ir import ValidationError


def tagged(operation: Operation, algorithm: str, **metadata: str | int | float) -> Operation:
    """Annotate an evolution operation with the ``unitary`` paradigm, registering the algorithm name and metadata whose correctness is pending.

    Args:
        operation: The evolution operation to annotate.
        algorithm: Algorithm name registered into the attributes.
        **metadata: Additional algorithm metadata to register; string or numeric values.

    Returns:
        Operation: The operation annotated with the ``unitary`` paradigm and correctness=pending.
    """
    return annotate(
        operation,
        "unitary",
        algorithm=algorithm,
        correctness="pending",
        validation_stage="paradigm",
        **metadata,
    )


def _lcu_dynamics(
    model: LinearODE,
    time: float,
    nodes: Sequence[float],
    weights: Sequence[complex],
    hamiltonian_function: Callable[[BlockEncoding, float], BlockEncoding],
    algorithm: str,
    **metadata: str,
) -> StateOracle:
    """Build the Hamiltonian branches K=H+kL node by node, combine them via an LCU, and apply the result to the initial state."""
    require_instance(model, LinearODE, algorithm + ".model")
    finite_real(time, algorithm + ".time", minimum=0)
    operator_state_contract(algorithm).check(
        generator=model.parts.hermitian, initial=model.initial
    ).require()
    operator_state_contract(algorithm).check(
        generator=model.parts.h, initial=model.initial
    ).require()
    if not callable(hamiltonian_function):
        raise ValidationError("hamiltonian_function must be callable")
    terms: list[tuple[complex, BlockEncoding]] = []
    for node, weight in zip(nodes, weights, strict=True):
        hk = lcu([(1, model.parts.h), (node, model.parts.hermitian)])
        encoded = hamiltonian_function(hk, time)
        if not isinstance(encoded, BlockEncoding):
            raise ValidationError(
                "The Hamiltonian-function protocol must return a BlockEncoding to preserve alpha"
            )
        terms.append((weight, encoded))
    evolution = lcu(terms)
    state = apply_be_to_state(evolution, model.initial)
    return StateOracle(
        tagged(
            state.operation,
            algorithm,
            evolution_alpha=evolution.alpha,
            branch_count=len(nodes),
            input_assumption="L>=0; autonomous; homogeneous",
            **metadata,
        )
    )
