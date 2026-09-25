"""Compatibility implementations of the early evolution factories; new code should use the dedicated entry points of each method."""

from __future__ import annotations

from collections.abc import Callable, Iterable

from oracq.algorithms.common.state_preparation import (
    apply_be_to_state,
    extend_initial,
    select_subspace,
)
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    annotate,
)
from oracq.infrastructure.builder import Operation
from oracq.infrastructure.ir import ValidationError


def make_lchs_qode(
    ham_sim: Callable[[BlockEncoding, float], Operation],
    times: Iterable[float],
    weights: Iterable[complex],
) -> Callable[[BlockEncoding, StatePreparation, float], StateOracle]:
    """Wrap a Hamiltonian simulation protocol into an LCHS weighted-sum QODE generation function.

    Args:
        ham_sim: Hamiltonian simulation protocol taking a block encoding and a duration and returning the implementing operation.
        times: Discrete quadrature time sequence.
        weights: Complex weight sequence of the same length as ``times``.

    Returns:
        callable: A generation function of the form ``(generator, initial, final_time) -> StateOracle``
        that combines the per-time evolutions via an LCU with the given weights and applies
        the result to the initial state.

    Raises:
        ValidationError: ``times`` and ``weights`` have different lengths.
    """
    times, weights = tuple(times), tuple(weights)
    if len(times) != len(weights):
        raise ValidationError("LCHS discrete times and weights have different lengths")

    def generate(
        generator: BlockEncoding, initial: StatePreparation, final_time: float
    ) -> StateOracle:
        """LCU-weight the per-time evolutions by the weights and apply the sum to the initial state."""
        terms: list[tuple[complex, BlockEncoding]] = []
        for weight, time in zip(weights, times, strict=True):
            op = ham_sim(generator, time * final_time)
            terms.append((weight, BlockEncoding(annotate(op, "block_encoding", be_alpha=1.0))))
        return apply_be_to_state(lcu(terms), initial)

    return generate


def make_schrodingerisation_qode(
    embedding: Callable[[BlockEncoding], BlockEncoding],
    ham_sim: Callable[[BlockEncoding, float], Operation],
    *,
    extra_width: int = 1,
) -> Callable[[BlockEncoding, StatePreparation, float], StateOracle]:
    """Keep the implementation boundary of the Hamiltonian lift explicit, then assemble the evolution and the physical channel.

    Args:
        embedding: Hamiltonian lift; maps the generator block encoding to a block encoding
            enlarged by ``extra_width`` bits.
        ham_sim: Hamiltonian simulation implementation of the form (BE, time) returning an Operation.
        extra_width: Auxiliary bits added by the lift; a positive integer, defaulting to 1.

    Returns:
        Callable[[BlockEncoding, StatePreparation, float], StateOracle]: A generation function
        taking the generator, the initial state, and the final time, and outputting the
        physical-channel readout state oracle.
    """

    def generate(
        generator: BlockEncoding, initial: StatePreparation, final_time: float
    ) -> StateOracle:
        """Assemble the post-lift evolution and read out the physical channel subspace."""
        hamiltonian = embedding(generator)
        if hamiltonian.width != generator.width + extra_width:
            raise ValidationError("Schrodingerisation lift register width mismatch")
        evolution = ham_sim(hamiltonian, final_time)
        encoded = BlockEncoding(annotate(evolution, "block_encoding", be_alpha=1.0))
        lifted_state = apply_be_to_state(encoded, extend_initial(initial, extra_width))
        return select_subspace(
            lifted_state, generator.width, 0, label="schrodingerisation_physical_channel"
        )

    return generate
