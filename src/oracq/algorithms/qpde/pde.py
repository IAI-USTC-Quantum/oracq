"""Assembly interface between spatial discretization results and replaceable QODE generators."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
)


@dataclass(frozen=True)
class DiscretePDE:
    """Linear PDE input object after spatial discretization.

    It stores only the discrete generator, the initial state preparation, and a label;
    application metadata such as the grid, boundaries, and dimensions is not carried by
    this type.

    Attributes:
        generator: Block encoding of the discrete generator G, corresponding to ``u' = Gu``.
        initial: Initial state preparation.
        label: Problem label.
    """

    generator: BlockEncoding
    initial: StatePreparation
    label: str = "linear_pde"


def make_qpde(
    qode: Callable[[BlockEncoding, StatePreparation, float], StateOracle],
    discretizer: Callable[[object], DiscretePDE] = lambda problem: cast("DiscretePDE", problem),
) -> Callable[[object, float], StateOracle]:
    """Wrap a three-argument linear QODE protocol into a QPDE generation function.

    Args:
        qode: A linear solve protocol of the form ``(generator, initial, time) -> StateOracle``.
        discretizer: A callable mapping the problem object to an object with ``generator`` and ``initial`` attributes (e.g. ``DiscretePDE``); returns its input unchanged by default.

    Returns:
        callable: A generation function of the form ``(problem, final_time) -> StateOracle``.
    """
    def generate(problem: object, final_time: float) -> StateOracle:
        """Discretize the problem object and hand it to the linear solve protocol for evolution."""
        discrete = discretizer(problem)
        result = qode(discrete.generator, discrete.initial, final_time)
        return result

    return generate


@dataclass(frozen=True)
class PDEInput:
    """Input model after spatial discretization; may directly be a PolynomialODE; no dense matrix is required."""

    model: object
    label: str = "open_spatial_discretization"


def qpde_solver(
    qode: Callable[[object, float], StateOracle],
    spatial_discretizer: Callable[[object], object] = lambda problem: cast("PDEInput", problem).model,
) -> Callable[[object, float], StateOracle]:
    """Wrap a ``(model, time)``-shaped solve protocol into a QPDE generation function.

    Unlike ``make_qpde``, the discretization result is passed to ``qode`` directly as a
    single model, which suits model-level protocols such as ``PolynomialODE``.

    Args:
        qode: A solve protocol of the form ``(model, time) -> StateOracle``.
        spatial_discretizer: A callable extracting the model from the problem object; defaults to ``problem.model``.

    Returns:
        callable: A generation function of the form ``(problem, time) -> StateOracle``.
    """
    def generate(problem: object, time: float) -> StateOracle:
        """Extract the problem model and hand it to the solve protocol for evolution."""
        return qode(spatial_discretizer(problem), time)

    return generate
