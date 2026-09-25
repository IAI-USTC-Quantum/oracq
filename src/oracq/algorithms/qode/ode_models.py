"""Shared input model of Hermitian decompositions and autonomous linear ODEs."""

from __future__ import annotations

from dataclasses import dataclass

from oracq.algorithms.input_model.block_encoding import adjoint_be, lcu
from oracq.algorithms.input_model.contracts import require_instance
from oracq.algorithms.input_model.interfaces import (
    as_block_encoding,
    as_state_preparation,
)
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
)
from oracq.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class HermitianParts:
    """A=L+iH; the Hermitian/positive-semidefinite properties are input-model declarations, not proved by the language."""

    hermitian: BlockEncoding
    h: BlockEncoding

    def __post_init__(self) -> None:
        """Normalize the two components into block encodings and validate matching widths."""
        object.__setattr__(self, "hermitian", as_block_encoding(self.hermitian))
        object.__setattr__(self, "h", as_block_encoding(self.h))
        require_instance(self.h, BlockEncoding, "HermitianParts.H")
        if self.hermitian.width != self.h.width:
            raise ValidationError("Hermitian parts widths do not match")

    @classmethod
    def from_operator(cls, a: BlockEncoding) -> HermitianParts:
        """Construct the Hermitian decomposition ``A = L + iH`` from the block encoding of operator A.

        Combines ``a`` and its adjoint via an LCU to build ``L = (A + A†)/2`` and
        ``H = (A - A†)/(2i)``.

        Args:
            a: Block encoding of the operator A.

        Returns:
            HermitianParts: The corresponding Hermitian decomposition; the Hermitian nature
            of each component remains caller-declared.
        """
        adj = adjoint_be(a)
        return cls(lcu([(0.5, a), (0.5, adj)]), lcu([(-0.5j, a), (0.5j, adj)]))


@dataclass(frozen=True)
class LinearODE:
    """Shared input model of the autonomous linear ODE ``u' = -Au`` (with ``A = L + iH``).

    Attributes:
        parts: The ``HermitianParts`` decomposition of A.
        initial: Initial state preparation; its target width must match ``parts``.
        label: Model label; the default spells out the equation form.
    """

    parts: HermitianParts
    initial: StatePreparation
    label: str = "du_dt_equals_minus_A_u"

    def __post_init__(self) -> None:
        """Validate the parts type, normalize the initial state, and check matching widths."""
        require_instance(self.parts, HermitianParts, "LinearODE.parts")
        object.__setattr__(self, "initial", as_state_preparation(self.initial))
        if self.parts.hermitian.width != self.initial.width:
            raise ValidationError("Linear ODE initial state and operator widths do not match")
