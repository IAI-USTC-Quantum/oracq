"""Assemble dissipative linear evolution as a finite weighted sum of Hermitian evolution branches."""

from __future__ import annotations

import json
import math
from collections.abc import Callable
from dataclasses import dataclass

from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import StateOracle
from oracq.algorithms.qode._dynamics import _lcu_dynamics
from oracq.algorithms.qode.ode_models import LinearODE
from oracq.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class QuadraturePlan:
    """Discrete quadrature plan for the LCHS integral.

    Attributes:
        nodes: Finite real quadrature nodes.
        weights: Finite complex weights of the same length as ``nodes``; must already include the integration kernel and not be all zero.
        kernel: Plan source marker, written verbatim into the result metadata.

    Raises:
        ValidationError: Node and weight lengths mismatch, non-finite values are present, or the weights are all zero.
    """

    nodes: tuple[float, ...]
    weights: tuple[complex, ...]
    kernel: str = "user_supplied"

    def __post_init__(self) -> None:
        """Normalize the nodes and weights into tuples, and validate matching lengths, finite values, and not-all-zero weights."""
        object.__setattr__(self, "nodes", tuple(self.nodes))
        object.__setattr__(self, "weights", tuple(self.weights))
        if not self.nodes or len(self.nodes) != len(self.weights):
            raise ValidationError("Discrete node and weight lengths do not match")
        for node in self.nodes:
            finite_real(node, "QuadraturePlan.node")
        if any(
            type(w) not in (int, float, complex)
            or not (math.isfinite(w.real) and math.isfinite(w.imag))
            for w in self.weights
        ):
            raise ValidationError("Integration weights must be finite complex numbers")
        if not any(self.weights):
            raise ValidationError("Integration weights cannot all be zero")

    @classmethod
    def cauchy(cls, cutoff: int = 2, spacing: float = 1.0) -> QuadraturePlan:
        """Construct a symmetric quadrature plan for the Cauchy kernel.

        Args:
            cutoff: Non-negative truncation; nodes are k*spacing with k in -cutoff..cutoff.
            spacing: Positive node spacing.

        Returns:
            QuadraturePlan: Weights are spacing/(pi*(1+k**2)), with kernel marked ``finite_cauchy``.

        Raises:
            ValidationError: The truncation is negative or the spacing is not positive.
        """
        positive_integer(cutoff, "QuadraturePlan.cutoff", minimum=0)
        finite_real(spacing, "QuadraturePlan.spacing", minimum=0, strict=True)
        if cutoff < 0 or spacing <= 0:
            raise ValidationError("Invalid Cauchy discretization parameters")
        nodes = tuple(k * spacing for k in range(-cutoff, cutoff + 1))
        return cls(nodes, tuple(spacing / (math.pi * (1 + k * k)) for k in nodes), "finite_cauchy")


def lchs_qode(
    model: LinearODE,
    time: float,
    *,
    plan: QuadraturePlan | None = None,
    hamiltonian_function: Callable[[BlockEncoding, float], BlockEncoding] = taylor_hamiltonian,
) -> StateOracle:
    """Assemble dissipative linear evolution as a finite weighted sum of Hermitian evolution branches per LCHS.

    Args:
        model: LinearODE input model; parts.hermitian is L, parts.h is H, and initial is the initial state preparation.
        time: Non-negative evolution time.
        plan: QuadraturePlan quadrature plan; the default Cauchy plan is used when omitted.
        hamiltonian_function: Replaceable protocol taking (K, time) and returning a BlockEncoding.

    Returns:
        StateOracle: Per-node K_j=H+k_j*L branches combined via LCU and applied to the initial state; the success subspace is all-zero signal.

    Raises:
        ValidationError: model or plan has the wrong type, time is invalid, or the input capability contract is not satisfied.

    Finite quadrature carries no tail-integral guarantee; the remainder is declared pending
    in the ``remainder`` metadata."""
    plan = plan or QuadraturePlan.cauchy()
    require_instance(plan, QuadraturePlan, "lchs.plan")
    return _lcu_dynamics(
        model,
        time,
        plan.nodes,
        plan.weights,
        hamiltonian_function,
        "lchs_qode",
        quadrature_kernel=plan.kernel,
        quadrature_nodes=json.dumps(plan.nodes),
        remainder="finite quadrature pending",
    )
