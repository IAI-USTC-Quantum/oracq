"""Contour-decomposition-based assembly of matrix functions and evolution, explicitly recording the finite series and omitted terms."""

from __future__ import annotations

import cmath
import json
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import StateOracle, annotate
from oracq.algorithms.qode._dynamics import _lcu_dynamics
from oracq.algorithms.qode.ode_models import HermitianParts, LinearODE
from oracq.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class ContourPlan:
    """QST Eq.12 main series; auxiliary poles and the finite truncation remainder are explicitly retained in the report."""

    a: float = 1.0
    cutoff: int = 2
    poles: tuple[complex, ...] = (2j, -1 + 1j, 1j, 1 + 1j)

    def __post_init__(self) -> None:
        """Validate a, the cutoff, and the finiteness, distinctness, and non-reality of the auxiliary poles."""
        finite_real(self.a, "ContourPlan.a", minimum=0, strict=True)
        positive_integer(self.cutoff, "ContourPlan.cutoff", minimum=0)
        object.__setattr__(self, "poles", tuple(self.poles))
        if any(
            type(p) not in (int, float, complex)
            or not (math.isfinite(p.real) and math.isfinite(p.imag))
            for p in self.poles
        ):
            raise ValidationError("CBMD poles must be finite")
        if not math.isfinite(self.a) or self.a <= 0 or self.cutoff < 0:
            raise ValidationError("CBMD a must be positive and the cutoff non-negative")
        if len(set(self.poles)) != len(self.poles) or any(
            p == -1j or p.imag == 0 for p in self.poles
        ):
            raise ValidationError("CBMD currently requires distinct non-real simple auxiliary poles avoiding -i")

    @property
    def nodes(self) -> tuple[float, ...]:
        """Real nodes of the main series; entry k is k/a with k in -cutoff..cutoff, 2*cutoff+1 in total."""
        return tuple(k / self.a for k in range(-self.cutoff, self.cutoff + 1))

    @property
    def weights(self) -> tuple[complex, ...]:
        """Complex weights of the main series nodes, in one-to-one correspondence with ``nodes``.

        Computed by the residue closed form of QST Eq.12; the denominator contains the (q+i)
        factor and all auxiliary poles."""
        numerator = math.expm1(-2 * math.pi * self.a)
        return tuple(
            numerator
            / (
                self.a
                * 2
                * math.pi
                * 1j
                * (q + 1j)
                * math.prod((q - p) / (-1j - p) for p in self.poles)
            )
            for q in self.nodes
        )

    @property
    def auxiliary_coefficients(self) -> tuple[complex, ...]:
        """Complex coefficients of the auxiliary poles, in one-to-one correspondence with ``poles``.

        The corresponding non-Hermitian evolution branches are not generated currently and
        are only declared in the omitted list of ``metadata``."""
        numerator = math.expm1(-2 * math.pi * self.a)
        return tuple(
            numerator
            / (
                (cmath.exp(-2 * math.pi * p * self.a * 1j) - 1)
                * math.prod((p - other) / (-1j - other) for other in self.poles if other != p)
            )
            for p in self.poles
        )

    def metadata(self) -> str:
        """Export the plan parameters and the omitted-terms declaration as JSON text.

        Returns:
            str: Contains a, cutoff, nodes, poles, weights, auxiliary coefficients, the omitted terms, the assumption preconditions, and the source reference."""
        def pair(z: complex) -> list[float]:
            """Expand a complex number into a two-element list of real and imaginary parts."""
            return [complex(z).real, complex(z).imag]

        return json.dumps(
            {
                "a": self.a,
                "cutoff": self.cutoff,
                "nodes": self.nodes,
                "poles": [pair(x) for x in self.poles],
                "weights": [pair(x) for x in self.weights],
                "auxiliary_coefficients": [pair(x) for x in self.auxiliary_coefficients],
                "omitted": ["auxiliary_nonhermitian_evolutions", "infinite_series_tail"],
                "assumption": "L>=0 and norm(integral L dt)<=2*pi*a",
                "source": "QST 11 035027 (2026), Eq.12-13; arXiv:2511.10267v3",
            },
            separators=(",", ":"),
        )


def cbmd_qode(
    model: LinearODE,
    time: float,
    *,
    plan: ContourPlan | None = None,
    hamiltonian_function: Callable[[BlockEncoding, float], BlockEncoding] = taylor_hamiltonian,
) -> StateOracle:
    """Assemble the quantum simulation program for u'=-Au via contour decomposition.

    Args:
        model: LinearODE input model; parts.hermitian is L, parts.h is H, and initial is the initial state preparation.
        time: Non-negative evolution time.
        plan: ContourPlan contour plan; the default plan is used when omitted.
        hamiltonian_function: Replaceable protocol taking (K, time) and returning a BlockEncoding.

    Returns:
        StateOracle: Per-node K_k=H+q_k*L branches combined via LCU and applied to the initial state; the success subspace is all-zero signal.

    Raises:
        ValidationError: model or plan has the wrong type, time is invalid, or the input capability contract is not satisfied.

    Auxiliary-pole branches and the infinite series tail are not generated and are only
    declared explicitly in the contour_plan metadata."""
    plan = plan or ContourPlan()
    require_instance(plan, ContourPlan, "cbmd.plan")
    return _lcu_dynamics(
        model,
        time,
        plan.nodes,
        plan.weights,
        hamiltonian_function,
        "cbmd_qode",
        contour_plan=plan.metadata(),
        remainder="auxiliary pole contribution and truncation pending",
    )


def cbmd_function(
    a: BlockEncoding,
    nodes: Sequence[float],
    residue_weights: Sequence[complex],
    hermitian_function: Callable[[BlockEncoding], BlockEncoding],
) -> BlockEncoding:
    """Generic f(A) assembly point: the Hermitian function protocol stays open and is not silently swapped for matrix inversion.

    Args:
        a: Block encoding of the target operator A.
        nodes: Residue pole locations, entering the H + q*L combination point by point.
        residue_weights: Complex residue weights in one-to-one correspondence with the poles; the sign convention is the caller's responsibility.
        hermitian_function: Hermitian function protocol implementation; the input is the
            block encoding after the pole combination, and it must return a BlockEncoding.

    Returns:
        BlockEncoding: The residue-weighted f(A) block encoding, with correctness marked pending.
    """
    parts = HermitianParts.from_operator(a)
    terms = [
        (weight, hermitian_function(lcu([(node, parts.h), (1, parts.hermitian)])))
        for node, weight in zip(nodes, residue_weights, strict=True)
    ]
    result = lcu(terms)
    return BlockEncoding(
        annotate(
            result.operation,
            "block_encoding",
            be_alpha=result.alpha,
            algorithm="cbmd_matrix_function",
            correctness="pending",
            residue_sign_convention="caller supplies target-side weights from contour identity",
        )
    )
