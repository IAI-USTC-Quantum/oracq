"""基于轮廓分解的矩阵函数与演化组装，显式记录有限级数与遗漏项。"""

from __future__ import annotations

import cmath
import json
import math
from dataclasses import dataclass

from pyqecclang.algorithms._dynamics import _lcu_dynamics
from pyqecclang.algorithms.block_encoding import lcu
from pyqecclang.algorithms.contracts import finite_real, positive_integer, require_instance
from pyqecclang.algorithms.hamiltonian import taylor_hamiltonian
from pyqecclang.algorithms.ode_models import HermitianParts
from pyqecclang.algorithms.operators import BlockEncoding
from pyqecclang.algorithms.oracles import (
    annotate,
)
from pyqecclang.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class ContourPlan:
    """QST Eq.12 主级数；辅助极点与有限截断余项在报告中明确保留。"""

    a: float = 1.0
    cutoff: int = 2
    poles: tuple[complex, ...] = (2j, -1 + 1j, 1j, 1 + 1j)

    def __post_init__(self):
        finite_real(self.a, "ContourPlan.a", minimum=0, strict=True)
        positive_integer(self.cutoff, "ContourPlan.cutoff", minimum=0)
        object.__setattr__(self, "poles", tuple(self.poles))
        if any(
            type(p) not in (int, float, complex)
            or not (math.isfinite(p.real) and math.isfinite(p.imag))
            for p in self.poles
        ):
            raise ValidationError("CBMD 极点必须有限")
        if not math.isfinite(self.a) or self.a <= 0 or self.cutoff < 0:
            raise ValidationError("CBMD a 必须为正且截断非负")
        if len(set(self.poles)) != len(self.poles) or any(
            p == -1j or p.imag == 0 for p in self.poles
        ):
            raise ValidationError("CBMD 当前要求非实互异简单辅助极点，且避开 -i")

    @property
    def nodes(self):
        return tuple(k / self.a for k in range(-self.cutoff, self.cutoff + 1))

    @property
    def weights(self):
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
    def auxiliary_coefficients(self):
        numerator = math.expm1(-2 * math.pi * self.a)
        return tuple(
            numerator
            / (
                (cmath.exp(-2 * math.pi * p * self.a * 1j) - 1)
                * math.prod((p - other) / (-1j - other) for other in self.poles if other != p)
            )
            for p in self.poles
        )

    def metadata(self):
        def pair(z):
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


def cbmd_qode(model, time, *, plan=None, hamiltonian_function=taylor_hamiltonian):
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


def cbmd_function(a, nodes, residue_weights, hermitian_function):
    """通用 f(A) 组装点：Hermitian function protocol 保持开放，不偷换为矩阵求逆。"""
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
