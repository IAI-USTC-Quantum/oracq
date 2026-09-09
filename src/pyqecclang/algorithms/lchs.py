"""将耗散线性演化组装为 Hermitian 演化分支的有限加权和。"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass

from pyqecclang.algorithms._dynamics import _lcu_dynamics
from pyqecclang.algorithms.contracts import finite_real, positive_integer, require_instance
from pyqecclang.algorithms.hamiltonian import taylor_hamiltonian
from pyqecclang.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class QuadraturePlan:
    nodes: tuple[float, ...]
    weights: tuple[complex, ...]
    kernel: str = "user_supplied"

    def __post_init__(self):
        object.__setattr__(self, "nodes", tuple(self.nodes))
        object.__setattr__(self, "weights", tuple(self.weights))
        if not self.nodes or len(self.nodes) != len(self.weights):
            raise ValidationError("离散节点和权重长度不符")
        for node in self.nodes:
            finite_real(node, "QuadraturePlan.node")
        if any(
            type(w) not in (int, float, complex)
            or not (math.isfinite(w.real) and math.isfinite(w.imag))
            for w in self.weights
        ):
            raise ValidationError("积分权重必须是有限复数")
        if not any(self.weights):
            raise ValidationError("积分权重不能全部为零")

    @classmethod
    def cauchy(cls, cutoff=2, spacing=1.0):
        positive_integer(cutoff, "QuadraturePlan.cutoff", minimum=0)
        finite_real(spacing, "QuadraturePlan.spacing", minimum=0, strict=True)
        if cutoff < 0 or spacing <= 0:
            raise ValidationError("Cauchy 离散参数无效")
        nodes = tuple(k * spacing for k in range(-cutoff, cutoff + 1))
        return cls(nodes, tuple(spacing / (math.pi * (1 + k * k)) for k in nodes), "finite_cauchy")


def lchs_qode(model, time, *, plan=None, hamiltonian_function=taylor_hamiltonian):
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
