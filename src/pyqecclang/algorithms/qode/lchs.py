"""将耗散线性演化组装为 Hermitian 演化分支的有限加权和。"""

from __future__ import annotations

import json
import math
from collections.abc import Callable
from dataclasses import dataclass

from pyqecclang.algorithms.common.hamiltonian import taylor_hamiltonian
from pyqecclang.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from pyqecclang.algorithms.input_model.operators import BlockEncoding
from pyqecclang.algorithms.input_model.oracles import StateOracle
from pyqecclang.algorithms.qode._dynamics import _lcu_dynamics
from pyqecclang.algorithms.qode.ode_models import LinearODE
from pyqecclang.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class QuadraturePlan:
    """LCHS 积分的离散求积计划。

    Attributes:
        nodes: 有限实数求积节点。
        weights: 与 ``nodes`` 等长的有限复权重，需已包含积分核且不全为零。
        kernel: 计划来源标记，原样写入结果元数据。

    Raises:
        ValidationError: 节点与权重长度不符、含非有限数值或权重全为零。
    """

    nodes: tuple[float, ...]
    weights: tuple[complex, ...]
    kernel: str = "user_supplied"

    def __post_init__(self) -> None:
        """把节点与权重规范化为元组，并校验长度一致、数值有限且权重非全零。"""
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
    def cauchy(cls, cutoff: int = 2, spacing: float = 1.0) -> QuadraturePlan:
        """构造 Cauchy 核的对称求积计划。

        Args:
            cutoff: 非负截断，节点取 k*spacing、k 在 -cutoff..cutoff 内。
            spacing: 正的节点间距。

        Returns:
            QuadraturePlan: 权重为 spacing/(pi*(1+k**2))，kernel 标记为 ``finite_cauchy``。

        Raises:
            ValidationError: 截断为负或间距非正。
        """
        positive_integer(cutoff, "QuadraturePlan.cutoff", minimum=0)
        finite_real(spacing, "QuadraturePlan.spacing", minimum=0, strict=True)
        if cutoff < 0 or spacing <= 0:
            raise ValidationError("Cauchy 离散参数无效")
        nodes = tuple(k * spacing for k in range(-cutoff, cutoff + 1))
        return cls(nodes, tuple(spacing / (math.pi * (1 + k * k)) for k in nodes), "finite_cauchy")


def lchs_qode(
    model: LinearODE,
    time: float,
    *,
    plan: QuadraturePlan | None = None,
    hamiltonian_function: Callable[[BlockEncoding, float], BlockEncoding] = taylor_hamiltonian,
) -> StateOracle:
    """按 LCHS 把耗散线性演化组装为 Hermitian 演化分支的有限加权和。

    Args:
        model: LinearODE 输入模型，parts.hermitian 为 L、parts.h 为 H，initial 为初态制备。
        time: 非负演化时间。
        plan: QuadraturePlan 求积计划；省略时使用 Cauchy 默认计划。
        hamiltonian_function: 接受 (K, time) 并返回 BlockEncoding 的可替换协议。

    Returns:
        StateOracle: 逐节点 K_j=H+k_j*L 分支经 LCU 组合后作用到初态，成功子空间为 signal 全零。

    Raises:
        ValidationError: model 或 plan 类型不符、time 非法或输入能力契约不满足。

    有限求积没有尾积分保证，余项以 ``remainder`` 元数据声明为 pending。"""
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
