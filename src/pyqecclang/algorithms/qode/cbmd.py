"""基于轮廓分解的矩阵函数与演化组装，显式记录有限级数与遗漏项。"""

from __future__ import annotations

import cmath
import json
import math
from dataclasses import dataclass

from pyqecclang.algorithms.common.hamiltonian import taylor_hamiltonian
from pyqecclang.algorithms.input_model.block_encoding import lcu
from pyqecclang.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from pyqecclang.algorithms.input_model.operators import BlockEncoding
from pyqecclang.algorithms.input_model.oracles import (
    annotate,
)
from pyqecclang.algorithms.qode._dynamics import _lcu_dynamics
from pyqecclang.algorithms.qode.ode_models import HermitianParts
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
        """主级数实节点，第 k 项为 k/a，k 取 -cutoff..cutoff，共 2*cutoff+1 个。"""
        return tuple(k / self.a for k in range(-self.cutoff, self.cutoff + 1))

    @property
    def weights(self):
        """主级数各节点的复权重，与 ``nodes`` 一一对应。

        按 QST Eq.12 的留数闭式计算，分母包含 (q+i) 因子与全部辅助极点。"""
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
        """各辅助极点的复系数，与 ``poles`` 一一对应。

        对应的非 Hermitian 演化分支当前不生成，仅在 ``metadata`` 的 omitted 列表中声明。"""
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
        """导出计划参数与遗漏项声明的 JSON 文本。

        Returns:
            str: 含 a、cutoff、节点、极点、权重、辅助系数、omitted 遗漏项、assumption 前提与 source 文献来源。"""
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
    """基于轮廓分解组装 u'=-Au 的量子模拟程序。

    Args:
        model: LinearODE 输入模型，parts.hermitian 为 L、parts.h 为 H，initial 为初态制备。
        time: 非负演化时间。
        plan: ContourPlan 轮廓计划；省略时使用默认计划。
        hamiltonian_function: 接受 (K, time) 并返回 BlockEncoding 的可替换协议。

    Returns:
        StateOracle: 逐节点 K_k=H+q_k*L 分支经 LCU 组合后作用到初态，成功子空间为 signal 全零。

    Raises:
        ValidationError: model 或 plan 类型不符、time 非法或输入能力契约不满足。

    辅助极点分支与无穷级数尾不生成，只在 contour_plan 元数据中显式声明。"""
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
