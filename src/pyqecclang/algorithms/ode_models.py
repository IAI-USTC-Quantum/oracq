"""Hermitian 分解与自治线性 ODE 的共享输入模型。"""

from __future__ import annotations

from dataclasses import dataclass

from pyqecclang.algorithms.block_encoding import adjoint_be, lcu
from pyqecclang.algorithms.contracts import require_instance
from pyqecclang.algorithms.interfaces import (
    as_block_encoding,
    as_state_preparation,
)
from pyqecclang.algorithms.operators import BlockEncoding
from pyqecclang.algorithms.oracles import (
    StatePreparation,
)
from pyqecclang.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class HermitianParts:
    """A=L+iH；Hermitian/半正定性质是输入模型声明，不由语言证明。"""

    hermitian: BlockEncoding
    h: BlockEncoding

    def __post_init__(self):
        object.__setattr__(self, "hermitian", as_block_encoding(self.hermitian))
        object.__setattr__(self, "h", as_block_encoding(self.h))
        require_instance(self.h, BlockEncoding, "HermitianParts.H")
        if self.hermitian.width != self.h.width:
            raise ValidationError("Hermitian parts 宽度不匹配")

    @classmethod
    def from_operator(cls, a):
        """从算子 A 的块编码构造 Hermitian 分解 ``A = L + iH``。

        以 ``a`` 与其伴随的 LCU 组合出 ``L = (A + A†)/2`` 与 ``H = (A - A†)/(2i)``。

        Args:
            a: 算子 A 的块编码。

        Returns:
            HermitianParts: 对应的 Hermitian 分解；各分量的 Hermitian 性质仍由调用方声明。
        """
        adj = adjoint_be(a)
        return cls(lcu([(0.5, a), (0.5, adj)]), lcu([(-0.5j, a), (0.5j, adj)]))


@dataclass(frozen=True)
class LinearODE:
    """自治线性 ODE ``u' = -Au``（其中 ``A = L + iH``）的共享输入模型。

    Attributes:
        parts: A 的 ``HermitianParts`` 分解。
        initial: 初态制备，目标宽度须与 ``parts`` 一致。
        label: 模型标签，默认标明方程形式。
    """

    parts: HermitianParts
    initial: StatePreparation
    label: str = "du_dt_equals_minus_A_u"

    def __post_init__(self):
        require_instance(self.parts, HermitianParts, "LinearODE.parts")
        object.__setattr__(self, "initial", as_state_preparation(self.initial))
        if self.parts.hermitian.width != self.initial.width:
            raise ValidationError("线性 ODE 初态和算子宽度不匹配")
