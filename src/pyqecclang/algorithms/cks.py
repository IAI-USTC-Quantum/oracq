"""CKS §4 的基本 Chebyshev/LCU 路线；不包含 VTAA，不宣称最优复杂度。"""

import math
from dataclasses import dataclass

from ..combinators import lcu
from ..ir import ValidationError
from ..oracles import StateOracle, annotate
from ..qlss import QLSSProtocol
from ..sparse_models import chebyshev_block, real_symmetric_sparse_encoding
from .solvers import apply_be_to_state


@dataclass(frozen=True)
class CKSConfig:
    order: int = 2
    terms: int | None = None

    def coefficients(self):
        if type(self.order) is not int or not 1 <= self.order <= 128:
            raise ValidationError("CKS 基础原型的 order 范围为 1..128")
        terms = self.order if self.terms is None else self.terms
        if not 1 <= terms <= self.order:
            raise ValidationError("CKS 截断项数无效")
        return tuple(
            4
            * (-1) ** j
            * sum(math.comb(2 * self.order, self.order + i) for i in range(j + 1, self.order + 1))
            / 2 ** (2 * self.order)
            for j in range(terms)
        )


def cks_chebyshev(system, config=None):
    config = config or CKSConfig()
    if not system.hermitian:
        raise ValidationError("CKS 稀疏输入需要 Hermitian 声明或显式 Hermitian dilation")
    a = real_symmetric_sparse_encoding(
        system.access,
        system.value_format,
        system.entry_bound,
        diagonal_nonnegative=system.diagonal_nonnegative,
    )
    coefficients = config.coefficients()
    inverse = lcu(
        [(weight, chebyshev_block(a, 2 * j + 1)) for j, weight in enumerate(coefficients)]
    )
    state = apply_be_to_state(inverse, system.rhs)
    return StateOracle(
        annotate(
            state.operation,
            "unitary",
            algorithm="cks_chebyshev_basic",
            input_model="sparse_location_inplace_and_entry_xor",
            input_alpha=a.alpha,
            encoded_inverse_bound=system.spectrum.inverse_bound(a.alpha),
            polynomial_order=config.order,
            polynomial_terms=len(coefficients),
            inverse_lcu_normalization=inverse.alpha,
            implementation_scope="CKS section 4 basic LCU; no VTAA",
            correctness="pending",
        )
    )


def make_cks_qlss(config=None):
    return QLSSProtocol(
        "cks_chebyshev_basic", "sparse", lambda system: cks_chebyshev(system, config)
    )
