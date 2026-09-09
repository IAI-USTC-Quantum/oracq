"""早期演化工厂的兼容实现；新代码使用各方法的专门入口。"""

from __future__ import annotations

from pyqecclang.algorithms.block_encoding import lcu
from pyqecclang.algorithms.operators import BlockEncoding
from pyqecclang.algorithms.oracles import (
    annotate,
)
from pyqecclang.algorithms.state_preparation import (
    apply_be_to_state,
    extend_initial,
    select_subspace,
)
from pyqecclang.infrastructure.ir import ValidationError


def make_lchs_qode(ham_sim, times, weights):
    times, weights = tuple(times), tuple(weights)
    if len(times) != len(weights):
        raise ValidationError("LCHS 离散时间与权重长度不同")

    def generate(generator, initial, final_time):
        terms = []
        for weight, time in zip(weights, times, strict=True):
            op = ham_sim(generator, time * final_time)
            terms.append((weight, BlockEncoding(annotate(op, "block_encoding", be_alpha=1.0))))
        return apply_be_to_state(lcu(terms), initial)

    return generate


def make_schrodingerisation_qode(embedding, ham_sim, *, extra_width=1):
    """显式保留 Hamiltonian lift 的实现边界，随后组装演化和物理通道。"""

    def generate(generator, initial, final_time):
        hamiltonian = embedding(generator)
        if hamiltonian.width != generator.width + extra_width:
            raise ValidationError("Schrodingerisation lift 的寄存器宽度不符")
        evolution = ham_sim(hamiltonian, final_time)
        encoded = BlockEncoding(annotate(evolution, "block_encoding", be_alpha=1.0))
        lifted_state = apply_be_to_state(encoded, extend_initial(initial, extra_width))
        return select_subspace(
            lifted_state, generator.width, 0, label="schrodingerisation_physical_channel"
        )

    return generate
