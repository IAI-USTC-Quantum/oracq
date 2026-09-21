"""早期演化工厂的兼容实现；新代码使用各方法的专门入口。"""

from __future__ import annotations

from collections.abc import Callable, Iterable

from pyqecclang.algorithms.common.state_preparation import (
    apply_be_to_state,
    extend_initial,
    select_subspace,
)
from pyqecclang.algorithms.input_model.block_encoding import lcu
from pyqecclang.algorithms.input_model.operators import BlockEncoding
from pyqecclang.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    annotate,
)
from pyqecclang.infrastructure.builder import Operation
from pyqecclang.infrastructure.ir import ValidationError


def make_lchs_qode(
    ham_sim: Callable[[BlockEncoding, float], Operation],
    times: Iterable[float],
    weights: Iterable[complex],
) -> Callable[[BlockEncoding, StatePreparation, float], StateOracle]:
    """把 Hamiltonian 模拟协议包装成 LCHS 加权和的 QODE 生成函数。

    Args:
        ham_sim: Hamiltonian 模拟协议，接受块编码与时长并返回实现操作。
        times: 离散求积时刻序列。
        weights: 与 ``times`` 等长的复权重序列。

    Returns:
        callable: 形如 ``(generator, initial, final_time) -> StateOracle`` 的生成函数，
        各时刻演化按权重经 LCU 组合后作用到初态。

    Raises:
        ValidationError: ``times`` 与 ``weights`` 长度不一致。
    """
    times, weights = tuple(times), tuple(weights)
    if len(times) != len(weights):
        raise ValidationError("LCHS 离散时间与权重长度不同")

    def generate(
        generator: BlockEncoding, initial: StatePreparation, final_time: float
    ) -> StateOracle:
        """按权重对各时刻演化做 LCU 加权求和并作用到初态。"""
        terms: list[tuple[complex, BlockEncoding]] = []
        for weight, time in zip(weights, times, strict=True):
            op = ham_sim(generator, time * final_time)
            terms.append((weight, BlockEncoding(annotate(op, "block_encoding", be_alpha=1.0))))
        return apply_be_to_state(lcu(terms), initial)

    return generate


def make_schrodingerisation_qode(
    embedding: Callable[[BlockEncoding], BlockEncoding],
    ham_sim: Callable[[BlockEncoding, float], Operation],
    *,
    extra_width: int = 1,
) -> Callable[[BlockEncoding, StatePreparation, float], StateOracle]:
    """显式保留 Hamiltonian lift 的实现边界，随后组装演化和物理通道。

    Args:
        embedding: Hamiltonian lift；把生成元块编码映射为扩大 ``extra_width`` 位
            的块编码。
        ham_sim: 形如 (BE, time) 返回 Operation 的哈密顿量模拟实现。
        extra_width: lift 增加的辅助位数，取正整数，缺省为 1。

    Returns:
        Callable[[BlockEncoding, StatePreparation, float], StateOracle]: 接受生成元、
        初态与末时刻，输出物理通道读出态 oracle 的生成函数。
    """

    def generate(
        generator: BlockEncoding, initial: StatePreparation, final_time: float
    ) -> StateOracle:
        """组装 lift 后的演化并读出物理通道子空间。"""
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
