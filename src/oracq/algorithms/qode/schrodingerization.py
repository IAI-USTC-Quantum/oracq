"""通过辅助坐标和 Fourier 变换构造线性非酉演化的量子表示。"""

from __future__ import annotations

import json
import math
from collections.abc import Callable
from dataclasses import dataclass

from oracq.algorithms.common.fourier import qft_with_work as qft
from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.common.state_preparation import apply_be_to_state, select_subspace
from oracq.algorithms.input_model.block_encoding import lcu, tensor
from oracq.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.interfaces import (
    as_block_encoding,
    as_state_preparation,
    operator_state_contract,
)
from oracq.algorithms.input_model.operators import BlockEncoding, _name, identity
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.algorithms.qode._dynamics import tagged
from oracq.algorithms.qode.ode_models import HermitianParts
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, ValidationError


def fourier_momentum(width: int, period: float) -> BlockEncoding:
    """频率对角 BE，按位的投影求和；避免构造稠密矩阵。

    Args:
        width: 频率寄存器位宽，投影按位求和时逐位扫描。
        period: 辅助网格周期，必须为正的有限实数，决定频率间隔 2π/period。

    Returns:
        BlockEncoding: 以 2πk/period 为对角频率的动量算符块编码。
    """
    terms: list[tuple[float, BlockEncoding]] = []
    for bit in range(width):
        b = Builder(_name("momentum_bit", width, bit), {"target": Bits(width), "signal": Bits(1)})
        b.x(b["signal"])
        with b.control(b["target"][bit]):
            b.x(b["signal"])
        be = BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))
        coefficient = (1 << bit) * 2 * math.pi / period * (-1 if bit == width - 1 else 1)
        terms.append((coefficient, be))
    return lcu(terms)


@dataclass(frozen=True)
class SchrodingerPlan:
    """Schrödingerization 辅助网格与读出通道的配置计划。

    Attributes:
        auxiliary_width: 辅助 p 寄存器位数，有效范围为 1..63。
        period: 辅助网格的周期，必须为正的有限实数。
        selected_index: 最终读出物理解的通道编号，范围为 0..2**auxiliary_width-1。

    Raises:
        ValidationError: 任一字段越界或不是有限数值。
    """

    auxiliary_width: int = 2
    period: float = 8.0
    selected_index: int = 1

    def __post_init__(self) -> None:
        """校验辅助位数、周期与读出通道编号的范围。"""
        positive_integer(self.auxiliary_width, "SchrodingerPlan.auxiliary_width", maximum=63)
        finite_real(self.period, "SchrodingerPlan.period", minimum=0, strict=True)
        positive_integer(
            self.selected_index,
            "SchrodingerPlan.selected_index",
            minimum=0,
            maximum=(1 << self.auxiliary_width) - 1,
        )
        if (
            self.auxiliary_width < 1
            or self.period <= 0
            or not 0 <= self.selected_index < (1 << self.auxiliary_width)
        ):
            raise ValidationError("Schrodingerization 辅助网格/通道无效")


def schrodinger_qode(
    generator: BlockEncoding,
    initial: StatePreparation,
    time: float,
    *,
    plan: SchrodingerPlan | None = None,
    hamiltonian_function: Callable[[BlockEncoding, float], BlockEncoding] = taylor_hamiltonian,
) -> StateOracle:
    """u'=Gu，G=H1+iH2；Fourier lift −P⊗H1-I⊗H2（正 QFT 约定下恢复正向流）。

    Args:
        generator: 演化生成元 G 的块编码。
        initial: 物理初态的制备句柄。
        time: 演化时长，取非负有限实数。
        plan: 辅助网格与读出通道配置；缺省为 2 位辅助位、周期 8、通道 1。
        hamiltonian_function: 形如 (BE, time) 返回 BlockEncoding 的哈密顿量模拟实现，
            缺省为截断 Taylor。

    Returns:
        StateOracle: 选定辅助通道上的物理解读出态 oracle，correctness 标记为 pending。
    """
    operator_state_contract("schrodingerization").check(
        generator=generator, initial=initial
    ).require()
    generator, initial = as_block_encoding(generator), as_state_preparation(initial)
    finite_real(time, "schrodingerization.time", minimum=0)
    plan = plan or SchrodingerPlan()
    require_instance(plan, SchrodingerPlan, "schrodingerization.plan")
    parts = HermitianParts.from_operator(generator)
    n, p = generator.width, plan.auxiliary_width
    momentum = fourier_momentum(p, plan.period)
    hamiltonian = lcu([(-1, tensor(momentum, parts.hermitian)), (-1, tensor(identity(p), parts.h))])
    evolution = hamiltonian_function(hamiltonian, time)
    require_instance(evolution, BlockEncoding, "schrodingerization.hamiltonian_function.output")
    grid = [
        (j if j < (1 << (p - 1)) else j - (1 << p)) * plan.period / (1 << p) for j in range(1 << p)
    ]
    warp = gate_state_prep([math.exp(-abs(x)) for x in grid])
    transform = qft(p)
    prep = Builder(
        _name("schrod_warp_initial", initial.operation, plan),
        {"target": Bits(n + p), "work": Bits(initial.work_width)},
        resources_for(("initial", initial.operation)),
    )
    invoke(prep, initial.operation, "initial", target=prep["target"][:n], work=prep["work"])
    invoke(prep, warp.operation, target=prep["target"][n:], work=prep["work"][:0])
    invoke(prep, transform, target=prep["target"][n:], work=prep["work"][:0])
    state = apply_be_to_state(
        evolution, StatePreparation(annotate(prep.finish(), "state_prep_isometry", zero_input=True))
    )
    out = Builder(
        _name("schrod_inverse_fourier", state.operation),
        {"target": Bits(n + p), "signal": Bits(state.signal_qubits)},
        resources_for(("state", state.operation)),
    )
    invoke(out, state.operation, "state", target=out["target"], signal=out["signal"])
    with out.adjoint():
        invoke(out, transform, target=out["target"][n:], work=out["signal"][:0])
    selected = select_subspace(
        StateOracle(out.finish()),
        n,
        plan.selected_index,
        label="schrodingerization_physical_channel",
    )
    return StateOracle(
        tagged(
            selected.operation,
            "schrodingerization_qode",
            auxiliary_grid=json.dumps(grid),
            selected_p=grid[plan.selected_index],
            recovery_scale=math.exp(grid[plan.selected_index]),
            recovery_assumption="selected p in valid warped region; periodic truncation pending",
        )
    )
