"""通过辅助坐标和 Fourier 变换构造线性非酉演化的量子表示。"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass

from pyqecclang.algorithms._dynamics import tagged
from pyqecclang.algorithms.block_encoding import lcu, tensor
from pyqecclang.algorithms.contracts import finite_real, positive_integer, require_instance
from pyqecclang.algorithms.fourier import qft_with_work as qft
from pyqecclang.algorithms.hamiltonian import taylor_hamiltonian
from pyqecclang.algorithms.interfaces import (
    as_block_encoding,
    as_state_preparation,
    operator_state_contract,
)
from pyqecclang.algorithms.ode_models import HermitianParts
from pyqecclang.algorithms.operators import BlockEncoding, _name, identity
from pyqecclang.algorithms.oracles import (
    StateOracle,
    StatePreparation,
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from pyqecclang.algorithms.state_preparation import apply_be_to_state, select_subspace
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError


def fourier_momentum(width, period):
    """频率对角 BE，按位的投影求和；避免构造稠密矩阵。"""
    terms = []
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
    auxiliary_width: int = 2
    period: float = 8.0
    selected_index: int = 1

    def __post_init__(self):
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
    generator, initial, time, *, plan=None, hamiltonian_function=taylor_hamiltonian
):
    """u'=Gu，G=H1+iH2；Fourier lift P⊗H1-I⊗H2。"""
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
    hamiltonian = lcu([(1, tensor(momentum, parts.hermitian)), (-1, tensor(identity(p), parts.h))])
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
