"""非酉演化算法共享的 LCU 分支组装工具。"""

from __future__ import annotations

from pyqecclang.algorithms.block_encoding import lcu
from pyqecclang.algorithms.contracts import finite_real, require_instance
from pyqecclang.algorithms.interfaces import (
    operator_state_contract,
)
from pyqecclang.algorithms.ode_models import LinearODE
from pyqecclang.algorithms.operators import BlockEncoding
from pyqecclang.algorithms.oracles import (
    StateOracle,
    annotate,
)
from pyqecclang.algorithms.state_preparation import apply_be_to_state
from pyqecclang.infrastructure.ir import ValidationError


def tagged(operation, algorithm, **metadata):
    return annotate(
        operation,
        "unitary",
        algorithm=algorithm,
        correctness="pending",
        validation_stage="paradigm",
        **metadata,
    )


def _lcu_dynamics(model, time, nodes, weights, hamiltonian_function, algorithm, **metadata):
    require_instance(model, LinearODE, algorithm + ".model")
    finite_real(time, algorithm + ".time", minimum=0)
    operator_state_contract(algorithm).check(
        generator=model.parts.hermitian, initial=model.initial
    ).require()
    operator_state_contract(algorithm).check(
        generator=model.parts.h, initial=model.initial
    ).require()
    if not callable(hamiltonian_function):
        raise ValidationError("hamiltonian_function 必须可调用")
    terms = []
    for node, weight in zip(nodes, weights, strict=True):
        hk = lcu([(1, model.parts.h), (node, model.parts.hermitian)])
        encoded = hamiltonian_function(hk, time)
        if not isinstance(encoded, BlockEncoding):
            raise ValidationError(
                "Hamiltonian-function protocol 必须返回 BlockEncoding，保留 alpha"
            )
        terms.append((weight, encoded))
    evolution = lcu(terms)
    state = apply_be_to_state(evolution, model.initial)
    return StateOracle(
        tagged(
            state.operation,
            algorithm,
            evolution_alpha=evolution.alpha,
            branch_count=len(nodes),
            input_assumption="L>=0; autonomous; homogeneous",
            **metadata,
        )
    )
