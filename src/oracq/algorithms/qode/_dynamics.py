"""非酉演化算法共享的 LCU 分支组装工具。"""

from __future__ import annotations

from collections.abc import Callable, Sequence

from oracq.algorithms.common.state_preparation import apply_be_to_state
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.contracts import finite_real, require_instance
from oracq.algorithms.input_model.interfaces import (
    operator_state_contract,
)
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    annotate,
)
from oracq.algorithms.qode.ode_models import LinearODE
from oracq.infrastructure.builder import Operation
from oracq.infrastructure.ir import ValidationError


def tagged(operation: Operation, algorithm: str, **metadata: str | int | float) -> Operation:
    """以 ``unitary`` 范式标注演化操作，登记算法名与正确性待验证的元数据。

    Args:
        operation: 待标注的演化操作。
        algorithm: 登记到属性中的算法名。
        **metadata: 追加登记的算法元数据，取值为字符串或数值。

    Returns:
        Operation: 带 ``unitary`` 范式与 correctness=pending 标注的操作。
    """
    return annotate(
        operation,
        "unitary",
        algorithm=algorithm,
        correctness="pending",
        validation_stage="paradigm",
        **metadata,
    )


def _lcu_dynamics(
    model: LinearODE,
    time: float,
    nodes: Sequence[float],
    weights: Sequence[complex],
    hamiltonian_function: Callable[[BlockEncoding, float], BlockEncoding],
    algorithm: str,
    **metadata: str,
) -> StateOracle:
    """逐节点构造 K=H+kL 的 Hamiltonian 分支，经 LCU 组合并作用到初态。"""
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
    terms: list[tuple[complex, BlockEncoding]] = []
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
