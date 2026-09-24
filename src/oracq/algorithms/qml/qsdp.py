"""量子半定规划（QSDP）框架：Gibbs 采样 + 迹估计 + 矩阵乘权外层循环。

按 Brandão–Svore 2017（"Quantum speed-ups for semidefinite programming"，FOCS）
与 van Apeldoorn–Gilyén 2019 的结构组织：可行性 SDP 的解经矩阵乘权（MMW）
迭代逼近，每轮的关键量子子程序是 (a) 惩罚 Hamiltonian 的 Gibbs 态制备
（复用 density.gibbs_purification 的 QSVT 纯化路线）与 (b) 对观测量 A_i 的
迹估计 ``Tr(A_i ρ)`` （纯化态 + BE 的 Hadamard 型探针电路）。

外层 MMW 驱动是经典的（与各 QSDP 论文一致——量子的加速正在于 Gibbs 制备与
迹估计两个内层原语）；驱动通过 estimator 回调消费 ``Tr(A_i ρ)``，缺省使用
density.gibbs_state 的经典参考实现以便小实例对拍，量子路径由
iteration_circuits 逐轮生成。
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.density import (
    ApproximatePurification,
    PurificationAccess,
    _check_hermitian,
    gibbs_purification,
    gibbs_state,
)
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import invoke, resources_for
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError


def trace_estimate_circuit(
    purification: PurificationAccess | ApproximatePurification,
    observable: BlockEncoding,
    *,
    component: str = "real",
    name: str | None = None,
) -> Operation:
    """``Tr(M ρ)/α`` 的探针估计电路：纯化态 + BE 的受控调用（Hadamard test 结构）。

    probe 的 Z 期望等于 ``Re/Im⟨ψ_ρ| (M/α ⊗ I) |ψ_ρ⟩``；``Tr(M ρ)`` 由
    trace_from_probe 解码。environment 寄存器只旁观（对其取偏迹即 ρ），
    BE 的 signal 从全零参与调用，电路不执行测量。
    purification 也可以是 ApproximatePurification（近似纯化）：此时纯化态
    只存在于其 signal == 0 分支，probe 读数必须对该分支条件化，解码用
    trace_from_joint。

    Args:
        purification: Gibbs 态的纯化访问句柄，system 宽度须与观测量一致。
        observable: 观测量 M 的块编码。
        component: 读取的分量，取 ``real`` 或 ``imag``。
        name: 生成的模块名；缺省自动生成。

    Returns:
        Operation: 由 probe 寄存器读出的迹估计探针电路，不执行测量。
    """
    approximate = isinstance(purification, ApproximatePurification)
    if not approximate and not isinstance(purification, PurificationAccess):
        raise ValidationError("trace_estimate_circuit 需要 PurificationAccess 或 ApproximatePurification")
    require_instance(observable, BlockEncoding, "trace_estimate_circuit.observable")
    if observable.width != purification.width:
        raise ValidationError("观测量 BE 的目标宽度必须等于纯化的 system 宽度")
    if component not in ("real", "imag"):
        raise ValidationError("component 必须是 real 或 imag")
    n, m, s = purification.width, purification.environment_width, observable.signal_qubits
    registers = {
        "system": Bits(n),
        "environment": Bits(m),
        "signal": Bits(s),
        "probe": Bits(1),
    }
    purif_args = {"system": "system", "environment": "environment"}
    if approximate:
        registers["purification_signal"] = Bits(cast("ApproximatePurification", purification).signal_qubits)
        purif_args["signal"] = "purification_signal"
    b = Builder(
        name or _name("trace_estimate", purification.operation, observable.operation, component),
        registers,
        resources_for(("purification", purification.operation), ("be", observable.operation)),
        attributes={
            "algorithm": "trace_estimate",
            "readout_register": "probe",
            "component": component,
            "be_alpha": observable.alpha,
            "decoder": "trace_from_joint" if approximate else "trace_from_probe",
        },
    )
    invoke(
        b,
        purification.operation,
        "purification",
        **{key: b[value] for key, value in purif_args.items()},
    )
    b.h(b["probe"])
    with b.control(b["probe"]):
        invoke(b, observable.operation, "be", target=b["system"], signal=b["signal"])
    if component == "imag":
        b.gate("phase", b["probe"], -math.pi / 2)
    b.h(b["probe"])
    return b.finish()


def trace_from_probe(probability_one: float, alpha: float) -> float:
    """由 probe 测得 1 的概率解码 ``Tr(M ρ)``：Z 期望 = 1 − 2p，乘上 BE 的 α。

    Args:
        probability_one: probe 测得 1 的概率，取 [0,1]。
        alpha: 观测量块编码的尺度 α，取正实数。

    Returns:
        float: ``Tr(M ρ)`` 的估计值。
    """
    finite_real(probability_one, "trace_from_probe.probability_one", minimum=0)
    if probability_one > 1:
        raise ValidationError("probe 概率必须在 [0,1] 内")
    finite_real(alpha, "trace_from_probe.alpha", minimum=0, strict=True)
    return alpha * (1.0 - 2.0 * probability_one)


def trace_from_joint(
    probe_expectation_joint: float, weight_signal_zero: float, alpha: float
) -> float:
    """近似纯化路径的解码：``Tr(M ρ) = α · E[Z_probe·1_{signal=0}] / P(signal=0)``。

    probe_expectation_joint 是 probe 的 Z 期望与纯化 signal == 0 指示的联合期望，
    weight_signal_zero 是纯化 signal == 0 分支的概率。

    Args:
        probe_expectation_joint: probe 的 Z 期望与 signal==0 指示的联合期望。
        weight_signal_zero: 纯化 signal==0 分支的概率，取正实数。
        alpha: 观测量块编码的尺度 α，取正实数。

    Returns:
        float: 对成功分支条件化后的 ``Tr(M ρ)`` 估计值。
    """
    finite_real(weight_signal_zero, "trace_from_joint.weight_signal_zero", minimum=0, strict=True)
    finite_real(alpha, "trace_from_joint.alpha", minimum=0, strict=True)
    return alpha * probe_expectation_joint / weight_signal_zero


@dataclass(frozen=True)
class SdpInstance:
    """可行性 SDP 输入模型：找 ``X ⪰ 0``、``Tr X = 1`` 使 ``|Tr(A_i X) − b_i| ≤ ε``。

    constraints 为 (A_i, b_i) 序列：A_i 是显式小 Hermitian 矩阵，b_i 为实数。
    显式矩阵路径面向小实例；大规模时 A_i 应替换为稀疏/BE 访问模型，
    惩罚 Hamiltonian 的 BE 组装随之替换（见 iteration_circuits 的参数）。
    """

    constraints: tuple[tuple[tuple[tuple[complex, ...], ...], float], ...]

    def __post_init__(self) -> None:
        """校验并规范化各约束为 ``(A_i, b_i)`` 二元组后存回字段。"""
        constraints = tuple(self.constraints)
        if not constraints:
            raise ValidationError("SdpInstance 至少需要一个约束")
        normalized = []
        for i, item in enumerate(constraints):
            if len(item) != 2:
                raise ValidationError("每个约束必须是 (A_i, b_i) 对")
            matrix, bound = item
            matrix = _check_hermitian(matrix, f"SdpInstance.constraints[{i}]")
            if len(matrix) != len(constraints[0][0]):
                raise ValidationError("各约束矩阵维度必须一致")
            finite_real(bound, f"SdpInstance.constraints[{i}].b")
            normalized.append((matrix, float(bound)))
        object.__setattr__(self, "constraints", tuple(normalized))

    @property
    def width(self) -> int:
        """系统量子位数。"""
        return (len(self.constraints[0][0]) - 1).bit_length()

    @property
    def num_constraints(self) -> int:
        """约束条数 m，即 ``constraints`` 中 (A_i, b_i) 对的数量。"""
        return len(self.constraints)


def penalty_hamiltonian(
    instance: SdpInstance, weights: Iterable[float]
) -> tuple[tuple[complex, ...], ...]:
    """MMW 惩罚 Hamiltonian ``H = Σ_i w_i (A_i − b_i I)`` （显式小矩阵）。

    Args:
        instance: 可行性 SDP 输入模型，提供各 (A_i, b_i) 约束。
        weights: 与约束一一对应的权重序列，长度等于约束数。

    Returns:
        tuple[tuple[complex, ...], ...]: 惩罚矩阵，按行嵌套的元组表示。
    """
    require_instance(instance, SdpInstance, "penalty_hamiltonian.instance")
    weights = tuple(weights)
    if len(weights) != instance.num_constraints:
        raise ValidationError("权重数量必须等于约束数")
    d = len(instance.constraints[0][0])
    return tuple(
        tuple(
            sum(w * (a[i][j] - (b if i == j else 0)) for w, (a, b) in zip(weights, instance.constraints, strict=True))
            for j in range(d)
        )
        for i in range(d)
    )


def classical_estimator(
    hamiltonian_matrix: Iterable[Iterable[complex]],
    beta: float,
    instance: SdpInstance,
) -> tuple[tuple[tuple[complex, ...], ...], tuple[float, ...]]:
    """经典参考估计：由显式 Gibbs 态计算全部 ``Tr(A_i ρ)`` （供对拍与驱动缺省）。

    返回 (rho, estimates)：rho 为密度矩阵（嵌套元组），estimates 为各约束的迹。

    Args:
        hamiltonian_matrix: 显式 Hermitian 惩罚矩阵。
        beta: Gibbs 分布的逆温度，取正实数。
        instance: 提供 (A_i, b_i) 约束序列的 SDP 输入模型。

    Returns:
        tuple[tuple[tuple[complex, ...], ...], tuple[float, ...]]: (Gibbs 密度矩阵
        ρ, 各约束的迹估计序列)。
    """
    rho = gibbs_state(hamiltonian_matrix, beta)
    estimates = tuple(
        sum((a[i][j] * rho[j][i]).real for i in range(len(rho)) for j in range(len(rho)))
        for a, _ in instance.constraints
    )
    return rho, estimates


def qsdp_gibbs_solve(
    instance: SdpInstance,
    *,
    epsilon: float,
    max_iterations: int | None = None,
    estimator: (
        Callable[
            [Iterable[Iterable[complex]], float, SdpInstance],
            tuple[tuple[tuple[complex, ...], ...], tuple[float, ...]],
        ]
        | None
    ) = None,
) -> dict[str, tuple[tuple[complex, ...], ...] | tuple[float, ...] | int | bool]:
    """矩阵乘权（MMW）可行性求解：对称零和博弈 ``min_X max_w Σ_i w_i v_i`` 的 Hedge 迭代。

    约束语义为单边不等式 ``Tr(A_i X) ≤ b_i``（违反量 ``v_i = Tr(A_i ρ) − b_i``，
    只有正违反被惩罚）；等式约束应拆成 ``(A_i, b_i)`` 与 ``(−A_i, −b_i)`` 两个
    不等式。每轮：惩罚 Hamiltonian ``H_t = Σ_i w_i (A_i − b_i I)`` 的 Gibbs 态
    ρ_t，约束违反量 ``v_i = Tr(A_i ρ_t) − b_i``，权重按
    ``w_i ∝ exp(η·v_i)`` 更新。平均迭代 ``ρ̄`` 的最大违反量随
    ``O(√(ln m / T))`` 下降（Brandão–Svore 2017 的结构；量子的加速在
    Gibbs 制备与迹估计两个内层原语，本驱动的经典循环与论文一致）。

    Args:
        instance: SdpInstance。
        epsilon: 目标违反量上界，学习率取 η = ε/4。
        max_iterations: 迭代上限，缺省 ``ceil(64·ln(m+1)/ε²)``；平均迭代的
            最大违反量 ≤ epsilon 时提前停止。
        estimator: 回调 (hamiltonian_matrix, beta, instance) -> (rho, estimates)，
            缺省为 classical_estimator；量子路径应传入基于
            iteration_circuits 的估计器。

    Returns:
        dict: rho（平均迭代，嵌套元组）、violations（最终违反量）、
        iterations、converged（最大违反量 ≤ epsilon）。
    """
    require_instance(instance, SdpInstance, "qsdp_gibbs_solve.instance")
    finite_real(epsilon, "qsdp_gibbs_solve.epsilon", minimum=0, strict=True)
    m = instance.num_constraints
    if max_iterations is None:
        max_iterations = math.ceil(64 * math.log(m + 1) / (epsilon * epsilon)) + 1
    positive_integer(max_iterations, "qsdp_gibbs_solve.max_iterations", minimum=1)
    estimator = estimator or classical_estimator
    eta = epsilon / 4.0
    log_weights = [0.0] * m
    d = len(instance.constraints[0][0])
    rho_average = [[0.0 + 0j] * d for _ in range(d)]
    iterations_run = 0

    def current_violations(count: int) -> tuple[float, ...]:
        """以 ``rho_average / count`` 为平均迭代时各约束的违反量。"""
        return tuple(
            sum((a[i][j] * rho_average[j][i]).real for i in range(d) for j in range(d)) / count
            - bound
            for a, bound in instance.constraints
        )

    converged = False
    for iteration in range(1, max_iterations + 1):
        shift = max(log_weights)
        weights = [math.exp(w - shift) for w in log_weights]
        total = sum(weights)
        weights = [w / total for w in weights]
        # H_t = η·Σ_{s<t} M_s 的增量形式：直接以 β = 1、矩阵为惩罚阵的 Gibbs 态。
        hamiltonian = penalty_hamiltonian(instance, weights)
        rho, estimates = estimator(hamiltonian, 1.0, instance)
        for i in range(d):
            for j in range(d):
                rho_average[i][j] += rho[i][j]
        iterations_run = iteration
        for k, (_, bound) in enumerate(instance.constraints):
            log_weights[k] += eta * (estimates[k] - bound)
        if iteration % 32 == 0 or iteration == max_iterations:
            if max(current_violations(iteration)) <= epsilon:
                converged = True
                break
    violations = current_violations(iterations_run)
    converged = converged or max(violations) <= epsilon
    return {
        "rho": tuple(
            tuple(value / iterations_run for value in row) for row in rho_average
        ),
        "violations": violations,
        "iterations": iterations_run,
        "converged": converged,
    }


def iteration_circuits(
    instance: SdpInstance,
    weights: Iterable[float],
    beta: float,
    *,
    hamiltonian_encoding: Callable[[tuple[tuple[complex, ...], ...]], BlockEncoding],
    error: float = 0.05,
) -> tuple[ApproximatePurification, list[Operation]]:
    """生成单轮 MMW 迭代的量子子程序：Gibbs 纯化 + 各约束的迹估计电路。

    hamiltonian_encoding 是把显式惩罚矩阵编成 BE 的 callable（小实例可用
    matrix_pauli_encoding；大规模换稀疏/低秩访问模型）。返回
    (ApproximatePurification, [trace_estimate Operation, ...])，供量子路径
    的估计器逐轮调用。

    Args:
        instance: 可行性 SDP 输入模型，提供各 (A_i, b_i) 约束。
        weights: 本轮迭代的约束权重序列，长度等于约束数。
        beta: Gibbs 纯化的逆温度，取正实数。
        hamiltonian_encoding: 把显式惩罚矩阵编成块编码的回调。
        error: 纯化近似误差，取正实数。

    Returns:
        tuple[ApproximatePurification, list[Operation]]: Gibbs 纯化句柄与各约束
        的迹估计电路列表。
    """
    require_instance(instance, SdpInstance, "iteration_circuits.instance")
    finite_real(beta, "iteration_circuits.beta", minimum=0, strict=True)
    finite_real(error, "iteration_circuits.error", minimum=0, strict=True)
    hamiltonian = penalty_hamiltonian(instance, weights)
    be = hamiltonian_encoding(hamiltonian)
    require_instance(be, BlockEncoding, "iteration_circuits.hamiltonian_encoding")
    purification = gibbs_purification(be, beta, error=error)
    from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding

    traces = [
        trace_estimate_circuit(purification, matrix_pauli_encoding(a))
        for a, _ in instance.constraints
    ]
    return purification, traces
