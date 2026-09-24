"""QMatrix sample-and-query 数据结构 + KP 量子推荐系统演示。"""

from __future__ import annotations

from oracq import estimate_resources, simulate
from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.input_model.qdata import QMatrix
from oracq.algorithms.qml.recommendation import KPRecommendationConfig, kp_recommendation


def preference_matrix() -> QMatrix:
    """秩一主导的偏好矩阵：三类商品中第 0 类最受欢迎，附小对角扰动。"""
    base = (0.5, 0.25, 0.25, 0.0)
    scales = (1.0, 0.5, 0.5, 1.0)
    eps = 0.0625
    return QMatrix(
        [[scales[i] * base[j] + (eps if i == j else 0.0) for j in range(4)] for i in range(4)],
        fmt=FixedFormat(8, 4),
        angle_width=12,
    )


def main() -> None:
    """演示 QMatrix 查询结构与 KP 推荐电路的读出和资源统计。"""
    matrix = preference_matrix()
    print("‖A‖_F =", round(matrix.frobenius, 6))
    print("用户分布 |Ã⟩ =", [round(a, 4) for a in matrix.user_amplitudes()])
    print("用户 1 的行向量 =", [round(a, 4) for a in matrix.row_amplitudes(1)])

    result = kp_recommendation(matrix, 1, KPRecommendationConfig(precision=5, sigma=0.45))
    state = simulate(result.operation.program(), result.memories())
    success, distribution = result.readout(state)
    print("flag 成功概率 =", round(success, 4))
    print("条件推荐分布 =", {k: round(p, 4) for k, p in sorted(distribution.items())})
    print("(经典 A_{≥σ}A⁺_{≥σ}Ā₁ 的主方向：商品 0)")

    estimate = estimate_resources(result.operation.program())
    print(
        "资源台账：qubits =", estimate.qubits,
        "toffoli =", estimate.toffoli,
        "qram_queries =", dict(estimate.qram_queries),
    )


if __name__ == "__main__":
    main()
