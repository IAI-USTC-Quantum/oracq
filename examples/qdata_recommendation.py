"""Demo of the QMatrix sample-and-query data structure + the KP quantum recommendation system."""

from __future__ import annotations

from oracq import estimate_resources, simulate
from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.input_model.qdata import QMatrix
from oracq.algorithms.qml.recommendation import KPRecommendationConfig, kp_recommendation


def preference_matrix() -> QMatrix:
    """A rank-one dominated preference matrix: category 0 of the three product categories is the most popular, plus a small diagonal perturbation."""
    base = (0.5, 0.25, 0.25, 0.0)
    scales = (1.0, 0.5, 0.5, 1.0)
    eps = 0.0625
    return QMatrix(
        [[scales[i] * base[j] + (eps if i == j else 0.0) for j in range(4)] for i in range(4)],
        fmt=FixedFormat(8, 4),
        angle_width=12,
    )


def main() -> None:
    """Demonstrate the readout and resource statistics of the QMatrix query structure and the KP recommendation circuit."""
    matrix = preference_matrix()
    print("‖A‖_F =", round(matrix.frobenius, 6))
    print("user distribution |Ã⟩ =", [round(a, 4) for a in matrix.user_amplitudes()])
    print("row vector of user 1 =", [round(a, 4) for a in matrix.row_amplitudes(1)])

    result = kp_recommendation(matrix, 1, KPRecommendationConfig(precision=5, sigma=0.45))
    state = simulate(result.operation.program(), result.memories())
    success, distribution = result.readout(state)
    print("flag success probability =", round(success, 4))
    print("conditional recommendation distribution =", {k: round(p, 4) for k, p in sorted(distribution.items())})
    print("(principal direction of the classical A_{≥σ}A⁺_{≥σ}Ā₁: product 0)")

    estimate = estimate_resources(result.operation.program())
    print(
        "resource ledger: qubits =", estimate.qubits,
        "toffoli =", estimate.toffoli,
        "qram_queries =", dict(estimate.qram_queries),
    )


if __name__ == "__main__":
    main()
