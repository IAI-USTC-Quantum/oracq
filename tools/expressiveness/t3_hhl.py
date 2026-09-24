"""表达力基准 T3：HHL 2×2 的 oracq QLSS 路径。

任务：A=[[1,-1/3],[-1/3,1]]，b=[1,0]；用 CKS Chebyshev 求解器
（oracq.algorithms.qlss.cks_chebyshev，稀疏访问输入模型）求条件解态，
见证为恢复向量方向对 numpy.linalg.solve 独立参考的误差
（Chebyshev 截断方法误差 + 定点量化误差，阈值见 benchmarks/t3/SPEC.md）。
"""

import math
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np

import oracq
from oracq import FixedFormat, simulate
from oracq.algorithms.input_model.oracles import (
    SparseAccess,
    basis_state,
    gate_database,
    sparse_entry,
    sparse_location_gate,
)
from oracq.algorithms.qlss.qlss import CKSConfig, SparseSystem, SpectralPromise, cks_chebyshev

ORDER = 40
FRACTION = 8


def package_version():
    try:
        return version("oracq")
    except PackageNotFoundError:
        pkg_info = Path(oracq.__file__).parent.parent / "oracq.egg-info" / "PKG-INFO"
        for line in pkg_info.read_text(encoding="utf-8").splitlines():
            if line.startswith("Version:"):
                return line.split(":", 1)[1].strip()
        return "unknown"


def cks_coefficients(order):
    """CKS 系数独立闭式重算（math.comb，不依赖库内实现）。"""
    return [
        4.0
        * (-1) ** j
        * sum(math.comb(2 * order, order + i) for i in range(j + 1, order + 1))
        / 2.0 ** (2 * order)
        for j in range(order)
    ]


def chebyshev_apply(matrix, alpha, coefficients, vector):
    """P(M)v，P = Σ c_j T_{2j+1}，M = matrix/alpha；独立矩阵递推。"""
    m = np.asarray(matrix, dtype=float) / alpha
    t2 = 2 * (m @ m) - np.eye(m.shape[0])
    u_prev, u_curr = m.copy(), 2 * (t2 @ m) - m
    out = coefficients[0] * (u_prev @ vector)
    if len(coefficients) > 1:
        out = out + coefficients[1] * (u_curr @ vector)
    for j in range(2, len(coefficients)):
        u_next = 2 * (t2 @ u_curr) - u_prev
        u_prev, u_curr = u_curr, u_next
        out = out + coefficients[j] * (u_curr @ vector)
    return out


def direction_error(vector, reference):
    v = np.asarray(vector) / np.linalg.norm(vector)
    r = np.asarray(reference) / np.linalg.norm(reference)
    return float(np.sqrt(max(0.0, 2.0 - 2.0 * abs(np.vdot(v, r)))))


def main():
    print(f"oracq {package_version()} / numpy {version('numpy')}")
    matrix = [[1.0, -1 / 3], [-1 / 3, 1.0]]
    b = np.array([1.0, 0.0])
    fmt = FixedFormat(FRACTION + 2, FRACTION)
    quantized = [[fmt.decode(fmt.encode(value)) for value in row] for row in matrix]
    entry = sparse_entry(
        gate_database(
            2,
            fmt.width,
            {r + (c << 1): fmt.encode(quantized[r][c]) for r in range(2) for c in range(2)},
        ),
        1,
    )
    system = SparseSystem(
        SparseAccess(
            sparse_location_gate(1, [[0, 1], [0, 1]], work_width=0), entry, 1, fmt.width, 2
        ),
        fmt,
        1.0,
        basis_state(1),
        SpectralPromise(1.35, 0.65),
        diagonal_nonnegative=True,
        hermitian=True,
    )
    state = cks_chebyshev(system, CKSConfig(order=ORDER))
    amplitudes = dict(
        simulate(state.operation.program(), max_steps=2_000_000_000, max_states=1 << 22).amplitudes
    )
    good = {key[0]: amp for key, amp in amplitudes.items() if key[1] == 0}
    quantum = np.array([good.get(0, 0j), good.get(1, 0j)])
    success = sum(abs(value) ** 2 for value in good.values())
    a_q = np.array(quantized)
    alpha = 2.0
    polynomial = chebyshev_apply(a_q, alpha, cks_coefficients(ORDER), b)
    impl = direction_error(quantum, polynomial)
    method = direction_error(polynomial, np.linalg.solve(a_q, b))
    task = direction_error(polynomial, np.linalg.solve(np.array(matrix), b))
    print(f"order={ORDER} alpha={alpha} lcu_normalization={sum(abs(c) for c in cks_coefficients(ORDER)):.6f}")
    print(f"success_probability={success:.6f}")
    print(f"quantum_vector={quantum.tolist()}")
    print(f"reference_solve={np.linalg.solve(np.array(matrix), b).tolist()}")
    print(f"impl_direction_error_vs_chebyshev_oracle={impl:.3e}")
    print(f"method_direction_error_vs_solve(A_quantized)={method:.3e}")
    print(f"task_direction_error_vs_solve(A)={task:.3e}")
    ok = impl < 1e-8 and task < 1e-2
    print(f"T3 oracq QLSS: {'PASS' if ok else 'FAIL'} (threshold: impl<1e-8, direction<1e-2)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
