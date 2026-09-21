"""pyqecclang 侧 T2：QSVT 2×2 矩阵求逆方向（A=[[1,-1/3],[-1/3,1]], |b>=[1,0]）。

规格与判定阈值见 ~/projects/pyqecclang-dev/benchmarks/t2/SPEC.md。
独立可运行：

    cd ~/projects/qcfd-dev/pyqecclang && \
    PYTHONPATH=src ~/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python \
    tools/expressiveness/t2_qsvt_inversion.py

组装路径为主仓文档化原语（tests/verification/verify_fourier.py 同款）：
matrix_pauli_encoding（α = Pauli 1-范数 = 4/3 = ‖A‖）+ qsp_phases（J-多项式求逆相位，
虚部补全同 qsvt_matrix_inversion）+ qsvt_sequence，参考执行器读零信号块。
"""

import math
from importlib.metadata import version

import numpy as np

import pyqecclang
from pyqecclang import simulate
from pyqecclang.algorithms.common.qsvt import qsp_phases
from pyqecclang.algorithms.common.transforms import qsvt_sequence
from pyqecclang.algorithms.input_model.block_encoding import matrix_pauli_encoding
from pyqecclang.infrastructure.ir import ValidationError

KAPPA, EPS, THRESHOLD = 8, 1e-2, 1e-2
A = np.array([[1.0, -1.0 / 3.0], [-1.0 / 3.0, 1.0]])
B = np.array([1.0, 0.0])


def inversion_poly(b):
    f = [0.0] * (2 * b)
    for m in range(1, b + 1):
        f[2 * m - 1] = ((-1.0) ** (m + 1)) * math.comb(b, m)
    grid = [math.cos(math.pi * j / 2048) for j in range(2049)]
    sup = max(abs(sum(c * x**i for i, c in enumerate(f))) for x in grid)
    scale = 1.0 / (3.0 * sup)
    return [v * scale for v in f], scale


def run(b):
    f, _ = inversion_poly(b)
    phases = qsp_phases(tuple(f), imag=(0.0, math.sqrt(1.0 - f[-1] ** 2)))
    be = matrix_pauli_encoding(A.tolist())
    amps = simulate(qsvt_sequence(be, phases).program()).amplitudes
    v = np.array([complex(amps[(0, 0)]), complex(amps[(1, 0)])])
    scaled = A / be.alpha
    h1 = math.sqrt(1.0 - f[-1] ** 2)
    coeffs = [f[i] + 1j * (h1 if i == 1 else 0.0) for i in range(len(f))]
    expected = sum(c * np.linalg.matrix_power(scaled, i)
                   for i, c in enumerate(coeffs)) @ B
    block_err = float(np.abs(v - expected).max())
    direction = v.real / np.linalg.norm(v.real)
    exact = np.linalg.inv(A) @ B
    exact_dir = exact / np.linalg.norm(exact)
    return float(np.linalg.norm(direction - exact_dir)), block_err, be.alpha


def main():
    print(f"pyqecclang {version('pyqecclang')} ({pyqecclang.__file__})")
    b_spec = math.ceil(math.log(1.0 / EPS) / -math.log(1.0 - 1.0 / KAPPA**2))
    try:
        run(b_spec)
        b = b_spec
        print(f"spec b={b_spec} (degree {2 * b_spec - 1}) accepted")
    except ValidationError as err:
        print(f"spec b={b_spec} (degree {2 * b_spec - 1}) rejected: {err}")
        b = 0
        while True:
            try:
                run(b + 1)
                b += 1
            except ValidationError:
                break
        print(f"max feasible b={b} (degree {2 * b - 1})")
    dir_err, block_err, alpha = run(b)
    print(f"block_encoding_alpha={alpha:.6f} phases={2 * b}")
    print(f"zero_block_vs_poly_err={block_err:.3e}")
    print(f"direction_err={dir_err:.3e} (threshold {THRESHOLD})")
    print("status:", "ok" if dir_err <= THRESHOLD else "failed")


if __name__ == "__main__":
    main()
