"""表达力基准 T4：LCHS 装配 u'=Gu 的 oracq 路径。

任务：G=[[-1,0.5],[-0.5,-1]]（耗散+反对称），t=0.05，u0=[1,0]；
linear_qode("lchs", QuadraturePlan.cauchy(cutoff=1, spacing=1.0),
每支路 degree-1 Taylor)（An–Childs–Lin–Ying 伪代码的有限求积装配）。
见证：后选择块对独立 numpy 仿真（实现误差）与 scipy expm 精确解（方法误差）。
"""

import math
from functools import partial
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import scipy.linalg

import oracq
from oracq import simulate
from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding
from oracq.algorithms.input_model.oracles import gate_state_prep
from oracq.algorithms.qode.lchs import QuadraturePlan
from oracq.algorithms.qode.ode import linear_qode

G = [[-1.0, 0.5], [-0.5, -1.0]]
TIME = 0.05
U0 = [1.0, 0.0]
PLAN = QuadraturePlan.cauchy(cutoff=1, spacing=1.0)
DEGREE = 1


def package_version():
    try:
        return version("oracq")
    except PackageNotFoundError:
        pkg_info = Path(oracq.__file__).parent.parent / "oracq.egg-info" / "PKG-INFO"
        for line in pkg_info.read_text(encoding="utf-8").splitlines():
            if line.startswith("Version:"):
                return line.split(":", 1)[1].strip()
        return "unknown"


def lchs_reference(alpha_v):
    """独立 numpy 仿真：V = Σ_k w_k Σ_{p≤degree} (-it)^p/p! (H+kL)^p，A=-G，L=sym(A)，H=skew(A)/i。"""
    a_mat = -np.array(G)
    l_mat = (a_mat + a_mat.conj().T) / 2
    h_mat = (a_mat - a_mat.conj().T) / 2j
    v_mat = np.zeros((2, 2), dtype=complex)
    for node, weight in zip(PLAN.nodes, PLAN.weights, strict=True):
        branch = h_mat + node * l_mat
        taylor = np.zeros((2, 2), dtype=complex)
        power = np.eye(2, dtype=complex)
        for order in range(DEGREE + 1):
            taylor = taylor + ((-1j * TIME) ** order / math.factorial(order)) * power
            power = branch @ power
        v_mat = v_mat + weight * taylor
    return v_mat @ np.array(U0) / alpha_v


def direction_error(vector, reference):
    v = np.asarray(vector) / np.linalg.norm(vector)
    r = np.asarray(reference) / np.linalg.norm(reference)
    return float(np.sqrt(max(0.0, 2.0 - 2.0 * abs(np.vdot(v, r)))))


def module_attribute(program, key):
    for module in program.modules:
        attributes = dict(module.attributes)
        if key in attributes:
            return attributes[key]
    raise KeyError(key)


def main():
    print(f"oracq {package_version()} / numpy {version('numpy')} / scipy {version('scipy')}")
    solver = linear_qode(
        "lchs", plan=PLAN, hamiltonian_function=partial(taylor_hamiltonian, degree=DEGREE)
    )
    state = solver(matrix_pauli_encoding(G), gate_state_prep(U0), TIME)
    program = state.operation.program()
    alpha_v = module_attribute(program, "evolution_alpha")
    amplitudes = dict(simulate(program, max_steps=10_000_000).amplitudes)
    good = {key[0]: amp for key, amp in amplitudes.items() if all(v == 0 for v in key[1:])}
    quantum = np.array([good.get(0, 0j), good.get(1, 0j)])
    expected = lchs_reference(alpha_v)
    exact = scipy.linalg.expm(np.array(G) * TIME) @ np.array(U0)
    impl = float(np.abs(quantum - expected).max())
    direction = direction_error(quantum, exact)
    method_amplitude = float(np.abs(quantum - exact / alpha_v).max())
    success = sum(abs(value) ** 2 for value in good.values())
    print(f"alpha_V={alpha_v} success_probability={success:.6f}")
    print(f"quantum_block={quantum.tolist()}")
    print(f"independent_emulation={expected.tolist()}")
    print(f"expm_exact={exact.tolist()}")
    print(f"impl_error={impl:.3e}")
    print(f"direction_error_vs_expm={direction:.3e}")
    print(f"method_amplitude_error={method_amplitude:.3e} (有限求积+Taylor 余项)")
    ok = impl < 1e-9 and direction < 1e-2
    print(f"T4 oracq LCHS: {'PASS' if ok else 'FAIL'} (threshold: impl<1e-9, direction<1e-2)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
