"""Expressiveness benchmark T5: Schrödingerization solving the 1D heat equation (oracq side).

Specification (benchmarks/t5/SPEC.md): u_t = u_xx, a 4-point periodic grid,
second-difference generator G = S + S^T - 2I (spectral radius 4), t = 0.05,
u0 = [1, 0.5, 0, -0.5]/‖·‖; the reference solution is scipy.linalg.expm(G t) u0.
Metric: L2 relative error of the post-selected physical channel after the
recovery factor is applied.

Run: cd ~/projects/qcfd-dev/oracq && PYTHONPATH=src <python> \
tools/expressiveness/t5_schrodingerization_heat.py
"""

import importlib.metadata
import json
import math
import sys
from functools import partial

import numpy as np
import scipy.linalg

from oracq import simulate
from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.input_model.block_encoding import lcu, matrix_pauli_encoding, tensor
from oracq.algorithms.input_model.operators import identity
from oracq.algorithms.input_model.oracles import gate_state_prep
from oracq.algorithms.qode.ode import linear_qode
from oracq.algorithms.qode.ode_models import HermitianParts
from oracq.algorithms.qode.schrodingerization import SchrodingerPlan, fourier_momentum

TIME, DEGREE = 0.05, 2


def main():
    print(f"oracq {importlib.metadata.version('oracq')} python {sys.version.split()[0]}")
    print(f"numpy {np.__version__} scipy {scipy.__version__}")

    shift = np.roll(np.eye(4), 1, axis=1)
    g_mat = shift + shift.T - 2 * np.eye(4)
    u0 = np.array([1.0, 0.5, 0.0, -0.5])
    u0 = u0 / np.linalg.norm(u0)
    exact = scipy.linalg.expm(g_mat * TIME) @ u0
    time_reversed = scipy.linalg.expm(-g_mat * TIME) @ u0

    plan = SchrodingerPlan(auxiliary_width=2, period=1.2, selected_index=1)
    generator = matrix_pauli_encoding(g_mat.tolist())
    solver = linear_qode(
        "schrodingerization",
        plan=plan,
        hamiltonian_function=partial(taylor_hamiltonian, degree=DEGREE),
    )
    state = solver(generator, gate_state_prep(list(u0)), TIME)

    # Recovery factor alpha_E: rebuild the Taylor block encoding of K'=-P⊗H1-I⊗H2 from public combinators and read off alpha
    parts = HermitianParts.from_operator(generator)
    momentum = fourier_momentum(plan.auxiliary_width, plan.period)
    k_be = lcu(
        [(1, tensor(momentum, parts.hermitian)), (-1, tensor(identity(plan.auxiliary_width), parts.h))]
    )
    alpha_e = taylor_hamiltonian(k_be, TIME, degree=DEGREE).alpha
    attrs = dict(state.operation.module.attributes)
    grid = json.loads(attrs["auxiliary_grid"])
    p_sel = attrs["selected_p"]
    z_norm = math.sqrt(sum(math.exp(-2 * abs(p)) for p in grid))
    recovery = alpha_e * z_norm * math.exp(p_sel)

    amplitudes = dict(
        simulate(
            state.operation.program(),
            max_steps=1 << 30,
            max_states=1 << 22,
        ).amplitudes
    )
    block = np.array([amplitudes.get((i, 0), 0j) for i in range(4)])
    recovered = block * recovery

    err = float(np.linalg.norm(recovered - exact) / np.linalg.norm(exact))
    reversed_fit = float(np.linalg.norm(recovered - time_reversed) / np.linalg.norm(time_reversed))
    print(f"alpha_E={alpha_e:.6f} p_sel={p_sel:.3f} Z={z_norm:.6f} recovery={recovery:.6f}")
    print(f"recovered={np.round(recovered, 6)}")
    print(f"exact     ={np.round(exact, 6)}")
    print(f"L2_rel_vs_expm={err:.3e} L2_rel_vs_time_reversed={reversed_fit:.3e}")
    print(f"direction={'forward' if err < reversed_fit else 'REVERSED(!)'}")


if __name__ == "__main__":
    main()
