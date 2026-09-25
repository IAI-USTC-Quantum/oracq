"""Expressiveness benchmark T6: QSP solving the 1D advection equation (oracq side, re-implementation of a qsp4pde example).

Specification (benchmarks/t6/SPEC.md, aligned with the qsp4pde README quick
start):
du/dt + r du/dx = 0, N = 16 points on a periodic grid, r = 1.0, t = 0.2,
initial value gaussian_shifted (mu = -0.25, Sigma = 0.1); the reference
solution is the exact Fourier spectral shift
psi(t) = ifft(exp(-i 2π k r t) fft(u0)).

Assembly: Fourier diagonal generator D = diag(2π·k_signed) (diagonal block
encoding via an angle database, alpha = ‖D‖ = 8); propagation of e^{i τ D/8}
(τ = +2π r t·8, which under the positive-sign QFT convention is the forward
shift) is realized via the QSVT phase polynomial - qsvt_hamiltonian_simulation
internally is qsp_phases phase synthesis + qsvt_sequence assembly (the
Jacobi-Anger branch). The library's QSP phase synthesis applies for
|τ| ≲ 2.5, so the product formula is split into m = 5 steps (each
|τ'| = 2.011), re-encoding the initial state after each step from the
signal post-selection block × sim_scale (a measure-and-prepare protocol).
QFT/inverse QFT are provided by fourier.qft_with_work.

Run: cd ~/projects/qcfd-dev/oracq && PYTHONPATH=src <python> \
tools/expressiveness/t6_qsp_advection.py
"""

import importlib.metadata
import math
import sys

import numpy as np

from oracq import simulate
from oracq.algorithms.common.fourier import inverse_qft
from oracq.algorithms.common.fourier import qft_with_work as qft
from oracq.algorithms.common.qsvt import qsvt_hamiltonian_simulation
from oracq.algorithms.common.state_preparation import apply_be_to_state
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    annotate,
    diagonal_block_encoding,
    gate_database,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits

N_QUBITS, N = 4, 16
SPEED, TIME, STEPS, WORD_WIDTH = 1.0, 0.2, 5, 14
TAU = 2 * math.pi * SPEED * TIME * 8.0


def generator_be():
    signed = [k if k < N // 2 else k - N for k in range(N)]
    step = 2 * math.pi / (1 << WORD_WIDTH)
    words = [
        min(round(2 * math.acos(max(-1.0, min(1.0, k / 8.0))) / step), (1 << WORD_WIDTH) - 1)
        for k in signed
    ]
    return diagonal_block_encoding(gate_database(N_QUBITS, WORD_WIDTH, words), alpha=8.0)


def apply_stage(evo, amplitudes, *, sandwich):
    prep = gate_state_prep(list(amplitudes))
    stage = Builder(
        "advect_stage",
        {"target": Bits(N_QUBITS), "work": Bits(prep.work_width)},
        resources_for(("initial", prep.operation)),
    )
    invoke(stage, prep.operation, "initial", target=stage["target"], work=stage["work"])
    if sandwich == "forward":
        invoke(stage, qft(N_QUBITS), target=stage["target"], work=stage["work"][:0])
    state = apply_be_to_state(
        evo, StatePreparation(annotate(stage.finish(), "state_prep_isometry", zero_input=True))
    )
    out = Builder(
        "advect_readout",
        {"target": Bits(N_QUBITS), "signal": Bits(state.signal_qubits)},
        resources_for(("state", state.operation)),
    )
    invoke(out, state.operation, "state", target=out["target"], signal=out["signal"])
    if sandwich == "inverse":
        invoke(out, inverse_qft(N_QUBITS), target=out["target"])
    amps = dict(
        simulate(out.finish().program(), max_steps=1 << 30, max_states=1 << 24).amplitudes
    )
    return np.array([amps.get((k, 0), 0j) for k in range(N)])


def main():
    print(f"oracq {importlib.metadata.version('oracq')} python {sys.version.split()[0]}")
    print(f"numpy {np.__version__}")

    x = np.arange(N) / N - 0.5 + 1 / (2 * N)
    u0 = np.exp(-((x - (-0.25)) ** 2) / (2 * 0.1**2))
    u0 = u0 / np.linalg.norm(u0)
    signed = np.array([k if k < N // 2 else k - N for k in range(N)])
    exact = np.fft.ifft(np.exp(-1j * 2 * np.pi * signed * SPEED * TIME) * np.fft.fft(u0))

    evo = qsvt_hamiltonian_simulation(generator_be(), TAU / STEPS, error=0.01)
    attrs = dict(evo.operation.module.attributes)
    sim_scale = attrs["sim_scale"]

    psi = u0
    for step in range(STEPS):
        psi = apply_stage(
            evo,
            psi,
            sandwich=("forward" if step == 0 else "inverse" if step == STEPS - 1 else "none"),
        ) * sim_scale

    err = float(np.linalg.norm(psi - exact) / np.linalg.norm(exact))
    reverse = np.fft.ifft(np.exp(+1j * 2 * np.pi * signed * SPEED * TIME) * np.fft.fft(u0))
    rev_err = float(np.linalg.norm(psi - reverse) / np.linalg.norm(reverse))
    print(f"qsp_degree={attrs['qsp_degree']} sim_scale={sim_scale:.4f} steps={STEPS}")
    print(f"psi   ={np.round(psi.real, 6) + 0j}")
    print(f"exact={np.round(exact.real, 6) + 0j}")
    print(f"L2_rel_vs_shift={err:.3e} L2_rel_vs_reverse_shift={rev_err:.3e}")
    print(f"direction={'forward' if err < rev_err else 'REVERSED(!)'}")


if __name__ == "__main__":
    main()
