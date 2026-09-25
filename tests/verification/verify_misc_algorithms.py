"""Publication-grade numerical validation of the misc_algorithms group.

Covers real-backend numerical experiments of six algorithm modules:
- qpca.density_matrix_exponentiation: LMR density-matrix exponentiation,
  trace distance against the exact numpy evolution e^{-i*rho*t}, verifying
  the first-order convergence scaling as copies double;
- qpca.qpca: phase-estimation readout of the eigenvalues of rho, against
  numpy eigh;
- density: gate_purification recovering rho by partial trace,
  gibbs_purification post-selected reduced state against the numpy Gibbs
  state e^{-beta*H}/Z (with an error convergence sweep);
- qsdp: trace_estimate_circuit probe readouts against the numpy trace (both
  the exact-purification and the Gibbs-approximate-purification paths), MMW
  driver convergence, and estimator agreement of a single quantum round;
- dqi: syndrome distribution against the Krawtchouk closed form, expected
  satisfied count against classical brute force, abstract-decoder binding
  consistency, and Dicke-state uniformity;
- variational: hardware-efficient ansatz states against a gate-by-gate numpy
  oracle, QAOA distributions against exact numpy simulation with optimal-cut
  probability compared to the random baseline, and VQE Pauli expectations
  and total energy against numpy;
- error_correction: logical amplitude recovery and syndrome values of the
  3-bit repetition code after a single injected error, plus the effective
  block of the encode-recover composite unitary.

All classical oracles are independent (numpy/math closed forms) and do not
reuse quantum-side helpers of the implementations under test.

Run: PYTHONPATH=src <python with pysparq+uniqc> tests/verification/verify_misc_algorithms.py
"""

from __future__ import annotations

import cmath
import math
from math import comb

import numpy as np
from harness import (
    Report,
    adapter_pysparq,
    amplitude_error,
    amplitudes_to_statevector,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    statevector_error,
    tvd,
)

from oracq import Bits, Builder, bind
from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding
from oracq.algorithms.input_model.density import gate_purification, gibbs_purification
from oracq.algorithms.input_model.oracles import StatePreparation, gate_state_prep
from oracq.algorithms.optimization.dqi import (
    XorSatInstance,
    abstract_decoder,
    bruteforce_decoder,
    dicke_state,
    dqi,
)
from oracq.algorithms.optimization.variational import (
    hardware_efficient_ansatz,
    qaoa_maxcut,
    vqe_measurements,
)
from oracq.algorithms.qec.error_correction import repetition_encode, repetition_recover
from oracq.algorithms.qml.qpca import (
    density_matrix_exponentiation,
    eigenvalue_from_phase,
    qpca,
)
from oracq.algorithms.qml.qsdp import (
    SdpInstance,
    classical_estimator,
    iteration_circuits,
    penalty_hamiltonian,
    qsdp_gibbs_solve,
    trace_estimate_circuit,
    trace_from_joint,
    trace_from_probe,
)
from oracq.infrastructure.layout import workspace_table

# ---------------------------------------------------------------------------
# Independent classical utilities (numpy dense linear algebra and closed forms)
# ---------------------------------------------------------------------------

SQRT2 = math.sqrt(2)
PLUS = [1 / SQRT2, 1 / SQRT2]
# Mixed state and Hermitian observables used by the validation
RHO_MIXED = ((0.75 + 0j, 0j), (0j, 0.25 + 0j))
RHO_GENERAL = ((0.7 + 0j, 0.1 - 0.05j), (0.1 + 0.05j, 0.3 + 0j))
PAULI_Z = ((1.0 + 0j, 0j), (0j, -1.0 + 0j))
PAULI_X = ((0j, 1.0 + 0j), (1.0 + 0j, 0j))
PAULI_Y = ((0j, -1.0j), (1.0j, 0j))
# DQI planted instance: all 7 non-zero rows, n = 3, right-hand side planted from x* = 0b101
PLANTED_ROWS = ((0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2))
PLANTED_RHS = (1, 0, 1, 1, 0, 1, 0)


def widths_of(program):
    return [r.type.width for r in program.main.registers]


def full_dense(amplitudes, widths):
    """Sparse amplitude dictionary -> dense state vector (registers occupy rising high bits in index order)."""
    vector = np.zeros(1 << sum(widths), dtype=complex)
    for key, amplitude in amplitudes.items():
        index = 0
        shift = 0
        for value, w in zip(key, widths, strict=True):
            index |= value << shift
            shift += w
        vector[index] += amplitude
    return vector


def partial_trace_np(vector, system_width, env_width):
    """Partial trace over the high-position env: index = system | (env << system_width).

    rho[i,j] = sum_e v[i + e*2^s]*conj(v[j + e*2^s]), i.e. arr.T @ conj(arr).
    """
    arr = np.asarray(vector, dtype=complex).reshape(1 << env_width, 1 << system_width)
    return arr.T @ arr.conj()


def trace_distance_np(rho, sigma):
    values = np.linalg.eigvalsh(np.asarray(rho, dtype=complex) - np.asarray(sigma, dtype=complex))
    return 0.5 * float(np.abs(values).sum())


def exact_evolved(rho, time, sigma):
    """e^{-i*rho*t} sigma e^{i*rho*t}: an independent exact numpy oracle via eigendecomposition."""
    rho = np.asarray(rho, dtype=complex)
    sigma = np.asarray(sigma, dtype=complex)
    values, vectors = np.linalg.eigh(rho)
    unitary = vectors @ np.diag(np.exp(-1j * time * values)) @ vectors.conj().T
    return unitary @ sigma @ unitary.conj().T


def gibbs_np(hamiltonian, beta):
    """Independent numpy Gibbs state e^{-beta*H}/Z."""
    values, vectors = np.linalg.eigh(np.asarray(hamiltonian, dtype=complex))
    weights = np.exp(-beta * values)
    return (vectors * weights) @ vectors.conj().T / weights.sum()


def ry_matrix(theta):
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def rz_matrix(phi):
    return np.diag([np.exp(-0.5j * phi), np.exp(0.5j * phi)]).astype(complex)


def apply_single(state, qubit, matrix):
    out = np.zeros_like(state)
    for i in range(len(state)):
        bit = (i >> qubit) & 1
        base = i & ~(1 << qubit)
        for b in (0, 1):
            out[base | (b << qubit)] += matrix[b, bit] * state[i]
    return out


def apply_cnot(state, control, target):
    out = np.zeros_like(state)
    for i in range(len(state)):
        out[i ^ (1 << target) if (i >> control) & 1 else i] += state[i]
    return out


def ansatz_numpy(width, layers):
    """Gate-by-gate independent numpy oracle for the hardware-efficient ansatz."""
    state = np.zeros(1 << width, dtype=complex)
    state[0] = 1.0
    for layer in layers:
        for qubit, (ry, rz) in enumerate(layer):
            state = apply_single(state, qubit, ry_matrix(ry))
            state = apply_single(state, qubit, rz_matrix(rz))
        for qubit in range(width - 1):
            state = apply_cnot(state, qubit, qubit + 1)
    return state


def qaoa_numpy(width, edges, gammas, betas):
    """Exact numpy simulation of QAOA: diagonal cost phase + per-qubit mixer e^{-i*beta*X}."""
    dim = 1 << width
    state = np.ones(dim, dtype=complex) / math.sqrt(dim)
    for gamma, beta in zip(gammas, betas, strict=True):
        phase = np.ones(dim, dtype=complex)
        for z in range(dim):
            cost = 0.0
            for u, v, w in edges:
                zu, zv = 1 - 2 * ((z >> u) & 1), 1 - 2 * ((z >> v) & 1)
                cost += w * (1 - zu * zv) / 2
            phase[z] = np.exp(-1j * gamma * cost)
        state *= phase
        rx = np.array(
            [[math.cos(beta), -1j * math.sin(beta)], [-1j * math.sin(beta), math.cos(beta)]]
        )
        for qubit in range(width):
            state = apply_single(state, qubit, rx)
    return state


def krawtchouk(m, weight, j):
    """Krawtchouk polynomial K_l(j) (math.comb closed form; the distribution oracle of the DQI paper)."""
    return sum(((-1) ** t) * comb(j, t) * comb(m - j, weight - t) for t in range(weight + 1))


def satisfied_count(rows, rhs, assignment):
    """Count satisfied rows directly from the row table / right-hand side (not reusing XorSatInstance methods)."""
    return sum(
        (sum((assignment >> j) & 1 for j in row) & 1) == bit
        for row, bit in zip(rows, rhs, strict=True)
    )


def phase_distribution(amplitudes):
    """qpca program readout: phase is the last register."""
    dist = {}
    for key, amplitude in amplitudes.items():
        dist[key[-1]] = dist.get(key[-1], 0.0) + abs(amplitude) ** 2
    return dist


def circular_eigenvalue(dist, precision, step):
    """Circular-mean decoding of the phase distribution (a lambda estimate robust to broadened peaks)."""
    z = sum(p * cmath.exp(2j * math.pi * v / (1 << precision)) for v, p in dist.items())
    return -2 * math.pi * (cmath.phase(z) / (2 * math.pi)) / step


def originir_fits(program, budget=24):
    """OriginIR state-vector budget precheck (entry registers + workspace)."""
    width = sum(r.type.width for r in program.main.registers)
    return width + workspace_table(program)[program.entry] <= budget


# ---------------------------------------------------------------------------
# Density-matrix exponentiation (LMR)
# ---------------------------------------------------------------------------


def verify_dm_exponentiation(report):
    sigma0 = np.array([[1, 0], [0, 0]], dtype=complex)  # system initial state |0><0|
    configs = (
        ("pure", gate_state_prep(PLUS), np.array([[0.5, 0.5], [0.5, 0.5]]), 0.4, None),
        (
            "mixed",
            gate_purification(RHO_MIXED).as_state_preparation(),
            np.diag([0.75, 0.25]),
            0.5,
            1,
        ),
    )
    for label, prep, rho, time, swap_width in configs:
        distances = []
        cross = 0.0
        copies_seen = []
        for copies in (1, 2, 4):
            kwargs = {} if swap_width is None else {"swap_width": swap_width}
            operation = density_matrix_exponentiation(
                prep, time=time, copies=copies, **kwargs
            )
            program = operation.program()
            widths = widths_of(program)
            env_width = sum(widths) - 1
            exact = exact_evolved(rho, time, sigma0)
            ref = reference(program)
            reduced = partial_trace_np(full_dense(ref, widths), 1, env_width)
            distances.append(trace_distance_np(reduced, exact))
            # Cross-backend: pysparq RIR / adapter vs the reference executor, amplitude by amplitude
            cross = max(cross, amplitude_error(rir_pysparq(program), ref))
            cross = max(cross, amplitude_error(adapter_pysparq(program), ref))
            if originir_fits(program):
                origin = originir_ext(program)
                cross = max(
                    cross,
                    statevector_error(list(origin), list(full_dense(ref, widths))),
                )
            copies_seen.append(copies)
        ratios = [distances[i + 1] / distances[i] for i in range(len(distances) - 1)]
        report.case(
            f"dm-exponentiation-{label}-convergence",
            paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
            parameters={
                "rho": label,
                "time": time,
                "copies": copies_seen,
                "system_qubits": 1,
                "max_total_qubits": 1 + 4 * (1 if label == "pure" else 2),
            },
            metrics={
                "trace_distance": [round(d, 10) for d in distances],
                "error_ratio": [round(r, 6) for r in ratios],
                "cross_backend_max_error": cross,
            },
            criterion=(
                "trace distance roughly halves as copies double (ratio < 0.65, LMR first-order scaling)"
                " and cross-backend amplitude deviation < 1e-9"
            ),
            passed=all(r < 0.65 for r in ratios)
            and distances[-1] > 0
            and cross < 1e-9,
        )


# ---------------------------------------------------------------------------
# QPCA eigenvalue readout
# ---------------------------------------------------------------------------


def verify_qpca(report):
    # Pure state rho = |+><+|: system input |+> (symmetric SWAP eigenstate), deterministic readout
    precision, step = 3, math.pi / 4
    operation = qpca(
        gate_state_prep(PLUS), precision=precision, step_time=step,
        system=gate_state_prep(PLUS),
    )
    program = operation.program()
    ref = reference(program)
    dist = phase_distribution(ref)
    mode = max(dist, key=dist.get)
    lam = eigenvalue_from_phase(mode, precision, step)
    # Cross-backend phase-distribution TVD
    widths = widths_of(program)
    tvd_cross = 0.0
    for runner in (rir_pysparq, adapter_pysparq):
        tvd_cross = max(tvd_cross, tvd(phase_distribution(runner(program)), dist))
    if originir_fits(program):
        origin = originir_ext(program)
        shift = sum(widths) - precision
        dist_o = {}
        for index, amplitude in enumerate(origin):
            dist_o[index >> shift] = dist_o.get(index >> shift, 0.0) + abs(amplitude) ** 2
        tvd_cross = max(tvd_cross, tvd(dist_o, dist))
    report.case(
        "qpca-pure-plus-deterministic",
        paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        parameters={
            "rho": "|+><+|",
            "precision": precision,
            "step_time": step,
            "total_qubits": sum(widths),
        },
        metrics={
            "peak_probability": dist[mode],
            "eigenvalue_error": abs(lam - 1.0),
            "cross_backend_tvd": tvd_cross,
        },
        criterion="deterministic readout (peak >= 1-1e-9), decoded lambda = 1 (error < 1e-12), cross-backend TVD < 1e-9",
        passed=dist[mode] >= 1 - 1e-9 and abs(lam - 1.0) < 1e-12 and tvd_cross < 1e-9,
    )

    # Mixed state rho = diag(0.75, 0.25): numpy eigh gives eigenvalues [0.25, 0.75]
    eigenvalues = np.linalg.eigvalsh(np.asarray(RHO_MIXED))
    prep = gate_purification(RHO_MIXED).as_state_preparation()
    step = 2 * math.pi / 6  # for lambda=0.75, phi=7/8 lands exactly on the 3-bit grid
    for state, target_lam, label in (([1, 0], 0.75, "primary"), ([0, 1], 0.25, "secondary")):
        assert abs(eigenvalues[1] - 0.75) < 1e-12 and abs(eigenvalues[0] - 0.25) < 1e-12
        operation = qpca(
            prep, precision=3, step_time=step,
            system=gate_state_prep(state), swap_width=1,
        )
        program = operation.program()
        widths = widths_of(program)
        ref = reference(program)
        dist = phase_distribution(ref)
        mode = max(dist, key=dist.get)
        lam_mode = eigenvalue_from_phase(mode, 3, step)
        lam_circ = circular_eigenvalue(dist, 3, step)
        tvd_cross = max(
            tvd(phase_distribution(rir_pysparq(program)), dist),
            tvd(phase_distribution(adapter_pysparq(program)), dist),
        )
        if originir_fits(program):
            origin = originir_ext(program)
            shift = sum(widths) - 3
            dist_o = {}
            for index, amplitude in enumerate(origin):
                dist_o[index >> shift] = dist_o.get(index >> shift, 0.0) + abs(amplitude) ** 2
            tvd_cross = max(tvd_cross, tvd(dist_o, dist))
        if label == "primary":
            # Peak (mode) decoding is exactly 0.75; peak height >= 0.4; circular-mean tolerance 0.15
            ok = (
                abs(lam_mode - target_lam) < 1e-12
                and dist[mode] >= 0.4
                and abs(lam_circ - target_lam) <= 0.15
                and tvd_cross < 1e-9
            )
            criterion = (
                "mode decoding is exactly lambda = 0.75 (error < 1e-12), peak height >= 0.4, "
                "circular-mean tolerance 0.15, cross-backend TVD < 1e-9"
            )
        else:
            # lambda=0.25 is not on the 3-bit grid for this dt: compare via the circular mean, report the mode as information
            ok = abs(lam_circ - target_lam) <= 0.15 and tvd_cross < 1e-9
            criterion = "circular-mean lambda approx 0.25 (tolerance 0.15; peak broadening is the LMR first-order error), cross-backend TVD < 1e-9"
        report.case(
            f"qpca-mixed-{label}",
            paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
            parameters={
                "rho": "diag(0.75, 0.25)",
                "system_input": f"|{label == 'secondary'}>",
                "precision": 3,
                "step_time": step,
                "total_qubits": sum(widths),
                "numpy_eigvalsh": [round(float(v), 12) for v in eigenvalues],
            },
            metrics={
                "mode_eigenvalue": lam_mode,
                "mode_eigenvalue_error": abs(lam_mode - target_lam),
                "peak_probability": dist[mode],
                "circular_mean_eigenvalue": lam_circ,
                "circular_mean_error": abs(lam_circ - target_lam),
                "cross_backend_tvd": tvd_cross,
            },
            criterion=criterion,
            passed=ok,
        )


# ---------------------------------------------------------------------------
# density: purification witnesses and Gibbs-state preparation
# ---------------------------------------------------------------------------


def verify_purification(report):
    access = gate_purification(RHO_GENERAL)
    program = access.operation.program()
    widths = widths_of(program)
    ref = reference(program)
    reduced = partial_trace_np(full_dense(ref, widths), widths[0], widths[1])
    expected = np.asarray(RHO_GENERAL, dtype=complex)
    max_error = float(np.abs(reduced - expected).max())
    cross = max(
        amplitude_error(rir_pysparq(program), ref),
        amplitude_error(adapter_pysparq(program), ref),
    )
    if originir_fits(program):
        cross = max(
            cross, statevector_error(list(originir_ext(program)), list(full_dense(ref, widths)))
        )
    report.case(
        "purification-partial-trace",
        paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        parameters={"rho": "2x2 complex Hermitian (off-diagonal 0.1+/-0.05j)", "qubits": sum(widths)},
        metrics={"max_element_error": max_error, "cross_backend_max_error": cross},
        criterion="partial trace recovers rho (matrix-element error < 1e-9) and cross-backend < 1e-9",
        passed=max_error < 1e-9 and cross < 1e-9,
    )


def verify_gibbs(report):
    for hamiltonian, beta, error, label, bound in (
        (((1.0 + 0j, 0j), (0j, -0.5 + 0j)), 0.6, 0.05, "diagonal", 0.05),
        (((0.3 + 0j, 0.2 - 0.1j), (0.2 + 0.1j, -0.4 + 0j)), 0.8, 0.05, "nondiagonal", 0.05),
    ):
        access = gibbs_purification(matrix_pauli_encoding(hamiltonian), beta, error=error)
        program = access.operation.program()
        widths = widths_of(program)
        ref = reference(program)
        full = full_dense(ref, widths)
        selector = (np.arange(len(full)) >> (widths[0] + widths[1])) == 0
        p_success = float((np.abs(full[selector]) ** 2).sum())
        kept = full[selector].reshape(1 << widths[1], 1 << widths[0])
        reduced = kept.T @ kept.conj()
        reduced = reduced / np.trace(reduced)
        distance = trace_distance_np(reduced, gibbs_np(hamiltonian, beta))
        cross = amplitude_error(rir_pysparq(program), ref)
        if originir_fits(program):
            cross = max(
                cross, statevector_error(list(originir_ext(program)), list(full))
            )
        report.case(
            f"gibbs-purification-{label}",
            paths=["reference", "rir-pysparq", "originir-ext"],
            parameters={
                "beta": beta,
                "error": error,
                "qsp_degree": access.attributes["qsp_degree"],
                "gibbs_scale": access.attributes["gibbs_scale"],
                "qubits": sum(widths),
            },
            metrics={
                "trace_distance": distance,
                "success_probability": p_success,
                "cross_backend_max_error": cross,
            },
            criterion=(
                f"trace distance between the post-selected reduced state and the numpy Gibbs state <= {bound}"
                " (error is the uniform-approximation upper-bound convention) and cross-backend < 1e-9"
            ),
            passed=distance <= bound and cross < 1e-9,
        )
    # error convergence sweep: same H/beta, distance <= error and monotonically non-increasing in error
    hamiltonian, beta = ((1.0 + 0j, 0j), (0j, -0.5 + 0j)), 0.8
    exact = gibbs_np(hamiltonian, beta)
    errors = (0.4, 0.2, 0.1)
    distances = []
    for error in errors:
        access = gibbs_purification(matrix_pauli_encoding(hamiltonian), beta, error=error)
        program = access.operation.program()
        widths = widths_of(program)
        full = full_dense(reference(program), widths)
        kept = full[(np.arange(len(full)) >> (widths[0] + widths[1])) == 0]
        kept = kept.reshape(1 << widths[1], 1 << widths[0])
        reduced = kept.T @ kept.conj()
        reduced = reduced / np.trace(reduced)
        distances.append(trace_distance_np(reduced, exact))
    report.case(
        "gibbs-error-scaling",
        paths=["reference"],
        parameters={"beta": beta, "errors": list(errors)},
        metrics={"trace_distance": [round(d, 10) for d in distances]},
        criterion="per-level trace distance <= error and monotonically non-increasing in error (polynomial truncation converging)",
        passed=all(d <= e for d, e in zip(distances, errors, strict=True))
        and distances[1] <= distances[0]
        and distances[2] <= distances[1],
    )


# ---------------------------------------------------------------------------
# qsdp: trace-estimation circuits and the MMW driver
# ---------------------------------------------------------------------------


def verify_trace_estimate(report):
    purification = gate_purification(RHO_GENERAL)
    rho = np.asarray(RHO_GENERAL, dtype=complex)
    worst, cross = 0.0, 0.0
    estimates = {}
    for word, matrix in (("Z", PAULI_Z), ("X", PAULI_X), ("Y", PAULI_Y)):
        be = matrix_pauli_encoding(matrix)
        circuit = trace_estimate_circuit(purification, be)
        program = circuit.program()
        names = [r.name for r in program.main.registers]
        probe = names.index("probe")
        ref = reference(program)
        p_one = sum(abs(a) ** 2 for k, a in ref.items() if k[probe] == 1)
        estimate = trace_from_probe(p_one, be.alpha)
        exact = float(np.trace(np.asarray(matrix) @ rho).real)
        worst = max(worst, abs(estimate - exact))
        estimates[word] = estimate
        cross = max(cross, amplitude_error(rir_pysparq(program), ref))
        if originir_fits(program):
            widths = widths_of(program)
            cross = max(
                cross,
                statevector_error(list(originir_ext(program)), list(full_dense(ref, widths))),
            )
    report.case(
        "trace-estimate-pauli-xyz",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"rho": "2x2 complex Hermitian", "observables": ["Z", "X", "Y"]},
        metrics={
            "max_error": worst,
            "estimates": {k: round(v, 12) for k, v in estimates.items()},
            "cross_backend_max_error": cross,
        },
        criterion="Tr(P rho) probe estimate matches the numpy trace (max_error < 1e-9)",
        passed=worst < 1e-9 and cross < 1e-9,
    )

    # Approximate-purification path: Gibbs purification + trace_from_joint joint decoding
    hamiltonian, beta = ((1.0 + 0j, 0j), (0j, -0.5 + 0j)), 0.6
    gibbs = gibbs_purification(matrix_pauli_encoding(hamiltonian), beta, error=0.05)
    be_z = matrix_pauli_encoding(PAULI_Z)
    circuit = trace_estimate_circuit(gibbs, be_z)
    program = circuit.program()
    names = [r.name for r in program.main.registers]
    ps, probe = names.index("purification_signal"), names.index("probe")
    ref = reference(program)
    e_joint = w_zero = 0.0
    for key, amplitude in ref.items():
        prob = abs(amplitude) ** 2
        if key[ps] == 0:
            w_zero += prob
            e_joint += prob * (1 if key[probe] == 0 else -1)
    estimate = trace_from_joint(e_joint, w_zero, be_z.alpha)
    exact = float(np.trace(np.asarray(PAULI_Z) @ gibbs_np(hamiltonian, beta)).real)
    cross = amplitude_error(rir_pysparq(program), ref)
    report.case(
        "trace-estimate-joint-gibbs",
        paths=["reference", "rir-pysparq"],
        parameters={
            "beta": beta,
            "gibbs_error": 0.05,
            "weight_signal_zero": w_zero,
            "qubits": sum(widths_of(program)),
        },
        metrics={
            "estimate": estimate,
            "exact": exact,
            "abs_error": abs(estimate - exact),
            "cross_backend_max_error": cross,
        },
        criterion="joint decoding of the approximate purification matches numpy for Tr(Z rho_gibbs) (error < 1e-2)",
        passed=abs(estimate - exact) < 1e-2 and cross < 1e-9,
    )


def verify_mmW_driver(report):
    neg_z = tuple(tuple(-v for v in row) for row in PAULI_Z)
    neg_x = tuple(tuple(-v for v in row) for row in PAULI_X)
    instance = SdpInstance(((PAULI_Z, 0.2), (neg_z, -0.2), (PAULI_X, 0.1), (neg_x, -0.1)))
    epsilon = 0.08
    result = qsdp_gibbs_solve(instance, epsilon=epsilon)
    rho = np.asarray(result["rho"])
    # Independent feasibility check: PSD, trace 1, constraint violations (computed directly with numpy)
    min_eig = float(np.linalg.eigvalsh(rho).min())
    trace = float(np.trace(rho).real)
    violations = [
        float(np.trace(np.asarray(a) @ rho).real) - b for a, b in instance.constraints
    ]
    bloch_x = float((2 * rho[0, 1]).real)
    bloch_z = float((rho[0, 0] - rho[1, 1]).real)
    report.case(
        "qsdp-mmw-driver-feasibility",
        paths=["classical-driver (estimator=classical_estimator)"],
        parameters={"constraints": 4, "epsilon": epsilon, "target_bloch": [0.1, 0.2]},
        metrics={
            "converged": result["converged"],
            "iterations": result["iterations"],
            "max_violation": max(violations),
            "min_eigenvalue": min_eig,
            "trace": trace,
            "bloch": [bloch_x, bloch_z],
        },
        criterion=(
            "converged and the average iterate rho-bar satisfies |Tr(A_i rho_bar) - b_i| <= epsilon (independently checked with numpy), "
            "PSD (min_eig >= -1e-9), trace 1, Bloch vector within 0.12 of the target"
        ),
        passed=result["converged"]
        and max(violations) <= epsilon
        and min_eig >= -1e-9
        and abs(trace - 1.0) < 1e-9
        and abs(bloch_x - 0.1) <= 0.12
        and abs(bloch_z - 0.2) <= 0.12,
    )

    # Single quantum round: the quantum estimates of iteration_circuits match the classical estimator
    inst = SdpInstance(((PAULI_Z, 0.2), (neg_z, -0.2)))
    weights = (0.6, 0.4)
    purification, traces = iteration_circuits(
        inst, weights, 1.0, hamiltonian_encoding=matrix_pauli_encoding, error=0.05
    )
    _, classical = classical_estimator(penalty_hamiltonian(inst, weights), 1.0, inst)
    diffs, cross = [], 0.0
    for circuit, classical_value in zip(traces, classical, strict=True):
        program = circuit.program()
        names = [r.name for r in program.main.registers]
        ps, probe = names.index("purification_signal"), names.index("probe")
        ref = reference(program)
        e_joint = w_zero = 0.0
        for key, amplitude in ref.items():
            prob = abs(amplitude) ** 2
            if key[ps] == 0:
                w_zero += prob
                e_joint += prob * (1 if key[probe] == 0 else -1)
        alpha = dict(circuit.module.attributes)["be_alpha"]
        diffs.append(abs(trace_from_joint(e_joint, w_zero, alpha) - classical_value))
        cross = max(cross, amplitude_error(rir_pysparq(program), ref))
    report.case(
        "qsdp-mmw-quantum-round",
        paths=["reference", "rir-pysparq"],
        parameters={"weights": list(weights), "beta": 1.0, "gibbs_error": 0.05},
        metrics={"max_estimate_diff": max(diffs), "cross_backend_max_error": cross},
        criterion="quantum trace estimates match the classical estimator constraint by constraint (max diff < 1e-2)",
        passed=max(diffs) < 1e-2 and cross < 1e-9,
    )


# ---------------------------------------------------------------------------
# DQI
# ---------------------------------------------------------------------------


def _dqi_distributions(program):
    """Returns (error distribution, syndrome distribution, cross-backend TVD)."""
    ref = reference(program)

    def split(amplitudes):
        error_dist, syndrome_dist = {}, {}
        for key, amplitude in amplitudes.items():
            prob = abs(amplitude) ** 2
            error_dist[key[0]] = error_dist.get(key[0], 0.0) + prob
            syndrome_dist[key[1]] = syndrome_dist.get(key[1], 0.0) + prob
        return error_dist, syndrome_dist

    error_dist, syndrome_dist = split(ref)
    cross = 0.0
    for runner in (rir_pysparq, adapter_pysparq):
        e2, s2 = split(runner(program))
        cross = max(cross, tvd(e2, error_dist), tvd(s2, syndrome_dist))
    if originir_fits(program):
        widths = widths_of(program)
        origin = originir_ext(program)
        e3, s3 = {}, {}
        for index, amplitude in enumerate(origin):
            prob = abs(amplitude) ** 2
            e_v = index & ((1 << widths[0]) - 1)
            s_v = index >> widths[0]
            e3[e_v] = e3.get(e_v, 0.0) + prob
            s3[s_v] = s3.get(s_v, 0.0) + prob
        cross = max(cross, tvd(e3, error_dist), tvd(s3, syndrome_dist))
    return error_dist, syndrome_dist, cross


def verify_dqi(report):
    instance = XorSatInstance(PLANTED_ROWS, PLANTED_RHS, 3)
    m, n = instance.num_constraints, instance.num_variables
    operation = dqi(instance, bruteforce_decoder(instance, max_weight=1), weight=1)
    error_dist, syndrome_dist, cross = _dqi_distributions(operation.program())
    # Krawtchouk closed-form expected distribution (the paper's oracle)
    weights = {}
    for x in range(1 << n):
        unsatisfied = m - satisfied_count(PLANTED_ROWS, PLANTED_RHS, x)
        weights[x] = krawtchouk(m, 1, unsatisfied) ** 2
    total = sum(weights.values())
    expected = {x: w / total for x, w in weights.items()}
    distribution_tvd = tvd(syndrome_dist, expected)
    max_point = max(abs(syndrome_dist.get(x, 0) - p) for x, p in expected.items())
    expected_satisfied = sum(
        p * satisfied_count(PLANTED_ROWS, PLANTED_RHS, x) for x, p in syndrome_dist.items()
    )
    # Classical brute force: best assignment and satisfied count
    brute = {
        x: satisfied_count(PLANTED_ROWS, PLANTED_RHS, x) for x in range(1 << n)
    }
    best = max(brute, key=brute.get)
    quantum_argmax = max(syndrome_dist, key=syndrome_dist.get)
    report.case(
        "dqi-planted-krawtchouk",
        paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        parameters={
            "num_constraints": m,
            "num_variables": n,
            "dicke_weight": 1,
            "planted_assignment": best,
            "random_baseline": m / 2,
        },
        metrics={
            "distribution_tvd": distribution_tvd,
            "max_pointwise_error": max_point,
            "p_error_register_zero": error_dist.get(0, 0.0),
            "expected_satisfied": expected_satisfied,
            "quantum_argmax": quantum_argmax,
            "bruteforce_best_satisfied": brute[best],
            "cross_backend_tvd": cross,
        },
        criterion=(
            "syndrome distribution vs the Krawtchouk closed form TVD < 1e-9, error register deterministically clean, "
            "expected satisfied count = 6.5 (> random baseline 3.5), probability peak is the brute-force best assignment"
        ),
        passed=distribution_tvd < 1e-9
        and max_point < 1e-9
        and abs(error_dist.get(0, 0.0) - 1.0) < 1e-12
        and abs(expected_satisfied - 6.5) < 1e-9
        and quantum_argmax == best
        and cross < 1e-9,
    )

    # Abstract-decoder binding path matches the direct witness
    decoder = abstract_decoder("DqiSyndromeDecoder", n, m)
    program = dqi(instance, decoder, weight=1).program()
    bound = bind(
        program, {"DqiSyndromeDecoder": bruteforce_decoder(instance, max_weight=1).operation}
    )
    deviation = amplitude_error(rir_pysparq(bound), rir_pysparq(operation.program()))
    report.case(
        "dqi-abstract-decoder-bind",
        paths=["rir-pysparq"],
        parameters={"slot": "DqiSyndromeDecoder"},
        metrics={"bind_vs_witness_max_error": deviation},
        criterion="the bound program matches the direct witness amplitude by amplitude (< 1e-12)",
        passed=deviation < 1e-12,
    )

    # Dicke state |D_2^5>: C(5,2)=10 equal-amplitude components
    prep = dicke_state(5, 2)
    program = prep.operation.program()
    ref = reference(program)
    target_amp = 1 / math.sqrt(comb(5, 2))
    support = [k for k, a in ref.items() if abs(a) > 1e-12]
    amp_error = max(abs(abs(a) - target_amp) for a in ref.values())
    weight_ok = all(k[0].bit_count() == 2 and k[1] == 0 for k in support)
    cross = 0.0
    if originir_fits(program):
        widths = widths_of(program)
        cross = statevector_error(
            list(originir_ext(program)), amplitudes_to_statevector(ref, widths)
        )
    report.case(
        "dicke-state-uniformity",
        paths=["reference", "originir-ext"],
        parameters={"m": 5, "weight": 2},
        metrics={
            "support_size": len(support),
            "amplitude_error": amp_error,
            "cross_backend_max_error": cross,
        },
        criterion="support is exactly the 10 weight-2 basis states, equal amplitude 1/sqrt(10) (error < 1e-12)",
        passed=len(support) == 10 and amp_error < 1e-12 and weight_ok and cross < 1e-12,
    )


# ---------------------------------------------------------------------------
# variational: ansatz, QAOA, VQE measurements
# ---------------------------------------------------------------------------


def verify_ansatz(report):
    width = 3
    layers = (
        ((0.7, 0.3), (-0.4, 0.9), (1.1, -0.2)),
        ((0.2, -0.6), (0.5, 0.8), (-0.9, 0.1)),
    )
    program = hardware_efficient_ansatz(width, layers).program()
    expected = ansatz_numpy(width, layers)
    ref = reference(program)
    actual = np.zeros(1 << width, dtype=complex)
    for key, amplitude in ref.items():
        actual[key[0]] = amplitude
    max_error = float(np.abs(actual - expected).max())
    fidelity = float(abs(np.vdot(actual, expected)) ** 2)
    cross = max(
        amplitude_error(rir_pysparq(program), ref),
        amplitude_error(adapter_pysparq(program), ref),
    )
    if originir_fits(program):
        cross = max(cross, statevector_error(list(originir_ext(program)), list(expected)))
    report.case(
        "ansatz-parameter-intent",
        paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        parameters={"width": width, "layers": 2, "angles_per_qubit": 2},
        metrics={
            "max_amplitude_error": max_error,
            "fidelity": fidelity,
            "cross_backend_max_error": cross,
        },
        criterion="prepared state matches the gate-by-gate numpy oracle (max_error < 1e-12, fidelity >= 1-1e-12)",
        passed=max_error < 1e-12 and fidelity >= 1 - 1e-12 and cross < 1e-12,
    )


def verify_qaoa(report):
    # Single-edge analytic anchor: gamma=pi/2, beta=pi/8 reaches the optimal cut in one layer
    program = qaoa_maxcut(2, ((0, 1, 1.0),), (math.pi / 2,), (math.pi / 8,)).program()
    ref = reference(program)
    p_optimal = sum(abs(a) ** 2 for k, a in ref.items() if k[0] in (1, 2))
    report.case(
        "qaoa-single-edge-optimal",
        paths=["reference"],
        parameters={"vertices": 2, "gamma": math.pi / 2, "beta": math.pi / 8},
        metrics={"p_optimal_cut": p_optimal},
        criterion="optimal-cut probability is exactly 1 (analytic result, tolerance 1e-9)",
        passed=abs(p_optimal - 1.0) < 1e-9,
    )

    # 4-vertex cycle C4: angles chosen by grid search over the numpy oracle, quantum distribution compared pointwise
    edges = ((0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0))
    best = max(
        (
            (
                float(np.abs(qaoa_numpy(4, edges, (g,), (b,)))[0b0101] ** 2)
                + float(np.abs(qaoa_numpy(4, edges, (g,), (b,)))[0b1010] ** 2),
                g,
                b,
            )
            for g in np.linspace(0.05, 1.5, 30)
            for b in np.linspace(0.05, 1.5, 30)
        ),
        key=lambda item: item[0],
    )
    _, gamma, beta = best
    expected = np.abs(qaoa_numpy(4, edges, (gamma,), (beta,))) ** 2
    program = qaoa_maxcut(4, edges, (float(gamma),), (float(beta),)).program()
    ref = reference(program)
    measured = np.zeros(16)
    for key, amplitude in ref.items():
        measured[key[0]] = abs(amplitude) ** 2
    distribution_tvd = 0.5 * float(np.abs(measured - expected).sum())
    p_optimal = float(measured[0b0101] + measured[0b1010])
    baseline = 2 / 16
    cross = 0.0
    if originir_fits(program):
        origin = originir_ext(program)
        cross = 0.5 * sum(
            abs(abs(origin[z]) ** 2 - expected[z]) for z in range(16)
        )
    report.case(
        "qaoa-c4-distribution",
        paths=["reference", "originir-ext"],
        parameters={
            "graph": "C4 (4-vertex cycle, unit weights)",
            "gamma": float(gamma),
            "beta": float(beta),
            "layers": 1,
            "random_baseline_p_optimal": baseline,
        },
        metrics={
            "distribution_tvd": distribution_tvd,
            "p_optimal_cut": p_optimal,
            "advantage_over_baseline": p_optimal / baseline,
            "originir_tvd": cross,
        },
        criterion="distribution vs exact numpy simulation TVD < 1e-12 and optimal-cut probability >= 4x the random baseline",
        passed=distribution_tvd < 1e-12 and p_optimal >= 4 * baseline and cross < 1e-12,
    )


def verify_vqe(report):
    ansatz = hardware_efficient_ansatz(2, (((0.9, 0.4), (-0.7, 1.2)),))
    preparation = StatePreparation.from_unitary(ansatz)
    terms = ((0.5, "ZI"), (-0.3, "IZ"), (0.7, "XX"), (0.2, "YY"), (-0.1, "ZZ"))
    circuits = vqe_measurements(preparation, terms)
    # Independent numpy expectations: compute the prepared state first, then <psi|P|psi> per Pauli-word tensor
    psi = ansatz_numpy(2, (((0.9, 0.4), (-0.7, 1.2)),))
    paulis = {
        "I": np.eye(2),
        "X": np.array([[0, 1], [1, 0]], dtype=complex),
        "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
        "Z": np.diag([1, -1]).astype(complex),
    }
    worst, energy_measured, energy_exact = 0.0, 0.0, 0.0
    for (coefficient, word), (_, circuit) in zip(terms, circuits, strict=True):
        ref = reference(circuit.program())
        bits = [i for i, letter in enumerate(word) if letter != "I"]
        expectation = 0.0
        for key, amplitude in ref.items():
            parity = (key[0] & sum(1 << b for b in bits)).bit_count() % 2
            expectation += abs(amplitude) ** 2 * (1 - 2 * parity)
        operator = np.kron(paulis[word[1]], paulis[word[0]])  # the word's first character is the least significant bit
        exact = float((psi.conj() @ operator @ psi).real)
        worst = max(worst, abs(expectation - exact))
        energy_measured += coefficient * expectation
        energy_exact += coefficient * exact
    report.case(
        "vqe-pauli-expectations",
        paths=["reference"],
        parameters={"terms": [w for _, w in terms], "ansatz_width": 2},
        metrics={
            "max_term_error": worst,
            "energy_measured": energy_measured,
            "energy_exact": energy_exact,
            "energy_error": abs(energy_measured - energy_exact),
        },
        criterion="each Pauli expectation matches numpy (< 1e-9) and the weighted total energy agrees (< 1e-9)",
        passed=worst < 1e-9 and abs(energy_measured - energy_exact) < 1e-9,
    )


# ---------------------------------------------------------------------------
# error_correction: the 3-bit repetition code
# ---------------------------------------------------------------------------


def _repetition_program(error_kind, position, logical):
    encode = repetition_encode(error=error_kind)
    recover = repetition_recover(error=error_kind)
    b = Builder(f"rep_{error_kind}_{position}", {"target": Bits(1), "syndrome": Bits(2)})
    ry, rz = logical
    b.ry(b["target"][0], ry)
    b.rz(b["target"][0], rz)
    b.call(encode, target=b["target"], syndrome=b["syndrome"])
    qubits = {"target": b["target"][0], "s0": b["syndrome"][0], "s1": b["syndrome"][1]}
    if position != "none":
        if error_kind == "bit":
            b.x(qubits[position])
        else:
            b.z(qubits[position])
    b.call(recover, target=b["target"], syndrome=b["syndrome"])
    return b.finish().program()


def verify_repetition(report):
    logical = (0.73, 0.29)
    # Exact amplitudes of an arbitrary logical state (independent numpy oracle)
    alpha = math.cos(logical[0] / 2) * np.exp(-0.5j * logical[1])
    beta = math.sin(logical[0] / 2) * np.exp(0.5j * logical[1])
    expected_syndrome = {"none": 0, "target": 3, "s0": 1, "s1": 2}
    for kind in ("bit", "phase"):
        amp_worst, cross = 0.0, 0.0
        syndrome_ok = True
        for position in ("none", "target", "s0", "s1"):
            program = _repetition_program(kind, position, logical)
            ref = reference(program)
            syndromes = set()
            for key, amplitude in ref.items():
                syndromes.add(key[1])
                amp_worst = max(
                    amp_worst, abs(amplitude - (alpha if key[0] == 0 else beta))
                )
            syndrome_ok = syndrome_ok and syndromes == {expected_syndrome[position]}
            cross = max(
                cross,
                amplitude_error(rir_pysparq(program), ref),
                amplitude_error(adapter_pysparq(program), ref),
            )
            if originir_fits(program):
                cross = max(
                    cross,
                    statevector_error(
                        list(originir_ext(program)), amplitudes_to_statevector(ref, [1, 2])
                    ),
                )
        report.case(
            f"repetition-{kind}-flip-injection",
            paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
            parameters={
                "logical_angles": list(logical),
                "positions": ["none", "target", "s0", "s1"],
                "expected_syndrome": expected_syndrome,
            },
            metrics={
                "max_amplitude_error": amp_worst,
                "syndrome_exact": syndrome_ok,
                "cross_backend_max_error": cross,
            },
            criterion=(
                "after a single injected error the logical amplitudes are recovered pointwise (< 1e-12) and the syndrome deterministically identifies the error position"
            ),
            passed=amp_worst < 1e-12 and syndrome_ok and cross < 1e-12,
        )
    # Error-free composite unitary: the effective block on the syndrome=0 subspace is the identity, with no leakage
    for kind in ("bit", "phase"):
        encode = repetition_encode(error=kind)
        recover = repetition_recover(error=kind)
        b = Builder(f"rep_unitary_{kind}", {"target": Bits(1), "syndrome": Bits(2)})
        b.call(encode, target=b["target"], syndrome=b["syndrome"])
        b.call(recover, target=b["target"], syndrome=b["syndrome"])
        unitary = originir_unitary(b.finish().program())
        block = unitary[:2, :2]
        leakage = float(np.abs(unitary[2:, :2]).max())
        block_error = float(np.abs(block - np.eye(2)).max())
        report.case(
            f"repetition-{kind}-encode-recover-unitary",
            paths=["originir-ext+to_matrix"],
            parameters={"error_kind": kind},
            metrics={"block_error": block_error, "leakage": leakage},
            criterion="the encode-recover composite is exactly the identity on the syndrome=0 input block with no leakage (< 1e-12)",
            passed=block_error < 1e-12 and leakage < 1e-12,
        )


def run():
    report = Report(
        "misc_algorithms",
        "Real-backend numerical validation of LMR exponentiation/QPCA, purification and Gibbs preparation, QSDP trace estimation and MMW, DQI, "
        "variational circuits, and the repetition code (classical oracles are independent numpy/closed forms).",
    )
    verify_dm_exponentiation(report)
    verify_qpca(report)
    verify_purification(report)
    verify_gibbs(report)
    verify_trace_estimate(report)
    verify_mmW_driver(report)
    verify_dqi(report)
    verify_ansatz(report)
    verify_qaoa(report)
    verify_vqe(report)
    verify_repetition(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
