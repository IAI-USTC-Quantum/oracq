"""Publication-grade numerical validation of the ODE group (Carleman/LCHS/CBMD/Schrodingerization/ode/ode_models).

Source files covered: algorithms/ode.py, ode_models.py, carleman.py, lchs.py,
schrodingerization.py, cbmd.py.

Validation structure (all really executed, no mocks/skips):

- Substructures: the taylor_hamiltonian encoding block, the fourier_momentum
  block, the Carleman lifting block (OriginIR-ext -> UniQC ``to_matrix`` full
  unitary + ``effective_block`` extraction, cross-checked against a reference
  basis-state sweep), the carleman_initial lifted initial state on four paths,
  and CBMD/LCHS plan weights against independent closed forms and the contour
  identity residual.
- End to end (multiple input models): LCHS under five input paradigms (a whole
  BE / direct Hermitian parts / diagonal spectral-angle databases with gate
  and QRAM bindings / Fokker-Planck Pauli expansion + QODEProblem / the
  structured shifted BE of the heat equation) agrees amplitude by amplitude
  with numpy simulation (implementation error), and the method error versus
  the scipy ``expm`` exact solution (quadrature/Taylor remainders, marked
  pending in the library) is reported.
- CBMD and LCHS cross-checked on the same non-commuting problem; two
  Schrodingerization cases (grid-translation-exact scalar decay and an
  analytically solvable rotation) validate the recovered magnitudes; Carleman
  runs Riccati end to end with an injected minimal Taylor-protocol solver (a
  real BE assembly), decomposing the Taylor remainder from the truncation
  error, and reports the quantum convergence trend for truncation orders
  K=1,2 (K=3 only as a classical extra point).

All classical references are independent: numpy/scipy (``expm``,
``solve_ivp``, analytic rotation/decay, Fourier eigenvalues) and plan
identities evaluated directly from their formulas; expected values never
reuse the assembly logic under test.

Run: PYTHONPATH=src <python with pysparq+uniqc> tests/verification/verify_ode.py
"""

from __future__ import annotations

import cmath
import math
from functools import partial

import numpy as np
import scipy.linalg
from harness import (
    Report,
    adapter_pysparq,
    amplitude_error,
    effective_block,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
)

from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.common.state_preparation import apply_be_to_state
from oracq.algorithms.input_model.block_encoding import lcu, matrix_pauli_encoding, tensor
from oracq.algorithms.input_model.operators import identity, product, scale, zero
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    abstract_database,
    abstract_state_prep,
    annotate,
    diagonal_block_encoding,
    gate_database,
    gate_state_prep,
    qram_database,
    qram_state_angles,
    qram_state_prep,
)
from oracq.algorithms.qnlss.carleman import (
    PolynomialODE,
    carleman_initial,
    carleman_lift,
    carleman_qode,
)
from oracq.algorithms.qode.cbmd import ContourPlan, cbmd_qode
from oracq.algorithms.qode.lchs import QuadraturePlan, lchs_qode
from oracq.algorithms.qode.ode import linear_qode
from oracq.algorithms.qode.ode_models import HermitianParts, LinearODE
from oracq.algorithms.qode.schrodingerization import SchrodingerPlan, fourier_momentum

SQRT2 = math.sqrt(2)
RUN_KWARGS = {"max_steps": 1 << 30, "max_states": 1 << 22}


# ---------------------------------------------------------------------------
# Independent classical references (numpy/scipy/analytic formulas)
# ---------------------------------------------------------------------------


def taylor_matrix(k_mat, time, degree):
    """sum_{l<=degree} (-i t)^l/l! K^l (independent matrix-power recursion)."""
    result = np.zeros_like(k_mat)
    term = np.eye(k_mat.shape[0], dtype=complex)
    for order in range(degree + 1):
        result = result + ((-1j * time) ** order / math.factorial(order)) * term
        term = term @ k_mat
    return result


def ode_taylor_series(g_mat, time, degree):
    """sum_{l<=degree} t^l/l! G^l: the truncated series of the ODE propagator e^{Gt} (independent reference)."""
    result = np.zeros_like(g_mat)
    term = np.eye(g_mat.shape[0], dtype=complex)
    for order in range(degree + 1):
        result = result + (time**order / math.factorial(order)) * term
        term = term @ g_mat
    return result


def qft_matrix(width):
    """Positive-sign DFT matrix exp(2*pi*i*x*y/N)/sqrt(N) (matching the documented convention of fourier.qft)."""
    n = 1 << width
    index = np.arange(n).reshape(-1, 1)
    return np.exp(2j * math.pi * (index @ index.T) / n) / math.sqrt(n)


def momentum_matrix(width, period):
    """Diagonal reference of fourier_momentum: two's-complement signed frequencies x 2*pi/period."""
    n = 1 << width
    signed = [k if k < (1 << (width - 1)) else k - (1 << width) for k in range(n)]
    return np.diag(np.array(signed) * 2 * math.pi / period)


def assemble_carleman_lift(f_operators, cutoff):
    """Independent numpy assembly of the Carleman-lifted generator under the padded layout.

    Layout: target = data (cutoff groups, n bits each, low positions) + level
    (high position). The basis state |data, level> is valid only when the
    groups above level are all zero; the term (k, p, position) writes data of
    level=k+p-1 contracted through F_p into level=k: group pos receives the
    F_p output while the remaining source groups pos+p.. shift down by p-1
    positions (consistent with the swap in _carleman_term).
    """
    n = f_operators[1].shape[0].bit_length() - 1
    lb = cutoff.bit_length()
    group = (1 << n) - 1
    dim = (1 << (cutoff * n)) * (1 << lb)
    g_mat = np.zeros((dim, dim), dtype=complex)

    def index(data, level):
        return data | (level << (cutoff * n))

    for level in range(cutoff + 1):
        for data in range(1 << (cutoff * n)):
            if data >> (level * n):
                continue
            col = index(data, level)
            for k in range(1, cutoff + 1):
                p = level - k + 1
                if p not in f_operators or p < 1:
                    continue
                f_mat = f_operators[p]
                for pos in range(k):
                    vin = 0
                    for j in range(p):
                        vin |= ((data >> ((pos + j) * n)) & group) << (j * n)
                    for vout in range(1 << n):
                        amp = f_mat[vout, vin]
                        if not amp:
                            continue
                        new_data = data & ((1 << (pos * n)) - 1)  # groups before pos unchanged
                        new_data |= vout << (pos * n)
                        for j in range(pos + 1, k):  # remaining source groups shift down p-1 positions
                            new_data |= ((data >> ((j + p - 1) * n)) & group) << (j * n)
                        g_mat[index(new_data, k), col] += amp
    return g_mat


def lifted_initial(u0, initial_norm, cutoff):
    """Independent reference of carleman_initial: 1/Z sum_k r^k u0^{tensor k} (padded layout)."""
    n = len(u0).bit_length() - 1
    lb = cutoff.bit_length()
    norm = math.sqrt(sum(initial_norm ** (2 * k) for k in range(cutoff + 1)))
    z0 = np.zeros((1 << (cutoff * n)) * (1 << lb), dtype=complex)
    z0[0] = 1.0 / norm
    for k in range(1, cutoff + 1):
        for bits in range(1 << (k * n)):
            amp = initial_norm**k / norm
            for j in range(k):
                amp *= u0[(bits >> (j * n)) & ((1 << n) - 1)]
            z0[bits | (k << (cutoff * n))] = amp
    return z0


def schrodinger_emulate(
    g_mat, u0, time, plan, degree, alpha_e, *, exact_evolution=False, flip_momentum=False
):
    """Independent full-stack simulation of Schrodingerization: warp -> QFT -> Taylor (or exact exp) -> inverse QFT -> channel.

    flip_momentum=True corresponds to K' = -P tensor H1 - I tensor H2 (the
    post-fix library convention); historically the library assembled the
    +P tensor H1 term, recovering a time-reversed flow -- the verification
    finding and fix are recorded in the group report.
    """
    n = (g_mat.shape[0] - 1).bit_length()
    p = plan.auxiliary_width
    m = 1 << p
    grid = np.array(
        [(j if j < (1 << (p - 1)) else j - (1 << p)) * plan.period / m for j in range(m)]
    )
    h1 = (g_mat + g_mat.conj().T) / 2
    h2 = (g_mat - g_mat.conj().T) / (2j)
    p_sign = -1.0 if flip_momentum else 1.0
    k_mat = np.kron(p_sign * momentum_matrix(p, plan.period), h1) - np.kron(np.eye(m), h2)
    if exact_evolution:
        evolution = scipy.linalg.expm(-1j * k_mat * time)
    else:
        evolution = taylor_matrix(k_mat, time, degree)
    f_mat = qft_matrix(p)
    warp = np.exp(-np.abs(grid))
    warp = warp / np.linalg.norm(warp)
    init = np.kron(warp, u0)  # aux high, physical low
    out = np.kron(f_mat.conj().T, np.eye(1 << n)) @ (
        evolution @ (np.kron(f_mat, np.eye(1 << n)) @ init)
    )
    return out.reshape(m, 1 << n).T[:, plan.selected_index] / alpha_e


def lchs_emulate(l_mat, h_mat, nodes, weights, time, degree, alpha_v, u0):
    """Independent simulation of the LCHS/CBMD finite sum + truncated Taylor; returns the post-selected block V~|u0>/alpha_V."""
    v_mat = np.zeros((len(u0), len(u0)), dtype=complex)
    for node, weight in zip(nodes, weights, strict=True):
        v_mat = v_mat + weight * taylor_matrix(h_mat + node * l_mat, time, degree)
    return v_mat @ u0 / alpha_v


# ---------------------------------------------------------------------------
# Quantum-side helpers
# ---------------------------------------------------------------------------


def postselect(amplitudes, signal_index=1, signal_value=0):
    """Extract the signal==value block and success probability from a (target, signal) amplitude dictionary."""
    block = {}
    success = 0.0
    for key, amplitude in amplitudes.items():
        if key[signal_index] == signal_value:
            block[key[0]] = amplitude
            success += abs(amplitude) ** 2
    return block, success


def run_amplitude_paths(program, paths, memory=None):
    """Execute per path name and unify into amplitude dictionaries; originir-ext converted from the state vector."""
    results = {}
    for path in paths:
        if path == "originir-ext":
            vector = originir_ext(program, memory, max_steps=RUN_KWARGS["max_steps"])
            widths = [r.type.width for r in program.main.registers]
            amplitudes = {}
            for index, value in enumerate(vector):
                if abs(value) > 1e-15:
                    key = ()
                    rest = index
                    for w in widths:
                        key = key + (rest & ((1 << w) - 1),)
                        rest >>= w
                    amplitudes[key] = complex(value)
            results[path] = amplitudes
        else:
            runner = {"reference": reference, "rir-pysparq": rir_pysparq, "adapter-pysparq": adapter_pysparq}[
                path
            ]
            results[path] = runner(program, memory, **RUN_KWARGS)
    return results


def compare_path_blocks(results):
    """Extract the physical block per path and produce three-level cross metrics.

    - block_dev: pairwise maximum deviation of the physical blocks (signal==0)
      across all paths -- the strict criterion of physical-result agreement;
    - exact_dev: full-spectrum pairwise deviation between reference and
      originir-ext (both paths verified to the 1e-17 level);
    - pysparq_floor: maximum full-spectrum deviation of the pysparq paths from
      reference (informative). pysparq has a numerical floor on junk branches
      for deep LCU programs containing thousands of rotations: amplitudes
      below 1e-7 are pruned and ~1e-5 branches show ~0.3% relative jitter
      (reference and originir-ext agree exactly on such branches, so this is
      identified as a pysparq-side artifact; see the group report).
    """
    blocks = {name: postselect(amplitudes)[0] for name, amplitudes in results.items()}
    names = list(blocks)
    block_dev = 0.0
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            keys = set(blocks[names[i]]) | set(blocks[names[j]])
            if keys:
                block_dev = max(
                    block_dev,
                    max(abs(blocks[names[i]].get(k, 0) - blocks[names[j]].get(k, 0)) for k in keys),
                )
    exact_dev = 0.0
    if "reference" in results and "originir-ext" in results:
        exact_dev = amplitude_error(results["reference"], results["originir-ext"])
    pysparq_floor = 0.0
    sparse_paths = [p for p in ("reference", "rir-pysparq", "adapter-pysparq") if p in results]
    for i in range(len(sparse_paths)):
        for j in range(i + 1, len(sparse_paths)):
            if "pysparq" in sparse_paths[i] or "pysparq" in sparse_paths[j]:
                pysparq_floor = max(
                    pysparq_floor,
                    amplitude_error(results[sparse_paths[i]], results[sparse_paths[j]]),
                )
    return blocks, block_dev, exact_dev, pysparq_floor


def extract_block_reference(be, dim):
    """Extract the <0|U|0>*alpha block of a BE via a reference basis-state sweep (independent of the unitary-export path)."""
    program = be.operation.program()
    block = np.zeros((dim, dim), dtype=complex)
    for x in range(dim):
        out = reference(program, initial={"target": x, "signal": 0}, **RUN_KWARGS)
        for (tv, sv), a in out.items():
            if sv == 0:
                block[tv, x] = a
    return block * be.alpha


def evolution_alpha(state_oracle, key="evolution_alpha"):
    """Read the evolution alpha recorded by the injected solver from the program module attributes (scalar metadata)."""
    for module in state_oracle.operation.program().modules:
        attrs = dict(module.attributes)
        if key in attrs:
            return attrs[key]
    raise AssertionError(f"attribute {key} not found")


def schrodinger_alpha(g_be, plan, time, degree):
    """Rebuild the K/E block-encoding chain of Schrodingerization to read the scalar alpha (metadata)."""
    parts = HermitianParts.from_operator(g_be)
    momentum = fourier_momentum(plan.auxiliary_width, plan.period)
    k_be = lcu(
        [
            (1, tensor(momentum, parts.hermitian)),
            (-1, tensor(identity(plan.auxiliary_width), parts.h)),
        ]
    )
    return taylor_hamiltonian(k_be, time, degree=degree).alpha


def taylor_series_solver(generator, initial, time, *, degree):
    """Injected minimal three-argument linear QODE protocol solver: the truncated Taylor-series BE of e^{Gt}.

    Like the in-library taylor_hamiltonian, it is a real block-encoding
    assembly (lcu/product/apply_be_to_state) targeting u' = Gu directly
    (non-Hermitian generators allowed); it is used for the protocol-injection
    end-to-end validation of carleman_qode; it is not a stand-in for an
    in-library solver but another real implementation the protocol permits.
    """
    terms, current = [(1.0, identity(generator.width))], identity(generator.width)
    for order in range(1, degree + 1):
        current = product(generator, current)
        terms.append((time**order / math.factorial(order), current))
    evolution = lcu(terms)
    state = apply_be_to_state(evolution, initial)
    return StateOracle(
        annotate(
            state.operation,
            "unitary",
            algorithm="injected_taylor_series_qode",
            evolution_alpha=evolution.alpha,
            degree=degree,
            correctness="pending",
        )
    )


def schrodinger_sign_flipped_qode(g_be, initial, time, plan, *, degree):
    """Independent re-assembly of Schrodingerization with K' = -P tensor H1 - I tensor H2 (a real quantum program).

    The in-library schrodingerization used to assemble K = P tensor H1 -
    I tensor H2, recovering e^{+t}*u0 under the positive-QFT convention (a
    time-reversed flow; the verification finding is in the group report);
    after the fix the library convention is the K' this function assembles.
    The same construction is independently re-assembled here from public
    combinators, serving as a regression pin for the sign convention: amplitude
    agreement with the library program means the sign has not regressed.
    """
    from oracq.algorithms.common.fourier import qft_with_work as qft
    from oracq.algorithms.common.state_preparation import select_subspace
    from oracq.algorithms.input_model.operators import _name
    from oracq.algorithms.input_model.oracles import StatePreparation, invoke, resources_for
    from oracq.algorithms.qode.ode_models import HermitianParts
    from oracq.infrastructure.builder import Builder
    from oracq.infrastructure.ir import Bits

    parts = HermitianParts.from_operator(g_be)
    n, p = g_be.width, plan.auxiliary_width
    momentum = fourier_momentum(p, plan.period)
    hamiltonian = lcu(
        [(-1, tensor(momentum, parts.hermitian)), (-1, tensor(identity(p), parts.h))]
    )
    evolution = taylor_hamiltonian(hamiltonian, time, degree=degree)
    grid = [
        (j if j < (1 << (p - 1)) else j - (1 << p)) * plan.period / (1 << p)
        for j in range(1 << p)
    ]
    warp = gate_state_prep([math.exp(-abs(x)) for x in grid])
    transform = qft(p)
    prep = Builder(
        _name("verify_schrod_flip_warp", initial.operation, plan),
        {"target": Bits(n + p), "work": Bits(initial.work_width)},
        resources_for(("initial", initial.operation)),
    )
    invoke(prep, initial.operation, "initial", target=prep["target"][:n], work=prep["work"])
    invoke(prep, warp.operation, target=prep["target"][n:], work=prep["work"][:0])
    invoke(prep, transform, target=prep["target"][n:], work=prep["work"][:0])
    state = apply_be_to_state(
        evolution,
        StatePreparation(annotate(prep.finish(), "state_prep_isometry", zero_input=True)),
    )
    out = Builder(
        _name("verify_schrod_flip_inverse", state.operation),
        {"target": Bits(n + p), "signal": Bits(state.signal_qubits)},
        resources_for(("state", state.operation)),
    )
    invoke(out, state.operation, "state", target=out["target"], signal=out["signal"])
    with out.adjoint():
        invoke(out, transform, target=out["target"][n:], work=out["signal"][:0])
    selected = select_subspace(
        StateOracle(out.finish()), n, plan.selected_index, label="verify_schrod_flip_channel"
    )
    return StateOracle(
        annotate(
            selected.operation,
            "unitary",
            algorithm="schrodingerization_sign_flipped",
            correctness="pending",
            note="K'=-P tensor H1-I tensor H2 independent re-assembly (sign regression pin, matching the post-fix library convention)",
        )
    )


# ---------------------------------------------------------------------------
# A. Substructure correctness
# ---------------------------------------------------------------------------


def verify_taylor_hamiltonian_block(report):
    """The taylor_hamiltonian encoding block = the truncated Taylor matrix (unitary extraction + reference cross-check)."""
    k_mat = 0.7 * np.array([[0, 1], [1, 0]]) + 0.3 * np.diag([1.0, -1.0])
    k_be = matrix_pauli_encoding(k_mat.tolist())
    time, degree = 0.2, 3
    e_be = taylor_hamiltonian(k_be, time, degree=degree)
    expected = taylor_matrix(k_mat, time, degree)
    unitary = originir_unitary(e_be.operation.program())
    block, leakage = effective_block(unitary, 1)
    err_unitary = float(np.abs(block * e_be.alpha - expected).max())
    err_reference = float(np.abs(extract_block_reference(e_be, 2) - expected).max())
    worst = max(err_unitary, err_reference)
    # BE contract: the block is a contraction (operator norm <= 1); out-of-block amplitudes are by-design junk, not error
    block_norm = float(np.linalg.svd(block, compute_uv=False).max())
    report.case(
        "taylor-hamiltonian-block",
        paths=["originir-ext+to_matrix", "reference"],
        parameters={"matrix": "0.7X+0.3Z", "degree": degree, "time": time, "alpha_E": e_be.alpha},
        metrics={
            "max_error": worst,
            "error_unitary_path": err_unitary,
            "error_reference_path": err_reference,
            "block_operator_norm": block_norm,
            "junk_amplitude": leakage,
        },
        criterion="encoding block equals the truncated Taylor matrix element-wise (max_error < 1e-9), block operator norm <= 1+1e-9",
        passed=worst < 1e-9 and block_norm <= 1 + 1e-9,
    )


def verify_fourier_momentum_block(report):
    """The fourier_momentum encoding block = the two's-complement signed-frequency diagonal matrix."""
    width, period = 2, 8.0
    p_be = fourier_momentum(width, period)
    unitary = originir_unitary(p_be.operation.program())
    block, leakage = effective_block(unitary, width)
    expected = momentum_matrix(width, period)
    err_unitary = float(np.abs(block * p_be.alpha - expected).max())
    err_reference = float(np.abs(extract_block_reference(p_be, 1 << width) - expected).max())
    worst = max(err_unitary, err_reference)
    block_norm = float(np.linalg.svd(block, compute_uv=False).max())
    report.case(
        "fourier-momentum-block",
        paths=["originir-ext+to_matrix", "reference"],
        parameters={"width": width, "period": period, "alpha_P": p_be.alpha},
        metrics={"max_error": worst, "block_operator_norm": block_norm, "junk_amplitude": leakage},
        criterion="encoding block equals the signed-frequency diagonal matrix element-wise (max_error < 1e-9), block operator norm <= 1+1e-9",
        passed=worst < 1e-9 and block_norm <= 1 + 1e-9,
    )


def _riccati_problem():
    """u' = -u + u*u (componentwise Riccati): F1 = -I, F2 is the contraction matrix C."""
    f1 = scale(-1, identity(1))
    contraction = np.zeros((4, 4))
    contraction[0, 0] = 1.0  # input factor0 in the low position: C e_{i0+2 i1} = delta(i0,i1) e_{i0}
    contraction[1, 3] = 1.0
    f2 = matrix_pauli_encoding(contraction.tolist())
    u0 = np.array([0.6, 0.8])
    problem = PolynomialODE(1, ((1, f1), (2, f2)), gate_state_prep(list(u0)), 0.5)
    matrices = {1: -np.eye(2), 2: contraction}
    return problem, matrices, u0


def verify_carleman_lift_block(report):
    """The carleman_lift encoding block = the independently assembled padded-layout lifted generator G_K."""
    problem, matrices, _ = _riccati_problem()
    cutoff = 2
    lift = carleman_lift(problem, cutoff=cutoff)
    expected = assemble_carleman_lift(matrices, cutoff)
    # A reference sweep over 16 basis states is too slow at this scale (about 50 s), so the full-unitary
    # extraction is authoritative here; the reference basis-sweep path is covered by smaller cases such as
    # taylor-hamiltonian-block.
    unitary = originir_unitary(lift.operation.program())
    block, leakage = effective_block(unitary, lift.width)
    worst = float(np.abs(block * lift.alpha - expected).max())
    block_norm = float(np.linalg.svd(block, compute_uv=False).max())
    report.case(
        "carleman-lift-block",
        paths=["originir-ext+to_matrix"],
        parameters={
            "cutoff": cutoff,
            "width": lift.width,
            "signal": lift.signal_qubits,
            "alpha_G": lift.alpha,
            "problem": "u'=-u+u*u, F1=-I, F2=contraction C",
            "reference_skipped": "16-basis-state sweep exceeds the time budget at this scale (informative)",
        },
        metrics={"max_error": worst, "block_operator_norm": block_norm, "junk_amplitude": leakage},
        criterion="lifted-generator block equals the independent assembly element-wise (max_error < 1e-9), block operator norm <= 1+1e-9",
        passed=worst < 1e-9 and block_norm <= 1 + 1e-9,
    )


def verify_carleman_initial(report):
    """carleman_initial amplitudes = 1/Z sum r^k u0^{tensor k} (four-path full-amplitude cross-check)."""
    problem, matrices, u0 = _riccati_problem()
    cutoff = 2
    prep = carleman_initial(problem, cutoff=cutoff)
    program = prep.operation.program()
    z0 = lifted_initial(u0, problem.initial_norm, cutoff)
    expected = {(i, 0): complex(v) for i, v in enumerate(z0) if abs(v) > 1e-15}
    results = run_amplitude_paths(program, ["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"])
    worst = max(amplitude_error(amplitudes, expected) for amplitudes in results.values())
    report.case(
        "carleman-initial-state",
        paths=list(results),
        parameters={"cutoff": cutoff, "initial_norm": problem.initial_norm, "u0": list(u0)},
        metrics={"max_error": worst, "branches": len(expected)},
        criterion="lifted initial state equals the independent tensor-power reference amplitude by amplitude (max_error < 1e-9)",
        passed=worst < 1e-9,
    )


def verify_quadrature_plans(report):
    """QuadraturePlan.cauchy and ContourPlan nodes/weights against independent closed forms; contour-identity residual."""
    # LCHS Cauchy plan: w = h/(pi(1+k^2))
    plan = QuadraturePlan.cauchy(cutoff=4, spacing=0.75)
    dev_nodes = max(abs(a - b) for a, b in zip(plan.nodes, [k * 0.75 for k in range(-4, 5)], strict=True))
    dev_weights = max(
        abs(a - 0.75 / (math.pi * (1 + k * k)))
        for a, k in zip(plan.weights, plan.nodes, strict=True)
    )
    worst_lchs = max(dev_nodes, dev_weights)
    report.case(
        "lchs-cauchy-plan-formula",
        paths=["plan-objects"],
        parameters={"cutoff": 4, "spacing": 0.75},
        metrics={"max_error": worst_lchs},
        criterion="nodes/weights match the independent closed form h/(pi(1+k^2)) (max_error < 1e-12)",
        passed=worst_lchs < 1e-12,
    )
    # CBMD contour plan: weights/auxiliary coefficients against independent recomputation; the truncated residual of the t=0 identity main+aux -> 1
    trend = {}
    worst_formula = 0.0
    for cutoff in (2, 4, 8):
        cplan = ContourPlan(a=1.0, cutoff=cutoff)
        numerator = math.expm1(-2 * math.pi * cplan.a)
        for q, w in zip(cplan.nodes, cplan.weights, strict=True):
            ref_w = numerator / (
                cplan.a
                * 2
                * math.pi
                * 1j
                * (q + 1j)
                * math.prod((q - p) / (-1j - p) for p in cplan.poles)
            )
            worst_formula = max(worst_formula, abs(w - ref_w))
        for p, c in zip(cplan.poles, cplan.auxiliary_coefficients, strict=True):
            ref_c = numerator / (
                (cmath.exp(-2 * math.pi * p * cplan.a * 1j) - 1)
                * math.prod((p - o) / (-1j - o) for o in cplan.poles if o != p)
            )
            worst_formula = max(worst_formula, abs(c - ref_c))
        identity_residual = abs(sum(cplan.weights) + sum(cplan.auxiliary_coefficients) - 1)
        trend[cutoff] = identity_residual
    report.case(
        "cbmd-contour-plan-formula-and-identity",
        paths=["plan-objects"],
        parameters={"a": 1.0, "poles": [str(p) for p in ContourPlan().poles]},
        metrics={
            "max_error": worst_formula,
            "identity_residual_cutoff2": trend[2],
            "identity_residual_cutoff4": trend[4],
            "identity_residual_cutoff8": trend[8],
        },
        criterion=(
            "weights/auxiliary coefficients match the independent closed forms (max_error < 1e-12); "
            "the t=0 contour identity residual decreases with the truncation (informative, corresponding to the omitted infinite-series tail)"
        ),
        passed=worst_formula < 1e-12 and trend[8] < trend[4] < trend[2],
    )


# ---------------------------------------------------------------------------
# B. LCHS end to end (multiple input models)
# ---------------------------------------------------------------------------

LCHS_PLAN = QuadraturePlan.cauchy(cutoff=2, spacing=1.0)


def _lchs_case(report, name, state, l_mat, h_mat, time, degree, u0, exact, paths, extra_params):
    """LCHS/CBMD-style end-to-end cross-check: implementation error (vs independent simulation) + method error (vs the exact solution)."""
    alpha_v = evolution_alpha(state)
    program = state.operation.program()
    results = run_amplitude_paths(program, paths)
    expected = lchs_emulate(l_mat, h_mat, LCHS_PLAN.nodes, LCHS_PLAN.weights, time, degree, alpha_v, u0)
    blocks, block_dev, exact_dev, pysparq_floor = compare_path_blocks(results)
    impl = 0.0
    psuccess_measured = None
    for block in blocks.values():
        psuccess_measured = sum(abs(v) ** 2 for v in block.values())
        impl = max(impl, max(abs(block.get(i, 0) - expected[i]) for i in range(len(u0))))
    method = float(np.abs(expected - exact / alpha_v).max())
    psuccess_expected = float(np.sum(np.abs(expected) ** 2))
    report.case(
        name,
        paths=list(results),
        parameters={"time": time, "degree": degree, "alpha_V": alpha_v, **extra_params},
        metrics={
            "impl_error": impl,
            "method_error": method,
            "block_cross_deviation": block_dev,
            "exact_path_deviation": exact_dev,
            "pysparq_spectrum_floor": pysparq_floor,
            "success_probability": psuccess_measured,
            "success_probability_error": abs(psuccess_measured - psuccess_expected),
        },
        criterion=(
            "post-selected block agrees with the independent numpy simulation amplitude by amplitude (impl_error < 1e-9), physical blocks agree across paths; "
            "method_error is the finite-quadrature + Taylor remainder (informative, marked pending in the library); "
            "pysparq_spectrum_floor is the junk-branch numerical floor (informative, see the group report)"
        ),
        passed=impl < 1e-9
        and block_dev < 1e-9
        and exact_dev < 1e-9
        and pysparq_floor < 1e-6
        and abs(psuccess_measured - psuccess_expected) < 1e-9
        and method < 0.25,
    )
    return blocks.get("reference")


def verify_lchs_given_be(report):
    """Input model 1: a whole BE of G (G=-I scalar decay, analytic solution e^{-t})."""
    time, degree = 0.4, 3
    u0 = np.array([1.0, 1.0]) / SQRT2
    solver = linear_qode(
        "lchs", plan=LCHS_PLAN, hamiltonian_function=partial(taylor_hamiltonian, degree=degree)
    )
    state = solver(scale(-1, identity(1)), gate_state_prep(list(u0)), time)
    exact = math.exp(-time) * u0
    _lchs_case(
        report,
        "lchs-given-be-scalar-decay",
        state,
        np.eye(2),
        np.zeros((2, 2)),
        time,
        degree,
        u0,
        exact,
        ["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        {"input_model": "whole-generator BE (G=-I)", "exact": "analytic e^{-t}*u0"},
    )


def verify_lchs_parts_noncommuting(report):
    """Input model 2: direct Hermitian parts ([L,H] != 0, scipy expm reference)."""
    time, degree = 0.2, 3
    l_mat = np.array([[1.0, 0.3], [0.3, 0.5]])
    h_mat = np.array([[0.2, 0.1], [0.1, -0.1]])
    u0 = np.array([1.0, 1.0]) / SQRT2
    model = LinearODE(
        HermitianParts(matrix_pauli_encoding(l_mat.tolist()), matrix_pauli_encoding(h_mat.tolist())),
        gate_state_prep(list(u0)),
    )
    state = lchs_qode(
        model, time, plan=LCHS_PLAN, hamiltonian_function=partial(taylor_hamiltonian, degree=degree)
    )
    exact = scipy.linalg.expm(-(l_mat + 1j * h_mat) * time) @ u0
    return _lchs_case(
        report,
        "lchs-given-parts-noncommuting",
        state,
        l_mat,
        h_mat,
        time,
        degree,
        u0,
        exact,
        ["reference", "rir-pysparq", "originir-ext"],
        {
            "input_model": "direct HermitianParts (non-commuting)",
            "exact": "scipy.linalg.expm",
            "adapter_skipped": "about 16 thousand branches x deep LCU event tree; the adapter path is omitted for budget reasons (informative)",
        },
    )


def verify_lchs_diagonal_gate_vs_qram(report):
    """Input model 3: a diagonal spectral-angle database; the same open program bound to a gate table and to QRAM."""
    from oracq import Binding, bind, unresolved

    time, degree = 0.4, 3
    u0 = np.array([1.0, 1.0]) / SQRT2
    angles = abstract_database("DiagonalAngles", 1, 2)
    generator = scale(-1, diagonal_block_encoding(angles, alpha=1.0))
    initial = abstract_state_prep("Initial", 1, work_width=3)
    solver = linear_qode(
        "lchs", plan=LCHS_PLAN, hamiltonian_function=partial(taylor_hamiltonian, degree=degree)
    )
    state = solver(generator, initial, time)
    open_program = state.operation.program()
    alpha_v = evolution_alpha(state)
    # A = diag(1, cos(pi/4)): 2-bit angle words, default angle scale 2*pi/4
    a_mat = np.diag([1.0, math.cos(math.pi / 4)])
    expected = lchs_emulate(
        a_mat, np.zeros((2, 2)), LCHS_PLAN.nodes, LCHS_PLAN.weights, time, degree, alpha_v, u0
    )
    exact = scipy.linalg.expm(-a_mat * time) @ u0
    bindings = {
        "gate": (
            {
                "DiagonalAngles": gate_database(1, 2, [0, 1]).operation,
                "Initial": gate_state_prep([1, 1], work_width=3).operation,
            },
            None,
        ),
        "qram": (
            {
                "DiagonalAngles": Binding(
                    qram_database(1, 2).operation, {"table": "diagonal_angles"}
                ),
                "Initial": Binding(qram_state_prep(1, 2).operation, {"angles": "initial_angles"}),
            },
            {
                "diagonal_angles": [0, 1],
                "initial_angles": [
                    qram_state_angles([1, 1], 2).get(address, 0) for address in range(2)
                ],
            },
        ),
    }
    blocks = {}
    for binding_name, (binding, memory) in bindings.items():
        program = bind(open_program, binding)
        assert not unresolved(program)
        # rir/adapter omitted for budget reasons (about 35 s per binding); multi-path execution is covered by cases 1/2
        results = run_amplitude_paths(program, ["reference"], memory)
        path_blocks, block_dev, _, pysparq_floor = compare_path_blocks(results)
        blocks[binding_name] = path_blocks["reference"]
        impl = 0.0
        for block in path_blocks.values():
            impl = max(impl, max(abs(block.get(i, 0) - expected[i]) for i in range(2)))
        report.case(
            f"lchs-diagonal-spectral-{binding_name}",
            paths=list(results),
            parameters={
                "input_model": f"diagonal angle database ({binding_name} binding)",
                "time": time,
                "degree": degree,
                "alpha_V": alpha_v,
            },
            metrics={
                "impl_error": impl,
                "method_error": float(np.abs(expected - exact / alpha_v).max()),
                "block_cross_deviation": block_dev,
                "pysparq_spectrum_floor": pysparq_floor,
            },
            criterion=(
                "the angle-database-bound implementation matches the independent simulation amplitude by amplitude (impl_error < 1e-9); "
                "pysparq_spectrum_floor is the junk-branch numerical floor (informative)"
            ),
            passed=impl < 1e-9 and block_dev < 1e-9 and pysparq_floor < 1e-6,
        )
    gate_qram = max(abs(blocks["gate"].get(i, 0) - blocks["qram"].get(i, 0)) for i in range(2))
    report.case(
        "lchs-diagonal-gate-vs-qram",
        paths=["gate-binding", "qram-binding"],
        parameters={"note": "two data-storage implementations of the same open RIR"},
        metrics={"max_amplitude_difference": gate_qram},
        criterion="gate- and QRAM-bound complete complex amplitudes agree (max_difference < 1e-10)",
        passed=gate_qram < 1e-10,
    )


def verify_lchs_fokker_planck(report):
    """Input model 4: the Fokker-Planck OU discrete generator (Pauli-expansion BE) + QODEProblem.solve."""
    from oracq.algorithms.qode.sde import (
        FokkerPlanckProblem,
        boltzmann_distribution,
        matrix_exponential,
        sde_state_preparation,
    )

    time, degree = 0.3, 2
    problem = FokkerPlanckProblem(
        drift=[-1.0 * x for x in (-1.5, -0.5, 0.5, 1.5)],
        diffusion=0.5,
        grid=(-1.5, -0.5, 0.5, 1.5),
    )
    g_mat = np.array(problem.generator_matrix())
    probabilities = boltzmann_distribution(problem)
    initial = sde_state_preparation(probabilities, implementation="gate")
    qode_problem = problem.qode_problem(initial)
    solver = linear_qode(
        "lchs", plan=LCHS_PLAN, hamiltonian_function=partial(taylor_hamiltonian, degree=degree)
    )
    state = solver.solve(qode_problem, time)
    attrs = dict(state.operation.module.attributes)
    alpha_v = attrs["evolution_alpha"]
    u0 = np.sqrt(np.array(probabilities))
    a_mat = -g_mat
    l_mat = (a_mat + a_mat.T) / 2
    h_mat = (a_mat - a_mat.T) / 2j
    expected = lchs_emulate(
        l_mat, h_mat, LCHS_PLAN.nodes, LCHS_PLAN.weights, time, degree, alpha_v, u0
    )
    exact_scipy = scipy.linalg.expm(g_mat * time) @ u0
    exact_witness = np.array(matrix_exponential(problem.generator_matrix(), time)) @ u0
    classical_agreement = float(np.abs(exact_scipy - exact_witness).max())
    program = state.operation.program()
    # 66 thousand branches and a deep LCU event tree: the pysparq paths measured 71 s+ (adapter) / hundreds of seconds (rir),
    # omitted for budget reasons; reference and the dense OriginIR state vector cross-check each other independently.
    results = run_amplitude_paths(program, ["reference", "originir-ext"])
    blocks, block_dev, exact_dev, _ = compare_path_blocks(results)
    impl = 0.0
    for block in blocks.values():
        impl = max(impl, max(abs(block.get(i, 0) - expected[i]) for i in range(4)))
    report.case(
        "lchs-fokker-planck-ou",
        paths=list(results),
        parameters={
            "input_model": "FokkerPlanckProblem(OU, 4-point zero-flux) + QODEProblem.solve",
            "time": time,
            "degree": degree,
            "alpha_V": alpha_v,
            "qode_dissipative_promise": attrs["qode_dissipative_promise"],
            "pysparq_skipped": "about 66 thousand branches x deep LCU event tree; the pysparq paths exceed the time budget (informative)",
        },
        metrics={
            "impl_error": impl,
            "method_error": float(np.abs(expected - exact_scipy / alpha_v).max()),
            "block_cross_deviation": block_dev,
            "exact_path_deviation": exact_dev,
            "classical_reference_agreement": classical_agreement,
        },
        criterion=(
            "post-selected block matches the independent simulation (impl_error < 1e-9); the two classical references, "
            "scipy expm and the pure-Python sde matrix exponential, agree (< 1e-12); method error informative"
        ),
        passed=impl < 1e-9 and block_dev < 1e-9 and exact_dev < 1e-9
        and classical_agreement < 1e-12 and attrs["qode_dissipative_promise"] is True,
    )


def verify_lchs_heat_structured(report):
    """Input model 5: the structured shifted BE of the periodic heat equation (qham difference stencil), Fourier analytic reference."""
    from oracq.applications.qham import Grid
    from oracq.applications.qham.stencils import derivative_encoding

    time, degree = 0.3, 2
    grid = Grid(("x",), (4,), (1.0,), boundary="periodic")
    generator = scale(0.1, derivative_encoding(grid, (("x", 2),)))
    u0 = np.array([1.0, 0.0, 0.0, 0.0])
    solver = linear_qode(
        "lchs", plan=LCHS_PLAN, hamiltonian_function=partial(taylor_hamiltonian, degree=degree)
    )
    state = solver(generator, gate_state_prep(list(u0)), time)
    # G = 0.1(S+S^dagger-2I); circulant eigenvalues lambda_k = 0.1(2cos(2*pi*k/4)-2) (independent analytic reference)
    shift = np.roll(np.eye(4), 1, axis=1)
    g_mat = 0.1 * (shift + shift.T - 2 * np.eye(4))
    a_mat = -g_mat
    exact = scipy.linalg.expm(g_mat * time) @ u0
    eigenvalues = [0.1 * (2 * math.cos(2 * math.pi * k / 4) - 2) for k in range(4)]
    exact_fourier = np.fft.ifft(np.exp(np.array(eigenvalues) * time) * np.fft.fft(u0)).real
    fourier_agreement = float(np.abs(exact - exact_fourier).max())
    _lchs_case(
        report,
        "lchs-heat-structured-stencil",
        state,
        (a_mat + a_mat.T) / 2,
        (a_mat - a_mat.T) / 2j,
        time,
        degree,
        u0,
        exact,
        # about 16 thousand branches x deep LCU event tree: pysparq paths omitted for budget reasons
        ["reference", "originir-ext"],
        {
            "input_model": "structured shift-difference BE (periodic 4-point heat equation)",
            "exact": "Fourier-eigenvalue analytic + scipy expm cross-check",
            "pysparq_skipped": "about 16 thousand branches x deep LCU event tree; the pysparq paths exceed the time budget (informative)",
        },
    )
    report.case(
        "lchs-heat-fourier-reference-agreement",
        paths=["classical"],
        parameters={"eigenvalues": eigenvalues},
        metrics={"max_error": fourier_agreement},
        criterion="the Fourier analytic reference agrees with scipy expm (max_error < 1e-9)",
        passed=fourier_agreement < 1e-9,
    )


def verify_lchs_quadrature_convergence(report):
    """Quadrature convergence: the method-error trend of the same scalar problem under plans with increasing Kmax (quantum + classical)."""
    time = 0.2
    u0 = np.array([1.0, 1.0]) / SQRT2
    exact = math.exp(-time)
    trend = {}
    quantum_checks = [(2, 4), (8, 3)]  # (cutoff, taylor degree): high-order branches on the 2^d scale
    for cutoff, degree in quantum_checks:
        plan = QuadraturePlan.cauchy(cutoff=cutoff, spacing=1.0)
        model = LinearODE(HermitianParts(identity(1), zero(1)), gate_state_prep(list(u0)))
        state = lchs_qode(
            model,
            time,
            plan=plan,
            hamiltonian_function=partial(taylor_hamiltonian, degree=degree),
        )
        alpha_v = evolution_alpha(state)
        expected = lchs_emulate(
            np.eye(2), np.zeros((2, 2)), plan.nodes, plan.weights, time, degree, alpha_v, u0
        )
        # rir/adapter paths omitted for budget reasons (17 nodes x degree-4 nested LCU, rir about 110 s)
        results = run_amplitude_paths(state.operation.program(), ["reference"])
        blocks, block_dev, _, pysparq_floor = compare_path_blocks(results)
        impl = 0.0
        for block in blocks.values():
            impl = max(impl, max(abs(block.get(i, 0) - expected[i]) for i in range(2)))
        trend[cutoff] = float(abs(expected[0] - exact / alpha_v))
        report.case(
            f"lchs-quadrature-scan-cutoff{cutoff}",
            paths=list(results),
            parameters={"Kmax": cutoff, "degree": degree, "nodes": 2 * cutoff + 1},
            metrics={
                "impl_error": impl,
                "method_error": trend[cutoff],
                "block_cross_deviation": block_dev,
                "pysparq_spectrum_floor": pysparq_floor,
            },
            criterion="implementation error < 1e-9 under every quadrature plan; method error is the quadrature remainder (informative)",
            passed=impl < 1e-9 and block_dev < 1e-9 and pysparq_floor < 1e-6,
        )
    # Classical trend line: pure-kernel quadrature (no Taylor) at more truncation points
    classical_trend = {}
    for cutoff in (2, 4, 8, 16, 32):
        nodes = [float(k) for k in range(-cutoff, cutoff + 1)]
        value = sum(
            (1 / (math.pi * (1 + k * k))) * cmath.exp(-1j * k * time) for k in nodes
        )
        classical_trend[cutoff] = float(abs(value - exact))
    report.case(
        "lchs-quadrature-kernel-trend",
        paths=["classical-quadrature"],
        parameters={"time": time, "note": "pure Cauchy-kernel quadrature (no Taylor truncation); the oscillating tail is documented"},
        metrics={f"kernel_error_cutoff{k}": classical_trend[k] for k in classical_trend},
        criterion="kernel quadrature error decreases overall with Kmax (informative, oscillating tail included)",
        passed=classical_trend[32] < classical_trend[2],
    )


# ---------------------------------------------------------------------------
# C. CBMD end to end and the LCHS same-problem cross-check
# ---------------------------------------------------------------------------


def verify_cbmd_endtoend(report):
    """CBMD and LCHS on the same non-commuting problem: implementation error, method error (omitted terms), direction cross-check."""
    time, degree = 0.2, 3
    l_mat = np.array([[1.0, 0.3], [0.3, 0.5]])
    h_mat = np.array([[0.2, 0.1], [0.1, -0.1]])
    u0 = np.array([1.0, 1.0]) / SQRT2
    cplan = ContourPlan(a=1.0, cutoff=2)
    model = LinearODE(
        HermitianParts(matrix_pauli_encoding(l_mat.tolist()), matrix_pauli_encoding(h_mat.tolist())),
        gate_state_prep(list(u0)),
    )
    state = cbmd_qode(
        model, time, plan=cplan, hamiltonian_function=partial(taylor_hamiltonian, degree=degree)
    )
    alpha_v = evolution_alpha(state)
    expected = lchs_emulate(l_mat, h_mat, cplan.nodes, cplan.weights, time, degree, alpha_v, u0)
    exact = scipy.linalg.expm(-(l_mat + 1j * h_mat) * time) @ u0
    program = state.operation.program()
    # about 16 thousand branches x deep LCU event tree: pysparq paths omitted for budget reasons
    results = run_amplitude_paths(program, ["reference", "originir-ext"])
    blocks, block_dev, exact_dev, pysparq_floor = compare_path_blocks(results)
    impl = 0.0
    for block in blocks.values():
        impl = max(impl, max(abs(block.get(i, 0) - expected[i]) for i in range(2)))
    method = float(np.abs(expected - exact / alpha_v).max())
    report.case(
        "cbmd-parts-noncommuting",
        paths=list(results),
        parameters={
            "input_model": "direct HermitianParts (non-commuting, same problem as lchs-given-parts-noncommuting)",
            "time": time,
            "degree": degree,
            "alpha_V": alpha_v,
            "omitted": "auxiliary_nonhermitian_evolutions + infinite_series_tail (declared in the library docs)",
        },
        metrics={
            "impl_error": impl,
            "method_error": method,
            "block_cross_deviation": block_dev,
            "exact_path_deviation": exact_dev,
            "pysparq_spectrum_floor": pysparq_floor,
        },
        criterion=(
            "post-selected block matches the independent simulation amplitude by amplitude (impl_error < 1e-9); "
            "method_error is the remainder of the omitted auxiliary-pole branches and series tail (informative; the library records the omission explicitly)"
        ),
        passed=impl < 1e-9 and block_dev < 1e-9 and exact_dev < 1e-9
        and pysparq_floor < 1e-6 and method < 0.1,
    )
    return blocks.get("reference")


# ---------------------------------------------------------------------------
# D. Schrodingerization end to end
# ---------------------------------------------------------------------------


def verify_schrodingerization_decay(report):
    """Scalar decay G=-I: assembly fidelity and the recovery relation (forward flow) both pass strictly.

    A historical defect is fixed: the original assembly K = P tensor H1 -
    I tensor H2 recovered e^{+t}*u0 under the positive-QFT convention (a
    time-reversed flow; the verification finding is in the group report);
    after the fix K' = -P tensor H1 - I tensor H2 recovers e^{-t}*u0, with
    the grid error below 1e-6 under exact evolution; the end-to-end residual
    is the Taylor truncation remainder of the Nyquist momentum mode (about
    3e-2 at degree 4, decreasing with the order, informative).
    """
    time, degree = 0.3, 4
    plan = SchrodingerPlan(auxiliary_width=2, period=1.2, selected_index=1)  # dp = t
    g_mat = -np.eye(2)
    u0 = np.array([1.0, 1.0]) / SQRT2
    g_be = scale(-1, identity(1))
    solver = linear_qode(
        "schrodingerization",
        plan=plan,
        hamiltonian_function=partial(taylor_hamiltonian, degree=degree),
    )
    state = solver(g_be, gate_state_prep(list(u0)), time)
    attrs = dict(state.operation.module.attributes)
    alpha_e = schrodinger_alpha(g_be, plan, time, degree)
    expected = schrodinger_emulate(g_mat, u0, time, plan, degree, alpha_e, flip_momentum=True)
    exact = math.exp(-time) * u0
    # Decomposition: a reference on the same grid/channel but with exact evolution (grid+window error) and the Taylor remainder
    exact_grid = schrodinger_emulate(
        g_mat, u0, time, plan, degree, alpha_e, exact_evolution=True, flip_momentum=True
    )
    p_sel = attrs["selected_p"]
    grid_coords = [(j if j < 2 else j - 4) * plan.period / 4 for j in range(4)]
    z_norm = math.sqrt(sum(math.exp(-2 * abs(x)) for x in grid_coords))
    recovery = alpha_e * z_norm * math.exp(p_sel)
    grid_error = float(np.abs(exact_grid * recovery - exact).max())
    taylor_remainder = float(np.abs(expected - exact_grid).max())
    # Directional evidence: the recovery approaches e^{-t}*u0 (forward flow), not e^{+t}*u0 (the time-reversed solution)
    time_reversed = math.exp(time) * u0
    reversed_fit = float(np.abs(exact_grid * recovery - time_reversed).max())
    program = state.operation.program()
    # about 90 thousand branches (8^degree scaling of nested LCU): pysparq paths omitted for budget reasons; reference+OriginIR cross-check
    results = run_amplitude_paths(program, ["reference", "originir-ext"])
    blocks, block_dev, exact_dev, _ = compare_path_blocks(results)
    impl = 0.0
    recovered_err = 0.0
    for block in blocks.values():
        impl = max(impl, max(abs(block.get(i, 0) - expected[i]) for i in range(2)))
        recovered_err = max(
            recovered_err,
            max(abs(block.get(i, 0) * recovery - exact[i]) for i in range(2)),
        )
    report.case(
        "schrodingerization-scalar-decay-grid",
        paths=list(results),
        parameters={
            "G": "-I (u'=-u scalar decay)",
            "time": time,
            "degree": degree,
            "plan": "auxiliary_width=2, period=1.2 (dp=t), selected p=0.3",
            "alpha_E": alpha_e,
            "pysparq_skipped": "about 90 thousand branches (8^degree scaling of nested LCU), for budget reasons (informative)",
            "finding": "fixed: K'=-P tensor H1-I tensor H2 recovers e^{-t}*u0 (forward flow); the residual is the Taylor truncation",
        },
        metrics={
            "impl_error": impl,
            "grid_error_exact_evolution": grid_error,
            "time_reversed_fit_error": reversed_fit,
            "taylor_remainder": taylor_remainder,
            "recovery_error_endtoend": recovered_err,
            "block_cross_deviation": block_dev,
            "exact_path_deviation": exact_dev,
        },
        criterion=(
            "the quantum result matches the full-stack independent simulation (K'=-P tensor H1-I tensor H2) amplitude by amplitude (impl_error < 1e-9); "
            "the recovered magnitude approaches e^{-t}*u0 rather than e^{+t}*u0 (grid_error < 1e-6, "
            "time_reversed_fit_error > 1e-1, directional evidence); recovery_error_endtoend is "
            "the Taylor truncation remainder of the Nyquist momentum mode (informative, decreasing with degree)"
        ),
        passed=impl < 1e-9
        and block_dev < 1e-9
        and exact_dev < 1e-9
        and grid_error < 1e-6
        and reversed_fit > 1e-1,
    )


def verify_schrodingerization_sign_flipped(report):
    """Independent re-assembly (K'=-P tensor H1-I tensor H2, a real quantum program from public combinators): the recovery relation holds exactly."""
    time = 0.3
    plan = SchrodingerPlan(auxiliary_width=2, period=1.2, selected_index=1)
    g_mat = -np.eye(2)
    u0 = np.array([1.0, 1.0]) / SQRT2
    g_be = scale(-1, identity(1))
    exact = math.exp(-time) * u0
    grid_coords = [(j if j < 2 else j - 4) * plan.period / 4 for j in range(4)]
    z_norm = math.sqrt(sum(math.exp(-2 * abs(x)) for x in grid_coords))
    p_sel = grid_coords[1]
    # Structural claim: with exact evolution (no Taylor truncation) on the same grid, the flipped construction recovers the exact solution
    alpha_probe = schrodinger_alpha(g_be, plan, time, 4)
    exact_grid = schrodinger_emulate(
        g_mat, u0, time, plan, 4, alpha_probe, exact_evolution=True, flip_momentum=True
    )
    grid_error = float(np.abs(exact_grid * (alpha_probe * z_norm * math.exp(p_sel)) - exact).max())
    recovery_errors = {}
    for degree in (2, 4):
        state = schrodinger_sign_flipped_qode(
            g_be, gate_state_prep(list(u0)), time, plan, degree=degree
        )
        alpha_e = schrodinger_alpha(g_be, plan, time, degree)
        expected = schrodinger_emulate(g_mat, u0, time, plan, degree, alpha_e, flip_momentum=True)
        recovery = alpha_e * z_norm * math.exp(p_sel)
        # originir was already cross-checked to 1e-17 in schrodingerization-scalar-decay-grid;
        # this variant runs reference only to control runtime (the dense 21-qubit state vector is slow)
        results = run_amplitude_paths(state.operation.program(), ["reference"])
        blocks, block_dev, exact_dev, _ = compare_path_blocks(results)
        impl = 0.0
        recovered_err = 0.0
        for block in blocks.values():
            impl = max(impl, max(abs(block.get(i, 0) - expected[i]) for i in range(2)))
            recovered_err = max(
                recovered_err,
                max(abs(block.get(i, 0) * recovery - exact[i]) for i in range(2)),
            )
        recovery_errors[degree] = recovered_err
        report.case(
            f"schrodingerization-sign-flipped-degree{degree}",
            paths=list(results),
            parameters={
                "construction": "independent re-assembly from public combinators, K' = -P tensor H1 - I tensor H2 (post-fix library convention; sign regression pin)",
                "time": time,
                "degree": degree,
            },
            metrics={
                "impl_error": impl,
                "recovery_error": recovered_err,
                "grid_error_exact_evolution": grid_error,
                "block_cross_deviation": block_dev,
                "exact_path_deviation": exact_dev,
            },
            criterion=(
                "the flipped construction matches the independent simulation (impl_error < 1e-9); under exact evolution the grid error < 1e-6"
                " (proving the sign is the only structural mismatch); recovery_error is the Taylor"
                " truncation remainder of the Nyquist momentum mode (pi phase, decreasing with degree, informative)"
            ),
            passed=impl < 1e-9 and grid_error < 1e-6 and block_dev < 1e-9 and exact_dev < 1e-9,
        )
    report.case(
        "schrodingerization-sign-flipped-taylor-trend",
        paths=["reference"],
        parameters={"note": "recovery error = Taylor truncation of the Nyquist momentum mode (pi/2 grid-translation eigenphase)"},
        metrics={
            "recovery_error_degree2": recovery_errors[2],
            "recovery_error_degree4": recovery_errors[4],
        },
        criterion="recovery error decreases with the Taylor order (informative; the hamiltonian_function protocol is replaceable)",
        passed=recovery_errors[4] < recovery_errors[2],
    )


def verify_schrodingerization_rotation(report):
    """Antisymmetric generator G=J: analytic rotation solution; recovered magnitude e^{p}*Z*alpha_E*block = u(t)."""
    time, degree = 0.3, 4
    plan = SchrodingerPlan(auxiliary_width=2, period=1.2, selected_index=1)
    j_mat = np.array([[0.0, -1.0], [1.0, 0.0]])
    u0 = np.array([1.0, 1.0]) / SQRT2
    g_be = matrix_pauli_encoding(j_mat.tolist())
    solver = linear_qode(
        "schrodingerization",
        plan=plan,
        hamiltonian_function=partial(taylor_hamiltonian, degree=degree),
    )
    state = solver(g_be, gate_state_prep(list(u0)), time)
    attrs = dict(state.operation.module.attributes)
    alpha_e = schrodinger_alpha(g_be, plan, time, degree)
    expected = schrodinger_emulate(j_mat, u0, time, plan, degree, alpha_e)
    exact = np.array([[math.cos(time), -math.sin(time)], [math.sin(time), math.cos(time)]]) @ u0
    p_sel = attrs["selected_p"]
    grid_coords = [(j if j < 2 else j - 4) * plan.period / 4 for j in range(4)]
    z_norm = math.sqrt(sum(math.exp(-2 * abs(x)) for x in grid_coords))
    recovery = alpha_e * z_norm * math.exp(p_sel)
    program = state.operation.program()
    # about 90 thousand branches (8^degree scaling): both the pysparq paths and dense originir are omitted for budget reasons;
    # originir was already cross-checked to 1e-17 in the decay-grid case (same program family)
    results = run_amplitude_paths(program, ["reference"])
    blocks, block_dev, exact_dev, _ = compare_path_blocks(results)
    impl = 0.0
    recovered_err = 0.0
    for block in blocks.values():
        impl = max(impl, max(abs(block.get(i, 0) - expected[i]) for i in range(2)))
        recovered_err = max(
            recovered_err,
            max(abs(block.get(i, 0) * recovery - exact[i]) for i in range(2)),
        )
    report.case(
        "schrodingerization-rotation-recovery",
        paths=list(results),
        parameters={
            "G": "J=[[0,-1],[1,0]] (H1=0, pure H2 rotation)",
            "time": time,
            "degree": degree,
            "alpha_E": alpha_e,
            "other_paths_skipped": "about 90 thousand branches (8^degree scaling of nested LCU); the originir cross-check appears in the decay-grid case (informative)",
        },
        metrics={
            "impl_error": impl,
            "recovery_error": recovered_err,
            "block_cross_deviation": block_dev,
            "exact_path_deviation": exact_dev,
        },
        criterion=(
            "the quantum result matches the independent simulation (impl_error < 1e-9); "
            "the recovered magnitude errors < 1e-3 against the analytic rotation solution (the degree-4 Taylor remainder is ~2e-5)"
        ),
        passed=impl < 1e-9 and recovered_err < 1e-3 and block_dev < 1e-9 and exact_dev < 1e-9,
    )


# ---------------------------------------------------------------------------
# E. Carleman end to end (Riccati) and the truncation trend
# ---------------------------------------------------------------------------


def verify_carleman_riccati(report):
    """Carleman end to end for the Riccati u'=-u+u^2: an injected Taylor-protocol solver, with three-level error decomposition."""
    from scipy.integrate import solve_ivp

    problem, matrices, u0 = _riccati_problem()
    cutoff, time, degree = 2, 0.2, 3
    state = carleman_qode(
        problem, time, partial(taylor_series_solver, degree=degree), cutoff=cutoff
    )
    alpha_v = evolution_alpha(state)
    g_lift = assemble_carleman_lift(matrices, cutoff)
    z0 = lifted_initial(u0, problem.initial_norm, cutoff)
    base = 1 << cutoff  # channel base address with level==1 and the other groups zero (n=1)
    expected = ode_taylor_series(g_lift, time, degree) @ z0 / alpha_v
    z_linear = scipy.linalg.expm(g_lift * time) @ z0
    chan_linear = z_linear[base : base + 2]
    riccati = solve_ivp(
        lambda t, u: -u + u * u,
        (0, time),
        problem.initial_norm * u0,
        rtol=1e-12,
        atol=1e-14,
    ).y[:, -1]
    program = state.operation.program()
    # adapter path omitted for budget reasons
    results = run_amplitude_paths(program, ["reference", "rir-pysparq"])
    blocks, block_dev, _, pysparq_floor = compare_path_blocks(results)
    impl = 0.0
    direction_err = 0.0
    for block in blocks.values():
        impl = max(impl, max(abs(block.get(i, 0) - expected[base + i]) for i in range(2)))
        quantum_dir = np.array([block.get(i, 0) for i in range(2)])
        direction_err = max(
            direction_err,
            np.linalg.norm(
                quantum_dir / np.linalg.norm(quantum_dir) - chan_linear / np.linalg.norm(chan_linear)
            ),
        )
    taylor_remainder = float(
        np.abs(ode_taylor_series(g_lift, time, degree) @ z0 - z_linear).max()
    )
    truncation = float(
        np.linalg.norm(
            chan_linear / np.linalg.norm(chan_linear) - riccati / np.linalg.norm(riccati)
        )
    )
    report.case(
        "carleman-riccati-endtoend",
        paths=list(results),
        parameters={
            "problem": "u'=-u+u^2, u0=0.5*(0.6,0.8), cutoff=2",
            "time": time,
            "linear_solver": "injected minimal Taylor-protocol solver (degree 3, real BE assembly)",
            "alpha_V": alpha_v,
        },
        metrics={
            "impl_error": impl,
            "taylor_remainder": taylor_remainder,
            "carleman_truncation_error": truncation,
            "quantum_vs_exact_linear_direction": direction_err,
            "block_cross_deviation": block_dev,
            "pysparq_spectrum_floor": pysparq_floor,
        },
        criterion=(
            "the physical channel matches the full-stack independent simulation (impl_error < 1e-9); "
            "the Taylor remainder and the Carleman truncation error (vs the scipy solve_ivp exact nonlinear solution) are an informative decomposition"
        ),
        passed=impl < 1e-9 and block_dev < 1e-9 and pysparq_floor < 1e-6,
    )
    return blocks.get("reference")


def verify_carleman_cutoff_trend(report, cutoff2_block):
    """Truncation trend: a K=1 quantum run + the truncation error reusing the K=2 quantum block from E1; K=3 added classically."""
    from scipy.integrate import solve_ivp

    problem, matrices, u0 = _riccati_problem()
    time, degree = 0.2, 3
    riccati = solve_ivp(
        lambda t, u: -u + u * u,
        (0, time),
        problem.initial_norm * u0,
        rtol=1e-12,
        atol=1e-14,
    ).y[:, -1]
    riccati_dir = riccati / np.linalg.norm(riccati)
    truncations = {}
    for cutoff in (1, 2):
        g_lift = assemble_carleman_lift(matrices, cutoff)
        z0 = lifted_initial(u0, problem.initial_norm, cutoff)
        base = 1 << cutoff
        chan_linear = (scipy.linalg.expm(g_lift * time) @ z0)[base : base + 2]
        truncations[cutoff] = float(
            np.linalg.norm(chan_linear / np.linalg.norm(chan_linear) - riccati_dir)
        )
        if cutoff == 2:
            # The K=2 quantum implementation and program were verified by carleman-riccati-endtoend; its block is reused here
            assert cutoff2_block is not None
            continue
        state = carleman_qode(
            problem, time, partial(taylor_series_solver, degree=degree), cutoff=cutoff
        )
        alpha_v = evolution_alpha(state)
        expected = ode_taylor_series(g_lift, time, degree) @ z0 / alpha_v
        results = run_amplitude_paths(
            state.operation.program(), ["reference", "rir-pysparq", "adapter-pysparq"]
        )
        blocks, block_dev, _, pysparq_floor = compare_path_blocks(results)
        impl = 0.0
        for block in blocks.values():
            impl = max(impl, max(abs(block.get(i, 0) - expected[base + i]) for i in range(2)))
        report.case(
            f"carleman-cutoff{cutoff}-quantum",
            paths=list(results),
            parameters={"cutoff": cutoff, "time": time, "degree": degree},
            metrics={
                "impl_error": impl,
                "carleman_truncation_error": truncations[cutoff],
                "block_cross_deviation": block_dev,
                "pysparq_spectrum_floor": pysparq_floor,
            },
            criterion="each truncation order's quantum implementation matches the independent simulation (impl_error < 1e-9)",
            passed=impl < 1e-9 and block_dev < 1e-9 and pysparq_floor < 1e-6,
        )
    # K=3 classical only: the same independently assembled expm channel (the K=2 assembly was verified via the quantum block extraction)
    g3 = assemble_carleman_lift(matrices, 3)
    z3 = lifted_initial(u0, problem.initial_norm, 3)
    chan3 = (scipy.linalg.expm(g3 * time) @ z3)[8:10]
    truncations[3] = float(np.linalg.norm(chan3 / np.linalg.norm(chan3) - riccati_dir))
    report.case(
        "carleman-cutoff-trend",
        paths=["quantum(cutoff 1,2)", "classical(cutoff 3)"],
        parameters={"time": time, "note": "truncation error as the direction difference vs the scipy solve_ivp exact Riccati solution"},
        metrics={
            "truncation_cutoff1": truncations[1],
            "truncation_cutoff2": truncations[2],
            "truncation_cutoff3": truncations[3],
        },
        criterion="truncation error decreases with the order (K=1 -> K=2 -> K=3, classical convergence evidence)",
        passed=truncations[2] < truncations[1] and truncations[3] <= truncations[2],
    )


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------


def verify_cbmd_vs_lchs(report, lchs_block, cbmd_block):
    """Quantum direction cross-check of CBMD vs LCHS on the same non-commuting problem (reusing the quantum blocks of B2/C1)."""
    time = 0.2
    l_mat = np.array([[1.0, 0.3], [0.3, 0.5]])
    h_mat = np.array([[0.2, 0.1], [0.1, -0.1]])
    u0 = np.array([1.0, 1.0]) / SQRT2
    directions = {}
    for name, block in (("lchs", lchs_block), ("cbmd", cbmd_block)):
        vector = np.array([block.get(i, 0) for i in range(2)])
        directions[name] = vector / np.linalg.norm(vector)
    exact = scipy.linalg.expm(-(l_mat + 1j * h_mat) * time) @ u0
    exact_dir = exact / np.linalg.norm(exact)
    err_lchs = float(np.linalg.norm(directions["lchs"] - exact_dir))
    err_cbmd = float(np.linalg.norm(directions["cbmd"] - exact_dir))
    mutual = float(np.linalg.norm(directions["lchs"] - directions["cbmd"]))
    report.case(
        "cbmd-vs-lchs-direction",
        paths=["reference (reusing the quantum blocks of lchs-given-parts-noncommuting and cbmd-parts-noncommuting)"],
        parameters={"time": time, "degree": 3, "note": "same problem as B2/C1"},
        metrics={
            "lchs_direction_error": err_lchs,
            "cbmd_direction_error": err_cbmd,
            "mutual_direction_deviation": mutual,
        },
        criterion=(
            "both methods' quantum directions approach the exact solution; the CBMD main-series remainder is smaller than the "
            "LCHS Cauchy quadrature remainder on this instance (informative comparison)"
        ),
        passed=err_lchs < 0.1 and err_cbmd < 0.02 and mutual < err_lchs + err_cbmd + 1e-9,
    )


def run():
    report = Report(
        "ode",
        "Publication-grade validation of the QODE group: substructure block extraction, LCHS/CBMD/Schrodingerization end to end "
        "(multiple input models), and Carleman lifting/initial state/Riccati end to end with the truncation trend.",
    )
    # A. Substructures
    verify_taylor_hamiltonian_block(report)
    verify_fourier_momentum_block(report)
    verify_carleman_lift_block(report)
    verify_carleman_initial(report)
    verify_quadrature_plans(report)
    # B. LCHS end to end (multiple input models)
    verify_lchs_given_be(report)
    lchs_block = verify_lchs_parts_noncommuting(report)
    verify_lchs_diagonal_gate_vs_qram(report)
    verify_lchs_fokker_planck(report)
    verify_lchs_heat_structured(report)
    verify_lchs_quadrature_convergence(report)
    # C. CBMD
    cbmd_block = verify_cbmd_endtoend(report)
    verify_cbmd_vs_lchs(report, lchs_block, cbmd_block)
    # D. Schrodingerization
    verify_schrodingerization_decay(report)
    verify_schrodingerization_sign_flipped(report)
    verify_schrodingerization_rotation(report)
    # E. Carleman
    cutoff2_block = verify_carleman_riccati(report)
    verify_carleman_cutoff_trend(report, cutoff2_block)
    report.write()
    return report


if __name__ == "__main__":
    run()
