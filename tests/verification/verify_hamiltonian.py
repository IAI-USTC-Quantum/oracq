"""Publication-grade numerical validation of Hamiltonian simulation and standard QSVT transforms.

Covers src/oracq/algorithms/hamiltonian.py and src/oracq/algorithms/qsvt.py:

- trotter_hamsim / hamiltonian_simulation: circuit unitaries (UniQC
  ``Circuit.to_matrix``) compared element-wise against independent classical
  matrices from the product formula (single-term exact, identity term global
  phase, 1-4 qubits); the convergence order of first-order Lie-Trotter (error
  proportional to 1/steps) is fitted with scipy.linalg.expm as the independent
  oracle.
- qsp_phases / qsvt_sequence / qubitization_walk: phase-synthesis round trip
  (QSP response independently implemented with numpy per the documented
  convention), circuit zero-signal blocks vs the target polynomial pointwise
  over four backend paths, and the W^n zero-signal block of qubitization =
  T_n(x) recurrence witness.
- qsvt_hamiltonian_simulation / qsvt_matrix_inversion / eigenstate_filter /
  gibbs_purification (density.py, QSVT consumers): zero-signal blocks compared
  against scipy expm / numpy.linalg.inv / the analytic filter polynomial /
  the classical Gibbs state; reporting implementation vs method error,
  recovered-inverse condition-number comparison, pass/stop-band suppression,
  trace distance, and success probability.
- oblivious_amplification: the library operator W = U*[R U^dagger R U] (a
  historical defect is fixed: the original iterate [R U^dagger R U] lacked the
  trailing U, and its zero-signal block degenerated to 2B^dagger B - I,
  performing no amplification). After the fix the zero-signal block satisfies
  the Chebyshev amplification identity Pi W Pi = B(4B^dagger B - 3I) (verified
  exactly); the script independently assembles the literature sequence
  U R U^dagger R U and asserts the library operator matches it amplitude by
  amplitude (a regression pin); the V/2 -> -V amplification semantics appear
  in the numerical-validation section of
  docs/manual/algorithms/oblivious-amplification.md.

All classical oracles are independent: numpy/scipy/math closed forms and
matrix routines are constructed directly in this script without reusing
internal helpers of the modules under test.

Run: PYTHONPATH=src <python with pysparq+uniqc> tests/verification/verify_hamiltonian.py
"""

from __future__ import annotations

import cmath
import math
import random

import numpy as np
from harness import (
    Report,
    adapter_pysparq,
    amplitudes_to_statevector,
    basis_program,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    statevector_error,
    superposition_program,
)
from scipy.linalg import expm as scipy_expm
from scipy.special import jv as bessel_j

from oracq.algorithms.common.hamiltonian import (
    EncodedOperator,
    PauliHamiltonian,
    hamiltonian_simulation,
    trotter_hamsim,
)
from oracq.algorithms.common.qsvt import (
    eigenstate_filter,
    qsp_phases,
    qsvt_hamiltonian_simulation,
    qsvt_matrix_inversion,
)
from oracq.algorithms.common.transforms import (
    oblivious_amplification,
    qsvt_sequence,
    qubitization_walk,
)
from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding, reflect_zero
from oracq.algorithms.input_model.density import gibbs_purification
from oracq.algorithms.input_model.operators import block_encoding, linear_combination, zero
from oracq.algorithms.input_model.oracles import invoke, resources_for
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits
from oracq.infrastructure.layout import workspace_table

# ---------------------------------------------------------------------------
# Independent classical oracles: Pauli matrices, product formula, QSP response, target polynomials.
# ---------------------------------------------------------------------------

_PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_matrix(word):
    """Dense matrix of a Pauli word; word[i] acts on qubit i (weight 2^i, matching the RIR register)."""
    result = _PAULI[word[0]]
    for letter in word[1:]:
        result = np.kron(_PAULI[letter], result)
    return result


def hamiltonian_matrix(terms):
    return sum(coefficient * pauli_matrix(word) for coefficient, word in terms)


def term_unitary(coefficient, word, time):
    """Closed form of exp(-i*c*t*P): since P^2 = I it is cos(ct)*I - i*sin(ct)*P."""
    theta = coefficient * time
    dim = 1 << len(word)
    return math.cos(theta) * np.eye(dim) - 1j * math.sin(theta) * pauli_matrix(word)


def product_formula(terms, final_time, steps):
    """Classical matrix of the first-order Lie-Trotter product formula (term order matches the circuit application order)."""
    dim = 1 << len(terms[0][1])
    step = np.eye(dim, dtype=complex)
    for coefficient, word in terms:
        step = term_unitary(coefficient, word, final_time / steps) @ step
    return np.linalg.matrix_power(step, steps)


def qsp_response_ref(x, phases):
    """QSP response independently implemented with numpy, following the reflection convention of the module docstring.

    p(x) = [S(phi_0) W(x) S(phi_1) W(x) ... S(phi_d)]_00, S(phi) = diag(e^{i*phi}, e^{-i*phi}),
    W(x) = [[x, s], [s, -x]]. Independent of the qsvt.qsp_response under test.
    """
    s = math.sqrt(max(0.0, 1.0 - x * x))
    walk = np.array([[x, s], [s, -x]])
    total = np.eye(2, dtype=complex)
    for i, phi in enumerate(phases):
        if i:
            total = walk @ total
        total = np.diag([cmath.exp(1j * phi), cmath.exp(-1j * phi)]) @ total
    return total[0, 0]


def chebyshev_coeffs(n):
    """Ascending monomial coefficients of T_n (triple recursion, independently implemented)."""
    if n == 0:
        return (1.0,)
    if n == 1:
        return (0.0, 1.0)
    a, b = (1.0,), (0.0, 1.0)
    for _ in range(2, n + 1):
        a, b = b, tuple(
            2.0 * (b[i - 1] if i else 0.0) - (a[i] if i < len(a) else 0.0)
            for i in range(len(b) + 1)
        )
    return b


def poly_eval(coeffs, x):
    result = 0.0
    for c in reversed(coeffs):
        result = result * x + c
    return result


def matrix_polynomial(coeffs, mat):
    """Matrix Horner evaluation of a polynomial with ascending coefficients (independent oracle)."""
    result = np.zeros_like(mat)
    for c in reversed(coeffs):
        result = result @ mat + c * np.eye(mat.shape[0])
    return result


def chebyshev_value(n, x):
    """T_n(x), using the hyperbolic continuation when |x| > 1."""
    if abs(x) <= 1.0:
        return math.cos(n * math.acos(x))
    return math.cosh(n * math.acosh(abs(x))) * (1.0 if x > 0 or n % 2 == 0 else -1.0)


def filter_polynomial(gap, degree, x):
    """Independent evaluation of the Lin-Tong filter polynomial f(x) = T_d(g(x^2))/T_d(r)."""
    r = (1.0 + gap**2) / (1.0 - gap**2)
    g = 2.0 * (x * x - gap**2) / (1.0 - gap**2) - 1.0
    return chebyshev_value(degree, g) / math.cosh(degree * math.acosh(r))


def jacobi_anger_tail(t, kmax):
    """Truncated Jacobi-Anger tail bound 2*sum_{j>K} |J_j(t)| of e^{itx} (independently evaluated with scipy)."""
    return float(2.0 * sum(abs(bessel_j(j, t)) for j in range(kmax + 1, kmax + 200)))


def gibbs_reference(matrix, beta):
    """Classical Gibbs state e^{-beta*H}/Tr (scipy expm, independent oracle)."""
    weights = scipy_expm(-beta * np.asarray(matrix, dtype=complex))
    return weights / np.trace(weights)


def trace_distance(rho, sigma):
    return 0.5 * float(np.abs(np.linalg.eigvalsh(rho - sigma)).sum())


# ---------------------------------------------------------------------------
# Backend execution helpers: column-readout zero-signal block, budget precheck, unified entry decoding.
# ---------------------------------------------------------------------------


def budget_qubits(program):
    """Registers + workspace total (OriginIR state-vector budget is 24 qubits)."""
    return sum(r.type.width for r in program.main.registers) + workspace_table(program)[
        program.entry
    ]


def embed_unitary(circuit, dim):
    """Embed the UniQC to_matrix unitary into dimension dim.

    UniQC sizes the matrix by the highest referenced qubit and truncates
    trailing qubits that carry no gates (e.g. an I trailing a Pauli word); the
    truncated part acts as the identity, so a right Kronecker identity factor
    is applied. In this validation all used qubits of every program form a
    prefix, so the embedding is exact.
    """
    current = circuit.shape[0]
    if current == dim:
        return circuit
    if current > dim or dim % current:
        raise AssertionError(f"unitary dimension {current} cannot be embedded into {dim}")
    return np.kron(np.eye(dim // current, dtype=complex), circuit)


def _register_widths(operation):
    return [r.type.width for r in operation.module.registers]


def _decode_index(index, widths):
    values = []
    for w in widths:
        values.append(index & ((1 << w) - 1))
        index >>= w
    return tuple(values)


def _iter_entries(state, widths):
    """Unify dict sparse states and dense originir state vectors into (register-value tuple, amplitude) iteration."""
    if isinstance(state, dict):
        yield from state.items()
    else:
        for index, amplitude in enumerate(state):
            if amplitude:
                yield _decode_index(index, widths), amplitude


def run_column(operation, runner, column):
    """Prepare the target basis-state column with X gates, then execute (all four paths use the driver uniformly, starting from |0>)."""
    return runner(basis_program(operation, {"target": column}, name=f"column_{column}"))


def zero_signal_block(operation, runner, dim):
    """Read out the target block of the signal == 0 branch column by column, plus the per-column success probability."""
    widths = _register_widths(operation)
    block = np.zeros((dim, dim), dtype=complex)
    success = []
    for column in range(dim):
        state = run_column(operation, runner, column)
        total = 0.0
        for key, amplitude in _iter_entries(state, widths):
            if key[1] == 0:
                block[key[0], column] = amplitude
                total += abs(amplitude) ** 2
        success.append(total)
    return block, success


# ---------------------------------------------------------------------------
# A. Trotter circuit semantics: unitary vs independent product-formula matrices.
# ---------------------------------------------------------------------------


def verify_trotter_single_term(report):
    """Single-term Pauli evolution (steps=1, no Trotter error): unitary matches the closed form element-wise."""
    instances = [
        (0.7, "X", 0.4),
        (-0.3, "Y", 1.1),
        (0.55, "Z", -0.8),
        (0.25, "I", 0.9),
        (0.9, "ZZ", 0.6),
        (-0.4, "XY", 1.2),
        (0.35, "YZI", -0.5),
        (-0.8, "ZXI", 0.7),
        (0.45, "XYZI", 0.3),
        (-0.6, "IYXZ", -1.0),
    ]
    worst = 0.0
    per_word = {}
    for coefficient, word, time in instances:
        circuit = embed_unitary(
            originir_unitary(trotter_hamsim([(coefficient, word)], time, steps=1).program()),
            1 << len(word),
        )
        error = float(np.abs(circuit - term_unitary(coefficient, word, time)).max())
        per_word[f"{coefficient}·{word}·t={time}"] = error
        worst = max(worst, error)
    report.case(
        "trotter-single-term-exact",
        paths=["originir-ext+to_matrix"],
        parameters={"instances": len(instances)},
        metrics={"max_error": worst, "per_word": per_word},
        criterion="single-term evolution unitary equals cos(ct)I - i*sin(ct)P (max_error < 1e-12)",
        passed=worst < 1e-12,
    )


def verify_trotter_product_formula(report):
    """Multi-term non-commuting decompositions: unitary matches the product-formula classical matrix element-wise over step/time sweeps."""
    instances = [
        (
            "1q",
            [(0.3, "I"), (0.7, "X"), (-0.2, "Z")],
            [1, 3],
            [0.4, -0.7],
        ),
        (
            "2q",
            [(0.9, "ZZ"), (0.6, "XI"), (0.45, "IX"), (0.25, "II")],
            [1, 2, 5],
            [0.4, 1.3],
        ),
        (
            "3q",
            [(0.9, "ZZI"), (0.9, "IZZ"), (0.6, "XII"), (0.6, "IXI"), (0.6, "IIX")],
            [2, 3],
            [0.4, -0.7],
        ),
        (
            "4q",
            [(0.4, "ZIZI"), (-0.7, "XIII"), (0.2, "IXII"), (0.35, "IIXI"), (-0.15, "IIIX")],
            [3],
            [0.9],
        ),
    ]
    for label, terms, steps_list, times in instances:
        worst = 0.0
        gap = 0.0
        exact = scipy_expm(-1j * hamiltonian_matrix(terms) * times[0])
        for steps in steps_list:
            for time in times:
                program = trotter_hamsim(terms, time, steps=steps).program()
                circuit = embed_unitary(originir_unitary(program), 1 << len(terms[0][1]))
                expected = product_formula(terms, time, steps)
                worst = max(worst, float(np.abs(circuit - expected).max()))
                if time == times[0]:
                    gap = max(gap, float(np.linalg.norm(circuit - exact, 2)))
        report.case(
            f"trotter-product-formula-{label}",
            paths=["originir-ext+to_matrix"],
            parameters={
                "terms": len(terms),
                "steps": steps_list,
                "times": times,
                "qubits": budget_qubits(program),
            },
            metrics={"max_error": worst, "trotter_gap_vs_expm": gap},
            criterion="circuit unitary equals the product-formula matrix (max_error < 1e-12)",
            passed=worst < 1e-12,
        )


def verify_trotter_superposition_cross(report):
    """3-qubit TFIM on uniform-superposition input, four-path cross-check (one run covers all basis states)."""
    terms = [(0.9, "ZZI"), (0.9, "IZZ"), (0.6, "XII"), (0.6, "IXI"), (0.6, "IIX")]
    operation = trotter_hamsim(terms, 0.9, steps=3)
    program = superposition_program(operation, ["target"])
    expected = product_formula(terms, 0.9, 3) @ np.full(8, 1.0 / math.sqrt(8))
    widths = [3, 0]
    deviations = {}
    for name, runner in (
        ("reference", reference),
        ("rir-pysparq", rir_pysparq),
        ("adapter-pysparq", adapter_pysparq),
        ("originir-ext", originir_ext),
    ):
        state = runner(program)
        vector = (
            amplitudes_to_statevector(state, widths) if isinstance(state, dict) else list(state)
        )
        deviations[name] = float(statevector_error(vector, list(expected)))
    worst = max(deviations.values())
    report.case(
        "trotter-superposition-cross-3q",
        paths=list(deviations),
        parameters={"steps": 3, "time": 0.9, "qubits": budget_qubits(program)},
        metrics={"max_deviation": worst, "per_path": deviations},
        criterion="superposition state agrees with the product formula amplitude by amplitude on all four paths (max_deviation < 1e-9)",
        passed=worst < 1e-9,
    )


def _convergence_scan(terms, time, grid):
    """Sweep the step count r for fixed t; returns the spectral-norm error sequence (unitary path, vs scipy expm)."""
    exact = scipy_expm(-1j * hamiltonian_matrix(terms) * time)
    errors = []
    for steps in grid:
        circuit = embed_unitary(
            originir_unitary(trotter_hamsim(terms, time, steps=steps).program()),
            1 << len(terms[0][1]),
        )
        errors.append(float(np.linalg.norm(circuit - exact, 2)))
    return errors


def verify_trotter_convergence(report):
    """First-order Lie-Trotter convergence-order fit: the log-log slope should be close to -1 (error proportional to t^2/r)."""
    instances = [
        ("2q", [(0.9, "ZZ"), (0.6, "XI"), (0.45, "IX")], 1.0, [1, 2, 4, 8, 16, 32, 64]),
        (
            "3q",
            [(0.9, "ZZI"), (0.9, "IZZ"), (0.6, "XII"), (0.6, "IXI"), (0.6, "IIX")],
            1.0,
            [1, 2, 4, 8, 16, 32],
        ),
    ]
    for label, terms, time, grid in instances:
        errors = _convergence_scan(terms, time, grid)
        fit_r = np.log([r for r in grid if r >= 2])
        fit_e = np.log([e for r, e in zip(grid, errors, strict=True) if r >= 2])
        slope = float(np.polyfit(fit_r, fit_e, 1)[0])
        report.case(
            f"trotter-convergence-{label}",
            paths=["originir-ext+to_matrix"],
            parameters={"time": time, "steps_grid": grid, "terms": len(terms)},
            metrics={
                "errors": {str(r): e for r, e in zip(grid, errors, strict=True)},
                "fitted_order": -slope,
                "error_at_max_steps": errors[-1],
            },
            criterion="fitted convergence order within [0.8, 1.3] (theoretical value 1 for a first-order product formula)",
            passed=0.8 <= -slope <= 1.3,
        )


# ---------------------------------------------------------------------------
# B. hamiltonian_simulation protocol layer (Trotter routing and QSP injection).
# ---------------------------------------------------------------------------


def verify_protocol_trotter(report):
    """PauliHamiltonian routed through the protocol to Trotter: unitary exactly equals the product formula."""
    terms = [(0.9, "ZZ"), (0.6, "XI"), (0.45, "IX")]
    operator = PauliHamiltonian(terms)
    exact = scipy_expm(-1j * hamiltonian_matrix(terms))
    worst_formula = 0.0
    gaps = {}
    for steps in (2, 16):
        result = hamiltonian_simulation(operator, 1.0, steps=steps)
        circuit = embed_unitary(originir_unitary(result.operation.program()), 4)
        worst_formula = max(
            worst_formula, float(np.abs(circuit - product_formula(terms, 1.0, steps)).max())
        )
        gaps[str(steps)] = float(np.linalg.norm(circuit - exact, 2))
    report.case(
        "hamsim-protocol-trotter-2q",
        paths=["originir-ext+to_matrix"],
        parameters={"time": 1.0, "steps": [2, 16], "alpha": 1.0},
        metrics={
            "max_error_vs_product_formula": worst_formula,
            "trotter_gap_vs_expm": gaps,
        },
        criterion="protocol-assembled unitary equals the product formula (< 1e-12) and the expm deviation at r=16 is smaller than at r=2",
        passed=worst_formula < 1e-12 and gaps["16"] < gaps["2"],
    )


def verify_protocol_qsp_injection(report):
    """EncodedOperator routed via auto to an injected QSVT kernel: zero-signal block vs scipy expm."""
    terms = [(0.6, "Z"), (0.4, "X")]
    ham = PauliHamiltonian(terms)
    operator = EncodedOperator(ham.block_encoding(), True)
    time = 0.5

    def injected_qsp(be, evolution_time):
        # taking tau = -t*alpha in e^{i*tau*A/alpha} yields e^{-i*t*H}
        return qsvt_hamiltonian_simulation(be, -evolution_time * be.alpha, error=0.02)

    result = hamiltonian_simulation(operator, time, method="auto", qsp=injected_qsp)
    scale = dict(result.operation.module.attributes)["sim_scale"]
    expected = scipy_expm(-1j * hamiltonian_matrix(terms) * time) / scale
    worst = 0.0
    success = []
    for runner in (reference, rir_pysparq, originir_ext):
        block, probs = zero_signal_block(result.operation, runner, 2)
        worst = max(worst, float(np.abs(block - expected).max()))
        success = probs
    report.case(
        "hamsim-protocol-qsp-injection",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"time": time, "qsp_error": 0.02, "sim_scale": scale},
        metrics={"max_error": worst, "success_probability": success},
        criterion="zero-signal block approximates e^{-itH}/sim_scale (max_error < 3e-2, QSVT approximation budget included)",
        passed=worst < 3e-2,
    )


# ---------------------------------------------------------------------------
# C. QSP phase synthesis, QSVT sequence convention, and qubitization.
# ---------------------------------------------------------------------------


def verify_phase_synthesis(report):
    """qsp_phases round trip: synthesized phases fed through the independent response recover the target polynomial (401-point grid)."""
    grid = [-1.0 + 2.0 * i / 400 for i in range(401)]
    targets = {f"chebyshev-T{n}": (chebyshev_coeffs(n), None) for n in range(1, 7)}
    targets["explicit-imag"] = ((0.0, 0.5), (0.0, math.sqrt(0.75)))
    errors = {}
    for name, (coeffs, imag) in targets.items():
        phases = qsp_phases(coeffs, imag=imag)
        errors[name] = float(
            max(
                abs(
                    qsp_response_ref(x, phases)
                    - poly_eval(coeffs, x)
                    - (1j * poly_eval(imag, x) if imag else 0.0)
                )
                for x in grid
            )
        )
    worst = max(errors.values())
    report.case(
        "qsp-phase-synthesis-roundtrip",
        paths=["numpy-independent-response"],
        parameters={"grid_points": len(grid), "targets": len(targets)},
        metrics={"max_error": worst, "per_target": errors},
        criterion="synthesized-phase round-trip error < 1e-8 (module claims a typical 1e-9 order)",
        passed=worst < 1e-8,
    )


def verify_qsvt_sequence_convention(report):
    """qsvt_sequence circuit zero-signal block vs the independent response under random phases, pointwise (four paths)."""
    rng = random.Random(20240901)
    phases = tuple(rng.uniform(-math.pi, math.pi) for _ in range(6))
    matrix = [[0.85, 0.0], [0.0, -0.55]]
    be = matrix_pauli_encoding(matrix)
    operation = qsvt_sequence(be, phases)
    xs = [matrix[0][0] / be.alpha, matrix[1][1] / be.alpha]
    deviations = {}
    for name, runner in (
        ("reference", reference),
        ("rir-pysparq", rir_pysparq),
        ("adapter-pysparq", adapter_pysparq),
        ("originir-ext", originir_ext),
    ):
        block, _ = zero_signal_block(operation, runner, 2)
        deviations[name] = float(
            max(abs(block[col, col] - qsp_response_ref(xs[col], phases)) for col in (0, 1))
        )
    worst = max(deviations.values())
    report.case(
        "qsvt-sequence-convention-1q",
        paths=list(deviations),
        parameters={"phases": len(phases), "spectral_points": xs},
        metrics={"max_deviation": worst, "per_path": deviations},
        criterion="zero-signal block on all four paths matches the independent QSP response (max_deviation < 1e-9)",
        passed=worst < 1e-9,
    )


def verify_qsvt_sequence_matrix_block(report):
    """2-qubit non-diagonal BE: the full zero-signal block of qsvt_sequence equals the matrix polynomial T_4(A/alpha)."""
    terms = [(0.45, "ZZ"), (0.25, "XI"), (-0.15, "IZ"), (0.10, "II")]
    matrix = hamiltonian_matrix(terms)
    be = matrix_pauli_encoding(matrix)
    coeffs = chebyshev_coeffs(4)
    operation = qsvt_sequence(be, qsp_phases(coeffs))
    expected = matrix_polynomial(coeffs, matrix / be.alpha)
    worst = 0.0
    used = []
    for name, runner in (
        ("reference", reference),
        ("rir-pysparq", rir_pysparq),
        ("originir-ext", originir_ext),
    ):
        block, _ = zero_signal_block(operation, runner, 4)
        worst = max(worst, float(np.abs(block - expected).max()))
        used.append(name)
    report.case(
        "qsvt-sequence-matrix-block-2q",
        paths=used,
        parameters={"degree": 4, "alpha": be.alpha, "qubits": 4},
        metrics={"max_error": worst},
        criterion="zero-signal block equals the matrix polynomial T_4(A/alpha) (max_error < 1e-9)",
        passed=worst < 1e-9,
    )


def verify_qubitization(report):
    """qubitization_walk: the zero-signal block of W^n equals T_n(x) (n = 1..5, two eigenstates)."""
    matrix = [[0.85, 0.0], [0.0, -0.55]]
    be = matrix_pauli_encoding(matrix)
    walk = qubitization_walk(be)
    xs = [matrix[0][0] / be.alpha, matrix[1][1] / be.alpha]
    worst = 0.0
    for n in (1, 2, 3, 5):
        for col in (0, 1):
            b = Builder(
                f"walk_pow_{n}_{col}",
                {r.name: r.type for r in walk.module.registers},
                resources_for(("w", walk)),
            )
            for bit in range(b["target"].width):
                if (col >> bit) & 1:
                    b.x(b["target"][bit])
            with b.repeat(n):
                invoke(b, walk, "w", target=b["target"], signal=b["signal"])
            program = b.finish().program()
            widths = [r.type.width for r in program.main.registers]
            expected = chebyshev_value(n, max(-1.0, min(1.0, xs[col])))
            for runner in (reference, rir_pysparq, originir_ext):
                state = runner(program)
                amplitude = next(
                    (a for k, a in _iter_entries(state, widths) if k[0] == col and k[1] == 0),
                    0j,
                )
                worst = max(worst, abs(amplitude - expected))
    report.case(
        "qubitization-chebyshev-recurrence",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"powers": [1, 2, 3, 5], "spectral_points": xs},
        metrics={"max_error": worst},
        criterion="zero-signal block of W^n equals T_n(x) (max_error < 1e-12)",
        passed=worst < 1e-12,
    )


# ---------------------------------------------------------------------------
# D. Standard QSVT transforms: HamSim, matrix inversion, eigenstate filtering.
# ---------------------------------------------------------------------------


def verify_qsvt_hamsim(report):
    """qsvt_hamiltonian_simulation: zero-signal block vs scipy expm; method error and implementation error listed separately."""
    matrix = np.array([[0.5, 0.2], [0.2, -0.3]])
    be = matrix_pauli_encoding(matrix)
    for time in (0.7, 2.0):
        result = qsvt_hamiltonian_simulation(be, time, error=0.01)
        attrs = dict(result.operation.module.attributes)
        scale, degree = attrs["sim_scale"], attrs["qsp_degree"]
        expected = scipy_expm(1j * time * matrix / be.alpha) / scale
        worst = 0.0
        for runner in (reference, rir_pysparq, originir_ext):
            block, _ = zero_signal_block(result.operation, runner, 2)
            worst = max(worst, float(np.abs(block - expected).max()))
        tail = jacobi_anger_tail(time, degree)
        report.case(
            f"qsvt-hamsim-t{time}",
            paths=["reference", "rir-pysparq", "originir-ext"],
            parameters={"time": time, "error": 0.01, "qsp_degree": degree, "sim_scale": scale},
            metrics={
                "impl_error": worst,
                "method_tail_bound": tail,
                "total_error_bound": worst + tail,
            },
            criterion="total deviation of block*sim_scale from e^{itA/alpha} < 2x error (0.02)",
            passed=worst + tail < 0.02,
        )


def _inversion_case(report, name, matrix, kappa, error, runners):
    be = matrix_pauli_encoding(matrix)
    result = qsvt_matrix_inversion(be, kappa, error=error)
    attrs = dict(result.operation.module.attributes)
    scale, degree = attrs["inverse_scale"], attrs["qsp_degree"]
    matrix = np.asarray(matrix, dtype=complex)
    exact_inverse = np.linalg.inv(matrix)
    expected = scale * exact_inverse
    dim = matrix.shape[0]
    worst = 0.0
    success = []
    for runner in runners:
        block, probs = zero_signal_block(result.operation, runner, dim)
        worst = max(worst, float(np.abs(block - expected).max()))
        success = probs
    recovered = np.linalg.svd(block / scale, compute_uv=False)
    kappa_recovered = float(recovered[0] / recovered[-1])
    relative = float(
        np.linalg.norm(block / scale - exact_inverse, 2) / np.linalg.norm(exact_inverse, 2)
    )
    # Each singular value carries a relative deviation <= error, so the recovered condition number lies in kappa*[(1-e)/(1+e), (1+e)/(1-e)]
    kappa_exact = float(np.linalg.cond(matrix))
    band = (kappa_exact * (1 - error) / (1 + error), kappa_exact * (1 + error) / (1 - error))
    report.case(
        name,
        paths=[r.__name__ for r in runners],
        parameters={"kappa": kappa, "error": error, "qsp_degree": degree, "inverse_scale": scale},
        metrics={
            "max_error": worst,
            "relative_error": relative,
            "kappa_recovered": kappa_recovered,
            "kappa_exact": kappa_exact,
            "kappa_band": list(band),
            "success_probability": success,
        },
        criterion="relative spectral error of the recovered inverse <= 1.05x error and the recovered condition number lies in the kappa*(1+/-e)/(1-/+e) band",
        passed=relative <= 1.05 * error and band[0] <= kappa_recovered <= band[1],
    )


def verify_qsvt_inversion(report):
    """qsvt_matrix_inversion: block/scale vs the numpy exact inverse; condition-number comparison; error sweep."""
    matrix_k2 = [[0.6, -0.2], [-0.2, 0.6]]  # eigenvalues 0.4/0.8, kappa = 2
    for error in (0.15, 0.10):
        _inversion_case(
            report,
            f"qsvt-inversion-kappa2-error{error}",
            matrix_k2,
            2.0,
            error,
            (reference, rir_pysparq),
        )
    _inversion_case(
        report,
        "qsvt-inversion-kappa3-1q",
        [[0.4, -0.2], [-0.2, 0.4]],  # eigenvalues 0.2/0.6, kappa = 3
        3.0,
        0.4,
        (reference, rir_pysparq, originir_ext),
    )
    _inversion_case(
        report,
        "qsvt-inversion-kappa3-2q",
        [[0.3, 0, 0, 0], [0, 0.4, 0, 0], [0, 0, 0.6, 0], [0, 0, 0, 0.9]],
        3.0,
        0.4,
        (reference, rir_pysparq, originir_ext),
    )


def verify_eigenstate_filter(report):
    """eigenstate_filter: zero-signal-block diagonal magnitudes vs the analytic filter polynomial; pass/stop-band suppression."""
    # 1-qubit band-center shift (docs case): after shifting the spectral variable, one pass and one stop.
    be = matrix_pauli_encoding([[0.05, 0.0], [0.0, 0.5]])
    gap, degree, center = 0.2, 8, 0.1
    flt = eigenstate_filter(be, gap, degree, center=center)
    attrs = dict(flt.operation.module.attributes)
    alpha_shifted = be.alpha + abs(center)
    xs = [(0.05 - center) / alpha_shifted, (0.5 - center) / alpha_shifted]
    expected = [filter_polynomial(gap, degree, x) for x in xs]
    worst = 0.0
    for runner in (reference, rir_pysparq, originir_ext):
        block, _ = zero_signal_block(flt.operation, runner, 2)
        worst = max(worst, float(max(abs(block[c, c] - expected[c]) for c in (0, 1))))
    report.case(
        "eigenstate-filter-centered-1q",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"gap": gap, "degree": degree, "center": center, "shifted_x": xs},
        metrics={
            "max_error": worst,
            "passband_amplitude": expected[0],
            "stopband_amplitude": expected[1],
            "suppression_attribute": attrs["suppression"],
        },
        criterion="block diagonal magnitudes equal the analytic filter polynomial (max_error < 1e-5, the stripping-noise order of a synthesis degree 16)"
        " and stopband <= suppression",
        passed=worst < 1e-5 and abs(expected[1]) <= attrs["suppression"] + 1e-9,
    )
    # 2-qubit four eigenvalues: two passbands and two stopbands.
    matrix = [[0.02, 0, 0, 0], [0, -0.03, 0, 0], [0, 0, 0.4, 0], [0, 0, 0, 0.55]]
    be2 = matrix_pauli_encoding(matrix)
    gap2, degree2 = 0.2, 6
    flt2 = eigenstate_filter(be2, gap2, degree2)
    attrs2 = dict(flt2.operation.module.attributes)
    eigenvalues = [0.02, -0.03, 0.4, 0.55]
    xs2 = [v / be2.alpha for v in eigenvalues]
    expected2 = [filter_polynomial(gap2, degree2, x) for x in xs2]
    worst2 = 0.0
    for runner in (reference, rir_pysparq, originir_ext):
        block, _ = zero_signal_block(flt2.operation, runner, 4)
        worst2 = max(worst2, float(max(abs(block[c, c] - expected2[c]) for c in range(4))))
    stopband = max(abs(v) for v in expected2[2:])
    report.case(
        "eigenstate-filter-2q",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"gap": gap2, "degree": degree2, "alpha": be2.alpha, "spectral_x": xs2},
        metrics={
            "max_error": worst2,
            "passband_min": min(abs(v) for v in expected2[:2]),
            "stopband_max": stopband,
            "suppression_attribute": attrs2["suppression"],
        },
        criterion="all four eigenstate magnitudes match the analytic polynomial (max_error < 1e-5, stripping self-check tolerance)"
        " and stopband <= suppression",
        passed=worst2 < 1e-5 and stopband <= attrs2["suppression"] + 1e-9,
    )


# ---------------------------------------------------------------------------
# E. Oblivious amplitude amplification and Gibbs-state preparation.
# ---------------------------------------------------------------------------


def verify_oaa_iterate_identity(report):
    """Algebraic identity of the library operator W = U*[R U^dagger R U]: Pi W Pi = B(4B^dagger B - 3I) (exact for any BE).

    A historical defect is fixed: the zero-signal block of the original iterate
    [R U^dagger R U] degenerated to 2B^dagger B - I, performing none of the OAA
    amplification from the literature (see the group report and the
    oaa-standard-sequence case); the fixed operator adds the trailing U and its
    zero-signal block is the Chebyshev amplification B(4B^dagger B - 3I) --
    which is exactly -V for an input whose zero-signal block is V/2. This case
    pins the semantics the library function actually implements to machine
    precision.
    """
    matrix = np.array([[0.6, -0.2], [-0.2, 0.6]])
    be = matrix_pauli_encoding(matrix)
    source_block = matrix / be.alpha
    expected = source_block @ (4.0 * source_block.conj().T @ source_block - 3.0 * np.eye(2))
    oaa = oblivious_amplification(be, iterations=1)
    worst = 0.0
    for runner in (reference, rir_pysparq, originir_ext):
        block, _ = zero_signal_block(oaa, runner, 2)
        worst = max(worst, float(np.abs(block - expected).max()))
    report.case(
        "oaa-iterate-block-identity",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"alpha": be.alpha, "iterations": 1},
        metrics={"max_error": worst},
        criterion="Pi W Pi equals B(4B^dagger B - 3I) (max_error < 1e-12)",
        passed=worst < 1e-12,
    )


def verify_oaa_standard_sequence(report):
    """Standard OAA sequence U R U^dagger R U: a V/2 block encoding recovers -V after one iteration (literature semantics).

    The sequence is assembled in-script from public library components
    (invoke + reflect_zero + adjoint); the fixed library
    oblivious_amplification matches it amplitude by amplitude -- this case
    asserts both equalities at once, serving as the regression pin that
    "library == the standard three-query sequence".
    """
    vb = Builder("oaa_v", {"target": Bits(1), "signal": Bits(0)})
    vb.ry(vb["target"][0], 0.9)
    vb.rz(vb["target"][0], 0.4)
    vbe = block_encoding(vb.finish())
    theta_y, theta_z = 0.9, 0.4
    unitary_v = np.array(
        [[cmath.exp(-1j * theta_z / 2), 0], [0, cmath.exp(1j * theta_z / 2)]]
    ) @ np.array(
        [
            [math.cos(theta_y / 2), -math.sin(theta_y / 2)],
            [math.sin(theta_y / 2), math.cos(theta_y / 2)],
        ]
    )
    be = linear_combination(1.0, vbe, 1.0, zero(1))  # zero-signal block = V/2, alpha = 2
    b = Builder(
        "oaa_standard",
        {"target": Bits(1), "signal": Bits(be.signal_qubits)},
        resources_for(("a", be.operation)),
    )
    # Literature standard sequence U R U^dagger R U: matches the fixed library operator amplitude by amplitude
    invoke(b, be.operation, "a", target=b["target"], signal=b["signal"])
    reflect_zero(b, b["signal"])
    with b.adjoint():
        invoke(b, be.operation, "a", target=b["target"], signal=b["signal"])
    reflect_zero(b, b["signal"])
    invoke(b, be.operation, "a", target=b["target"], signal=b["signal"])
    operation = b.finish()
    library = oblivious_amplification(be, iterations=1)
    before, _ = zero_signal_block(be.operation, reference, 2)
    worst = 0.0
    lib_vs_script = 0.0
    measured = None
    for runner in (reference, rir_pysparq, originir_ext):
        block, _ = zero_signal_block(operation, runner, 2)
        lib_block, _ = zero_signal_block(library, runner, 2)
        measured = block
        worst = max(worst, float(np.abs(block + unitary_v).max()))
        lib_vs_script = max(lib_vs_script, float(np.abs(block - lib_block).max()))
    amplification = float(abs(measured[0, 0]) / abs(before[0, 0]))
    report.case(
        "oaa-standard-sequence-amplification",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"iterations": 1, "input_block": "V/2 (V unitary)"},
        metrics={
            "max_error_vs_minus_V": worst,
            "library_vs_script_error": lib_vs_script,
            "amplitude_before": float(abs(before[0, 0])),
            "amplitude_after": float(abs(measured[0, 0])),
            "amplification": amplification,
        },
        criterion=(
            "zero-signal block equals -V (max_error < 1e-12), amplitude amplified 0.5 -> 1.0; "
            "the library operator matches the script-assembled standard sequence amplitude by amplitude (library_vs_script_error < 1e-12)"
        ),
        passed=worst < 1e-12 and lib_vs_script < 1e-12,
    )


def _gibbs_case(report, name, matrix, beta, error, runners):
    matrix = np.asarray(matrix, dtype=complex)
    be = matrix_pauli_encoding(matrix)
    gibbs = gibbs_purification(be, beta, error=error)
    operation = gibbs.operation
    widths = _register_widths(operation)
    n_sys, n_env = widths[0], widths[1]
    expected = gibbs_reference(matrix, beta)
    worst = 0.0
    success = 0.0
    for runner in runners:
        state = runner(basis_program(operation, {}, name="gibbs_run"))
        psi = np.zeros((1 << n_sys, 1 << n_env), dtype=complex)
        total = 0.0
        for key, amplitude in _iter_entries(state, widths):
            if key[2] == 0:
                psi[key[0], key[1]] = amplitude
                total += abs(amplitude) ** 2
        rho = psi @ psi.conj().T
        rho = rho / np.trace(rho)
        worst = max(worst, trace_distance(rho, expected))
        success = total
    report.case(
        name,
        paths=[r.__name__ for r in runners],
        parameters={"beta": beta, "error": error, "qubits": sum(widths)},
        metrics={"trace_distance": worst, "success_probability": success},
        criterion=f"trace distance between the reduced density matrix and the classical Gibbs state < 3x error ({3 * error:.2f})",
        passed=worst < 3 * error,
    )


def verify_gibbs(report):
    """gibbs_purification: post-selected partial trace vs the classical Gibbs state (trace distance; diagonal/non-diagonal/beta=0)."""
    _gibbs_case(
        report,
        "gibbs-trace-distance-nondiag",
        [[0.5, 0.2], [0.2, -0.3]],
        1.2,
        0.02,
        (reference, rir_pysparq, originir_ext),
    )
    _gibbs_case(
        report,
        "gibbs-trace-distance-diag",
        [[0.8, 0.0], [0.0, -0.4]],
        0.6,
        0.02,
        (reference, rir_pysparq),
    )
    _gibbs_case(
        report,
        "gibbs-beta0-maximally-mixed",
        [[0.5, 0.2], [0.2, -0.3]],
        0.0,
        0.02,
        (reference, rir_pysparq, originir_ext),
    )


def run():
    report = Report(
        "hamiltonian",
        "Unitary-vs-product-formula cross-checks for Trotter/the protocol layer/standard QSVT transforms, convergence-order fits, "
        "pointwise polynomial comparisons, inversion condition-number checks, and Gibbs trace distances, covering 1-4 qubits.",
    )
    verify_trotter_single_term(report)
    verify_trotter_product_formula(report)
    verify_trotter_superposition_cross(report)
    verify_trotter_convergence(report)
    verify_protocol_trotter(report)
    verify_protocol_qsp_injection(report)
    verify_phase_synthesis(report)
    verify_qsvt_sequence_convention(report)
    verify_qsvt_sequence_matrix_block(report)
    verify_qubitization(report)
    verify_qsvt_hamsim(report)
    verify_qsvt_inversion(report)
    verify_eigenstate_filter(report)
    verify_oaa_iterate_identity(report)
    verify_oaa_standard_sequence(report)
    verify_gibbs(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
