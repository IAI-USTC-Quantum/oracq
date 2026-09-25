"""Publication-grade numerical validation of estimation algorithms (estimation.py / gradient.py).

Five readout algorithms are covered, each validated against an independent
classical closed form using distributions / expectation values from real
backends:

- QPE: unitaries with known eigenphases (phase gate, Fourier eigenstates of
  add_const); the phase-register histogram is compared against the closed-form
  Dirichlet-kernel distribution with a multi-precision bit sweep; additionally
  the full OriginIR unitary matrix is compared against a QPE unitary
  independently assembled with numpy.
- QAE: amplitude estimation under uniform / non-uniform preparations; the
  phase distribution is compared against the closed-form Grover bimodal
  distribution and the peak-decode error against the Brassard bound.
- Hadamard test: the Z expectation of probe is compared against Re/Im
  <psi|U|psi> computed independently from the angle table / state vectors.
- Swap test: P(probe=0) compared against (1+|<a|c>|^2)/2 with the overlap
  computed independently with numpy.
- Jordan gradient: linear functions (gate phase table and mathfunc phase
  oracle) read out exactly the gradient, cross-checked against central finite
  differences; the failure-probability decay rate of a perturbed linear
  function decreases with the grid bit count.

Run: PYTHONPATH=src <python with pysparq+uniqc> tests/verification/verify_estimation.py
"""

from __future__ import annotations

import math

from harness import (
    Report,
    adapter_pysparq,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    tvd,
)

from oracq import Bits, Builder, UInt
from oracq.algorithms.common.estimation import (
    amplitude_estimation,
    hadamard_test,
    phase_estimation,
    swap_test,
)
from oracq.algorithms.common.fourier import qft
from oracq.algorithms.input_model.oracles import gate_state_prep, uniform_state
from oracq.algorithms.optimization.gradient import (
    function_phase_oracle,
    gate_phase_oracle,
    gradient_estimation,
    gradient_from_readout,
)
from oracq.infrastructure.layout import workspace_table

# ---------------------------------------------------------------------------
# Common helpers
# ---------------------------------------------------------------------------


def _marginal_from_amplitudes(amplitudes, registers, name):
    """Marginal distribution over the specified register of a dictionary sparse state (register tuple -> amplitude)."""
    index = [r.name for r in registers].index(name)
    result = {}
    for key, amplitude in amplitudes.items():
        result[key[index]] = result.get(key[index], 0.0) + abs(amplitude) ** 2
    return result


def _marginal_from_statevector(vector, registers, name):
    """Register marginal of an OriginIR full-amplitude state vector; registers occupy the low qubits in declaration order."""
    widths = [r.type.width for r in registers]
    index = [r.name for r in registers].index(name)
    offset = sum(widths[:index])
    mask = (1 << widths[index]) - 1
    result = {}
    for basis, amplitude in enumerate(vector):
        value = (basis >> offset) & mask
        result[value] = result.get(value, 0.0) + abs(amplitude) ** 2
    return result


def _readout_marginals(program, readout, *, use_originir=True):
    """Readout-register marginal distributions across backend paths; returns {path name: distribution}."""
    marginals = {}
    for label, runner in (
        ("reference", reference),
        ("rir-pysparq", rir_pysparq),
        ("adapter-pysparq", adapter_pysparq),
    ):
        marginals[label] = _marginal_from_amplitudes(
            runner(program), program.main.registers, readout
        )
    if use_originir:
        marginals["originir-ext"] = _marginal_from_statevector(
            originir_ext(program), program.main.registers, readout
        )
    return marginals


def _dirichlet_kernel(precision, delta):
    """QPE single-peak closed form: D(delta) = sin^2(pi*2^p*delta) / (4^p*sin^2(pi*delta)), equal to 1 when delta == 0 (mod 1)."""
    points = 1 << precision
    numerator = math.sin(math.pi * points * delta)
    denominator = math.sin(math.pi * delta)
    if abs(denominator) < 1e-12:
        return 1.0
    return (numerator / (points * denominator)) ** 2


def _qpe_closed_form(phase, precision):
    """Closed-form QPE phase-register distribution for eigenphase phi."""
    points = 1 << precision
    return {y: _dirichlet_kernel(precision, y / points - phase) for y in range(points)}


def _qae_closed_form(amplitude, precision):
    """Closed-form QAE phase distribution for good-state probability a: an equal-weight mixture of two Dirichlet peaks at +/-theta/pi."""
    points = 1 << precision
    theta = math.asin(math.sqrt(amplitude))
    return {
        y: 0.5
        * (
            _dirichlet_kernel(precision, y / points - theta / math.pi)
            + _dirichlet_kernel(precision, y / points + theta / math.pi)
        )
        for y in range(points)
    }


def _distribution_checks(marginals, closed_form, extra_metrics):
    """TVD and pairwise cross-checks of the readout distribution on every path; returns the metrics dictionary and the closed-form peak."""
    worst_tvd = 0.0
    for distribution in marginals.values():
        worst_tvd = max(worst_tvd, tvd(distribution, closed_form))
    paths = list(marginals)
    pairwise = 0.0
    for i, left in enumerate(paths):
        for right in paths[i + 1 :]:
            pairwise = max(pairwise, tvd(marginals[left], marginals[right]))
    peak = max(closed_form, key=closed_form.get)
    metrics = {
        "max_tvd_vs_closed_form": worst_tvd,
        "max_pairwise_tvd": pairwise,
        "peak_value": peak,
        "peak_probability": max(d[peak] for d in marginals.values()),
        **extra_metrics,
    }
    return metrics, peak


# ---------------------------------------------------------------------------
# QPE
# ---------------------------------------------------------------------------


def _phase_unitary(phase):
    """Single-qubit diagonal phase gate U|1> = exp(2*pi*i*phase)|1>."""
    label = str(phase).replace(".", "p").replace("-", "m")
    b = Builder(f"qpe_witness_{label}", {"target": Bits(1)})
    b.gate("phase", b["target"], 2 * math.pi * phase)
    return b.finish()


def _qpe_driver(operation, precision, prepare, name):
    """Driver that prepares the eigenstate and then invokes QPE."""
    qpe = phase_estimation(operation, precision=precision)
    b = Builder(name, {r.name: r.type for r in qpe.module.registers})
    prepare(b)
    b.call(qpe, **{r.name: b[r.name] for r in qpe.module.registers})
    return b.finish().program()


def verify_qpe_ongrid(report):
    """Eigenphases on the grid: the readout must be deterministically 2^p*phi, swept over multiple precisions."""
    for precision in range(2, 7):
        phase = 0.625  # phi = 5/8, exactly on the grid for p >= 3; for p = 2 take the phi' = 0.5 branch
        if precision == 2:
            phase = 0.5
        program = _qpe_driver(
            _phase_unitary(phase), precision, lambda b: b.x(b["target"]), f"qpe_grid_{precision}"
        )
        marginals = _readout_marginals(program, "phase")
        closed = _qpe_closed_form(phase, precision)
        expected = round((1 << precision) * phase)
        metrics, peak = _distribution_checks(marginals, closed, {"expected_peak": expected})
        report.case(
            f"qpe-phase-ongrid-w{precision}",
            paths=list(marginals),
            parameters={"phase": phase, "precision": precision},
            metrics=metrics,
            criterion="peak == 2^p*phi with probability 1, distribution TVD < 1e-9",
            passed=peak == expected
            and metrics["peak_probability"] > 1 - 1e-9
            and metrics["max_tvd_vs_closed_form"] < 1e-9,
        )


def verify_qpe_offgrid(report):
    """Off-grid phase phi = 0.3: peak at the nearest grid point, probability lower bound 4/pi^2, distribution against the Dirichlet kernel."""
    phase = 0.3
    for precision in (3, 4, 5, 6):
        program = _qpe_driver(
            _phase_unitary(phase), precision, lambda b: b.x(b["target"]), f"qpe_off_{precision}"
        )
        marginals = _readout_marginals(program, "phase")
        closed = _qpe_closed_form(phase, precision)
        expected = round((1 << precision) * phase) % (1 << precision)
        metrics, peak = _distribution_checks(marginals, closed, {"expected_peak": expected})
        report.case(
            f"qpe-phase-offgrid-w{precision}",
            paths=list(marginals),
            parameters={"phase": phase, "precision": precision},
            metrics=metrics,
            criterion="peak at the nearest grid point with probability >= 4/pi^2, distribution TVD < 1e-9",
            passed=peak == expected
            and metrics["peak_probability"] >= 4 / math.pi**2 - 1e-9
            and metrics["max_tvd_vs_closed_form"] < 1e-9,
        )


def verify_qpe_fourier_eigenstate(report):
    """Non-diagonal unitary: the Fourier eigenstate |phi~_j> of add_const(1), with eigenphase (-j/4) mod 1.

    The classical reference independently builds the permutation matrix and the
    Fourier state with numpy to obtain the eigenvalue, without reusing any
    in-repo implementation.
    """
    import numpy as np

    permutation = np.zeros((4, 4), dtype=complex)
    for x in range(4):
        permutation[(x + 1) % 4, x] = 1.0
    add = Builder("qpe_add_w2", {"w": UInt(2)})
    add.add_const(add["w"], 1)
    add_operation = add.finish()
    for j, precision in ((1, 3), (3, 3), (1, 4)):
        fourier_state = np.exp(2j * math.pi * j * np.arange(4) / 4) / 2
        eigenvalue = fourier_state.conjugate() @ permutation @ fourier_state
        phase = (math.atan2(eigenvalue.imag, eigenvalue.real) / (2 * math.pi)) % 1.0

        def prepare(b, j=j):
            for bit in range(2):
                if (j >> bit) & 1:
                    b.x(b["w"][bit])
            b.call(qft(2), target=b["w"].reinterpret("bits"))

        program = _qpe_driver(add_operation, precision, prepare, f"qpe_fourier_{j}_{precision}")
        marginals = _readout_marginals(program, "phase")
        closed = _qpe_closed_form(phase, precision)
        expected = round((1 << precision) * phase) % (1 << precision)
        metrics, peak = _distribution_checks(
            marginals, closed, {"expected_peak": expected, "eigenvalue_phase": phase}
        )
        report.case(
            f"qpe-addconst-fourier-j{j}-w{precision}",
            paths=list(marginals),
            parameters={"fourier_index": j, "precision": precision},
            metrics=metrics,
            criterion="deterministic readout of the numpy eigenphase (peak probability 1, TVD < 1e-9)",
            passed=peak == expected
            and metrics["peak_probability"] > 1 - 1e-9
            and metrics["max_tvd_vs_closed_form"] < 1e-9,
        )


def verify_qpe_unitary_matrix(report):
    """Unitary level: full OriginIR unitary compared against a QPE unitary independently assembled with numpy (preparation X included)."""
    import numpy as np

    phase, precision = 0.625, 3
    points = 1 << precision
    program = _qpe_driver(
        _phase_unitary(phase), precision, lambda b: b.x(b["target"]), "qpe_matrix"
    )
    actual = originir_unitary(program)
    # Independent numpy assembly: H^tensor(p) (phase register) -> controlled U^z -> inverse DFT; qubit 0 is the target
    u = np.diag([1.0, np.exp(2j * math.pi * phase)])
    identity2 = np.eye(2)
    hadamard_full = np.empty((points, points))
    for y in range(points):
        for z in range(points):
            hadamard_full[y, z] = (-1) ** bin(y & z).count("1") / math.sqrt(points)
    controlled = np.zeros((2 * points, 2 * points), dtype=complex)
    power = np.eye(2, dtype=complex)
    for z in range(points):
        controlled[2 * z : 2 * z + 2, 2 * z : 2 * z + 2] = power
        power = power @ u
    inverse_dft = np.empty((points, points), dtype=complex)
    for y in range(points):
        for z in range(points):
            inverse_dft[y, z] = np.exp(-2j * math.pi * y * z / points) / math.sqrt(points)
    expected = (
        np.kron(inverse_dft, identity2)
        @ controlled
        @ np.kron(hadamard_full, identity2)
        @ np.kron(np.eye(points), np.array([[0, 1], [1, 0]]))
    )
    error = float(np.abs(actual - expected).max())
    report.case(
        "qpe-unitary-matrix-w3",
        paths=["originir-ext+to_matrix"],
        parameters={"phase": phase, "precision": precision, "matrix_dim": actual.shape[0]},
        metrics={"max_error": error},
        criterion="circuit unitary agrees element-wise with the independent numpy assembly (max_error < 1e-12)",
        passed=error < 1e-12,
    )


# ---------------------------------------------------------------------------
# QAE
# ---------------------------------------------------------------------------


def _brassard_bound(amplitude, precision):
    """Peak-decode error bound of Brassard et al.: |a_hat - a| <= 2*pi*sqrt(a(1-a))/2^p + pi^2/4^p."""
    points = 1 << precision
    return 2 * math.pi * math.sqrt(amplitude * (1 - amplitude)) / points + math.pi**2 / points**2


def _verify_qae(report, name, preparation, marked, amplitude, precision, exact):
    operation = amplitude_estimation(preparation, marked, precision=precision)
    program = operation.program()
    marginals = _readout_marginals(program, "phase")
    closed = _qae_closed_form(amplitude, precision)
    points = 1 << precision
    peak = max(max(m.values(), default=0.0) for m in marginals.values())
    # Peak decoding: take each path's own argmax and keep the worst error (the mirror peaks y and 2^p-y decode identically)
    worst_decode = 0.0
    for distribution in marginals.values():
        mode = max(distribution, key=distribution.get)
        estimate = math.sin(math.pi * mode / points) ** 2
        worst_decode = max(worst_decode, abs(estimate - amplitude))
    metrics, _ = _distribution_checks(
        marginals,
        closed,
        {
            "decode_error": worst_decode,
            "brassard_bound": _brassard_bound(amplitude, precision),
            "mode_peak_probability": peak,
        },
    )
    if exact:
        criterion = "a lies exactly on the decoding grid: decode_error == 0 and distribution TVD < 1e-9"
        passed = worst_decode < 1e-12 and metrics["max_tvd_vs_closed_form"] < 1e-9
    else:
        criterion = "distribution TVD < 1e-9 and peak-decode error <= Brassard bound"
        passed = (
            metrics["max_tvd_vs_closed_form"] < 1e-9
            and worst_decode <= metrics["brassard_bound"]
        )
    report.case(
        name,
        paths=list(marginals),
        parameters={
            "amplitude": amplitude,
            "precision": precision,
            "marked": list(marked),
        },
        metrics=metrics,
        criterion=criterion,
        passed=passed,
    )


def verify_qae(report):
    # a = 1/2 lies exactly on the grid (theta/pi = 1/4): deterministic bimodal, exact decoding
    for precision in (3, 4, 5):
        _verify_qae(
            report,
            f"qae-uniform-half-w{precision}",
            uniform_state(1),
            (1,),
            0.5,
            precision,
            exact=True,
        )
    # a = 3/8 off the grid: compare against the closed-form Grover bimodal and the Brassard bound
    for precision in (4, 5, 6):
        _verify_qae(
            report,
            f"qae-uniform-3over8-w{precision}",
            uniform_state(3),
            (0, 5, 7),
            3 / 8,
            precision,
            exact=False,
        )
    # Non-uniform preparation with a = 0.3: gate_state_prep amplitudes are known, closed form is independent
    _verify_qae(
        report,
        "qae-nonuniform-0p3-w5",
        gate_state_prep([math.sqrt(0.7), math.sqrt(0.3)]),
        (1,),
        0.3,
        5,
        exact=False,
    )


# ---------------------------------------------------------------------------
# Hadamard test and swap test
# ---------------------------------------------------------------------------


def _diagonal_unitary(angles, name):
    """Diagonal unitary U|x> = exp(i*angles[x])|x> (controlled-global-phase realization)."""
    width = (len(angles) - 1).bit_length()
    b = Builder(name, {"target": Bits(width)})
    for value, angle in enumerate(angles):
        if angle:
            with b.control(b["target"], value):
                b.global_phase(angle)
    return b.finish()


def verify_hadamard(report):
    # Single-qubit phase gate acting on |1>: expectation e^{i*theta}
    for angle in (0.6, -1.1):
        unitary = _phase_unitary(angle / (2 * math.pi))
        from oracq.algorithms.input_model.oracles import basis_state

        for component, expected in (
            ("real", math.cos(angle)),
            ("imag", math.sin(angle)),
        ):
            program = hadamard_test(
                unitary, basis_state(1, 1), component=component
            ).program()
            marginals = _readout_marginals(program, "probe")
            worst = max(
                abs(dist.get(0, 0.0) - dist.get(1, 0.0) - expected)
                for dist in marginals.values()
            )
            report.case(
                f"hadamard-phase-{component}-{angle}",
                paths=list(marginals),
                parameters={"angle": angle, "component": component},
                metrics={"z_expectation_error": worst, "expected": expected},
                criterion="probe Z expectation == " + ("cos theta" if component == "real" else "sin theta") + " (error < 1e-9)",
                passed=worst < 1e-9,
            )
    # Two-qubit diagonal unitary + complex-amplitude preparation: expectation sum_x |psi_x|^2 e^{i*theta_x}
    angles = (0.35, -0.9, 1.7, 0.55)
    amplitudes = (0.5, 0.5j, 0.5, -0.5)
    unitary = _diagonal_unitary(angles, "hadamard_diag_2q")
    preparation = gate_state_prep(amplitudes)
    expected_value = sum(
        abs(a) ** 2 * math.cos(t) + 1j * abs(a) ** 2 * math.sin(t)
        for a, t in zip(amplitudes, angles, strict=True)
    )
    for component, expected in (
        ("real", expected_value.real),
        ("imag", expected_value.imag),
    ):
        program = hadamard_test(unitary, preparation, component=component).program()
        marginals = _readout_marginals(program, "probe")
        worst = max(
            abs(dist.get(0, 0.0) - dist.get(1, 0.0) - expected)
            for dist in marginals.values()
        )
        report.case(
            f"hadamard-diagonal-2q-{component}",
            paths=list(marginals),
            parameters={"angles": list(angles), "component": component},
            metrics={"z_expectation_error": worst, "expected": expected},
            criterion="probe Z expectation == component of sum_x |psi_x|^2 e^{i*theta_x} (error < 1e-9)",
            passed=worst < 1e-9,
        )


def verify_swap(report):
    # (amplitude vector 1, amplitude vector 2), overlap F computed independently with numpy
    pairs = [
        ("same", (0.6, 0.8), (0.6, 0.8)),
        ("orthogonal", (0.6, 0.8), (0.8, -0.6)),
        ("partial", (1.0, 0.0), (0.6, 0.8)),
        ("two-qubit", (1 / math.sqrt(2), 1 / math.sqrt(2), 0.0, 0.0), (0.0, 0.5, 0.5, 1 / math.sqrt(2))),
    ]
    for name, first, second in pairs:
        import numpy as np

        overlap = abs(np.vdot(np.asarray(first, dtype=complex), np.asarray(second, dtype=complex))) ** 2
        expected = float((1 + overlap) / 2)
        program = swap_test(gate_state_prep(first), gate_state_prep(second)).program()
        marginals = _readout_marginals(program, "probe")
        worst = max(abs(dist.get(0, 0.0) - expected) for dist in marginals.values())
        report.case(
            f"swap-test-{name}",
            paths=list(marginals),
            parameters={
                "first": [str(complex(a)) for a in first],
                "second": [str(complex(a)) for a in second],
            },
            metrics={
                "p0_error": worst,
                "expected_p0": expected,
                "fidelity_overlap": float(overlap),
            },
            criterion="P(probe=0) == (1+F)/2 (error < 1e-9)",
            passed=worst < 1e-9,
        )


# ---------------------------------------------------------------------------
# Jordan gradient
# ---------------------------------------------------------------------------


def _linear_phase_oracle(coefficients, grid_bits):
    """Explicit phase table for the linear function f(x) = sum_i c_i*x_i (Jordan scaling convention)."""
    dimension = len(coefficients)
    points = 1 << grid_bits
    angles = []
    for value in range(1 << (dimension * grid_bits)):
        phase = 0.0
        for i, coefficient in enumerate(coefficients):
            chunk = (value >> (i * grid_bits)) & (points - 1)
            phase += coefficient * chunk / points
        angles.append(2 * math.pi * points * phase)
    return gate_phase_oracle(dimension * grid_bits, angles, phase_scale=points)


def _finite_difference(function, point, step):
    """Central finite-difference gradient (independent classical reference)."""
    return tuple(
        (
            function(*[x + (step if i == j else 0) for i, x in enumerate(point)])
            - function(*[x - (step if i == j else 0) for i, x in enumerate(point)])
        )
        / (2 * step)
        for j in range(len(point))
    )


def _verify_jordan_exact(
    report, name, oracle, dimension, grid_bits, gradient, paths_note=None, use_originir=True
):
    """Linear function: the readout distribution must land deterministically on the encoded value, decoding exactly the gradient."""
    operation = gradient_estimation(oracle, dimension=dimension, grid_bits=grid_bits)
    program = operation.program()
    marginals = _readout_marginals(program, "target", use_originir=use_originir)
    worst_decode = 0.0
    worst_peak = 1.0
    for distribution in marginals.values():
        mode = max(distribution, key=distribution.get)
        decoded = gradient_from_readout(mode, dimension=dimension, grid_bits=grid_bits)
        worst_decode = max(
            worst_decode, max(abs(g - e) for g, e in zip(decoded, gradient, strict=True))
        )
        worst_peak = min(worst_peak, distribution[mode])
    pairwise = 0.0
    labels = list(marginals)
    for i, left in enumerate(labels):
        for right in labels[i + 1 :]:
            pairwise = max(pairwise, tvd(marginals[left], marginals[right]))
    parameters = {"dimension": dimension, "grid_bits": grid_bits, "gradient": list(gradient)}
    if paths_note:
        parameters["paths_note"] = paths_note
    report.case(
        name,
        paths=labels,
        parameters=parameters,
        metrics={
            "max_decode_error": worst_decode,
            "min_peak_probability": worst_peak,
            "max_pairwise_tvd": pairwise,
        },
        criterion="deterministic readout (peak probability 1) with decoding exactly equal to the gradient (error < 1e-12)",
        passed=worst_decode < 1e-12 and worst_peak > 1 - 1e-9,
    )


def verify_jordan_gate_linear(report):
    # One dimension: negative components and varied grids, gradient components taken at exact grid values
    for grid_bits, coefficient in ((3, 3 / 8), (4, -5 / 16), (5, 9 / 32)):
        _verify_jordan_exact(
            report,
            f"jordan-linear-gate-d1-w{grid_bits}",
            _linear_phase_oracle((coefficient,), grid_bits),
            1,
            grid_bits,
            (coefficient,),
        )
    # Multiple dimensions: d = 2 and d = 3
    _verify_jordan_exact(
        report,
        "jordan-linear-gate-d2-w4",
        _linear_phase_oracle((3 / 16, -2 / 16), 4),
        2,
        4,
        (3 / 16, -2 / 16),
    )
    _verify_jordan_exact(
        report,
        "jordan-linear-gate-d3-w3",
        _linear_phase_oracle((1 / 8, -3 / 8, 1 / 4), 3),
        3,
        3,
        (1 / 8, -3 / 8, 1 / 4),
    )


def verify_jordan_function_oracle(report):
    """mathfunc phase-oracle path: exact linear readout after fixed-point arithmetic + phase kickback.

    The work bits of the default fixed-point format push the total qubit count
    well beyond the 24-qubit OriginIR budget, so only the register-level paths
    are exercised.
    """
    # d = 1: f(x) = 0.25x
    oracle1 = function_phase_oracle(
        "def f(x0):\n return 0.25 * x0", dimension=1, grid_bits=4
    )
    program1 = gradient_estimation(oracle1, dimension=1, grid_bits=4).program()
    budget1 = sum(r.type.width for r in program1.main.registers) + workspace_table(program1)[
        program1.entry
    ]
    _verify_jordan_exact(
        report,
        "jordan-function-oracle-d1-w4",
        oracle1,
        1,
        4,
        (0.25,),
        paths_note=f"originir-ext skipped: total qubits {budget1} exceed the 24-qubit budget",
        use_originir=False,
    )
    # d = 2: f(x0, x1) = 0.25*x0 - 0.125*x1, cross-checked against central finite differences
    source_gradient = (0.25, -0.125)
    oracle2 = function_phase_oracle(
        "def f(x0, x1):\n return 0.25 * x0 - 0.125 * x1", dimension=2, grid_bits=3
    )
    program2 = gradient_estimation(oracle2, dimension=2, grid_bits=3).program()
    budget2 = sum(r.type.width for r in program2.main.registers) + workspace_table(program2)[
        program2.entry
    ]
    finite_difference = _finite_difference(
        lambda x0, x1: 0.25 * x0 - 0.125 * x1, (0.5, 0.5), 1 / 8
    )
    fd_gap = max(
        abs(a - b) for a, b in zip(finite_difference, source_gradient, strict=True)
    )
    _verify_jordan_exact(
        report,
        "jordan-function-oracle-d2-w3",
        oracle2,
        2,
        3,
        source_gradient,
        paths_note=f"originir-ext skipped: total qubits {budget2} exceed the 24-qubit budget; "
        f"central finite difference vs truth gap {fd_gap:.3e} (0 for a linear function)",
        use_originir=False,
    )


def verify_jordan_perturbed(report):
    """Perturbed linear f(x) = a*x + x^2/N^2: peak stays at the truth, failure probability decays approximately quadratically.

    Informative metrics additionally report the central finite difference (which
    differs from the linear coefficient by O(1/N^2) when f is nonlinear).
    """
    a_numerator, a_denominator = 3, 8
    success = []
    for grid_bits in (3, 4, 5):
        points = 1 << grid_bits
        angles = [
            2
            * math.pi
            * points
            * (
                a_numerator / a_denominator * value / points
                + (value / points) ** 2 / points**2
            )
            for value in range(points)
        ]
        oracle = gate_phase_oracle(grid_bits, angles, phase_scale=points)
        program = gradient_estimation(oracle, dimension=1, grid_bits=grid_bits).program()
        marginals = _readout_marginals(program, "target")
        exact = a_numerator * points // a_denominator
        peak_probability = min(d[exact] for d in marginals.values())
        mode_ok = all(max(d, key=d.get) == exact for d in marginals.values())
        success.append(peak_probability)
        # Central finite difference: f'(1/2) = a + 1/N^2 (informative)
        fd = _finite_difference(
            lambda x, points=points: a_numerator / a_denominator * x + x**2 / points**2,
            (0.5,),
            1 / points,
        )[0]
        report.case(
            f"jordan-perturbed-w{grid_bits}",
            paths=list(marginals),
            parameters={"dimension": 1, "grid_bits": grid_bits, "linear_coefficient": 3 / 8},
            metrics={
                "success_probability": peak_probability,
                "failure_probability": 1 - peak_probability,
                "peak_matches_truth": mode_ok,
                "finite_difference_at_center": fd,
                "peak_gradient": exact / points,
            },
            criterion="peak == round(N*a), success probability monotonically increasing with the grid bit count",
            passed=mode_ok
            and (not success[:-1] or peak_probability > success[-2]),
        )
    ratios = [
        (1 - success[i + 1]) / (1 - success[i]) for i in range(len(success) - 1)
    ]
    report.case(
        "jordan-perturbed-decay-rate",
        paths=["derived"],
        parameters={"grid_bits_sequence": [3, 4, 5]},
        metrics={
            "success_probabilities": success,
            "failure_decay_ratios": ratios,
        },
        criterion="failure-probability decay ratio q_{m+1}/q_m < 0.34 (approximately quadratic convergence)",
        passed=all(ratio < 0.34 for ratio in ratios),
    )


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def run():
    report = Report(
        "estimation",
        "Distribution-level and expectation-level numerical validation of QPE/QAE/Hadamard/swap/Jordan gradient, "
        "against Dirichlet-kernel closed forms, independent numpy references, and finite differences.",
    )
    verify_qpe_ongrid(report)
    verify_qpe_offgrid(report)
    verify_qpe_fourier_eigenstate(report)
    verify_qpe_unitary_matrix(report)
    verify_qae(report)
    verify_hadamard(report)
    verify_swap(report)
    verify_jordan_gate_linear(report)
    verify_jordan_function_oracle(report)
    verify_jordan_perturbed(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
