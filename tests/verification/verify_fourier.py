"""Publication-grade numerical validation of the Fourier and transforms modules.

Two source files are covered:

- ``algorithms/fourier.py`` (qft / qft_with_work / inverse_qft / fourier_add):
  * n <= 4 full unitaries compared element-wise against independent DFT
    matrices (OriginIR-ext -> UniQC ``to_matrix`` plus ``effective_block``
    extraction);
  * n = 6/8/12 pointwise phase comparison on superposition states and Fourier
    mode focusing (rir-pysparq);
  * fourier_add against the integer-addition truth table (unitary and
    state-vector exhaustive for n <= 4, superposition exhaustive for n = 8,
    deterministic sampled basis pairs for n = 12);
  * inverse-QFT round-trip identity (unitary level for n = 3/4,
    register level for n = 8/12).
- ``algorithms/transforms.py`` (qubitization_walk / qsvt_sequence /
  oblivious_amplification): the inputs are the open block-encoding interfaces;
  per the task convention, end-to-end validation uses small gate-level bound
  instances (a 1-qubit target Hermitian matrix from matrix_pauli_encoding,
  a gate-level X/2 fixture, and a zero-width-signal Pauli word):
  * the full walk unitary = (2*Pi - I)U with spectral rotation angles
    +/-arccos(lambda/alpha);
  * the full QSVT-sequence unitary = the alternating S(phi)*U/U^dagger product;
    the zero-signal block = p(A/alpha), with the classical reference being an
    independently implemented 2x2 product formula and Chebyshev matrix
    recursion (not reusing the assembly logic under test);
  * the full OAA-iterate unitary = U*[R U^dagger R U]^it (a historical defect
    is now fixed: the original implementation was missing the trailing U and
    only assembled [R U^dagger R U]^it); on the minimal instance whose zero
    block is X/2, the zero-signal block = (-1)^it*sin((2it+1)theta)*V, with
    amplitudes matching the standard three-query OAA narrative step by step
    (sin theta -> sin 3theta).

All classical references are constructed independently: DFT and permutation
matrices from numpy/cmath, bit-by-bit phase-gate preparation (bypassing the
qft assembly under test), an explicit Pauli expansion for alpha, and the 2x2
QSP product and Chebyshev recursion.

Run: PYTHONPATH=src <python with pysparq+uniqc> tests/verification/verify_fourier.py
"""

from __future__ import annotations

import cmath
import math
import random

import numpy as np
from harness import (
    Report,
    adapter_pysparq,
    amplitude_error,
    amplitudes_to_statevector,
    basis_program,
    effective_block,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    sampled_inputs,
    statevector_error,
)

from oracq import Bits, Builder
from oracq.algorithms.common.fourier import fourier_add, inverse_qft, qft, qft_with_work
from oracq.algorithms.common.qsvt import qsp_phases
from oracq.algorithms.common.transforms import (
    oblivious_amplification,
    qsvt_sequence,
    qubitization_walk,
)
from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding, pauli_word
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import annotate

# ---------------------------------------------------------------------------
# Independent classical references
# ---------------------------------------------------------------------------


def _dft_matrix(n):
    """Positive-sign DFT matrix: F[y, x] = exp(2*pi*i*x*y / 2^n) / sqrt(2^n)."""
    size = 1 << n
    scale = 1.0 / math.sqrt(size)
    return np.array(
        [
            [cmath.exp(2j * math.pi * x * y / size) * scale for x in range(size)]
            for y in range(size)
        ]
    )


def _add_permutation(n):
    """Classical truth-table permutation of fourier_add: |a, b> -> |a, (a+b) mod 2^n> (a in the low bits)."""
    size = 1 << n
    perm = np.zeros((size * size, size * size), dtype=complex)
    for a in range(size):
        for b in range(size):
            perm[a | (((a + b) % size) << n), a | (b << n)] = 1.0
    return perm


def _pauli_alpha(matrix):
    """Block-encoding normalization alpha = sum_P |Tr(P M)| / d via an explicit Pauli expansion (independent of the library implementation)."""
    m = np.asarray(matrix, dtype=complex)
    d = m.shape[0]
    n = (d - 1).bit_length()
    paulis = [
        np.eye(2),
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]),
        np.array([[1, 0], [0, -1]]),
    ]
    alpha = 0.0
    for index in range(4**n):
        word = np.array([[1.0 + 0j]])
        rest = index
        for _ in range(n):
            word = np.kron(paulis[rest % 4], word)
            rest //= 4
        alpha += abs(np.trace(word @ m)) / d
    return alpha


def _signal_projector(dim, data_dim):
    """Pi = |0><0|_signal tensor I (in the entry layout signal is the high part; the zero block is the low data_dim dimensions)."""
    proj = np.zeros((dim, dim))
    proj[:data_dim, :data_dim] = np.eye(data_dim)
    return proj


def _signal_phase(phi, dim, data_dim):
    """QSVT phase operator S(phi): e^{+i*phi} on the zero-signal branch, e^{-i*phi} elsewhere."""
    diag = np.full(dim, cmath.exp(-1j * phi), dtype=complex)
    diag[:data_dim] = cmath.exp(1j * phi)
    return np.diag(diag)


def _qsp_response_independent(x, phases):
    """2x2 product reference per the documented convention: p(x) = [S(phi_d) W ... W S(phi_0)]_00, rightmost factor applied first."""
    s = math.sqrt(max(0.0, 1.0 - x * x))
    w = np.array([[x, s], [s, -x]], dtype=complex)
    total = np.diag([cmath.exp(1j * phases[0]), cmath.exp(-1j * phases[0])])
    for phi in phases[1:]:
        total = np.diag([cmath.exp(1j * phi), cmath.exp(-1j * phi)]) @ w @ total
    return total[0, 0]


def _chebyshev_matrix(degree, m):
    """Independent recursion for T_d(M): T_0 = I, T_1 = M, T_{k+1} = 2 M T_k - T_{k-1}."""
    if degree == 0:
        return np.eye(m.shape[0], dtype=complex)
    t0, t1 = np.eye(m.shape[0], dtype=complex), np.asarray(m, dtype=complex).copy()
    for _ in range(2, degree + 1):
        t0, t1 = t1, 2 * m @ t1 - t0
    return t1


def _angle_distance(a, b):
    """Distance of two rotation angles on the circle (in [0, pi])."""
    d = abs(a - b) % (2 * math.pi)
    return min(d, 2 * math.pi - d)


# ---------------------------------------------------------------------------
# Drivers and fixtures
# ---------------------------------------------------------------------------


def _fourier_mode_program(n, k):
    """Bit-by-bit phase-gate preparation of the negative-sign Fourier mode followed by qft.

    The prepared state |phi_k> = sum_x exp(-2*pi*i*k*x / 2^n) |x> / sqrt(2^n)
    is a separable product state built only from H and single-qubit phase
    gates (bypassing the qft assembly under test); it is row -k of the
    positive-sign DFT, so after the QFT all amplitude must focus onto basis
    index k.
    """
    size = 1 << n
    b = Builder(f"fourier_mode_{n}_{k}", {"target": Bits(n)})
    for j in range(n):
        b.h(b["target"][j])
        b.gate("phase", b["target"][j], -2.0 * math.pi * k * (1 << j) / size)
    b.call(qft(n), target=b["target"])
    return b.finish().program()


def _roundtrip_unitary_program(n):
    """Bare composition of qft followed by inverse_qft (unitary-level round trip, no state preparation)."""
    b = Builder(f"qft_roundtrip_unitary_{n}", {"target": Bits(n)})
    b.call(qft(n), target=b["target"])
    b.call(inverse_qft(n), target=b["target"])
    return b.finish().program()


def _roundtrip_program(n, *, value=None):
    """Round-trip program of qft then inverse_qft; starts from full superposition when value is None."""
    b = Builder(f"qft_roundtrip_{n}_{value}", {"target": Bits(n)})
    if value is None:
        b.h(b["target"])
    else:
        for bit in range(n):
            if (value >> bit) & 1:
                b.x(b["target"][bit])
    b.call(qft(n), target=b["target"])
    b.call(inverse_qft(n), target=b["target"])
    return b.finish().program()


def _fourier_add_fixed_b(n, b0):
    """Invoke fourier_add with a in full superposition and b fixed to a known constant: branches (a, (a+b0) mod 2^n)."""
    b = Builder(f"fadd_superposed_a_{n}_{b0}", {"a": Bits(n), "b": Bits(n)})
    b.h(b["a"])
    for bit in range(n):
        if (b0 >> bit) & 1:
            b.x(b["b"][bit])
    b.call(fourier_add(n), a=b["a"], b=b["b"])
    return b.finish().program()


# Gate-level binding fixture for transforms: 1-qubit target Hermitian matrix (spectral values fall in (-1, 1) after normalization)
MATRIX = np.array([[0.5, 0.3], [0.3, -0.1]])


def _fixture_be():
    """matrix_pauli_encoding belongs to the block_encoding module; here it only serves as the source of the bound instance."""
    return matrix_pauli_encoding(MATRIX.tolist())


def _half_x_be():
    """Gate-level minimal OAA fixture: U = M tensor X, with the zero-signal block exactly X/2 (a V/2 instance with sin theta = 1/2)."""
    b = Builder("be_half_x", {"target": Bits(1), "signal": Bits(1)})
    b.x(b["target"])
    b.z(b["signal"])
    b.ry(b["signal"], 2.0 * math.pi / 3.0)
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))


def _unitary_from_runs(operation, target_width, signal_width):
    """Assemble the small-instance unitary column by column via the reference path (cross-checked against the OriginIR path)."""
    data = 1 << target_width
    dim = 1 << (target_width + signal_width)
    matrix = np.zeros((dim, dim), dtype=complex)
    for s in range(1 << signal_width):
        for t in range(data):
            state = reference(basis_program(operation, {"target": t, "signal": s}))
            for (tv, sv), amp in state.items():
                matrix[tv | (sv << target_width), t | (s << target_width)] = amp
    return matrix


# ---------------------------------------------------------------------------
# fourier.py: unitary-level validation (n <= 4, OriginIR-ext + UniQC to_matrix)
# ---------------------------------------------------------------------------


def verify_qft_unitary(report):
    for n in (1, 2, 3, 4):
        unitary = originir_unitary(qft(n).program())
        block, leakage = effective_block(unitary, n)
        error = float(np.abs(block - _dft_matrix(n)).max())
        report.case(
            f"qft-unitary-dft-n{n}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"width": n, "dimension": 1 << n},
            metrics={"max_error": error, "workspace_leakage": leakage},
            criterion="full unitary agrees element-wise with the positive-sign DFT matrix (max_error < 1e-12, zero leakage)",
            passed=error < 1e-12 and leakage < 1e-12,
        )


def verify_inverse_qft_unitary(report):
    for n in (1, 2, 3, 4):
        unitary = originir_unitary(inverse_qft(n).program())
        block, leakage = effective_block(unitary, n)
        error = float(np.abs(block - _dft_matrix(n).conj().T).max())
        report.case(
            f"inverse-qft-unitary-n{n}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"width": n, "dimension": 1 << n},
            metrics={"max_error": error, "workspace_leakage": leakage},
            criterion="full unitary agrees element-wise with the adjoint of the DFT matrix (max_error < 1e-12)",
            passed=error < 1e-12 and leakage < 1e-12,
        )
    for n in (3, 4):
        unitary = originir_unitary(_roundtrip_unitary_program(n))
        error = float(np.abs(unitary - np.eye(1 << n)).max())
        report.case(
            f"inverse-qft-roundtrip-unitary-n{n}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"width": n},
            metrics={"max_error": error},
            criterion="composite unitary of qft followed by inverse_qft is the identity (max_error < 1e-12)",
            passed=error < 1e-12,
        )


# ---------------------------------------------------------------------------
# fourier.py: register-level validation (rir-pysparq, n = 6/8/12)
# ---------------------------------------------------------------------------


def verify_qft_basis_rows(report):
    """Pointwise comparison of all 2^n amplitudes of a basis state |x> after QFT against the DFT row."""
    for n in (6, 8, 12):
        size = 1 << n
        uniform = 1.0 / math.sqrt(size)
        inputs = sorted({0, 1, size // 2, size - 1, 0x5A5A % size})
        worst = 0.0
        for x in inputs:
            state = rir_pysparq(basis_program(qft(n), {"target": x}), max_states=1 << (n + 1))
            for y in range(size):
                expected = cmath.exp(2j * math.pi * x * y / size) * uniform
                worst = max(worst, abs(state.get((y,), 0j) - expected))
        report.case(
            f"qft-basis-row-pointwise-n{n}",
            paths=["rir-pysparq"],
            parameters={"width": n, "inputs": inputs, "amplitudes_per_input": size},
            metrics={"max_error": worst},
            criterion="all 2^n amplitudes of every input agree pointwise with the DFT row (max_error < 1e-9)",
            passed=worst < 1e-9,
        )


def verify_qft_fourier_mode_focus(report):
    """A Fourier mode prepared by independent gate-level means focuses onto a single basis index after the QFT."""
    for n in (6, 8, 12):
        size = 1 << n
        modes = sorted({1, size // 3, size - 1})
        min_probability = 1.0
        worst_leak = 0.0
        for k in modes:
            state = rir_pysparq(_fourier_mode_program(n, k), max_states=1 << (n + 1))
            min_probability = min(min_probability, abs(state.get((k,), 0j)) ** 2)
            leak = max(
                (abs(v) for key, v in state.items() if key != (k,)),
                default=0.0,
            )
            worst_leak = max(worst_leak, leak)
        report.case(
            f"qft-fourier-mode-focus-n{n}",
            paths=["rir-pysparq"],
            parameters={"width": n, "modes": modes},
            metrics={
                "min_success_probability": min_probability,
                "max_leaked_amplitude": worst_leak,
            },
            criterion="focusing success probability >= 1 - 1e-12 with out-of-mode leaked amplitude < 1e-6",
            passed=min_probability >= 1 - 1e-12 and worst_leak < 1e-6,
        )


def verify_qft_cross_path(report):
    """n = 4 phase-rich state (QFT acting on a basis state) four-path cross-check, plus qft_with_work adapter equivalence."""
    program = basis_program(qft(4), {"target": 11})
    ref = reference(program)
    deviation = max(
        amplitude_error(ref, rir_pysparq(program)),
        amplitude_error(ref, adapter_pysparq(program)),
        statevector_error(originir_ext(program), amplitudes_to_statevector(ref, [4])),
    )
    legacy = reference(basis_program(qft_with_work(4), {"target": 11}))
    legacy_deviation = max(
        abs(legacy.get((y, 0), 0j) - ref.get((y,), 0j)) for y in range(16)
    )
    deviation = max(deviation, legacy_deviation)
    report.case(
        "qft-cross-path-n4",
        paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        parameters={"width": 4, "input": 11},
        metrics={"max_pairwise_deviation": deviation},
        criterion="pairwise deviation across the four paths and the qft_with_work adapter < 1e-9",
        passed=deviation < 1e-9,
    )


# ---------------------------------------------------------------------------
# fourier.py: Fourier addition truth table
# ---------------------------------------------------------------------------


def verify_fourier_add_unitary(report):
    for n in (1, 2, 3, 4):
        unitary = originir_unitary(fourier_add(n).program())
        block, leakage = effective_block(unitary, 2 * n)
        error = float(np.abs(block - _add_permutation(n)).max())
        report.case(
            f"fourier-add-unitary-n{n}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"width": n, "dimension": 1 << (2 * n)},
            metrics={"max_error": error, "workspace_leakage": leakage},
            criterion="full unitary agrees with the integer modular-addition permutation matrix (max_error < 1e-12)",
            passed=error < 1e-12 and leakage < 1e-12,
        )


def verify_fourier_add_truth_table(report):
    """n <= 4: fixing b0 and superposing a covers all 4^n (a, b) inputs across iterations (two-path cross-check)."""
    for n in (1, 2, 3, 4):
        size = 1 << n
        uniform = 1.0 / math.sqrt(size)
        worst = 0.0
        for b0 in range(size):
            program = _fourier_add_fixed_b(n, b0)
            expected = {(a, (a + b0) % size): uniform for a in range(size)}
            worst = max(worst, amplitude_error(reference(program), expected))
            vector = originir_ext(program)
            worst = max(
                worst,
                statevector_error(vector, amplitudes_to_statevector(expected, [n, n])),
            )
        report.case(
            f"fourier-add-truthtable-n{n}",
            paths=["reference", "originir-ext"],
            parameters={"width": n, "input_pairs": size * size},
            metrics={"max_error": worst},
            criterion="all (a, b) input pairs map to |a, (a+b) mod 2^n> with uniform amplitudes (max_error < 1e-9)",
            passed=worst < 1e-9,
        )


def verify_fourier_add_superposed_a_n8(report):
    """n = 8: a single rir-pysparq superposition run exhausts all 256 branches of a (intermediate-state peak 2^16)."""
    n, size = 8, 256
    uniform = 1.0 / math.sqrt(size)
    worst = 0.0
    for b0 in (0, 1, 0x55, 0x80, 0xFF):
        state = rir_pysparq(_fourier_add_fixed_b(n, b0), max_states=1 << 17)
        expected = {(a, (a + b0) % size): uniform for a in range(size)}
        worst = max(worst, amplitude_error(state, expected))
    report.case(
        "fourier-add-superposed-a-n8",
        paths=["rir-pysparq"],
        parameters={"width": n, "constants": [0, 1, 0x55, 0x80, 0xFF], "branches": size},
        metrics={"max_error": worst},
        criterion="all 256 branches of every b0 agree amplitude by amplitude (max_error < 1e-9)",
        passed=worst < 1e-9,
    )


def verify_fourier_add_sampled_n12(report):
    """n = 12: deterministic output check on sampled basis pairs (intermediate-state peak 2^12; superposing a would exceed the budget)."""
    n = 12
    values, _ = sampled_inputs(n)
    pairs = list(zip(values[:12], values[1:13], strict=True))
    operation = fourier_add(n)
    failures = 0
    for a, b0 in pairs:
        state = rir_pysparq(
            basis_program(operation, {"a": a, "b": b0}), max_states=1 << 13
        )
        if set(state) != {(a, (a + b0) % (1 << n))} or abs(
            next(iter(state.values())) - 1
        ) > 1e-12:
            failures += 1
    report.case(
        "fourier-add-sampled-basis-n12",
        paths=["rir-pysparq"],
        parameters={"width": n, "sampled_pairs": len(pairs)},
        metrics={"failures": failures},
        criterion="sampled inputs produce exactly |a, (a+b) mod 2^12> (failures == 0)",
        passed=failures == 0,
    )


# ---------------------------------------------------------------------------
# fourier.py: inverse-QFT round-trip identity (register level)
# ---------------------------------------------------------------------------


def verify_inverse_qft_roundtrip_wide(report):
    for n in (8, 12):
        size = 1 << n
        uniform = 1.0 / math.sqrt(size)
        state = rir_pysparq(_roundtrip_program(n), max_states=1 << (n + 1))
        expected = {(x,): uniform for x in range(size)}
        uniform_error = amplitude_error(state, expected)
        failures = 0
        for x in (0, 1, size // 2, size - 1, 0x5A5A % size):
            single = rir_pysparq(_roundtrip_program(n, value=x), max_states=1 << (n + 1))
            if set(single) != {(x,)} or abs(next(iter(single.values())) - 1) > 1e-12:
                failures += 1
        report.case(
            f"inverse-qft-roundtrip-wide-n{n}",
            paths=["rir-pysparq"],
            parameters={"width": n, "branches": size, "basis_samples": 5},
            metrics={"uniform_max_error": uniform_error, "basis_failures": failures},
            criterion="full-superposition round trip restores the uniform state (< 1e-9) and sampled basis states round-trip to themselves (failures == 0)",
            passed=uniform_error < 1e-9 and failures == 0,
        )


# ---------------------------------------------------------------------------
# transforms.py: end-to-end validation on small gate-level bound instances
# ---------------------------------------------------------------------------


def verify_qubitization_walk(report):
    be = _fixture_be()
    tw, sw = be.width, be.signal_qubits
    data, dim = 1 << tw, 1 << (tw + sw)
    u_runs = _unitary_from_runs(be.operation, tw, sw)
    u_origin = originir_unitary(be.operation.program())
    crosscheck = float(np.abs(u_runs - u_origin).max())
    proj = _signal_projector(dim, data)
    walk_expected = (2 * proj - np.eye(dim)) @ u_runs
    walk_actual = originir_unitary(qubitization_walk(be).program())
    unitary_error = float(np.abs(walk_actual - walk_expected).max())
    alpha = _pauli_alpha(MATRIX)
    zero_block_error = float(np.abs(walk_actual[:data, :data] - MATRIX / alpha).max())
    # Spectral property: every spectral value x = lambda/alpha corresponds to a pair of eigenangles +/-arccos(x); the rotation angle on the remaining space lies in {0, pi}
    eigenangles = [abs(float(a)) for a in np.angle(np.linalg.eigvals(walk_actual))]
    allowed = {0.0, math.pi}
    for x in np.linalg.eigvalsh(MATRIX / alpha):
        allowed.add(math.acos(float(np.clip(x, -1.0, 1.0))))
    spectrum_dev = max(
        min(_angle_distance(angle, expect) for expect in allowed) for angle in eigenangles
    )
    report.case(
        "qubitization-walk-unitary-spectrum",
        paths=["originir-ext+uniqc-to_matrix", "reference"],
        parameters={
            "matrix": MATRIX.tolist(),
            "alpha": alpha,
            "signal_qubits": sw,
            "dimension": dim,
        },
        metrics={
            "unitary_max_error": unitary_error,
            "be_path_crosscheck": crosscheck,
            "zero_block_max_error": zero_block_error,
            "spectrum_max_deviation": spectrum_dev,
        },
        criterion=(
            "walk unitary = (2*Pi - I)U (< 1e-12), zero-signal block = A/alpha, "
            "eigenangles fall in {+/-arccos(lambda/alpha)} U {0, pi} (angle deviation < 1e-9)"
        ),
        passed=unitary_error < 1e-12
        and zero_block_error < 1e-12
        and spectrum_dev < 1e-9,
    )


def verify_qsvt_sequence(report):
    be = _fixture_be()
    tw, sw = be.width, be.signal_qubits
    data, dim = 1 << tw, 1 << (tw + sw)
    u = originir_unitary(be.operation.program())
    rng = random.Random(20260916)
    phases = tuple(rng.uniform(-math.pi, math.pi) for _ in range(7))
    actual = originir_unitary(qsvt_sequence(be, phases).program())
    # Full-unitary reference: alternating S(phi)*U/U^dagger in chronological order (U taken as the measured unitary of the BE)
    expected = _signal_phase(phases[0], dim, data)
    for k in range(len(phases) - 1):
        expected = (u if k % 2 == 0 else u.conj().T) @ expected
        expected = _signal_phase(phases[k + 1], dim, data) @ expected
    unitary_error = float(np.abs(actual - expected).max())
    # Zero-signal-block spectral semantics: p(A/alpha) with p from the independent 2x2 product formula
    alpha = _pauli_alpha(MATRIX)
    values, vectors = np.linalg.eigh(MATRIX / alpha)
    response = [complex(_qsp_response_independent(float(x), phases)) for x in values]
    expected_block = (vectors * response) @ vectors.conj().T
    block_error = float(np.abs(actual[:data, :data] - expected_block).max())
    report.case(
        "qsvt-sequence-random-phases",
        paths=["originir-ext+uniqc-to_matrix"],
        parameters={"degree": len(phases) - 1, "seed": 20260916},
        metrics={
            "unitary_max_error": unitary_error,
            "zero_block_max_error": block_error,
        },
        criterion="sequence unitary = the alternating S(phi)*U/U^dagger product with zero block = p(A/alpha) (both < 1e-12)",
        passed=unitary_error < 1e-12 and block_error < 1e-12,
    )
    # Chebyshev polynomial end to end: qsp_phases only generates the phase input; the reference is the independent T_d recursion
    for degree in (4, 5):
        coeffs = tuple(float(c) for c in np.polynomial.chebyshev.cheb2poly([0] * degree + [1]))
        t_phases = qsp_phases(coeffs)
        actual = originir_unitary(qsvt_sequence(be, t_phases).program())
        expected_block = _chebyshev_matrix(degree, MATRIX / alpha)
        block_error = float(np.abs(actual[:data, :data] - expected_block).max())
        report.case(
            f"qsvt-sequence-chebyshev-t{degree}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"degree": degree, "alpha": alpha},
            metrics={"zero_block_max_error": block_error},
            criterion=f"zero-signal block = T_{degree}(A/alpha) (independent Chebyshev recursion, < 1e-12)",
            passed=block_error < 1e-12,
        )


def verify_qsvt_degenerate(report):
    """Degenerate branch with signal_qubits == 0: phases degenerate into unconditional global phases."""
    be = pauli_word("X")
    x_matrix = np.array([[0, 1], [1, 0]], dtype=complex)
    for phases in ((0.3, -0.7, 1.1), (0.3, -0.7, 1.1, 0.5)):
        degree = len(phases) - 1
        base = np.eye(2, dtype=complex) if degree % 2 == 0 else x_matrix
        expected = cmath.exp(1j * sum(phases)) * base
        actual = originir_unitary(qsvt_sequence(be, phases).program())
        error = float(np.abs(actual - expected).max())
        report.case(
            f"qsvt-sequence-degenerate-signal0-d{degree}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"degree": degree, "be": "pauli_word(X)"},
            metrics={"max_error": error},
            criterion="degenerate sequence = e^{i*sum(phi)}*X^d (max_error < 1e-12)",
            passed=error < 1e-12,
        )


def verify_oaa_unitary(report):
    be = _fixture_be()
    tw, sw = be.width, be.signal_qubits
    data, dim = 1 << tw, 1 << (tw + sw)
    u = originir_unitary(be.operation.program())
    proj = _signal_projector(dim, data)
    reflect = np.eye(dim) - 2 * proj
    iterate = reflect @ u.conj().T @ reflect @ u
    expected = u.copy()
    for iterations in (1, 2, 3):
        expected = expected @ iterate
        actual = originir_unitary(oblivious_amplification(be, iterations).program())
        error = float(np.abs(actual - expected).max())
        report.case(
            f"oaa-unitary-it{iterations}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"iterations": iterations, "dimension": dim},
            metrics={"max_error": error},
            criterion="iterate unitary = U*[R U^dagger R U]^it (R = I - 2*Pi, max_error < 1e-12)",
            passed=error < 1e-12,
        )


def verify_oaa_half_block(report):
    """Minimal instance with zero block X/2: validates the amplification behavior of the standard three-query OAA.

    On this fixture the fixed operator U*[R U^dagger R U]^it has zero-signal
    block = (-1)^it*sin((2it+1)theta)*V (V = X, theta = pi/6): one iteration
    gives amplitude sin 3theta = 1 (full signal recovery in amplitude terms),
    with amplitudes matching the standard OAA narrative step by step; the
    (-1)^it is a global phase of the observable unitary, recorded as an
    informative metric.
    """
    x_matrix = np.array([[0, 1], [1, 0]], dtype=complex)
    theta = math.pi / 6
    for iterations in (1, 2, 3):
        actual = originir_unitary(oblivious_amplification(_half_x_be(), iterations).program())
        zero_block = actual[:2, :2]
        sign = -1.0 if iterations % 2 else 1.0
        expected_block = sign * math.sin((2 * iterations + 1) * theta) * x_matrix
        error = float(np.abs(zero_block - expected_block).max())
        report.case(
            f"oaa-half-block-it{iterations}",
            paths=["originir-ext+uniqc-to_matrix"],
            parameters={"fixture": "zero block = X/2 (sin theta = 1/2)", "iterations": iterations},
            metrics={
                "amplified_zero_block_error": error,
                "amplitude": abs(math.sin((2 * iterations + 1) * theta)),
                "global_sign": int(sign),
            },
            criterion=(
                "zero-signal block = (-1)^it*sin((2it+1)theta)*X (max_error < 1e-12); "
                "amplitudes follow the standard OAA narrative sin((2it+1)theta) step by step"
            ),
            passed=error < 1e-12,
        )


def run():
    report = Report(
        "fourier",
        "Three-level (unitary / superposition / truth-table) validation of QFT/inverse QFT/Fourier addition (n <= 12), "
        "plus gate-level bound end-to-end validation of the transforms block-encoding operations (walk/QSVT/OAA).",
    )
    verify_qft_unitary(report)
    verify_inverse_qft_unitary(report)
    verify_qft_basis_rows(report)
    verify_qft_fourier_mode_focus(report)
    verify_qft_cross_path(report)
    verify_fourier_add_unitary(report)
    verify_fourier_add_truth_table(report)
    verify_fourier_add_superposed_a_n8(report)
    verify_fourier_add_sampled_n12(report)
    verify_inverse_qft_roundtrip_wide(report)
    verify_qubitization_walk(report)
    verify_qsvt_sequence(report)
    verify_qsvt_degenerate(report)
    verify_oaa_unitary(report)
    verify_oaa_half_block(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
