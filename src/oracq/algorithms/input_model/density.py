"""Density matrix (DM) input model and Gibbs state preparation.

The core abstraction of the DM input model is purification access: access to a
density matrix ρ is defined as a quantum operation U that prepares its
purified state, U ``|0>`` = ``|ψ_ρ⟩`` acting on the two registers system and
environment, with the partial trace over the environment giving Tr_env
``|ψ⟩⟨ψ|`` = ρ. This module provides three paradigm layers:
abstract_purification open declarations, gate_purification explicit
small-matrix witnesses (eigendecomposition followed by preparation via
multiplexed rotations), and the PurificationAccess.from_state_preparation pure
state adapter (a pure state is the trivial purification with the environment
restored to zero). This view is the stable interface relied on by B2 (quantum
SDP).

Gibbs state preparation follows the QSVT purification route (Chowdhury–Somma
2017, van Apeldoorn–Gilyén 2019, Gilyén et al. 2019, arXiv:1806.01838): first
prepare the purification of the maximally mixed state on system and
environment (n Bell pairs), then apply to system a QSVT block encoding whose
zero-signal block is proportional to g(H/α), where g(x) = exp(−βα(x+1)/2)
takes values in (0,1] on [−1,1]. g decomposes by parity into the two branches
e^{−c}·cosh(cx) and −e^{−c}·sinh(cx) (c = βα/2), each approximated by a
modified Bessel truncation, with phases synthesized via imaginary completion,
the real part extracted by (U_Φ + U_{−Φ})/2, and the branches finally summed
by LCU. After post-selecting signal == 0, the reduced density matrix of
system is proportional to g(H/α)² = e^{−βH} (up to normalization).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, replace
from typing import cast

from oracq.algorithms.common.qsvt import (
    _MAX_DEGREE,
    _chebyshev_t,
    _eval,
    _real_qsvt_be,
    _scale,
    _sup_norm,
    _synthesize_with_imag,
    _trim,
)
from oracq.algorithms.input_model.contracts import (
    OracleSpec,
    OracleView,
    fail,
    finite_real,
    positive_integer,
    require_instance,
    validate_signature,
)
from oracq.algorithms.input_model.operators import BlockEncoding, _name, linear_combination
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    annotate,
    declare,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError, fuse

__all__ = [
    "ApproximatePurification",
    "PurificationAccess",
    "abstract_purification",
    "gate_purification",
    "gibbs_purification",
    "gibbs_state",
    "maximally_mixed_purification",
    "partial_trace",
    "trace_distance",
]


# ---------------------------------------------------------------------------
# Classical small-matrix utilities (pure Python complex matrices, reused by
# witnesses and the classical side of downstream algorithms).
# ---------------------------------------------------------------------------


def _as_complex_matrix(
    matrix: Iterable[Iterable[complex]], path: str
) -> tuple[tuple[complex, ...], ...]:
    """Normalize the input into a nonempty square complex matrix with finite entries, as nested tuples."""
    try:
        result = tuple(tuple(complex(v) for v in row) for row in matrix)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{path} requires a numeric matrix in square form") from exc
    d = len(result)
    if d < 1 or any(len(row) != d for row in result):
        raise ValidationError(f"{path} requires a nonempty square matrix")
    if not all(math.isfinite(v.real) and math.isfinite(v.imag) for row in result for v in row):
        raise ValidationError(f"matrix elements of {path} must be finite")
    return result


def _check_hermitian(
    matrix: Iterable[Iterable[complex]], path: str, *, tol: float = 1e-9
) -> tuple[tuple[complex, ...], ...]:
    """Validate and return a Hermitian complex square matrix; raise ``ValidationError`` when not Hermitian."""
    matrix = _as_complex_matrix(matrix, path)
    d = len(matrix)
    scale = max(1.0, max(abs(v) for row in matrix for v in row))
    for i in range(d):
        for j in range(i + 1, d):
            if abs(matrix[i][j] - matrix[j][i].conjugate()) > tol * scale:
                raise ValidationError(f"{path} must be a Hermitian matrix")
    return matrix


def _check_density_matrix(
    rho: Iterable[Iterable[complex]], *, tol: float = 1e-7
) -> tuple[tuple[complex, ...], ...]:
    """Validate the dimension and trace of a density matrix, returning the normalized complex square matrix as nested tuples."""
    matrix = _check_hermitian(rho, "gate_purification.rho")
    d = len(matrix)
    if d < 2 or d & (d - 1):
        raise ValidationError("density matrix dimension must be a power of two of at least 2")
    if d > 16:
        raise ValidationError("explicit purification witnesses support only small density matrices of dimension at most 16")
    if abs(sum(matrix[i][i] for i in range(d)) - 1.0) > tol:
        raise ValidationError("the trace of a density matrix must be 1")
    return matrix


def _hermitian_eigendecomposition(
    matrix: Sequence[Sequence[complex]], *, tol: float = 1e-13, max_sweeps: int = 64
) -> tuple[tuple[float, ...], list[list[complex]]]:
    """Eigendecomposition of a small Hermitian matrix by the cyclic Jacobi method.

    Returns (eigenvalue tuple, eigenvector table); eigenvalues are sorted in
    descending order, and row i column j of the eigenvector table is component
    i of eigenvector j.
    """
    d = len(matrix)
    a = [[complex(matrix[i][j]) for j in range(d)] for i in range(d)]
    vectors = [[1.0 + 0j if i == j else 0.0 + 0j for j in range(d)] for i in range(d)]
    scale = max(1.0, max(abs(a[i][j]) for i in range(d) for j in range(d)))
    for _ in range(max_sweeps):
        off = max((abs(a[i][j]) for i in range(d) for j in range(i + 1, d)), default=0.0)
        if off <= tol * scale:
            break
        for p in range(d):
            for q in range(p + 1, d):
                if abs(a[p][q]) <= tol * scale:
                    continue
                # First apply a diagonal phase rotation to make the off-diagonal element a
                # positive real number, then eliminate it with a real Jacobi rotation.
                # The similarity transform D†AD does not change diagonal entries;
                # a[q][q] must be restored after the column scaling.
                phase = a[p][q] / abs(a[p][q])
                diagonal_qq = a[q][q].real
                for k in range(d):
                    a[k][q] *= phase
                    vectors[k][q] *= phase
                a[q][q] = diagonal_qq + 0j
                for k in range(d):
                    if k != q:
                        a[q][k] = a[k][q].conjugate()
                app, aqq, b = a[p][p].real, a[q][q].real, a[p][q].real
                tau = (aqq - app) / (2 * b)
                t = (1.0 if tau >= 0 else -1.0) / (abs(tau) + math.sqrt(1.0 + tau * tau))
                cos, sin = 1.0 / math.sqrt(1.0 + t * t), t / math.sqrt(1.0 + t * t)
                for k in range(d):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = cos * akp - sin * akq
                    a[k][q] = sin * akp + cos * akq
                for j in range(d):
                    apj, aqj = a[p][j], a[q][j]
                    a[p][j] = cos * apj - sin * aqj
                    a[q][j] = sin * apj + cos * aqj
                for k in range(d):
                    vkp, vkq = vectors[k][p], vectors[k][q]
                    vectors[k][p] = cos * vkp - sin * vkq
                    vectors[k][q] = sin * vkp + cos * vkq
    order = sorted(range(d), key=lambda k: -a[k][k].real)
    values = tuple(a[k][k].real for k in order)
    return values, [[vectors[i][k] for k in order] for i in range(d)]


def partial_trace(
    amplitudes: Iterable[complex], system_width: int, environment_width: int
) -> tuple[tuple[complex, ...], ...]:
    """Take the partial trace over environment, returning the reduced density matrix on system as a row-major nested tuple.

    amplitudes is the dense state vector of length 2^(system_width +
    environment_width), with the basis state index convention system |
    (environment << system_width).

    Args:
        amplitudes: Pure state amplitude sequence satisfying the length convention above.
        system_width: Bit width of the system register, range 0..20.
        environment_width: Bit width of the environment register to be partially traced, range 0..20.

    Returns:
        tuple[tuple[complex, ...], ...]: The 2^system_width-dimensional reduced
    density matrix as a row-major nested tuple.
    """
    positive_integer(system_width, "partial_trace.system_width", minimum=0, maximum=20)
    positive_integer(environment_width, "partial_trace.environment_width", minimum=0, maximum=20)
    values = [complex(v) for v in amplitudes]
    dim_s, dim_e = 1 << system_width, 1 << environment_width
    if len(values) != dim_s * dim_e:
        raise ValidationError("state vector length does not match the register widths")
    return tuple(
        tuple(
            sum(
                values[i + (e << system_width)] * values[j + (e << system_width)].conjugate()
                for e in range(dim_e)
            )
            for j in range(dim_s)
        )
        for i in range(dim_s)
    )


def gibbs_state(
    hamiltonian: Iterable[Iterable[complex]], beta: float
) -> tuple[tuple[complex, ...], ...]:
    """Classical reference Gibbs state e^{−βH}/Tr(e^{−βH}); used only as a classical witness for small matrices.

    Args:
        hamiltonian: Small Hermitian matrix given as a nested sequence of complex elements.
        beta: Inverse temperature, a finite real number in the same units as the spectrum of H.

    Returns:
        tuple[tuple[complex, ...], ...]: Normalized Gibbs state matrix as a row-major nested tuple.
    """
    matrix = _check_hermitian(hamiltonian, "gibbs_state.hamiltonian")
    finite_real(beta, "gibbs_state.beta")
    values, vectors = _hermitian_eigendecomposition(matrix)
    d = len(matrix)
    weights = [math.exp(-beta * v) for v in values]
    partition = sum(weights)
    return tuple(
        tuple(
            sum(weights[k] * vectors[i][k] * vectors[j][k].conjugate() for k in range(d))
            / partition
            for j in range(d)
        )
        for i in range(d)
    )


def trace_distance(
    rho: Iterable[Iterable[complex]], sigma: Iterable[Iterable[complex]]
) -> float:
    """Trace distance T(ρ,σ) = ‖ρ−σ‖₁/2, computed via the Hermitian eigendecomposition of the difference matrix.

    Args:
        rho: First Hermitian matrix, usually a density matrix.
        sigma: Second Hermitian matrix, dimension must match rho.

    Returns:
        float: Trace distance, with range [0,1].
    """
    a = _check_hermitian(rho, "trace_distance.rho")
    b = _check_hermitian(sigma, "trace_distance.sigma")
    if len(a) != len(b):
        raise ValidationError("trace distance requires both matrices to have the same dimension")
    d = len(a)
    values, _ = _hermitian_eigendecomposition(
        [[a[i][j] - b[i][j] for j in range(d)] for i in range(d)]
    )
    return 0.5 * sum(abs(v) for v in values)


# ---------------------------------------------------------------------------
# Purification access view: the stable interface of the DM input model
# (reused by B2 quantum SDP).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PurificationAccess(OracleView):
    """Purification access view of a density matrix: prepares ``|ψ_ρ⟩`` with Tr_env ``|ψ⟩⟨ψ|`` = ρ.

    The operation starts from the all-zero state and prepares the purification
    of ρ on the two registers system and environment; the environment width is
    at least the bits needed by the rank of system (taking it equal to the
    system width is usually fine).
    """

    oracle_kind = "purification_access"
    operation: Operation

    def purification_access(self) -> PurificationAccess:
        """Structured protocol accessor for purification access, returning ``self``.

        A host object implementing a method of the same name returning
        ``PurificationAccess`` can be adapted by protocol, isomorphic to
        accessors such as ``state_preparation`` and ``block_encoding``.

        Returns:
            PurificationAccess: This view itself.
        """
        return self

    def __post_init__(self) -> None:
        """Validate the system/environment register signature and bit widths of the wrapped operation."""
        validate_signature(self.operation, ("system", "environment"), "PurificationAccess")
        if self.width < 1:
            raise ValidationError("the system register of PurificationAccess cannot be empty")

    @property
    def width(self) -> int:
        """Bit width of the system register, per the RIR register signature."""
        return next(
            r.type.width for r in self.operation.module.registers if r.name == "system"
        )

    @property
    def environment_width(self) -> int:
        """Bit width of the environment register, per the RIR register signature."""
        return next(
            r.type.width for r in self.operation.module.registers if r.name == "environment"
        )

    def describe(self) -> OracleSpec:
        """Return the ``OracleSpec`` snapshot with type ``purification_access``.

        The system width is recorded in main_qubit and the environment width
        in anc_qubit.

        Returns:
            OracleSpec: The snapshot, with type ``purification_access``.
        """
        from oracq.algorithms.input_model.contracts import describe_oracle

        base = describe_oracle(self.operation)
        return replace(
            base,
            type="purification_access",
            main_qubit=self.width,
            anc_qubit=self.environment_width,
        )

    @classmethod
    def from_state_preparation(cls, preparation: StatePreparation) -> PurificationAccess:
        """A pure state is the trivial purification: adapt a work-restored StatePreparation into PurificationAccess.

        The work register plays the role of environment; the preparation
        contract promises that work is restored to zero, so after tracing out
        the environment the state on system is still the original pure state.

        Args:
            preparation: Pure state preparation view with work restored to
                zero; when the work width is greater than 0, the clean_work
                promise is required.

        Returns:
            PurificationAccess: The trivial purification view whose environment is the original work register.
        """
        require_instance(
            preparation, StatePreparation, "PurificationAccess.from_state_preparation"
        )
        attributes = dict(preparation.operation.module.attributes)
        if preparation.work_width and attributes.get("clean_work") is not True:
            fail(
                "INPUT_PROMISE",
                "PurificationAccess.from_state_preparation.clean_work",
                True,
                attributes.get("clean_work"),
                "must promise that work is restored to zero after preparation so the pure state can serve as the trivial purification",
            )
        n, m = preparation.width, preparation.work_width
        b = Builder(
            _name("pure_state_purification", preparation.operation),
            {"system": Bits(n), "environment": Bits(m)},
            resources_for(("prep", preparation.operation)),
        )
        invoke(b, preparation.operation, "prep", target=b["system"], work=b["environment"])
        return cls(
            annotate(
                b.finish(),
                "unitary",
                density_model="purification_access",
                implementation="pure_state_adapter",
                zero_input=True,
            )
        )

    def as_state_preparation(self) -> StatePreparation:
        """Treat the purification operation as a whole as a StatePreparation on system⊕environment (for B2 composition).

        Returns:
            StatePreparation: Preparation view whose target state is the purified state and whose work width is 0.
        """
        n, m = self.width, self.environment_width
        b = Builder(
            _name("purification_as_state_preparation", self.operation),
            {"target": Bits(n + m), "work": Bits(0)},
            resources_for(("purification", self.operation)),
        )
        invoke(
            b,
            self.operation,
            "purification",
            system=b["target"][:n],
            environment=b["target"][n:],
        )
        return StatePreparation(
            annotate(
                b.finish(),
                "state_prep_isometry",
                zero_input=True,
                clean_work=True,
                density_model="purification_access",
            )
        )


@dataclass(frozen=True)
class ApproximatePurification(OracleView):
    """Post-selected approximate purification: on the signal == 0 branch, system and environment carry the approximate purified state.

    It differs from PurificationAccess by having a signal register and a
    controllable approximation error; algorithm parameters (such as beta and
    error for Gibbs preparation) are stored in module attributes and read via
    attributes.
    """

    oracle_kind = "approximate_purification"
    operation: Operation

    def approximate_purification(self) -> ApproximatePurification:
        """Structured protocol accessor for approximate purification, returning ``self``.

        A host object implementing a method of the same name returning
        ``ApproximatePurification`` can be adapted by protocol, isomorphic to
        accessors such as ``state_preparation`` and ``block_encoding``.

        Returns:
            ApproximatePurification: This view itself.
        """
        return self

    def __post_init__(self) -> None:
        """Validate the system/environment/signal register signature and bit widths of the wrapped operation."""
        validate_signature(
            self.operation, ("system", "environment", "signal"), "ApproximatePurification"
        )
        if self.width < 1:
            raise ValidationError("the system register of ApproximatePurification cannot be empty")

    @property
    def width(self) -> int:
        """Bit width of the system register, per the RIR register signature."""
        return next(
            r.type.width for r in self.operation.module.registers if r.name == "system"
        )

    @property
    def environment_width(self) -> int:
        """Bit width of the environment register, per the RIR register signature."""
        return next(
            r.type.width for r in self.operation.module.registers if r.name == "environment"
        )

    @property
    def signal_qubits(self) -> int:
        """Bit width of the signal register; post-selection requires its readout to be all zeros."""
        return next(
            r.type.width for r in self.operation.module.registers if r.name == "signal"
        )

    @property
    def attributes(self) -> dict[str, str | int | float | bool]:
        """Copy of the module attribute dictionary; stores algorithm parameters such as beta and error."""
        return dict(self.operation.module.attributes)

    @property
    def beta(self) -> float | None:
        """Inverse temperature β read from the module attributes; ``None`` when not recorded."""
        return cast("float | None", self.attributes.get("beta"))

    @property
    def error(self) -> float | None:
        """Polynomial uniform approximation error parameter read from the module attributes; ``None`` when not recorded."""
        return cast("float | None", self.attributes.get("error"))

    def describe(self) -> OracleSpec:
        """Return the ``OracleSpec`` snapshot with type ``approximate_purification``.

        The system width is recorded in main_qubit, and the sum of the
        environment and signal widths in anc_qubit.

        Returns:
            OracleSpec: The snapshot, with type ``approximate_purification``.
        """
        from oracq.algorithms.input_model.contracts import describe_oracle

        base = describe_oracle(self.operation)
        return replace(
            base,
            type="approximate_purification",
            main_qubit=self.width,
            anc_qubit=self.environment_width + self.signal_qubits,
        )


# ---------------------------------------------------------------------------
# Three paradigm layers: abstract open declarations and gate witness
# implementations.
# ---------------------------------------------------------------------------


def abstract_purification(
    name: str, width: int, environment_width: int | None = None, *, reversible: bool = True
) -> PurificationAccess:
    """Open declaration of the DM input model: a purification access slot preparing ``|ψ_ρ⟩``, for batched binding.

    environment_width defaults to width (every density matrix has a
    purification with an equal-width environment); the declaration is closed
    by binding a witness implementation such as gate_purification via
    linking.bind.

    Args:
        name: Slot name, referenced when binding a witness implementation via linking.bind.
        width: Bit width of the system register, range 1..63.
        environment_width: Bit width of the environment register, range 0..63; defaults to width.
        reversible: Whether the declaration slot also promises adjoint and controlled capabilities.

    Returns:
        PurificationAccess: The open declaration slot awaiting a witness implementation.
    """
    positive_integer(width, "abstract_purification.width", maximum=63)
    environment_width = width if environment_width is None else environment_width
    positive_integer(
        environment_width, "abstract_purification.environment_width", minimum=0, maximum=63
    )
    return PurificationAccess(
        declare(
            name,
            {"system": Bits(width), "environment": Bits(environment_width)},
            paradigm="unitary",
            attributes={"density_model": "purification_access", "zero_input": True},
            supports_adjoint=reversible,
            supports_controlled=reversible,
        )
    )


def gate_purification(
    rho: Iterable[Iterable[complex]], *, name: str | None = None
) -> PurificationAccess:
    """Purification witness for an explicit small density matrix: eigendecomposition ρ = Σ_j p_j ``|v_j⟩⟨v_j|`` followed by controlled preparation.

    The purified state is ``|ψ_ρ⟩ = Σ_j √p_j |v_j⟩_s |j⟩_e``; its amplitude
    vector is prepared on the concatenated register of system and environment
    by the multiplexed rotation tree of gate_state_prep; taking the partial
    trace over the environment returns exactly ρ.

    Args:
        rho: Explicit density matrix, which must be a positive semidefinite
        power-of-two square matrix of trace 1 and dimension at most 16.
        name: Name of the generated operation; derived from the matrix content by default.

    Returns:
        PurificationAccess: Purification access view implemented by the eigendecomposition witness.
    """
    matrix = _check_density_matrix(rho)
    d = len(matrix)
    n = (d - 1).bit_length()
    values, vectors = _hermitian_eigendecomposition(matrix)
    if values[-1] < -1e-7:
        fail(
            "INPUT_PROMISE",
            "gate_purification.rho",
            "positive semidefinite",
            values[-1],
            "the density matrix must be positive semidefinite",
        )
    amplitudes = [0j] * (d * d)
    for j, p in enumerate(values):
        root = math.sqrt(max(p, 0.0))
        for i in range(d):
            amplitudes[i + (j << n)] = root * vectors[i][j]
    preparation = gate_state_prep(amplitudes)
    b = Builder(
        name or _name("purification", matrix),
        {"system": Bits(n), "environment": Bits(n)},
        resources_for(("prep", preparation.operation)),
    )
    invoke(
        b,
        preparation.operation,
        "prep",
        target=fuse(b["system"], b["environment"]),
        work=b["system"][:0],
    )
    return PurificationAccess(
        annotate(
            b.finish(),
            "unitary",
            density_model="purification_access",
            implementation="gate_eigendecomposition",
            zero_input=True,
        )
    )


def maximally_mixed_purification(width: int, *, name: str | None = None) -> PurificationAccess:
    """Purification generator of the maximally mixed state I/2^n: the tensor product of n Bell pairs ``|Φ+⟩``.

    Args:
        width: Bit width of each of the system and environment registers, range 1..32.
        name: Name of the generated operation; ``bell_purification_{width}`` by default.

    Returns:
        PurificationAccess: Purification view preparing the tensor product of n Bell pairs.
    """
    positive_integer(width, "maximally_mixed_purification.width", maximum=32)
    b = Builder(
        name or f"bell_purification_{width}",
        {"system": Bits(width), "environment": Bits(width)},
    )
    for bit in range(width):
        b.h(b["system"][bit])
        b.xor(b["system"][bit], b["environment"][bit])
    return PurificationAccess(
        annotate(
            b.finish(),
            "unitary",
            density_model="purification_access",
            implementation="bell_pairs",
            zero_input=True,
        )
    )


# ---------------------------------------------------------------------------
# Gibbs state preparation: the QSVT purification route.
# ---------------------------------------------------------------------------


def _bessel_i(n: int, x: float) -> float:
    """Modified Bessel function of the first kind I_n(x), a pure Python power-series implementation."""
    term = (x / 2) ** n / math.factorial(n)
    total = term
    m = 0
    while term > 1e-18 * max(1.0, total) and m < 100000:
        term *= ((x / 2) ** 2) / ((m + 1) * (m + n + 1))
        total += term
        m += 1
    return total


def _gibbs_branches(
    c: float, error: float
) -> tuple[tuple[float, ...], tuple[float, ...], int, int]:
    """Even/odd Chebyshev truncations of g(x) = e^{−c(x+1)}: e^{−c}cosh(cx) and −e^{−c}sinh(cx).

    The truncation tail of each branch is controlled by 2e^{−c}·Σ_{k>d} I_k(c)
    ≤ error/8, so the total uniform error is at most error/4; returns
    (even-branch ascending coefficients, odd-branch ascending coefficients,
    even-branch degree, odd-branch degree).
    """
    kmax = min(_MAX_DEGREE, int(math.ceil(c)) + 8 * int(math.ceil(math.log10(8 / error))) + 8)
    ivals = [_bessel_i(k, c) for k in range(kmax + 2)]
    suffix = [0.0] * (kmax + 3)
    for j in range(kmax + 1, -1, -1):
        suffix[j] = suffix[j + 1] + ivals[j]

    def tail_ok(k: int) -> bool:
        """Check whether the Bessel truncation tail of degree ``k`` has been pushed within error/8."""
        return 2.0 * math.exp(-c) * suffix[k + 1] <= error / 8

    d_even = next((k for k in range(2, kmax + 1, 2) if tail_ok(k)), None)
    d_odd = next((k for k in range(1, kmax + 1, 2) if tail_ok(k)), None)
    if d_even is None or d_odd is None:
        raise ValidationError(
            "beta times alpha is too large: the Gibbs polynomial degree exceeds the synthesis limit; reduce beta or first shrink the spectral scale of H"
        )
    shift = math.exp(-c)
    even = [0.0] * (d_even + 1)
    for k in range(d_even // 2 + 1):
        coef = shift * (1.0 if k == 0 else 2.0) * ivals[2 * k]
        for i, v in enumerate(_chebyshev_t(2 * k)):
            even[i] += coef * v
    odd = [0.0] * (d_odd + 1)
    for k in range((d_odd + 1) // 2):
        coef = -2.0 * shift * ivals[2 * k + 1]
        for i, v in enumerate(_chebyshev_t(2 * k + 1)):
            odd[i] += coef * v
    return cast("tuple[float, ...]", _trim(even)), cast("tuple[float, ...]", _trim(odd)), d_even, d_odd


def _even_imag_candidates(f: Sequence[float]) -> Iterator[tuple[float, ...]]:
    """Enumerate candidate polynomials for the even-branch imaginary completion: the constant term and the ``x^{2m}`` term each take the saturation amplitude."""
    d = len(f) - 1
    a0 = math.sqrt(max(0.0, 1.0 - _eval(f, 0.0).real ** 2))
    a1 = math.sqrt(max(0.0, 1.0 - _eval(f, 1.0).real ** 2))
    for m in range(1, d // 2 + 1):
        yield tuple(
            (a0 if i == 0 else 0.0) + ((a1 - a0) if i == 2 * m else 0.0)
            for i in range(d + 1)
        )


def _odd_imag_candidates(f: Sequence[float]) -> Iterator[tuple[float, ...]]:
    """Enumerate candidate polynomials for the odd-branch imaginary completion: only the ``x^{2m+1}`` term takes the saturation amplitude."""
    d = len(f) - 1
    a1 = math.sqrt(max(0.0, 1.0 - _eval(f, 1.0).real ** 2))
    for m in range(0, (d - 1) // 2 + 1):
        yield tuple(a1 * (1.0 if i == 2 * m + 1 else 0.0) for i in range(d + 1))


def gibbs_purification(
    hamiltonian: BlockEncoding, beta: float, *, error: float = 0.01
) -> ApproximatePurification:
    """Approximate purification preparation of the Gibbs state ρ = e^{−βH}/Z (QSVT purification route).

    hamiltonian is the BlockEncoding of H, with the convention that the
    spectrum lies in [−α,α] (α = be_alpha); the target function on the
    spectral variable x = λ/α ∈ [−1,1] is g(x) = exp(−βα(x+1)/2) ∈ (0,1]. The
    construction has two steps: (a) prepare n Bell pairs on system and
    environment (the purification ``|Φ⟩`` of the maximally mixed state); (b)
    apply to system a QSVT block encoding whose zero-signal block is
    g(H/α)/(2s). Since the even branch e^{−c}cosh(cx) and the odd branch
    −e^{−c}sinh(cx) of g are both convex, the endpoint-matched imaginary
    completion always satisfies the unit-disk constraint (f²(x) does not
    exceed the chord between the endpoints), so the phase synthesis is always
    feasible. After post-selecting signal == 0, the state is proportional to
    (g(H/α) ⊗ I)``|Φ⟩``, and the reduced density matrix of system is exactly
    e^{−βH}/Z; the normalization factor 2s is unrelated to the partition
    function and does not affect the reduced state.

    error controls the polynomial uniform approximation error (each branch's
    truncation tail ≤ error/8); returns an ApproximatePurification whose
    module attributes include algorithm="gibbs_purification", beta, error,
    qsp_degree, and gibbs_scale = 2s. β = 0 degenerates to the maximally
    mixed purification.

    Args:
        hamiltonian: Block encoding of H, with spectrum contained in [−α,α] (α = be_alpha).
        beta: Inverse temperature, a nonnegative finite real number; 0 degenerates to the maximally mixed purification.
        error: Polynomial uniform approximation error, with range (0,1).

    Returns:
        ApproximatePurification: The approximate purification view whose reduced
    state on system after post-selecting signal == 0 is e^{−βH}/Z.
    """
    require_instance(hamiltonian, BlockEncoding, "gibbs_purification.hamiltonian")
    finite_real(beta, "gibbs_purification.beta", minimum=0)
    finite_real(error, "gibbs_purification.error", minimum=0, strict=True)
    if error >= 1:
        raise ValidationError("the approximation error must be between 0 and 1")
    n = hamiltonian.width
    if beta == 0:
        bells = maximally_mixed_purification(n)
        b = Builder(
            _name("gibbs_purification", hamiltonian.operation, 0.0, error),
            {"system": Bits(n), "environment": Bits(n), "signal": Bits(0)},
            resources_for(("bells", bells.operation)),
        )
        invoke(b, bells.operation, "bells", system=b["system"], environment=b["environment"])
        return ApproximatePurification(
            annotate(
                b.finish(),
                "unitary",
                density_model="approximate_purification",
                algorithm="gibbs_purification",
                beta=0.0,
                error=float(error),
                qsp_degree=0,
                gibbs_scale=1.0,
                success_condition="signal == 0",
            )
        )
    c = float(beta) * hamiltonian.alpha / 2.0
    g_even, g_odd, d_even, d_odd = _gibbs_branches(c, error)
    s = 1.5 * max(_sup_norm(g_even), _sup_norm(g_odd), 1e-3)
    f_even, f_odd = (
        cast("tuple[float, ...]", _scale(1.0 / s, g_even)),
        cast("tuple[float, ...]", _scale(1.0 / s, g_odd)),
    )
    phases_even = _synthesize_with_imag(f_even, _even_imag_candidates(f_even), "Gibbs even branch")
    phases_odd = _synthesize_with_imag(f_odd, _odd_imag_candidates(f_odd), "Gibbs odd branch")
    be_even = _real_qsvt_be(hamiltonian, phases_even)
    be_odd = _real_qsvt_be(hamiltonian, phases_odd)
    gibbs_be = linear_combination(1.0, be_even, 1.0, be_odd)
    b = Builder(
        _name("gibbs_purification", hamiltonian.operation, beta, error),
        {"system": Bits(n), "environment": Bits(n), "signal": Bits(gibbs_be.signal_qubits)},
        resources_for(("gibbs_be", gibbs_be.operation)),
    )
    for bit in range(n):
        b.h(b["system"][bit])
        b.xor(b["system"][bit], b["environment"][bit])
    invoke(b, gibbs_be.operation, "gibbs_be", target=b["system"], signal=b["signal"])
    return ApproximatePurification(
        annotate(
            b.finish(),
            "unitary",
            density_model="approximate_purification",
            algorithm="gibbs_purification",
            beta=float(beta),
            error=float(error),
            qsp_degree=max(d_even, d_odd),
            gibbs_scale=2.0 * s,
            success_condition="signal == 0",
            spectral_variable="x = eigenvalue/alpha in [-1,1]",
            target_function="g(x) = exp(-beta*alpha*(x+1)/2)",
        )
    )
