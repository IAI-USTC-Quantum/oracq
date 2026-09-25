"""Quantum semidefinite programming (QSDP) framework: Gibbs sampling plus trace
estimation plus a matrix multiplicative weights outer loop.

The structure follows Brandão–Svore 2017 ("Quantum speed-ups for semidefinite
programming", FOCS) and van Apeldoorn–Gilyén 2019: the solution of a
feasibility SDP is approached by matrix multiplicative weights (MMW)
iterations, and the key quantum subroutines of each round are (a) Gibbs state
preparation of the penalty Hamiltonian (reusing the QSVT purification route of
density.gibbs_purification) and (b) trace estimation ``Tr(A_i ρ)`` for the
observables A_i (a Hadamard-type probe circuit over the purified state plus a
block encoding).

The outer MMW driver is classical (consistent with the QSDP papers; the
quantum speedup lies exactly in the two inner primitives, Gibbs preparation
and trace estimation); the driver consumes ``Tr(A_i ρ)`` through an estimator
callback, defaulting to the classical reference implementation of
density.gibbs_state so small instances can be cross-checked, while the
quantum path is generated round by round by iteration_circuits.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.density import (
    ApproximatePurification,
    PurificationAccess,
    _check_hermitian,
    gibbs_purification,
    gibbs_state,
)
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import invoke, resources_for
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError


def trace_estimate_circuit(
    purification: PurificationAccess | ApproximatePurification,
    observable: BlockEncoding,
    *,
    component: str = "real",
    name: str | None = None,
) -> Operation:
    """Probe estimation circuit for ``Tr(M ρ)/α``: a controlled invocation of the block encoding on the purified state (Hadamard test structure).

    The Z expectation of probe equals ``Re/Im⟨ψ_ρ| (M/α ⊗ I) |ψ_ρ⟩``;
    ``Tr(M ρ)`` is decoded by trace_from_probe. The environment register is a
    bystander only (tracing it out yields ρ), the block encoding signal joins
    the invocation from all zeros, and the circuit performs no measurement.
    purification may also be an ApproximatePurification (approximate
    purification): the purified state then lives only on its signal == 0
    branch, the probe readout must be conditioned on that branch, and decoding
    uses trace_from_joint.

    Args:
        purification: Purification access handle of the Gibbs state; its system
            width must match the observable.
        observable: Block encoding of the observable M.
        component: Component to read, ``real`` or ``imag``.
        name: Name of the generated module; generated automatically by default.

    Returns:
        Operation: Trace estimation probe circuit read out from the probe
        register, without performing any measurement.
    """
    approximate = isinstance(purification, ApproximatePurification)
    if not approximate and not isinstance(purification, PurificationAccess):
        raise ValidationError("trace_estimate_circuit requires a PurificationAccess or an ApproximatePurification")
    require_instance(observable, BlockEncoding, "trace_estimate_circuit.observable")
    if observable.width != purification.width:
        raise ValidationError("the target width of the observable block encoding must equal the purification system width")
    if component not in ("real", "imag"):
        raise ValidationError("component must be real or imag")
    n, m, s = purification.width, purification.environment_width, observable.signal_qubits
    registers = {
        "system": Bits(n),
        "environment": Bits(m),
        "signal": Bits(s),
        "probe": Bits(1),
    }
    purif_args = {"system": "system", "environment": "environment"}
    if approximate:
        registers["purification_signal"] = Bits(cast("ApproximatePurification", purification).signal_qubits)
        purif_args["signal"] = "purification_signal"
    b = Builder(
        name or _name("trace_estimate", purification.operation, observable.operation, component),
        registers,
        resources_for(("purification", purification.operation), ("be", observable.operation)),
        attributes={
            "algorithm": "trace_estimate",
            "readout_register": "probe",
            "component": component,
            "be_alpha": observable.alpha,
            "decoder": "trace_from_joint" if approximate else "trace_from_probe",
        },
    )
    invoke(
        b,
        purification.operation,
        "purification",
        **{key: b[value] for key, value in purif_args.items()},
    )
    b.h(b["probe"])
    with b.control(b["probe"]):
        invoke(b, observable.operation, "be", target=b["system"], signal=b["signal"])
    if component == "imag":
        b.gate("phase", b["probe"], -math.pi / 2)
    b.h(b["probe"])
    return b.finish()


def trace_from_probe(probability_one: float, alpha: float) -> float:
    """Decode ``Tr(M ρ)`` from the probability of measuring 1 on probe: the Z expectation is 1 − 2p, multiplied by the block encoding α.

    Args:
        probability_one: Probability of measuring 1 on probe, in [0,1].
        alpha: Scale α of the observable block encoding, a positive real.

    Returns:
        float: Estimate of ``Tr(M ρ)``.
    """
    finite_real(probability_one, "trace_from_probe.probability_one", minimum=0)
    if probability_one > 1:
        raise ValidationError("the probe probability must lie between 0 and 1")
    finite_real(alpha, "trace_from_probe.alpha", minimum=0, strict=True)
    return alpha * (1.0 - 2.0 * probability_one)


def trace_from_joint(
    probe_expectation_joint: float, weight_signal_zero: float, alpha: float
) -> float:
    """Decoding for the approximate purification path: ``Tr(M ρ) = α · E[Z_probe·1_{signal=0}] / P(signal=0)``.

    probe_expectation_joint is the joint expectation of the probe Z expectation
    and the purification signal == 0 indicator, and weight_signal_zero is the
    probability of the purification signal == 0 branch.

    Args:
        probe_expectation_joint: Joint expectation of the probe Z expectation
            and the signal==0 indicator.
        weight_signal_zero: Probability of the purification signal==0 branch, a
            positive real.
        alpha: Scale α of the observable block encoding, a positive real.

    Returns:
        float: Estimate of ``Tr(M ρ)`` conditioned on the successful branch.
    """
    finite_real(weight_signal_zero, "trace_from_joint.weight_signal_zero", minimum=0, strict=True)
    finite_real(alpha, "trace_from_joint.alpha", minimum=0, strict=True)
    return alpha * probe_expectation_joint / weight_signal_zero


@dataclass(frozen=True)
class SdpInstance:
    """Feasibility SDP input model: find ``X ⪰ 0`` with ``Tr X = 1`` such that ``|Tr(A_i X) − b_i| ≤ ε``.

    constraints is a sequence of (A_i, b_i) pairs: A_i is an explicit small
    Hermitian matrix and b_i is a real number. The explicit-matrix path targets
    small instances; at scale, A_i should be replaced by a sparse-access or
    block-encoding access model, and the block encoding assembly of the penalty
    Hamiltonian is replaced accordingly (see the parameters of
    iteration_circuits).
    """

    constraints: tuple[tuple[tuple[tuple[complex, ...], ...], float], ...]

    def __post_init__(self) -> None:
        """Validate and normalize each constraint into an ``(A_i, b_i)`` pair, then store it back into the field."""
        constraints = tuple(self.constraints)
        if not constraints:
            raise ValidationError("SdpInstance requires at least one constraint")
        normalized = []
        for i, item in enumerate(constraints):
            if len(item) != 2:
                raise ValidationError("each constraint must be a pair of A_i and b_i")
            matrix, bound = item
            matrix = _check_hermitian(matrix, f"SdpInstance.constraints[{i}]")
            if len(matrix) != len(constraints[0][0]):
                raise ValidationError("all constraint matrices must share the same dimension")
            finite_real(bound, f"SdpInstance.constraints[{i}].b")
            normalized.append((matrix, float(bound)))
        object.__setattr__(self, "constraints", tuple(normalized))

    @property
    def width(self) -> int:
        """Number of system qubits."""
        return (len(self.constraints[0][0]) - 1).bit_length()

    @property
    def num_constraints(self) -> int:
        """Number of constraints m, i.e. the number of (A_i, b_i) pairs in ``constraints``."""
        return len(self.constraints)


def penalty_hamiltonian(
    instance: SdpInstance, weights: Iterable[float]
) -> tuple[tuple[complex, ...], ...]:
    """MMW penalty Hamiltonian ``H = Σ_i w_i (A_i − b_i I)`` (explicit small matrix).

    Args:
        instance: Feasibility SDP input model providing the (A_i, b_i)
            constraints.
        weights: Weight sequence in one-to-one correspondence with the
            constraints, of length equal to the number of constraints.

    Returns:
        tuple[tuple[complex, ...], ...]: The penalty matrix as row-nested
        tuples.
    """
    require_instance(instance, SdpInstance, "penalty_hamiltonian.instance")
    weights = tuple(weights)
    if len(weights) != instance.num_constraints:
        raise ValidationError("the number of weights must equal the number of constraints")
    d = len(instance.constraints[0][0])
    return tuple(
        tuple(
            sum(w * (a[i][j] - (b if i == j else 0)) for w, (a, b) in zip(weights, instance.constraints, strict=True))
            for j in range(d)
        )
        for i in range(d)
    )


def classical_estimator(
    hamiltonian_matrix: Iterable[Iterable[complex]],
    beta: float,
    instance: SdpInstance,
) -> tuple[tuple[tuple[complex, ...], ...], tuple[float, ...]]:
    """Classical reference estimation: compute all ``Tr(A_i ρ)`` from an explicit Gibbs state (for cross-checks and as the driver default).

    Returns (rho, estimates): rho is the density matrix (nested tuples) and
    estimates are the traces of the constraints.

    Args:
        hamiltonian_matrix: Explicit Hermitian penalty matrix.
        beta: Inverse temperature of the Gibbs distribution, a positive real.
        instance: SDP input model providing the (A_i, b_i) constraint sequence.

    Returns:
        tuple[tuple[tuple[complex, ...], ...], tuple[float, ...]]: (the Gibbs
        density matrix ρ, the sequence of trace estimates per constraint).
    """
    rho = gibbs_state(hamiltonian_matrix, beta)
    estimates = tuple(
        sum((a[i][j] * rho[j][i]).real for i in range(len(rho)) for j in range(len(rho)))
        for a, _ in instance.constraints
    )
    return rho, estimates


def qsdp_gibbs_solve(
    instance: SdpInstance,
    *,
    epsilon: float,
    max_iterations: int | None = None,
    estimator: (
        Callable[
            [Iterable[Iterable[complex]], float, SdpInstance],
            tuple[tuple[tuple[complex, ...], ...], tuple[float, ...]],
        ]
        | None
    ) = None,
) -> dict[str, tuple[tuple[complex, ...], ...] | tuple[float, ...] | int | bool]:
    """Matrix multiplicative weights (MMW) feasibility solving: Hedge iterations for the symmetric zero-sum game ``min_X max_w Σ_i w_i v_i``.

    The constraint semantics are one-sided inequalities ``Tr(A_i X) ≤ b_i``
    (violation ``v_i = Tr(A_i ρ) − b_i``, only positive violations are
    penalized); an equality constraint should be split into the two
    inequalities ``(A_i, b_i)`` and ``(−A_i, −b_i)``. Each round: the Gibbs
    state ρ_t of the penalty Hamiltonian ``H_t = Σ_i w_i (A_i − b_i I)``, the
    constraint violations ``v_i = Tr(A_i ρ_t) − b_i``, and the weights updated
    as ``w_i ∝ exp(η·v_i)``. The maximum violation of the averaged iterate
    ``ρ̄`` decreases as ``O(√(ln m / T))`` (the structure of Brandão–Svore
    2017; the quantum speedup lies in the two inner primitives of Gibbs
    preparation and trace estimation, and the classical loop of this driver
    matches the papers).

    Args:
        instance: SdpInstance.
        epsilon: Target violation bound; the learning rate is η = ε/4.
        max_iterations: Iteration cap, defaulting to ``ceil(64·ln(m+1)/ε²)``;
            stops early once the maximum violation of the averaged iterate is
            ≤ epsilon.
        estimator: Callback (hamiltonian_matrix, beta, instance) -> (rho,
            estimates), defaulting to classical_estimator; the quantum path
            should pass an estimator built on iteration_circuits.

    Returns:
        dict: rho (the averaged iterate, nested tuples), violations (final
        violations), iterations, and converged (maximum violation ≤ epsilon).
    """
    require_instance(instance, SdpInstance, "qsdp_gibbs_solve.instance")
    finite_real(epsilon, "qsdp_gibbs_solve.epsilon", minimum=0, strict=True)
    m = instance.num_constraints
    if max_iterations is None:
        max_iterations = math.ceil(64 * math.log(m + 1) / (epsilon * epsilon)) + 1
    positive_integer(max_iterations, "qsdp_gibbs_solve.max_iterations", minimum=1)
    estimator = estimator or classical_estimator
    eta = epsilon / 4.0
    log_weights = [0.0] * m
    d = len(instance.constraints[0][0])
    rho_average = [[0.0 + 0j] * d for _ in range(d)]
    iterations_run = 0

    def current_violations(count: int) -> tuple[float, ...]:
        """Violations of each constraint when the averaged iterate is ``rho_average / count``."""
        return tuple(
            sum((a[i][j] * rho_average[j][i]).real for i in range(d) for j in range(d)) / count
            - bound
            for a, bound in instance.constraints
        )

    converged = False
    for iteration in range(1, max_iterations + 1):
        shift = max(log_weights)
        weights = [math.exp(w - shift) for w in log_weights]
        total = sum(weights)
        weights = [w / total for w in weights]
        # Incremental form of H_t = η·Σ_{s<t} M_s: use the Gibbs state directly with β = 1 and the penalty matrix.
        hamiltonian = penalty_hamiltonian(instance, weights)
        rho, estimates = estimator(hamiltonian, 1.0, instance)
        for i in range(d):
            for j in range(d):
                rho_average[i][j] += rho[i][j]
        iterations_run = iteration
        for k, (_, bound) in enumerate(instance.constraints):
            log_weights[k] += eta * (estimates[k] - bound)
        if iteration % 32 == 0 or iteration == max_iterations:
            if max(current_violations(iteration)) <= epsilon:
                converged = True
                break
    violations = current_violations(iterations_run)
    converged = converged or max(violations) <= epsilon
    return {
        "rho": tuple(
            tuple(value / iterations_run for value in row) for row in rho_average
        ),
        "violations": violations,
        "iterations": iterations_run,
        "converged": converged,
    }


def iteration_circuits(
    instance: SdpInstance,
    weights: Iterable[float],
    beta: float,
    *,
    hamiltonian_encoding: Callable[[tuple[tuple[complex, ...], ...]], BlockEncoding],
    error: float = 0.05,
) -> tuple[ApproximatePurification, list[Operation]]:
    """Generate the quantum subroutines of a single MMW round: Gibbs purification plus a trace estimation circuit per constraint.

    hamiltonian_encoding is a callable that encodes the explicit penalty matrix
    into a block encoding (matrix_pauli_encoding works for small instances;
    switch to a sparse-access or low-rank access model at scale). Returns
    (ApproximatePurification, [trace_estimate Operation, ...]) for the quantum
    path's estimator to invoke round by round.

    Args:
        instance: Feasibility SDP input model providing the (A_i, b_i)
            constraints.
        weights: Constraint weight sequence for this round, of length equal to
            the number of constraints.
        beta: Inverse temperature of the Gibbs purification, a positive real.
        hamiltonian_encoding: Callback that encodes the explicit penalty matrix
            into a block encoding.
        error: Purification approximation error, a positive real.

    Returns:
        tuple[ApproximatePurification, list[Operation]]: The Gibbs purification
        handle and the list of trace estimation circuits per constraint.
    """
    require_instance(instance, SdpInstance, "iteration_circuits.instance")
    finite_real(beta, "iteration_circuits.beta", minimum=0, strict=True)
    finite_real(error, "iteration_circuits.error", minimum=0, strict=True)
    hamiltonian = penalty_hamiltonian(instance, weights)
    be = hamiltonian_encoding(hamiltonian)
    require_instance(be, BlockEncoding, "iteration_circuits.hamiltonian_encoding")
    purification = gibbs_purification(be, beta, error=error)
    from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding

    traces = [
        trace_estimate_circuit(purification, matrix_pauli_encoding(a))
        for a, _ in instance.constraints
    ]
    return purification, traces
