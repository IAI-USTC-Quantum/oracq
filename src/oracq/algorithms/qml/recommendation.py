"""Kerenidis-Prakash quantum recommendation system (arXiv:1603.08675): QSVE plus threshold projection plus sampling.

The data plane is QMatrix's sample-and-query structure
(algorithms/qdata.py). Following Lemma 5.3 of the paper, build W = U·V with
U = Ũ R₁ Ũ⁻¹ and V = Ṽ R₀ Ṽ⁻¹: Ũ maps (i, 0) to (i, the normalized vector
of row i) via the row tree, Ṽ maps (0, j) to (the user distribution Ã, j)
via the row-norm root tree, and R₀/R₁ are reflections about the zero basis
states of the row/item registers. Phase estimation on W yields θ satisfying
cos(θ_i/2) = σ_i/‖A‖_F; σ̂ is estimated per phase word and the flag is
flipped by threshold (the deterministic projection version of §5.3 Alg 2),
and after inverse phase estimation (flag, item) is measured: the item
register distribution on the flag=1 branch is exactly the recommendation
sampling distribution (§6).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from oracq.algorithms.common.estimation import phase_estimation
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.qdata import QMatrix
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.execution import RegisterState
from oracq.infrastructure.ir import QRAM, Bits, Ref, ValidationError


def sigma_from_phase(value: int, precision: int, frobenius: float) -> float:
    """Phase readout t → singular value estimate σ̂ = ‖A‖_F times the absolute value of cos(π t/2^precision).

    The eigenvalues of W come in pairs e^{±iθ} (the two rotation directions of
    the same σ), and the mirror phase 2^precision−t must map to the same
    singular value, hence the absolute value.

    Args:
        value: Phase register readout, an integer in 0..2^precision−1.
        precision: Number of phase register bits.
        frobenius: Frobenius norm ‖A‖_F of the matrix.

    Returns:
        float: Singular value estimate σ̂.
    """
    if not 0 <= value < 1 << precision:
        raise ValidationError("the phase readout is out of range")
    return frobenius * abs(math.cos(math.pi * value / (1 << precision)))


def _reflect_zero(builder: Builder, register: Ref) -> None:
    """Reflection about the zero basis state (twice the projector minus the identity): flip all with X, multi-controlled Z, flip all with X; for a single bit this is just Z."""
    if register.width == 1:
        builder.z(register)
        return
    for bit in range(register.width):
        builder.x(register[bit])
    with builder.control(register[1:], (1 << (register.width - 1)) - 1):
        builder.z(register[0])
    for bit in range(register.width):
        builder.x(register[bit])


@dataclass(frozen=True)
class KPRecommendationConfig:
    """precision is the number of phase register bits; sigma is the singular value threshold, defaulting to 0.5·‖A‖_F."""

    precision: int = 4
    sigma: float | None = None


@dataclass(frozen=True)
class RecommendationResult:
    """KP recommendation sampling circuit and its readout contract.

    Attributes:
        operation: Recommendation sampling circuit; registers row, item,
            phase, and flag, with QRAM resources row_angles and root_angles.
        matrix: The QMatrix input used to build the circuit.
        user: Target user index.
        precision: Number of phase register bits.
        sigma: Singular value threshold actually applied; 0.5·frobenius when
            config leaves it unspecified.
        frobenius: Frobenius norm of the matrix, used when decoding with
            sigma_from_phase.
    """

    operation: Operation
    matrix: QMatrix
    user: int
    precision: int
    sigma: float
    frobenius: float

    def memories(self) -> dict[str, dict[int, int]]:
        """Extract the initial values of the two QRAM angle banks required by the circuit.

        Returns:
            dict: Keys row_angles (the row tree) and root_angles (the
            row-norm root tree); values map addresses to angle words, for the
            execution entry point to bind the QRAM resources declared by the
            circuit."""
        snapshot = self.matrix.snapshot()
        return {"row_angles": snapshot["row_angles"], "root_angles": snapshot["root_angles"]}

    def readout(self, state: RegisterState) -> tuple[float, dict[int, float]]:
        """Reduce a simulation or execution result to (success probability, conditional recommendation distribution).

        Args:
            state: Amplitude state obtained from register-level simulation or
                execution.

        Returns:
            tuple[float, dict[int, float]]: The total probability of the
            flag=1 branch and the conditional recommendation distribution
            over items on that branch.
        """
        registers = self.operation.module.registers
        index = {r.name: i for i, r in enumerate(registers)}
        items: dict[int, float]
        success, items = 0.0, {}
        for key, amplitude in state.amplitudes.items():
            probability = abs(amplitude) ** 2
            if key[index["flag"]]:
                success += probability
                item = key[index["item"]]
                items[item] = items.get(item, 0.0) + probability
        distribution = {item: p / success for item, p in items.items()} if success else {}
        return success, distribution


def _walk_unitary(matrix: QMatrix) -> Operation:
    """W = Ũ R₁ Ũ⁻¹ · Ṽ R₀ Ṽ⁻¹; registers row/item, resources share names with the QMatrix banks."""
    r, c, aw = matrix.rows, matrix.cols, matrix.angle_width
    b = Builder(
        _name("kp_walk", matrix.rows, matrix.cols, matrix.angle_width),
        {"row": Bits(r), "item": Bits(c)},
        {"row_angles": QRAM(r + c, aw), "root_angles": QRAM(r, aw)},
        attributes={"algorithm": "kp_recommendation_walk"},
    )
    amp = matrix.amplitude_preparation()
    row = matrix.row_preparation()
    amp_work = b.local("amp_work", Bits(aw))
    row_work = b.local("row_work", Bits(aw))

    def sandwich(
        prep: Operation, register: Ref, work: Ref, resource: dict[str, str]
    ) -> None:
        """Zero-basis reflection conjugated by prep: adjoint first, reflection in the middle, forward call to finish."""
        with b.adjoint():
            b.call(prep, row=b["row"], item=b["item"], work=work, resources=resource)
        _reflect_zero(b, register)
        b.call(prep, row=b["row"], item=b["item"], work=work, resources=resource)

    # First V = Ṽ R₀ Ṽ⁻¹, then U = Ũ R₁ Ũ⁻¹; the operator product is U·V.
    sandwich(amp, b["row"], amp_work, {"root_angles": "root_angles"})
    sandwich(row, b["item"], row_work, {"row_angles": "row_angles"})
    return b.finish()


def kp_recommendation(
    matrix: QMatrix, user: int, config: KPRecommendationConfig | None = None
) -> RecommendationResult:
    """Generate the recommendation sampling circuit for user user; returns a RecommendationResult carrying the readout contract.

    Args:
        matrix: Recommendation matrix input model with the sample-and-query
            structure.
        user: Target user index, in 0..number of rows − 1.
        config: Phase bit count and singular value threshold configuration;
            defaults to 4 phase bits and threshold 0.5·‖A‖_F.

    Returns:
        RecommendationResult: Result containing the recommendation sampling
        circuit, the initial QRAM angle bank values, and the readout
        contract.
    """
    config = config or KPRecommendationConfig()
    if not isinstance(matrix, QMatrix):
        raise ValidationError("kp_recommendation requires a QMatrix input")
    if not 1 <= config.precision <= 12:
        raise ValidationError("the recommendation phase register width must be between 1 and 12")
    r, c, aw = matrix.rows, matrix.cols, matrix.angle_width
    sigma = config.sigma
    if sigma is None:
        sigma = 0.5 * matrix.frobenius
    if not 0 < sigma <= matrix.frobenius:
        raise ValidationError("the singular value threshold must be positive and at most the Frobenius norm")
    walk = _walk_unitary(matrix)
    qpe = phase_estimation(walk, precision=config.precision)
    prep = matrix.row_state_prep(user)
    b = Builder(
        _name("kp_recommendation", matrix.rows, matrix.cols, aw, user, config.precision, sigma),
        {
            "row": Bits(r),
            "item": Bits(c),
            "phase": Bits(config.precision),
            "flag": Bits(1),
        },
        {"row_angles": QRAM(r + c, aw), "root_angles": QRAM(r, aw)},
        attributes={
            "algorithm": "kp_recommendation",
            "correctness": "pending",
            "user": user,
            "sigma": sigma,
            "precision": config.precision,
            "frobenius": matrix.frobenius,
            "reference": "arXiv:1603.08675",
            "assumptions": "; ".join(
                (
                    "nonnegative fixed-point entries",
                    "rotation-angle quantization at angle_width bits",
                    "deterministic flag projection; no amplitude amplification",
                )
            ),
        },
    )
    qpe_resources = {f"u__{name}": name for name in ("row_angles", "root_angles")}
    amp = matrix.amplitude_preparation()
    prep_work = b.local("prep_work", Bits(aw))
    amp_work = b.local("amp_work", Bits(aw))
    b.call(prep.operation, target=b["item"], work=prep_work, resources={"row_angles": "row_angles"})
    b.call(amp, row=b["row"], item=b["item"], work=amp_work, resources={"root_angles": "root_angles"})
    b.call(qpe, row=b["row"], item=b["item"], phase=b["phase"], resources=qpe_resources)
    for value in range(1 << config.precision):
        if sigma_from_phase(value, config.precision, matrix.frobenius) >= sigma:
            with b.control(b["phase"], value):
                b.x(b["flag"])
    with b.adjoint():
        b.call(qpe, row=b["row"], item=b["item"], phase=b["phase"], resources=qpe_resources)
    with b.adjoint():
        b.call(
            amp,
            row=b["row"],
            item=b["item"],
            work=amp_work,
            resources={"root_angles": "root_angles"},
        )
    return RecommendationResult(
        b.finish(),
        matrix,
        user,
        config.precision,
        sigma,
        matrix.frobenius,
    )
