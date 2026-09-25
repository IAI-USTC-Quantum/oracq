"""QPCA quantum principal component analysis: density matrix exponentiation (the LMR protocol) plus phase estimation.

Implements the core primitive of Lloyd–Mohseni–Rebentrost 2014 ("Quantum
principal component analysis", Nature Physics 10, 631). Density matrix
exponentiation exploits the eigenvalue structure of the SWAP operator:
perform a partial swap ``e^{-iΔt·SWAP}`` between one copy of ρ and the
system, and after discarding the copy the effective channel on the system is
the first-order approximation of ``e^{-iρΔt}`` (error ``O(Δt²)``); chaining
copies copies gives total time t = copies·Δt. QPCA runs phase estimation on
that unitary and reads out the spectrum of ρ.

Input model: copies of ρ are provided by a state preparation oracle (SP/QRAM
bindable at three layers); mixed states are adapted via
density.PurificationAccess.as_state_preparation (environment qubits stay out
of the swap, which has the same effect as taking the partial trace). Multiple
copies mean multiple invocations of the preparation oracle — the resource
precondition of QPCA — recorded into the attributes together with the QRAM
assumption.
"""

from __future__ import annotations

import math

from oracq.algorithms.common.fourier import inverse_qft
from oracq.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import StatePreparation, invoke, resources_for
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Ref, RegType, ValidationError


def _pauli_pair_rotation(b: Builder, first: Ref, second: Ref, angle: float, axis: str) -> None:
    """``e^{-i(angle/2)·P⊗P}`` with P given by axis ∈ {x, y, z}; an exact decomposition based on parity."""
    if axis == "x":
        b.h(first)
        b.h(second)
    elif axis == "y":
        # Rx(π/2) = H·Rz(π/2)·H rotates the Z basis to −Y (the two sign terms cancel).
        for ref in (first, second):
            b.h(ref)
            b.rz(ref, math.pi / 2)
            b.h(ref)
    b.xor(first, second)
    b.rz(second, angle)
    b.xor(first, second)
    if axis == "y":
        for ref in (first, second):
            b.h(ref)
            b.rz(ref, -math.pi / 2)
            b.h(ref)
    elif axis == "x":
        b.h(first)
        b.h(second)


def _partial_swap(b: Builder, first: Ref, second: Ref, angle: float) -> None:
    """``e^{-i·angle·SWAP}``: an exact decomposition of XX+YY+ZZ = 2·SWAP − I (the three axes commute)."""
    b.global_phase(-angle / 2)
    for axis in ("z", "x", "y"):
        _pauli_pair_rotation(b, first, second, angle, axis)


def _registers_of(operation: Operation) -> dict[str, RegType]:
    """Map from register names of the entry module to storage types."""
    return {r.name: r.type for r in operation.module.registers}


def density_matrix_exponentiation(
    preparation: StatePreparation,
    *,
    time: float,
    copies: int,
    swap_width: int | None = None,
    name: str | None = None,
) -> Operation:
    """LMR density matrix exponentiation: approximate ``e^{-iρ·time}`` on the system using copies copies of ρ.

    Args:
        preparation: State preparation of a ρ copy (a StatePreparation,
            possibly a purification adapter).
        time: Total evolution time t; each step is Δt = t/copies.
        copies: Number of copies; the first-order error O(t·Δt) decreases
            linearly with it.
        swap_width: Prefix bit width participating in the swap, defaulting to
            the preparation's entire target (purification scenarios should
            use ρ's system width, leaving the environment qubits in the copy
            out of the swap).
        name: Overrides the automatically generated module name.

    Returns:
        Operation: Registers system (swap_width bits) and copies
        (copies × prep.width bits). The system's input state is prepared by
        the caller; the copy registers are filled from all zeros by the
        preparation oracle."""
    require_instance(preparation, StatePreparation, "density_matrix_exponentiation.preparation")
    finite_real(time, "density_matrix_exponentiation.time", minimum=0, strict=True)
    positive_integer(copies, "density_matrix_exponentiation.copies", minimum=1, maximum=63)
    prep_width = preparation.width
    swap_width = prep_width if swap_width is None else swap_width
    positive_integer(
        swap_width, "density_matrix_exponentiation.swap_width", minimum=1, maximum=prep_width
    )
    if preparation.work_width:
        raise ValidationError("copy preparation requires a StatePreparation with zero work width; adapt purifications by concatenating registers first")
    step = float(time) / copies
    b = Builder(
        name or _name("dm_exponentiation", preparation.operation, time, copies),
        {"system": Bits(swap_width), "copies": Bits(copies * prep_width)},
        resources_for(("prep", preparation.operation)),
        attributes={
            "algorithm": "density_matrix_exponentiation",
            "copies": copies,
            "step_time": step,
            "error_scaling": "O(time * step_time)",
            "reference": "Lloyd-Mohseni-Rebentrost 2014, Nature Physics 10, 631",
        },
    )
    for k in range(copies):
        copy = b["copies"][k * prep_width : (k + 1) * prep_width]
        invoke(b, preparation.operation, "prep", target=copy, work=copy[:0])
        for bit in range(swap_width):
            _partial_swap(b, b["system"][bit], copy[bit], step)
    return b.finish()


def qpca(
    preparation: StatePreparation,
    *,
    precision: int,
    step_time: float,
    system: StatePreparation | None = None,
    swap_width: int | None = None,
    name: str | None = None,
) -> Operation:
    """QPCA principal component analysis: run phase estimation on ``e^{-iρ·step_time}`` and read out the spectrum of ρ.

    Args:
        preparation: State preparation of a ρ copy; ``2**precision − 1`` copies
            are consumed in total.
        precision: Number of phase register bits, in 1..6 (the copy count
            grows exponentially).
        step_time: Single-step evolution time Δt; it must satisfy λ·Δt ≪ 2π
            to avoid readout aliasing.
        system: Optional state preparation of the system input (default
            ``|0>``); different eigenstate inputs correspond to different
            eigenvalue peaks in the readout.
        swap_width: Prefix bit width participating in the swap, same semantics
            as density_matrix_exponentiation.
        name: Overrides the automatically generated module name.

    Returns:
        Operation: Registers system, copies, and phase. After reading out
        phase, decode the eigenvalue with eigenvalue_from_phase."""
    require_instance(preparation, StatePreparation, "qpca.preparation")
    positive_integer(precision, "qpca.precision", minimum=1, maximum=6)
    finite_real(step_time, "qpca.step_time", minimum=0, strict=True)
    prep_width = preparation.width
    swap_width = prep_width if swap_width is None else swap_width
    positive_integer(swap_width, "qpca.swap_width", minimum=1, maximum=prep_width)
    if preparation.work_width:
        raise ValidationError("copy preparation requires a StatePreparation with zero work width; adapt purifications by concatenating registers first")
    if system is not None:
        require_instance(system, StatePreparation, "qpca.system")
        if system.width != swap_width:
            raise ValidationError("the system input state width must equal swap_width")
    total_copies = (1 << precision) - 1
    resources = [("prep", preparation.operation)]
    if system is not None:
        resources.append(("system", system.operation))
    b = Builder(
        name or _name("qpca", preparation.operation, precision, step_time),
        {
            "system": Bits(swap_width),
            "copies": Bits(total_copies * prep_width),
            "phase": Bits(precision),
        },
        resources_for(*resources),
        attributes={
            "algorithm": "qpca",
            "readout_register": "phase",
            "decoder": "eigenvalue_from_phase",
            "copies": total_copies,
            "step_time": float(step_time),
            "reference": "Lloyd-Mohseni-Rebentrost 2014, Nature Physics 10, 631",
        },
    )
    if system is not None:
        invoke(b, system.operation, "system", target=b["system"], work=b["system"][:0])
    for k in range(total_copies):
        copy = b["copies"][k * prep_width : (k + 1) * prep_width]
        invoke(b, preparation.operation, "prep", target=copy, work=copy[:0])
    b.h(b["phase"])
    cursor = 0
    for bit in range(precision):
        with b.control(b["phase"][bit]):
            for _ in range(1 << bit):
                copy = b["copies"][cursor * prep_width : (cursor + 1) * prep_width]
                for qubit in range(swap_width):
                    _partial_swap(b, b["system"][qubit], copy[qubit], step_time)
                cursor += 1
    invoke(b, inverse_qft(precision), "iqft", target=b["phase"])
    return b.finish()


def eigenvalue_from_phase(value: int, precision: int, step_time: float) -> float:
    """Decode the phase readout of qpca into an eigenvalue estimate of ρ.

    The eigenphase of the unitary step ``e^{-iρΔt}`` is ``φ = -λΔt/(2π) (mod
    1)``; this function decodes on the λ ∈ [0, π/Δt) branch, and aliasing
    occurs when λ·Δt leaves that range (the caller's responsibility).

    Args:
        value: phase register readout, an integer in 0..2^precision−1.
        precision: Number of phase register bits, in 1..63.
        step_time: Duration Δt of a partial swap step, a positive real.

    Returns:
        float: Eigenvalue estimate λ of the density matrix.
    """
    positive_integer(precision, "eigenvalue_from_phase.precision", maximum=63)
    positive_integer(value, "eigenvalue_from_phase.value", minimum=0, maximum=(1 << precision) - 1)
    finite_real(step_time, "eigenvalue_from_phase.step_time", minimum=0, strict=True)
    if value > 1 << (precision - 1):
        value -= 1 << precision
    return -2.0 * math.pi * value / ((1 << precision) * step_time)
