"""Parameterized quantum circuits, MaxCut QAOA, and Pauli measurement circuits for VQE."""

from __future__ import annotations

from collections.abc import Sequence
from typing import cast

from oracq.algorithms.input_model.contracts import finite_real, positive_integer
from oracq.algorithms.input_model.interfaces import (
    StatePreparationProtocol,
    checked_state_preparation,
)
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import invoke, resources_for
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError


def hardware_efficient_ansatz(
    width: int, layers: Sequence[Sequence[Sequence[float]]]
) -> Operation:
    """Generate a parameterized circuit of Ry/Rz layers with nearest-neighbor CNOTs.

    Args:
        width: Width of target.
        layers: Parameters shaped [layer][qubit][Ry,Rz]; all angles are given
            in radians.

    Returns:
        Operation: target only. Parameter optimization and repeated execution
        are organized by the caller."""
    positive_integer(width, "ansatz.width", maximum=64)
    layers = tuple(tuple(tuple(pair) for pair in layer) for layer in layers)
    if not layers or any(
        len(layer) != width or any(len(pair) != 2 for pair in layer) for layer in layers
    ):
        raise ValidationError("ansatz parameters require a non-empty layout indexed by layer, qubit, and the Ry Rz pair")
    for layer in layers:
        for pair in layer:
            for angle in pair:
                finite_real(angle, "ansatz.angle")
    b = Builder(
        _name("hardware_ansatz", width, layers),
        {"target": Bits(width)},
        attributes={"algorithm": "hardware_efficient_ansatz", "layers": len(layers)},
    )
    for layer in layers:
        for bit, (ry, rz) in enumerate(layer):
            b.ry(b["target"][bit], ry)
            b.rz(b["target"][bit], rz)
        for bit in range(width - 1):
            b.xor(b["target"][bit], b["target"][bit + 1])
    return b.finish()


def qaoa_maxcut(
    width: int,
    edges: Sequence[tuple[int, int, float]],
    gammas: Sequence[float],
    betas: Sequence[float],
) -> Operation:
    """Generate the QAOA cost/mixer circuit for MaxCut.

    Args:
        width: Number of graph vertices, also the width of target.
        edges: (u,v,weight) triples with nonnegative weights; self-loops are
            not accepted.
        gammas: Cost evolution angles per layer.
        betas: Mixer evolution angles per layer, of the same length as
            gammas.

    Returns:
        Operation: The QAOA circuit starting from the uniform state. Reading
        target yields a candidate bit string of a cut.

    The cost is Σw(1-ZuZv)/2. This function does not run a classical
    optimizer."""
    positive_integer(width, "qaoa.width", maximum=64)
    edges = tuple(cast("tuple[int, int, float]", tuple(edge)) for edge in edges)
    gammas, betas = tuple(gammas), tuple(betas)
    if not gammas or len(gammas) != len(betas):
        raise ValidationError("QAOA requires non-empty gamma and beta lists of equal length")
    for angle in gammas + betas:
        finite_real(angle, "qaoa.angle")
    for edge in edges:
        if len(edge) != 3:
            raise ValidationError("each MaxCut edge must be a u v weight triple")
        u, v, weight = edge
        positive_integer(u, "qaoa.u", minimum=0, maximum=width - 1)
        positive_integer(v, "qaoa.v", minimum=0, maximum=width - 1)
        finite_real(weight, "qaoa.weight", minimum=0)
        if u == v:
            raise ValidationError("MaxCut does not accept self-loops")
    b = Builder(
        _name("qaoa_maxcut", width, edges, gammas, betas),
        {"target": Bits(width)},
        attributes={"algorithm": "qaoa_maxcut", "layers": len(gammas)},
    )
    b.h(b["target"])
    for gamma, beta in zip(gammas, betas, strict=True):
        for u, v, weight in edges:
            b.global_phase(-gamma * weight / 2)
            b.xor(b["target"][u], b["target"][v])
            b.rz(b["target"][v], -gamma * weight)
            b.xor(b["target"][u], b["target"][v])
        for bit in range(width):
            b.gate("rx", b["target"][bit], 2 * beta)
    return b.finish()


def pauli_measurement(preparation: StatePreparationProtocol, word: str) -> Operation:
    """Rotate the prepared state into the specified Pauli measurement basis.

    Args:
        preparation: A state preparation with zero input and a clean work
            space.
        word: An I/X/Y/Z string of the same width; the first character
            corresponds to the least significant bit.

    Returns:
        Operation: target/work interface. Reading the Z parity of the non-I
        positions estimates the expectation of that Pauli word."""
    prep = checked_state_preparation(preparation)
    if len(word) != prep.width or any(letter not in "IXYZ" for letter in word):
        raise ValidationError("the Pauli measurement word must have the same width as the state")
    b = Builder(
        _name("pauli_measurement", prep.operation, word),
        {"target": Bits(prep.width), "work": Bits(prep.work_width)},
        resources_for(("prep", prep.operation)),
        attributes={
            "algorithm": "pauli_measurement",
            "pauli_word": word,
            "readout_register": "target",
        },
    )
    invoke(b, prep.operation, "prep", target=b["target"], work=b["work"])
    for bit, letter in enumerate(word):
        if letter == "Y":
            b.gate("phase", b["target"][bit], -1.5707963267948966)
        if letter in "XY":
            b.h(b["target"][bit])
    return b.finish()


def vqe_measurements(
    preparation: StatePreparationProtocol, terms: Sequence[tuple[float, str]]
) -> tuple[tuple[float, Operation], ...]:
    """Generate measurement circuits for each term of the VQE Hamiltonian.

    Args:
        preparation: An ansatz with instantiated parameters or another state
            preparation.
        terms: List of (real coefficient, Pauli word) pairs.

    Returns:
        tuple: List of (coefficient, measurement Operation) pairs. Energy
        aggregation and parameter optimization happen on the classical
        side."""
    circuits: list[tuple[float, Operation]] = []
    for coefficient, word in terms:
        finite_real(coefficient, "vqe.coefficient")
        circuits.append((coefficient, pauli_measurement(preparation, word)))
    if not circuits:
        raise ValidationError("VQE requires a non-empty list of Hamiltonian terms")
    return tuple(circuits)
