"""Algorithm gallery: each entry provides a small runnable example and its readout instructions."""

from __future__ import annotations

import math
from dataclasses import dataclass

from oracq.algorithms.basics.number_theory import modular_multiply, order_finding
from oracq.algorithms.basics.oracle_algorithms import (
    affine_boolean_oracle,
    bernstein_vazirani,
    deutsch_jozsa,
    simon_sample,
)
from oracq.algorithms.common.estimation import (
    amplitude_estimation,
    hadamard_test,
    phase_estimation,
    swap_test,
)
from oracq.algorithms.common.fourier import fourier_add, qft
from oracq.algorithms.common.hamiltonian import PauliHamiltonian, hamiltonian_simulation
from oracq.algorithms.common.search import amplify_success, grover
from oracq.algorithms.common.walks import cycle_walk
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    abstract_database,
    basis_state,
    gate_database,
    phase_marks,
    uniform_state,
)
from oracq.algorithms.optimization.variational import (
    hardware_efficient_ansatz,
    pauli_measurement,
    qaoa_maxcut,
)
from oracq.algorithms.qec.error_correction import repetition_encode, repetition_recover
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, fuse
from oracq.infrastructure.linking import bind


@dataclass(frozen=True)
class GalleryCase:
    """A single demo entry in the algorithm gallery.

    Attributes:
        name: Entry name, unique within the gallery.
        family: Algorithm family, matching a module name of ``oracq.algorithms``.
        operation: The complete ``Operation`` used for the demo.
        readout: Readout instructions describing the expected measurement outcome or decision rule.
        opened: The open ``Operation`` before binding; provided only by entries demonstrating oracle binding, otherwise ``None``.
    """

    name: str
    family: str
    operation: Operation
    readout: str
    opened: Operation | None = None


def algorithm_gallery() -> tuple[GalleryCase, ...]:
    """Return small demo cases; no backend is run, no files are read or written, and no external data is requested.

    Returns:
        tuple[GalleryCase, ...]: The gallery of demo entries, each with a
        complete ``Operation``, its algorithm family and readout instructions.
    """
    cases: list[GalleryCase] = []

    def add(
        name: str,
        family: str,
        operation: Operation,
        readout: str,
        opened: Operation | None = None,
    ) -> None:
        """Append one demo entry to the gallery."""
        cases.append(GalleryCase(name, family, operation, readout, opened))

    f = affine_boolean_oracle(3, 5, bias=1)
    add("deutsch_jozsa", "oracle_algorithms", deutsch_jozsa(f), "input is non-zero, judged balanced")
    opened = bernstein_vazirani(abstract_database("SecretFunction", 3, 1))
    closed = Operation.from_program(bind(opened, {"SecretFunction": f.operation}))
    add("bernstein_vazirani", "oracle_algorithms", closed, "input=5", opened)
    add(
        "simon",
        "oracle_algorithms",
        simon_sample(gate_database(2, 1, [0, 1, 1, 0])),
        "input satisfies y*3=0; use repeated samples for GF(2) elimination",
    )
    add("qft", "fourier", qft(3), "the zero input becomes a uniform superposition")
    add("fourier_add", "fourier", fourier_add(3), "output is b=(a+b) mod 8, with a preserved")
    add("grover", "search", grover(phase_marks(2, [3]), 2).operation, "target=3, signal=0")
    b = Builder("GalleryQuarterSuccess", {"target": Bits(1), "signal": Bits(1)})
    b.ry(b["signal"], 2 * math.pi / 3)
    add(
        "amplitude_amplification",
        "search",
        amplify_success(StateOracle(b.finish())).operation,
        "after one amplification the probability of signal=0 is 1",
    )
    add(
        "amplitude_estimation",
        "estimation",
        amplitude_estimation(uniform_state(1), (1,), precision=3),
        "phase is 2 or 6, corresponding to probability 1/2",
    )
    add(
        "quantum_counting",
        "estimation",
        amplitude_estimation(uniform_state(2), (3,), precision=4),
        "multiply the amplitude estimate by 4 to estimate the number of marked items",
    )
    b = Builder("GalleryPhase", {"target": Bits(1)})
    b.gate("phase", b["target"], math.pi / 2)
    phase = b.finish()
    b = Builder("GalleryQPE", {"target": Bits(1), "phase": Bits(3)})
    b.x(b["target"])
    b.call(phase_estimation(phase, precision=3), target=b["target"], phase=b["phase"])
    add("phase_estimation", "estimation", b.finish(), "phase=2, representing phase 1/4")
    add(
        "hadamard_real",
        "estimation",
        hadamard_test(phase, basis_state(1, 1)),
        "the Z expectation of probe is 0",
    )
    add(
        "hadamard_imag",
        "estimation",
        hadamard_test(phase, basis_state(1, 1), component="imag"),
        "the Z expectation of probe is 1",
    )
    add(
        "swap_test",
        "estimation",
        swap_test(basis_state(1), basis_state(1, 1)),
        "two orthogonal states; the probability of probe=0 is 1/2",
    )
    add(
        "ansatz",
        "variational",
        hardware_efficient_ansatz(2, (((0.3, 0.1), (0.5, -0.2)),)),
        "a parameterized pure state, serving as the circuit generation step of the optimizer",
    )
    add(
        "qaoa_maxcut",
        "variational",
        qaoa_maxcut(2, ((0, 1, 1),), (math.pi / 2,), (math.pi / 8,)),
        "target=01 or 10, attaining MaxCut on the single edge",
    )
    add(
        "vqe_pauli_measurement",
        "variational",
        pauli_measurement(uniform_state(1), "X"),
        "target=0, corresponding to an X expectation of 1",
    )
    add("cycle_walk", "walks", cycle_walk(2, steps=2), "read position/coin; the coherent superposition is preserved")
    add(
        "modular_multiply",
        "number_theory",
        modular_multiply(2, 5),
        "maps to 2x mod 5 when x<5, unchanged otherwise",
    )
    add(
        "order_finding",
        "number_theory",
        order_finding(2, 3, precision=3),
        "phase is 0 or 4; a non-zero sample yields order 2",
    )
    for error in ("bit", "phase"):
        b = Builder("GalleryRepetition_" + error, {"target": Bits(1), "syndrome": Bits(2)})
        b.h(b["target"])
        b.call(repetition_encode(error=error), target=b["target"], syndrome=b["syndrome"])
        b.gate("x" if error == "bit" else "z", fuse(b["target"], b["syndrome"])[1])
        b.call(repetition_recover(error=error), target=b["target"], syndrome=b["syndrome"])
        add(
            "repetition_" + error,
            "error_correction",
            b.finish(),
            "recovers the logical |+>; the syndrome retains the error information",
        )
    h = PauliHamiltonian(((0.3, "I"), (0.7, "X")))
    add(
        "hamiltonian_trotter",
        "hamiltonian",
        hamiltonian_simulation(h, 0.4, steps=3).operation,
        "the full unitary evolution, including the global phase",
    )
    return tuple(cases)
