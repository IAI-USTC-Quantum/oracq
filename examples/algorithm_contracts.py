"""From a single gate to QLSS, LCU, and HamSim: algorithms as Python protocols."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Protocol, cast, runtime_checkable

from oracq import (
    Bits,
    BlockEncoding,
    BlockSystem,
    Builder,
    LinearSystem,
    QODEProblem,
    SpectralPromise,
    bind,
    dumps,
    export_originir,
    identity,
    requires,
    scale,
)
from oracq.algorithms.common.hamiltonian import (
    PauliOperator,
    TrotterizableProtocol,
    TrotterTerm,
    hamiltonian_simulation,
)
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.interfaces import (
    BlockEncodingProtocol,
    StateOracleProtocol,
    StatePreparationProtocol,
    UnitaryProtocol,
)
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    abstract_block_encoding,
    abstract_state_prep,
)
from oracq.algorithms.qlss.qlss import CostaConfig, make_costa_qlss
from oracq.algorithms.qode.lchs import QuadraturePlan
from oracq.algorithms.qode.ode import linear_qode


@runtime_checkable
class HasDiagonal(Protocol):
    """A protocol that belongs entirely to this application; the language does not need to know it."""

    def diagonal_values(self) -> tuple[float, ...]: ...


class GivenMatrix:
    def __init__(self) -> None:
        """Construct the host object, declaring matrix A as an open block encoding named GivenA."""
        self.encoding = abstract_block_encoding("GivenA", 1, 0, 1.0)

    def block_encoding(self) -> BlockEncoding:
        """Return the open block encoding declared at construction time."""
        return self.encoding

    def diagonal_values(self) -> tuple[float, float]:
        """Return the diagonal entries of matrix A."""
        return (1.0, 1.0)


class MyHamiltonian:
    hermitian = True

    def trotter_list(self) -> tuple[TrotterTerm, ...]:
        """Return the Trotter decomposition consisting of an I term and an X term."""
        return (TrotterTerm(0.3, PauliOperator("I")), TrotterTerm(0.7, PauliOperator("X")))


def main() -> None:
    """Demonstrate protocol checking and solving step by step, and write the QLSS, QODE, and HamSim artifacts to disk."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path, default=Path("out/algorithm-contracts"))
    root = parser.parse_args().output
    root.mkdir(parents=True, exist_ok=True)

    b = Builder("MyXGate", {"q": Bits(1)})
    b.x(b["q"])
    gate = b.finish()
    assert isinstance(gate, UnitaryProtocol)
    assert isinstance(gate, StatePreparationProtocol)
    assert isinstance(gate, BlockEncodingProtocol)
    # The gate satisfies the BlockEncoding/StatePreparation protocols at runtime (isinstance assertions above);
    # the builder annotates by concrete view type, so narrow to the known shape here.
    encoded_sum = lcu(((1, cast(BlockEncoding, gate)), (1, identity(1))))

    a = GivenMatrix()
    assert requires(a, HasDiagonal).diagonal_values() == (1, 1)
    info = a.block_encoding()
    print("A:", info.type, "main", info.main_qubit, "anc", info.anc_qubit, "alpha", info.alpha)
    qlss = make_costa_qlss(CostaConfig(steps=1))
    problem = LinearSystem(
        block=BlockSystem(
            cast(BlockEncoding, a), cast(StatePreparation, gate), SpectralPromise(1, 1)
        ),
        rhs_norm=1,
    )
    report = qlss.check(problem)
    report.require()
    result = qlss(problem)
    assert isinstance(result, StateOracleProtocol)
    assert all(isinstance(result, p) for p in qlss.provides)
    closed = bind(result.operation.program(), {"GivenA": identity(1).operation})
    (root / "qlss.open.rir.yaml").write_text(dumps(result.operation.program()), encoding="utf-8")
    (root / "qlss.closed.rir.yaml").write_text(dumps(closed), encoding="utf-8")
    (root / "qlss.originir").write_text(export_originir(closed).text, encoding="utf-8")

    limited = abstract_state_prep("ForwardOnly", 1, reversible=False)
    rejected = qlss.check(
        LinearSystem(block=BlockSystem(cast(BlockEncoding, a), limited, SpectralPromise(1, 1)))
    )
    assert not rejected.ok
    print("QLSS rejected:", [i.path for i in rejected.issues])

    model = QODEProblem(
        scale(-1, identity(1)), cast(StatePreparation, gate), dissipative=True, initial_norm=1
    )
    lchs = linear_qode("lchs", plan=QuadraturePlan.cauchy(cutoff=1))
    schrodinger = linear_qode("schrodingerization")
    for solver in (lchs, schrodinger):
        solver.check(model, time=0.1).require()
        state = solver.solve(model, 0.1)
        (root / (solver.name + ".rir.yaml")).write_text(
            dumps(state.operation.program()), encoding="utf-8"
        )

    h = MyHamiltonian()
    assert isinstance(h, TrotterizableProtocol)
    assert not isinstance(h, BlockEncodingProtocol)
    evolution = hamiltonian_simulation(h, 0.4, steps=3)
    (root / "trotter.rir.yaml").write_text(dumps(evolution.operation.program()), encoding="utf-8")
    (root / "trotter.originir").write_text(
        export_originir(evolution.operation.program()).text, encoding="utf-8"
    )
    (root / "unitary_lcu.rir.yaml").write_text(
        dumps(encoded_sum.operation.program()), encoding="utf-8"
    )
    (root / "contracts.json").write_text(
        json.dumps(
            {
                "A": info.spec.to_dict(),
                "qlss": qlss.contract.to_dict(),
                "accepted": report.to_dict(),
                "rejected": rejected.to_dict(),
                "lchs": lchs.contract.to_dict(),
                "schrodingerization": schrodinger.contract.to_dict(),
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print("Saved to", root)


if __name__ == "__main__":
    main()
