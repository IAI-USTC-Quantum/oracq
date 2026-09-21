"""从一个 gate 到 QLSS、LCU、HamSim：算法自己的 Python 协议。"""

import argparse
import json
from pathlib import Path
from typing import Protocol, runtime_checkable

from pyqecclang import (
    Bits,
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
from pyqecclang.algorithms.common.hamiltonian import (
    PauliOperator,
    TrotterizableProtocol,
    TrotterTerm,
    hamiltonian_simulation,
)
from pyqecclang.algorithms.input_model.block_encoding import lcu
from pyqecclang.algorithms.input_model.interfaces import (
    BlockEncodingProtocol,
    StateOracleProtocol,
    StatePreparationProtocol,
    UnitaryProtocol,
)
from pyqecclang.algorithms.input_model.oracles import abstract_block_encoding, abstract_state_prep
from pyqecclang.algorithms.qlss.qlss import CostaConfig, make_costa_qlss
from pyqecclang.algorithms.qode.lchs import QuadraturePlan
from pyqecclang.algorithms.qode.ode import linear_qode


@runtime_checkable
class HasDiagonal(Protocol):
    """完全属于本应用的协议，语言无需知道它。"""

    def diagonal_values(self) -> tuple[float, ...]: ...


class GivenMatrix:
    def __init__(self):
        self.encoding = abstract_block_encoding("GivenA", 1, 0, 1.0)

    def block_encoding(self):
        return self.encoding

    def diagonal_values(self):
        return (1.0, 1.0)


class MyHamiltonian:
    hermitian = True

    def trotter_list(self):
        return (TrotterTerm(0.3, PauliOperator("I")), TrotterTerm(0.7, PauliOperator("X")))


def main():
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
    encoded_sum = lcu(((1, gate), (1, identity(1))))

    a = GivenMatrix()
    assert requires(a, HasDiagonal).diagonal_values() == (1, 1)
    info = a.block_encoding()
    print("A:", info.type, "main", info.main_qubit, "anc", info.anc_qubit, "alpha", info.alpha)
    qlss = make_costa_qlss(CostaConfig(steps=1))
    problem = LinearSystem(block=BlockSystem(a, gate, SpectralPromise(1, 1)), rhs_norm=1)
    report = qlss.check(problem)
    report.require()
    result = qlss(problem)
    assert isinstance(result, StateOracleProtocol)
    assert all(isinstance(result, p) for p in qlss.provides)
    closed = bind(result.operation.program(), {"GivenA": identity(1).operation})
    (root / "qlss.open.rir.json").write_text(dumps(result.operation.program()), encoding="utf-8")
    (root / "qlss.closed.rir.json").write_text(dumps(closed), encoding="utf-8")
    (root / "qlss.originir").write_text(export_originir(closed).text, encoding="utf-8")

    limited = abstract_state_prep("ForwardOnly", 1, reversible=False)
    rejected = qlss.check(LinearSystem(block=BlockSystem(a, limited, SpectralPromise(1, 1))))
    assert not rejected.ok
    print("QLSS rejected:", [i.path for i in rejected.issues])

    model = QODEProblem(scale(-1, identity(1)), gate, dissipative=True, initial_norm=1)
    lchs = linear_qode("lchs", plan=QuadraturePlan.cauchy(cutoff=1))
    schrodinger = linear_qode("schrodingerization")
    for solver in (lchs, schrodinger):
        solver.check(model, time=0.1).require()
        state = solver.solve(model, 0.1)
        (root / (solver.name + ".rir.json")).write_text(
            dumps(state.operation.program()), encoding="utf-8"
        )

    h = MyHamiltonian()
    assert isinstance(h, TrotterizableProtocol)
    assert not isinstance(h, BlockEncodingProtocol)
    evolution = hamiltonian_simulation(h, 0.4, steps=3)
    (root / "trotter.rir.json").write_text(dumps(evolution.operation.program()), encoding="utf-8")
    (root / "trotter.originir").write_text(
        export_originir(evolution.operation.program()).text, encoding="utf-8"
    )
    (root / "unitary_lcu.rir.json").write_text(
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
