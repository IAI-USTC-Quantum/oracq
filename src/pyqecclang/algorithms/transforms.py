"""qubitization、显式 QSVT 相位序列和 oblivious amplification 的组装。"""

from pyqecclang.algorithms.block_encoding import reflect_zero
from pyqecclang.algorithms.operators import _name
from pyqecclang.algorithms.oracles import invoke, resources_for
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits


def qubitization_walk(a):
    b = Builder(
        _name("qubitization_walk", a.operation),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        resources_for(("a", a.operation)),
        attributes={"algorithm": "qubitization_walk"},
    )
    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
    reflect_zero(b, b["signal"], positive=True)
    return b.finish()


def qsvt_sequence(a, phases):
    phases = tuple(float(p) for p in phases)
    b = Builder(
        _name("qsvt", a.operation, phases),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        resources_for(("a", a.operation)),
        attributes={"algorithm": "qsvt_sequence", "validation_stage": "paradigm"},
    )
    for i, phase in enumerate(phases):
        b.global_phase(-phase)
        if b["signal"].width:
            with b.control(b["signal"], 0):
                b.global_phase(2 * phase)
        else:
            b.global_phase(2 * phase)
        if i + 1 < len(phases):
            if i % 2:
                with b.adjoint():
                    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
            else:
                invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
    return b.finish()


def oblivious_amplification(a, iterations=1):
    b = Builder(
        _name("oaa", a.operation, iterations),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        resources_for(("a", a.operation)),
        attributes={"algorithm": "oaa"},
    )
    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
    with b.repeat(iterations):
        reflect_zero(b, b["signal"])
        with b.adjoint():
            invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
        reflect_zero(b, b["signal"])
        invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
    return b.finish()
