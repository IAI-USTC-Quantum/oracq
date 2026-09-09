"""D-J、Grover、QPE 和 QSVT 的模块化生成形式。"""

import math

from ..builder import Builder
from ..combinators import reflect_zero
from ..ir import Bits, ValidationError
from ..library import _name
from ..oracles import StateOracle, XorDatabase, annotate, invoke, resources_for, uniform_state


def deutsch_jozsa(function: XorDatabase):
    if function.data_width != 1:
        raise ValidationError("D-J oracle 必须只有一个结果位")
    b = Builder(
        _name("deutsch_jozsa", function.operation),
        {"input": Bits(function.address_width), "answer": Bits(1)},
        resources_for(("function", function.operation)),
        attributes={
            "algorithm": "deutsch_jozsa",
            "readout_register": "input",
            "validation_stage": "paradigm",
        },
    )
    b.x(b["answer"])
    b.h(b["answer"])
    b.h(b["input"])
    invoke(b, function.operation, "function", address=b["input"], data=b["answer"])
    b.h(b["input"])
    return b.finish()


def grover(phase_oracle, width, *, iterations=1, preparation=None):
    prep = preparation or uniform_state(width)
    phase_work = next((r.type.width for r in phase_oracle.module.registers if r.name == "work"), 0)
    b = Builder(
        _name("grover", phase_oracle, prep.operation, iterations),
        {"target": Bits(width), "signal": Bits(phase_work + prep.work_width)},
        resources_for(("phase", phase_oracle), ("prep", prep.operation)),
        attributes={
            "algorithm": "grover",
            "iterations": iterations,
            "validation_stage": "paradigm",
        },
    )
    pw, bw = b["signal"][:phase_work], b["signal"][phase_work:]
    invoke(b, prep.operation, "prep", target=b["target"], work=bw)
    with b.repeat(iterations):
        arguments = {"target": b["target"]}
        if any(r.name == "work" for r in phase_oracle.module.registers):
            arguments["work"] = pw
        invoke(b, phase_oracle, "phase", **arguments)
        with b.adjoint():
            invoke(b, prep.operation, "prep", target=b["target"], work=bw)
        reflect_zero(b, b["target"], positive=True)
        invoke(b, prep.operation, "prep", target=b["target"], work=bw)
    return StateOracle(b.finish())


def phase_from_database(database: XorDatabase):
    if database.data_width != 1:
        raise ValidationError("谓词数据库必须有一个输出位")
    b = Builder(
        _name("phase_from_database", database.operation),
        {"target": Bits(database.address_width), "work": Bits(1)},
        resources_for(("db", database.operation)),
    )
    invoke(b, database.operation, "db", address=b["target"], data=b["work"])
    b.z(b["work"])
    with b.adjoint():
        invoke(b, database.operation, "db", address=b["target"], data=b["work"])
    return annotate(b.finish(), "phase_oracle")


def qft(width):
    b = Builder(f"qft_{width}", {"target": Bits(width)})
    for high in reversed(range(width)):
        b.h(b["target"][high])
        for low in reversed(range(high)):
            with b.control(b["target"][low]):
                b.gate("phase", b["target"][high], math.pi / (1 << (high - low)))
    for bit in range(width // 2):
        b.swap(b["target"][bit], b["target"][width - bit - 1])
    return b.finish()


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


def phase_estimation(operation, *, precision=2):
    registers = {r.name: r.type for r in operation.module.registers}
    if "phase" in registers:
        raise ValidationError("被调接口占用了 phase 参数名")
    b = Builder(
        _name("qpe", operation, precision),
        {**registers, "phase": Bits(precision)},
        resources_for(("u", operation)),
        attributes={"algorithm": "qpe"},
    )
    b.h(b["phase"])
    for bit in range(precision):
        with b.control(b["phase"][bit]):
            with b.repeat(1 << bit):
                invoke(b, operation, "u", **{name: b[name] for name in registers})
    with b.adjoint():
        invoke(b, qft(precision), target=b["phase"])
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
    with b.repeat(iterations):
        invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
        reflect_zero(b, b["signal"])
        with b.adjoint():
            invoke(b, a.operation, "a", target=b["target"], signal=b["signal"])
        reflect_zero(b, b["signal"])
    return b.finish()
