"""Grover 搜索、迭代算子与成功子空间的相干振幅放大。"""

from pyqecclang.algorithms.block_encoding import reflect_zero
from pyqecclang.algorithms.contracts import positive_integer, require_instance
from pyqecclang.algorithms.interfaces import checked_state_preparation
from pyqecclang.algorithms.operators import _name
from pyqecclang.algorithms.oracles import (
    StateOracle,
    XorDatabase,
    annotate,
    invoke,
    resources_for,
    uniform_state,
)
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse


def grover(phase_oracle, width, *, iterations=1, preparation=None):
    """生成带可替换初态的 Grover 搜索电路。

    Args:
        phase_oracle: 标记目标基态的相位 Operation，接口为 target，可另有 work。
        width: 搜索空间的目标位宽。
        iterations: 非负的 Grover 迭代次数。
        preparation: 初态制备；省略时使用均匀态。提供的 work 必须复净。

    Returns:
        StateOracle: target 为搜索结果，signal 保留相位查询和制备的工作空间。

    初态反射同时包括 target 与制备工作区。"""
    positive_integer(iterations, "grover.iterations", minimum=0)
    prep = (
        checked_state_preparation(preparation, adjoint=True)
        if preparation is not None
        else uniform_state(width)
    )
    if prep.width != width:
        raise ValidationError("Grover 初态宽度与搜索空间不一致")
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
        reflect_zero(b, fuse(b["target"], bw), positive=True)
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


def grover_iterate(preparation, marked):
    """返回 Q=A(2|0><0|-I)A†S_good；marked 是目标基态编号集合。"""
    from pyqecclang.algorithms.oracles import phase_marks

    prep = checked_state_preparation(preparation, adjoint=True)
    marker = phase_marks(prep.width, tuple(marked))
    b = Builder(
        _name("grover_iterate", prep.operation, marker),
        {"target": Bits(prep.width), "work": Bits(prep.work_width)},
        resources_for(("prep", prep.operation), ("marker", marker)),
    )
    invoke(b, marker, "marker", target=b["target"])
    with b.adjoint():
        invoke(b, prep.operation, "prep", target=b["target"], work=b["work"])
    reflect_zero(b, fuse(b["target"], b["work"]), positive=True)
    invoke(b, prep.operation, "prep", target=b["target"], work=b["work"])
    return b.finish()


def amplify_success(state, *, iterations=1):
    """相干放大态 oracle 的零信号成功子空间。

    Args:
        state: 支持伴随调用的 StateOracle。
        iterations: 非负放大次数。

    Returns:
        StateOracle: 公开接口与输入相同，成功条件仍为 signal==0。

    次数需要结合输入成功概率选择；过多迭代可能降低成功概率。"""
    require_instance(state, StateOracle, "amplify_success.state")
    positive_integer(iterations, "amplify_success.iterations", minimum=0)
    b = Builder(
        _name("amplify_success", state.operation, iterations),
        {"target": Bits(state.width), "signal": Bits(state.signal_qubits)},
        resources_for(("state", state.operation)),
        attributes={"algorithm": "amplitude_amplification", "success_condition": "signal == 0"},
    )
    invoke(b, state.operation, "state", target=b["target"], signal=b["signal"])
    with b.repeat(iterations):
        reflect_zero(b, b["signal"])
        with b.adjoint():
            invoke(b, state.operation, "state", target=b["target"], signal=b["signal"])
        reflect_zero(b, fuse(b["target"], b["signal"]), positive=True)
        invoke(b, state.operation, "state", target=b["target"], signal=b["signal"])
    return StateOracle(b.finish())
