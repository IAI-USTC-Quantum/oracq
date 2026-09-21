"""qubitization、显式 QSVT 相位序列和 oblivious amplification 的组装。"""

from pyqecclang.algorithms.input_model.block_encoding import reflect_zero
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.algorithms.input_model.oracles import invoke, resources_for
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits


def qubitization_walk(a):
    """组装块编码 ``a`` 的 qubitization walk 操作。

    先正向调用 ``a``，再对 signal 做关于零子空间的正反射，得到
    ``(2 P_0 - I) U`` 形式的 walk 算子；寄存器沿用 ``a`` 的 target 与
    signal 签名。

    Args:
        a: 输入 ``BlockEncoding``。

    Returns:
        Operation: 保留对 ``a`` 的模块调用的 walk 操作。
    """
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
    """按显式相位序列对块编码 ``a`` 组装 QSVT 线路。

    相位按时间顺序排列（``phases[0]`` 最先作用），``len(phases)`` 个相位
    伴随 ``len(phases) - 1`` 次 ``a`` 调用，正向与逆向交替；每个相位由
    全局相位与 signal 零子空间上的受控相位合成。相位约定与 ``qsvt``
    模块的相位合成一致。

    Args:
        a: 输入 ``BlockEncoding``。
        phases: 弧度相位序列，元素可为任何可转 ``float`` 的值。

    Returns:
        Operation: 对 ``a`` 的调用保留为模块调用的 QSVT 操作。
    """
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
    """对块编码 ``a`` 组装 oblivious amplitude amplification。

    先正向调用 ``a``；随后每轮依次执行 signal 零反射、逆向调用 ``a``、
    再次零反射和正向调用，用于放大零信号投影分量。

    Args:
        a: 输入 ``BlockEncoding``。
        iterations: 放大轮数，默认为一轮。

    Returns:
        Operation: 迭代以 RIR ``Repeat`` 保留的操作。
    """
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
