"""算法展示目录：每个条目提供可运行的小实例及其读出说明。"""

from __future__ import annotations

import math
from dataclasses import dataclass

from pyqecclang.algorithms.basics.number_theory import modular_multiply, order_finding
from pyqecclang.algorithms.basics.oracle_algorithms import (
    affine_boolean_oracle,
    bernstein_vazirani,
    deutsch_jozsa,
    simon_sample,
)
from pyqecclang.algorithms.common.estimation import (
    amplitude_estimation,
    hadamard_test,
    phase_estimation,
    swap_test,
)
from pyqecclang.algorithms.common.fourier import fourier_add, qft
from pyqecclang.algorithms.common.hamiltonian import PauliHamiltonian, hamiltonian_simulation
from pyqecclang.algorithms.common.search import amplify_success, grover
from pyqecclang.algorithms.common.walks import cycle_walk
from pyqecclang.algorithms.input_model.oracles import (
    StateOracle,
    abstract_database,
    basis_state,
    gate_database,
    phase_marks,
    uniform_state,
)
from pyqecclang.algorithms.optimization.variational import (
    hardware_efficient_ansatz,
    pauli_measurement,
    qaoa_maxcut,
)
from pyqecclang.algorithms.qec.error_correction import repetition_encode, repetition_recover
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, fuse
from pyqecclang.infrastructure.linking import bind


@dataclass(frozen=True)
class GalleryCase:
    """算法展示目录中的单个演示条目。

    Attributes:
        name: 条目名，在目录内唯一。
        family: 所属算法类别，与 ``pyqecclang.algorithms`` 的模块名对应。
        operation: 演示用完整 ``Operation``。
        readout: 读出说明，描述期望的测量结果或判定方式。
        opened: 绑定前的开放 ``Operation``；仅演示 oracle 绑定的条目提供，否则为 ``None``。
    """

    name: str
    family: str
    operation: Operation
    readout: str
    opened: Operation | None = None


def algorithm_gallery() -> tuple[GalleryCase, ...]:
    """返回小型案例；不运行后端，不读写文件，也不请求外部数据。

    Returns:
        tuple[GalleryCase, ...]: 演示条目目录，每条含完整 ``Operation``、
        算法类别与读出说明。
    """
    cases: list[GalleryCase] = []

    def add(
        name: str,
        family: str,
        operation: Operation,
        readout: str,
        opened: Operation | None = None,
    ) -> None:
        """把一条演示条目追加到目录。"""
        cases.append(GalleryCase(name, family, operation, readout, opened))

    f = affine_boolean_oracle(3, 5, bias=1)
    add("deutsch_jozsa", "oracle_algorithms", deutsch_jozsa(f), "input 非零，判为平衡函数")
    opened = bernstein_vazirani(abstract_database("SecretFunction", 3, 1))
    closed = Operation.from_program(bind(opened, {"SecretFunction": f.operation}))
    add("bernstein_vazirani", "oracle_algorithms", closed, "input=5", opened)
    add(
        "simon",
        "oracle_algorithms",
        simon_sample(gate_database(2, 1, [0, 1, 1, 0])),
        "input 满足 y·3=0；用多次采样做 GF(2) 消元",
    )
    add("qft", "fourier", qft(3), "零输入变为均匀叠加")
    add("fourier_add", "fourier", fourier_add(3), "输出 b=(a+b) mod 8，a 保留")
    add("grover", "search", grover(phase_marks(2, [3]), 2).operation, "target=3，signal=0")
    b = Builder("GalleryQuarterSuccess", {"target": Bits(1), "signal": Bits(1)})
    b.ry(b["signal"], 2 * math.pi / 3)
    add(
        "amplitude_amplification",
        "search",
        amplify_success(StateOracle(b.finish())).operation,
        "一次放大后 signal=0 的概率为 1",
    )
    add(
        "amplitude_estimation",
        "estimation",
        amplitude_estimation(uniform_state(1), (1,), precision=3),
        "phase 为 2 或 6，对应概率 1/2",
    )
    add(
        "quantum_counting",
        "estimation",
        amplitude_estimation(uniform_state(2), (3,), precision=4),
        "将振幅估计乘以 4，得到标记个数的估计",
    )
    b = Builder("GalleryPhase", {"target": Bits(1)})
    b.gate("phase", b["target"], math.pi / 2)
    phase = b.finish()
    b = Builder("GalleryQPE", {"target": Bits(1), "phase": Bits(3)})
    b.x(b["target"])
    b.call(phase_estimation(phase, precision=3), target=b["target"], phase=b["phase"])
    add("phase_estimation", "estimation", b.finish(), "phase=2，表示相位 1/4")
    add(
        "hadamard_real",
        "estimation",
        hadamard_test(phase, basis_state(1, 1)),
        "probe 的 Z 期望为 0",
    )
    add(
        "hadamard_imag",
        "estimation",
        hadamard_test(phase, basis_state(1, 1), component="imag"),
        "probe 的 Z 期望为 1",
    )
    add(
        "swap_test",
        "estimation",
        swap_test(basis_state(1), basis_state(1, 1)),
        "两个正交态，probe=0 的概率为 1/2",
    )
    add(
        "ansatz",
        "variational",
        hardware_efficient_ansatz(2, (((0.3, 0.1), (0.5, -0.2)),)),
        "参数化纯态，作为优化器的电路生成步骤",
    )
    add(
        "qaoa_maxcut",
        "variational",
        qaoa_maxcut(2, ((0, 1, 1),), (math.pi / 2,), (math.pi / 8,)),
        "target=01 或 10，达到单边 MaxCut",
    )
    add(
        "vqe_pauli_measurement",
        "variational",
        pauli_measurement(uniform_state(1), "X"),
        "target=0，对应 X 期望为 1",
    )
    add("cycle_walk", "walks", cycle_walk(2, steps=2), "读取 position/coin；保留相干叠加")
    add(
        "modular_multiply",
        "number_theory",
        modular_multiply(2, 5),
        "x<5 时映射到 2x mod 5，其余不变",
    )
    add(
        "order_finding",
        "number_theory",
        order_finding(2, 3, precision=3),
        "phase 为 0 或 4；非零样本给出阶 2",
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
            "恢复逻辑 |+>，syndrome 保留错误信息",
        )
    h = PauliHamiltonian(((0.3, "I"), (0.7, "X")))
    add(
        "hamiltonian_trotter",
        "hamiltonian",
        hamiltonian_simulation(h, 0.4, steps=3).operation,
        "完整酉演化，包含全局相位",
    )
    return tuple(cases)
