"可重复生成的范式案例目录：开放 IR、绑定计划与有限具体实例。"

from __future__ import annotations

from dataclasses import dataclass, field

from pyqecclang.algorithms.block_encoding import (
    direct_sum,
    kronecker_sum,
    lcu,
    matrix_pauli_encoding,
    pad_signal,
    tensor,
)
from pyqecclang.algorithms.estimation import phase_estimation
from pyqecclang.algorithms.legacy import make_lchs_qode
from pyqecclang.algorithms.ode import make_euler_history_qode
from pyqecclang.algorithms.operators import identity, pauli_x
from pyqecclang.algorithms.oracle_algorithms import deutsch_jozsa
from pyqecclang.algorithms.oracles import (
    abstract_block_encoding,
    abstract_database,
    abstract_sparse_access,
    abstract_state_prep,
    declare,
    diagonal_block_encoding,
    gate_database,
    gate_state_prep,
    phase_marks,
    qram_database,
    qram_state_angles,
    qram_state_prep,
    sparse_entry,
    sparse_location_gate,
    sparse_location_qram,
    uniform_state,
)
from pyqecclang.algorithms.pde import make_qpde
from pyqecclang.algorithms.qlss import CostaConfig, costa_qlss, make_costa_qlss
from pyqecclang.algorithms.search import grover
from pyqecclang.algorithms.sparse import (
    batch_lookup,
    reversible_lookup,
    sparse_block_encoding,
    word_rotation,
)
from pyqecclang.algorithms.transforms import (
    oblivious_amplification,
    qsvt_sequence,
    qubitization_walk,
)
from pyqecclang.applications.legacy import (
    qfvm_inputs,
    qfvm_step,
    qham_initial_vector,
    qham_lift_m1,
    qham_m1,
)
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, Program
from pyqecclang.infrastructure.linking import Binding, bind
from pyqecclang.infrastructure.readout import ReadoutAction


@dataclass(frozen=True)
class Case:
    """一条范式参考案例：开放程序、实现绑定与宿主内存数据。

    Attributes:
        name: 案例名，为 ``CASES`` 中的条目之一。
        program: 案例的 ``Program``，抽象槽位待 ``bindings`` 填充。
        bindings: 槽位名到 ``Binding`` 或 ``Operation`` 的映射；未列出的槽位保持开放。
        memory: QRAM 案例的宿主侧内存表，键为绑定资源映射所引用的名字。
        source: 出处标识元组（语言规范、规格测试或参考负载条目）。
        notes: 案例适用范围的限定说明元组。
        readout: 导出后按序执行的末端 ``ReadoutAction`` 读出动作元组。
    """

    name: str
    program: Program
    bindings: dict = field(default_factory=dict)
    memory: dict = field(default_factory=dict)
    source: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    readout: tuple[ReadoutAction, ...] = ()

    def closed(self):
        """按 ``bindings`` 绑定 ``program`` 的开放槽位，返回闭合后的 ``Program``。"""
        return bind(self.program, self.bindings)

    def artifact(self):
        """导出闭合程序的 OriginIR 产物。

        Returns:
            OriginIRArtifact: ``closed()`` 结果的 OriginIR 导出产物。

        OriginIR 导出器在本方法内延迟导入，核心代码不因此引入后端依赖。"""
        from pyqecclang.infrastructure.backends import export_originir

        return export_originir(self.closed())


CASES = (
    "bell",
    "ghz",
    "python_generators",
    "trotter_hamsim",
    "dj_gate",
    "dj_qram",
    "grover_gate",
    "grover_qram",
    "stateprep_gate",
    "stateprep_qram",
    "sparse_gate",
    "sparse_qram",
    "costa_gate",
    "costa_qram",
    "costa_sparse_qram",
    "qfvm_gate",
    "qfvm_qram",
    "qham_qode",
    "qham_qpde",
    "qpe",
    "qsvt",
    "oaa",
    "lchs",
    "heat_qode",
    "poisson_qlss",
    "carleman_step",
    "schrodingerisation",
    "be_algebra",
    "arithmetic",
    "batch_qram",
    "banked_qram",
    "register_views",
    "measure_reset",
)
"""``build_case`` 接受的全部案例名。

命名约定：多数案例以应用家族为前缀（如 ``dj_``、``grover_``、``qham_``）；
``_gate`` 后缀表示固定门级实现，``_qram`` 后缀表示绑定 QRAM 数据库并
携带宿主内存表的实现，``_qode``/``_qpde`` 后缀区分求解器封装层级；
其余无后缀的名字是原语或单一算法的演示案例。
"""


def _sparse_case(name):
    access = abstract_sparse_access("Sparse", 1, 2, 2)
    amplitude = declare(
        "SparseAmplitude", {"value": Bits(2), "amplitude": Bits(1)}, paradigm="reversible_function"
    )
    a = sparse_block_encoding(access, amplitude)
    binding = {"SparseAmplitude": word_rotation(2)}
    memory = {}
    if name.endswith("qram"):
        binding["Sparse_position"] = Binding(
            sparse_location_qram(1), {"forward": "positions", "inverse": "inverse_positions"}
        )
        binding["Sparse_entry"] = Binding(
            sparse_entry(qram_database(2, 2), 1), {"db__table": "entries"}
        )
        memory = {
            "positions": [0, 1, 1, 0],
            "inverse_positions": [0, 1, 1, 0],
            "entries": [0, 1, 1, 0],
        }
    else:
        binding["Sparse_position"] = sparse_location_gate(1, [[0, 1], [1, 0]])
        binding["Sparse_entry"] = sparse_entry(gate_database(2, 2, [0, 1, 1, 0]), 1)
    if name.startswith("costa"):
        rhs = abstract_state_prep("Rhs", 1)
        binding["Rhs"] = uniform_state(1).operation
        operation = costa_qlss(a, rhs, CostaConfig(steps=1)).operation
    else:
        operation = a.operation
    return Case(
        name,
        operation.program(),
        binding,
        memory,
        ("language-spec-v2 §E.4; reference-workloads 4/5; CKS §1.1",),
        ("T†SWAP T 稀疏候选构造；矩阵与归一化的验证留到下一阶段。",),
    )


def _qfvm_case(name):
    inputs = qfvm_inputs()
    state = qfvm_step(inputs, make_costa_qlss(CostaConfig(steps=1)))
    bindings, memory = {}, {}
    database_slots = [
        (inputs.position, "position", [0, 1, 1, 0]),
        (inputs.reverse_slot, "reverse_slot", [0, 0, 1, 1]),
        (inputs.column_position, "column_position", [0, 1, 1, 0]),
        (inputs.column_reverse_slot, "column_reverse_slot", [0, 0, 1, 1]),
        (inputs.flow, "flow", [0, 1]),
        (inputs.boundary, "boundary", [1, 0]),
    ]
    for slot, label, data in database_slots:
        impl = (
            qram_database(slot.address_width, slot.data_width)
            if name.endswith("qram")
            else gate_database(slot.address_width, slot.data_width, data)
        )
        if name.endswith("qram"):
            bindings[slot.operation.module.name] = Binding(impl.operation, {"table": label})
            memory[label] = data
        else:
            bindings[slot.operation.module.name] = impl.operation
    ins = {"row": 1, "slot": 1, "column": 1, "flow": 1, "boundary": 1}
    outs = {"value": 2, "status": 1}
    table = {i: (((i & 1) + ((i >> 3) & 1) + 2 * ((i >> 4) & 1)) % 4) | 4 for i in range(32)}
    physical = reversible_lookup(
        ins, outs, qram_database(5, 3) if name.endswith("qram") else gate_database(5, 3, table)
    )
    bindings[inputs.physical_entry.module.name] = (
        Binding(physical, {"db__table": "physical_entries"}) if name.endswith("qram") else physical
    )
    if name.endswith("qram"):
        memory["physical_entries"] = table
    bindings[inputs.amplitude.module.name] = word_rotation(2)
    bindings[inputs.residual.operation.module.name] = uniform_state(1).operation
    return Case(
        name,
        state.operation.program(),
        bindings,
        memory,
        ("reference-workloads 1: QFVM sparse inputs, T_L†SWAP T_R, QLSS/filter",),
        (
            "物理条目绑定为明确的玩具查表函数，不冒充完整 Roe 实现。",
            "本例验证双向访问、资源绑定和过滤求解器组装；流体与量子正确性待下一阶段。",
        ),
    )


def _qham_case(name):
    linear_impl = matrix_pauli_encoding([[-0.2, 0.1], [0.1, -0.2]])
    fold_impl = matrix_pauli_encoding([[0, 0.5, 0, 0], [0, 0, -0.5, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    linear = abstract_block_encoding("QhamLinear", 1, linear_impl.signal_qubits, linear_impl.alpha)
    fold = abstract_block_encoding("QhamFold", 2, fold_impl.signal_qubits, fold_impl.alpha)
    lifted = qham_lift_m1(linear, fold)
    initial = abstract_state_prep("QhamInitial", lifted.width)
    qode = make_euler_history_qode(make_costa_qlss(CostaConfig(steps=1)), steps=1)
    solver = make_qpde(qode) if name.endswith("qpde") else qode
    state = qham_m1(linear, fold, initial, solver, via_pde=name.endswith("qpde"))
    return Case(
        name,
        state.operation.program(),
        {
            "QhamLinear": linear_impl.operation,
            "QhamFold": fold_impl.operation,
            "QhamInitial": gate_state_prep(qham_initial_vector([0.3, 0.1])).operation,
        },
        {},
        ("reference-workloads 6; generator-protocols §8; QHAM m=1 lifted system",),
        (
            "m=1、两个空间点；提升维数 10，补齐到 16。",
            "通过 QODE 历史系统与 Costa/filter 生成；本阶段不验证 HAM 收敛与条件解态。",
        ),
    )


def build_case(name):
    """按名称构建一条范式参考案例。

    全部案例由固定常量生成，可重复构造；QRAM 案例同时给出宿主内存表。

    Args:
        name: 案例名，必须是 ``CASES`` 中的条目。

    Returns:
        Case: 对应的案例，含实现绑定、内存数据与读出动作（如适用）。

    Raises:
        KeyError: 案例名不在 ``CASES`` 中。
    """
    if name not in CASES:
        raise KeyError(name)
    if name in {"bell", "ghz"}:
        width = 2 if name == "bell" else 3
        b = Builder(name, {"q": Bits(width)})
        b.h(b["q"][0])
        for bit in range(1, width):
            b.xor(b["q"][0], b["q"][bit])
        return Case(name, b.finish().program(), source=("spec-tests 00-primitives",))
    if name == "python_generators":

        def count(n):
            return 1 if n < 2 else count(n - 1) + count(n - 2)

        width = count(3)

        def generate(depth):
            b = Builder(f"recursive_{depth}", {"q": Bits(width)}, attributes={"const_alpha": 2.5})
            if depth:
                with b.repeat(2):
                    b.call(generate(depth - 1), q=b["q"])
            else:
                for bit, angle in enumerate(tuple(0.1 * (i + 1) for i in range(width))):
                    b.rz(b["q"][bit], angle)
            return b.finish()

        return Case(name, generate(2).program(), source=("spec-tests 02-generics-const",))
    if name == "trotter_hamsim":
        from pyqecclang.algorithms.hamiltonian import trotter_hamsim

        return Case(
            name,
            trotter_hamsim([(0.5, "X"), (0.7, "Z")], 0.2).program(),
            source=("spec-tests 06-protocols/hamsim-trotter",),
        )
    if name == "carleman_step":
        implementation = matrix_pauli_encoding(
            [[0, 0, 0, 0], [0, -0.2, 0.1, 0], [0, 0, -0.4, 0.2], [0, 0, 0, -0.6]]
        )
        generator = abstract_block_encoding(
            "CarlemanGenerator", 2, implementation.signal_qubits, implementation.alpha
        )
        initial = abstract_state_prep("CarlemanInitial", 2)
        solver = make_euler_history_qode(make_costa_qlss(CostaConfig(steps=1)), steps=1)
        return Case(
            name,
            solver(generator, initial, 0.1).operation.program(),
            {
                "CarlemanGenerator": implementation.operation,
                "CarlemanInitial": gate_state_prep([1, 0.3, 0.09, 0.027]).operation,
            },
            {},
            ("spec-tests 07-scientific/carleman-step",),
            ("二次标量模型的三阶截断提升；只验收组装，不认证截断误差。",),
        )
    if name == "schrodingerisation":
        from pyqecclang.algorithms.legacy import make_schrodingerisation_qode

        hamiltonian = abstract_block_encoding("LiftedHamiltonian", 2, 1, 1.0)
        initial = abstract_state_prep("SchrodingerInitial", 1)
        solver = make_schrodingerisation_qode(
            lambda g: hamiltonian, lambda h, t: qsvt_sequence(h, [t, -t, t])
        )
        state = solver(identity(1), initial, 0.1)
        return Case(
            name,
            state.operation.program(),
            {
                "LiftedHamiltonian": pad_signal(pauli_x(2), 1).operation,
                "SchrodingerInitial": uniform_state(1).operation,
            },
            {},
            ("spec-tests 07-scientific/schrodingerization; QODE alternative",),
            ("Hamiltonian lift 是显式 oracle 边界；演示绑定只检查范式，不声称复现特定 PDE。",),
        )
    if name == "banked_qram":
        from pyqecclang.algorithms.oracles import banked_database

        slot = banked_database(2, 96, abstract=True)
        impl = banked_database(2, 96)
        return Case(
            name,
            slot.program(),
            {slot.module.name: Binding(impl, {"bank0": "low", "bank1": "high"})},
            {"low": {0: 5}, "high": {0: 7}},
            ("reference-workloads 1 packed flow words; register-level storage",),
            ("96 位逻辑数据拆成 64/32 位两个 bank；本例验收描述，不在密集模拟器执行。",),
        )
    if name == "register_views":
        from pyqecclang.infrastructure.ir import fuse

        b = Builder("views_demo", {"input": Bits(4), "output": Bits(4)})
        b.h(b["input"][:2])
        b.xor(fuse(b["input"][2:], b["input"][:2]), b["output"])
        return Case(name, b.finish().program(), source=("spec-tests 01-registers / 00-primitives",))
    if name == "measure_reset":
        b = Builder("readout_demo", {"q": Bits(2)})
        b.h(b["q"][0])
        b.xor(b["q"][0], b["q"][1])
        return Case(
            name,
            b.finish().program(),
            source=("spec-tests 00-primitives/measure-reset",),
            readout=(ReadoutAction("measure", "q"), ReadoutAction("reset", "q")),
        )
    if name.startswith("qfvm_"):
        return _qfvm_case(name)
    if name.startswith("qham_"):
        return _qham_case(name)
    if name.startswith("sparse_") or name == "costa_sparse_qram":
        return _sparse_case(name)
    if name.startswith("dj_"):
        slot = abstract_database("BooleanFunction", 2, 1)
        implementation = (
            qram_database(2, 1) if name.endswith("qram") else gate_database(2, 1, [0, 1, 1, 0])
        )
        return Case(
            name,
            deutsch_jozsa(slot).program(),
            {
                "BooleanFunction": Binding(implementation.operation, {"table": "truth"})
                if name.endswith("qram")
                else implementation.operation
            },
            {"truth": [0, 1, 1, 0]} if name.endswith("qram") else {},
            ("Deutsch–Jozsa 用户新增负载",),
        )
    if name.startswith("grover_"):
        from pyqecclang.algorithms.search import phase_from_database

        work = 1 if name.endswith("qram") else 0
        slot = declare(
            "Mark",
            {"target": Bits(2), **({"work": Bits(1)} if work else {})},
            paradigm="phase_oracle",
        )
        if work:
            impl = phase_from_database(qram_database(2, 1))
            binding = Binding(impl, {"db__table": "marks"})
            memory = {"marks": [0, 0, 0, 1]}
        else:
            binding, memory = phase_marks(2, [3]), {}
        return Case(
            name,
            grover(slot, 2).operation.program(),
            {"Mark": binding},
            memory,
            ("language-spec-v2 §15.1; spec-tests 03-oracle/grover2",),
        )
    if name.startswith("stateprep_"):
        amplitudes = (1, 2, 3, 4)
        if name.endswith("qram"):
            impl = qram_state_prep(2, 3)
            slot = abstract_state_prep("Prepare", 2, impl.work_width)
            binding = Binding(impl.operation, {"angles": "angles"})
            memory = {"angles": qram_state_angles(amplitudes, 3)}
        else:
            slot = abstract_state_prep("Prepare", 2)
            binding, memory = gate_state_prep(amplitudes).operation, {}
        return Case(
            name,
            slot.operation.program(),
            {"Prepare": binding},
            memory,
            ("spec-tests 03-oracle/state-prep-isometry; reference-workloads 4",),
        )
    if name in {"costa_gate", "costa_qram", "poisson_qlss"}:
        qram = name.endswith("qram")
        a = abstract_block_encoding("A", 1, 3, 1.0)
        bp = abstract_state_prep("B", 1, 4 if qram else 0)
        a_impl = diagonal_block_encoding(
            qram_database(1, 2) if qram else gate_database(1, 2, [0, 1])
        )
        b_impl = qram_state_prep(1, 3) if qram else gate_state_prep([1, 2])
        if name == "poisson_qlss":
            a_impl = matrix_pauli_encoding([[2, -1], [-1, 2]])
            a = abstract_block_encoding("A", 1, a_impl.signal_qubits, a_impl.alpha)
            b_impl = gate_state_prep([1, 0])
        bindings = {
            "A": Binding(a_impl.operation, {"db__table": "matrix_angles"})
            if qram
            else a_impl.operation,
            "B": Binding(b_impl.operation, {"angles": "rhs_angles"}) if qram else b_impl.operation,
        }
        memory = (
            {"matrix_angles": [0, 1], "rhs_angles": qram_state_angles([1, 2], 3)} if qram else {}
        )
        return Case(
            name,
            costa_qlss(a, bp).operation.program(),
            bindings,
            memory,
            ("reference-workloads 4/5; language-spec-v2 §E.10",),
            ("general walk 与实际 unary Dolph–Chebyshev LCU filtering；数值认证待后续。",),
        )
    if name in {"qpe", "qsvt", "oaa", "lchs", "heat_qode"}:
        impl = (
            matrix_pauli_encoding([[-0.2, 0.1], [0.1, -0.2]])
            if name == "heat_qode"
            else pad_signal(pauli_x(1), 1)
        )
        a = abstract_block_encoding("Operator", 1, impl.signal_qubits, impl.alpha)
        bindings = {"Operator": impl.operation}
        if name == "qpe":
            op = phase_estimation(qubitization_walk(a), precision=2)
        elif name == "qsvt":
            op = qsvt_sequence(a, [0.1, -0.2, 0.3])
        elif name == "oaa":
            op = oblivious_amplification(a)
        else:
            initial = abstract_state_prep("Initial", 1)
            bindings["Initial"] = uniform_state(1).operation
            if name == "heat_qode":
                solver = make_euler_history_qode(make_costa_qlss(CostaConfig(steps=1)), steps=1)
            else:
                solver = make_lchs_qode(
                    lambda g, t: qsvt_sequence(g, [t, -t, t]), [0.5, 1.0], [0.4, 0.6]
                )
            op = solver(a, initial, 0.1).operation
        return Case(
            name,
            op.program(),
            bindings,
            {},
            ("language-spec-v2 §E.7–E.11; spec-tests 05-qsvt-qpe/07-scientific",),
            ("此例验证调用范式；给定相位或生成元绑定不构成 PDE/指数函数精度保证。",),
        )
    if name == "be_algebra":
        a = abstract_block_encoding("A", 1, 0, 1.0)
        b = abstract_block_encoding("B", 1, 1, 2.0)
        combined = lcu([(1, direct_sum(a, b)), (0.5, tensor(a, b)), (-0.25, kronecker_sum(a))])
        return Case(
            name,
            combined.operation.program(),
            {"A": identity(1).operation, "B": _one_signal_identity(2.0)},
            {},
            ("reference-workloads 2; operation tuple/array heterogeneous LCU",),
        )
    if name == "arithmetic":
        slot = declare(
            "Multiply",
            {"left": Bits(2), "right": Bits(2), "output": Bits(4)},
            paradigm="reversible_function",
            attributes={"fraction_bits": 1},
        )
        table = {i: (i & 3) * (i >> 2) for i in range(16)}
        impl = reversible_lookup({"left": 2, "right": 2}, {"output": 4}, gate_database(4, 4, table))
        b = Builder("arithmetic_demo", {"a": Bits(2), "b": Bits(2), "out": Bits(4)})
        b.h(b["a"])
        b.call(slot, left=b["a"], right=b["b"], output=b["out"])
        return Case(
            name,
            b.finish().program(),
            {"Multiply": impl},
            {},
            ("reference-workloads 3; reversible arithmetic and fixed-point interpretation",),
        )
    if name == "batch_qram":
        slot = abstract_database("BatchData", 4, 8)
        return Case(
            name,
            batch_lookup(slot, 8).program(),
            {"BatchData": Binding(qram_database(4, 8).operation, {"table": "batch"})},
            {"batch": {0: 7, 3: 12}},
            ("language-spec-v2 §E.13; register arrays / paged QRAM",),
        )
    raise AssertionError(name)


def _one_signal_identity(alpha):
    from pyqecclang.algorithms.operators import scale

    return pad_signal(scale(alpha, identity(1)), 1).operation
