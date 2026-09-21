"""多项式 ODE 的有限阶张量提升、初态与物理通道选择。"""

from __future__ import annotations

from dataclasses import dataclass

from pyqecclang.algorithms.common.state_preparation import select_subspace
from pyqecclang.algorithms.input_model.block_encoding import lcu
from pyqecclang.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from pyqecclang.algorithms.input_model.interfaces import (
    as_state_preparation,
)
from pyqecclang.algorithms.input_model.operators import BlockEncoding, _name
from pyqecclang.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from pyqecclang.algorithms.qode._dynamics import tagged
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError


@dataclass(frozen=True)
class PolynomialODE:
    """多项式 ODE 输入模型与 Carleman 组装所需的输入声明。

    Attributes:
        width: 单层状态向量宽度 d（寄存器位数）。
        coefficients: (次数 p, 块编码 F_p) 二元组序列，次数不重复且非空；
            F_p 作用在 d^p 维张量幂上，宽度为 ``max(1, p) * width``。
        initial: 初态 u(0) 的制备。
        initial_norm: 初态范数，用于各张量层的加权，默认 1.0。

    Raises:
        ValidationError: 宽度/次数/范数取值非法、系数表为空或次数重复、
            F_p 宽度不匹配或初态布局不符时在构造时抛出。
    """

    width: int
    coefficients: tuple[tuple[int, BlockEncoding], ...]
    initial: StatePreparation
    initial_norm: float = 1.0

    def __post_init__(self):
        positive_integer(self.width, "PolynomialODE.width", maximum=64)
        object.__setattr__(self, "initial", as_state_preparation(self.initial))
        finite_real(self.initial_norm, "PolynomialODE.initial_norm", minimum=0)
        object.__setattr__(self, "coefficients", tuple(self.coefficients))
        if not self.coefficients or len({p for p, _ in self.coefficients}) != len(
            self.coefficients
        ):
            raise ValidationError("PolynomialODE 系数需要非空且次数不重复")
        if self.initial.width != self.width or self.initial_norm < 0:
            raise ValidationError("Carleman 初始数据布局无效")
        for order, coefficient in self.coefficients:
            positive_integer(order, "PolynomialODE.order", minimum=0)
            require_instance(coefficient, BlockEncoding, "PolynomialODE.coefficient")
            if order < 0 or coefficient.width != max(1, order) * self.width:
                raise ValidationError("F_p 必须为 d×d^p、补齐到 max(d,d^p) 的 BE")


def _carleman_term(coefficient, n, cutoff, output_level, order, position):
    level_bits = cutoff.bit_length()
    data_width, source_level = cutoff * n, output_level + order - 1
    b = Builder(
        _name("carleman_placement", coefficient.operation, cutoff, output_level, order, position),
        {"target": Bits(data_width + level_bits), "signal": Bits(coefficient.signal_qubits + 1)},
        resources_for(("f", coefficient.operation)),
        attributes={
            "carleman_row_level": output_level,
            "carleman_column_level": source_level,
            "tensor_position": position,
            "polynomial_order": order,
            "correctness": "pending",
        },
    )
    data, level, flag = (
        b["target"][:data_width],
        b["target"][data_width:],
        b["signal"][coefficient.signal_qubits],
    )
    b.x(flag)
    with b.control(level, source_level):
        if source_level * n < data.width:
            with b.control(data[source_level * n :], 0):
                b.x(flag)
        else:
            b.x(flag)
    if order == 0:
        for j in reversed(range(position, source_level)):
            b.swap(data[j * n : (j + 1) * n], data[(j + 1) * n : (j + 2) * n])
        group = data[position * n : (position + 1) * n]
    else:
        group = data[position * n : (position + order) * n]
    invoke(
        b, coefficient.operation, "f", target=group, signal=b["signal"][: coefficient.signal_qubits]
    )
    if order > 1:
        # F_p 输出保存在组内最低 n 位；其余零行移至数据区高端。
        for j in range(position + 1, output_level):
            b.swap(data[j * n : (j + 1) * n], data[(j + order - 1) * n : (j + order) * n])
    delta = source_level ^ output_level
    for bit in range(level_bits):
        if (delta >> bit) & 1:
            b.x(level[bit])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=coefficient.alpha))


def carleman_lift(problem, *, cutoff=2):
    """组装截断 Carleman 线性嵌入的块编码。

    对每个目标张量层与多项式次数枚举放置项，以等权 LCU 组合；输出属性
    记录 ``cutoff`` 与填充行零假设（首 d 行之外的补齐行、d^p 之外的列均为零）。

    Args:
        problem: ``PolynomialODE`` 输入模型。
        cutoff: Carleman 截断阶，必须为正。

    Returns:
        BlockEncoding: 作用在 ``cutoff * width`` 位数据加层级寄存器上的
        嵌入算子。

    Raises:
        ValidationError: problem 不是 ``PolynomialODE`` 或 cutoff 无效。
    """
    require_instance(problem, PolynomialODE, "carleman.problem")
    positive_integer(cutoff, "carleman.cutoff")
    if cutoff < 1:
        raise ValidationError("Carleman 截断阶必须为正")
    terms = []
    for k in range(1, cutoff + 1):
        for order, coefficient in problem.coefficients:
            source = k + order - 1
            if 0 <= source <= cutoff:
                for position in range(k):
                    terms.append(
                        (1, _carleman_term(coefficient, problem.width, cutoff, k, order, position))
                    )
    out = lcu(terms)
    return BlockEncoding(
        annotate(
            out.operation,
            "block_encoding",
            be_alpha=out.alpha,
            algorithm="carleman_lift",
            cutoff=cutoff,
            correctness="pending",
            coefficient_assumption="F_p padded rows outside first d and columns outside d^p are zero",
        )
    )


def carleman_initial(problem, *, cutoff=2):
    """构造张量初态的零输入制备。

    层级寄存器按 ``initial_norm`` 的幂加权制备，再受控地把初态复制到
    前 k 层，形成 u(0) 的张量幂向量。

    Args:
        problem: ``PolynomialODE`` 输入模型。
        cutoff: Carleman 截断阶，必须为正。

    Returns:
        StatePreparation: 目标为 ``cutoff * width`` 位数据加层级寄存器，
        工作区按 ``initial`` 的宽度逐层复用。

    Raises:
        ValidationError: problem 不是 ``PolynomialODE`` 或 cutoff 无效。
    """
    require_instance(problem, PolynomialODE, "carleman.problem")
    positive_integer(cutoff, "carleman.cutoff")
    n, lb = problem.width, cutoff.bit_length()
    weights = [problem.initial_norm**k for k in range(cutoff + 1)]
    weights += [0] * ((1 << lb) - len(weights))
    levels = gate_state_prep(weights)
    b = Builder(
        _name("carleman_initial", problem.initial.operation, problem.initial_norm, cutoff),
        {"target": Bits(n * cutoff + lb), "work": Bits(problem.initial.work_width * cutoff)},
        resources_for(("initial", problem.initial.operation)),
    )
    level = b["target"][n * cutoff :]
    invoke(b, levels.operation, target=level, work=level[:0])
    for k in range(1, cutoff + 1):
        with b.control(level, k):
            for j in range(k):
                invoke(
                    b,
                    problem.initial.operation,
                    "initial",
                    target=b["target"][j * n : (j + 1) * n],
                    work=b["work"][
                        j * problem.initial.work_width : (j + 1) * problem.initial.work_width
                    ],
                )
    return StatePreparation(
        annotate(
            b.finish(),
            "state_prep_isometry",
            zero_input=True,
            algorithm="carleman_tensor_initial",
            correctness="pending",
        )
    )


def carleman_qode(problem, time, linear_solver, *, cutoff=2):
    """把 Carleman 嵌入交给线性求解器并投影回第一层。

    先组装 ``carleman_lift`` 与 ``carleman_initial``，再由 ``linear_solver``
    演化到时刻 ``time``，最后选取第一张量层作为解态；截断尾项的影响
    以 ``truncation_assumption`` 记录为待验证。

    Args:
        problem: ``PolynomialODE`` 输入模型。
        time: 目标演化时刻，非负。
        linear_solver: 形如 ``(generator, initial, time) -> StateOracle``
            的可调用线性求解器。
        cutoff: Carleman 截断阶，必须为正。

    Returns:
        StateOracle: 宽度为 ``problem.width`` 的第一层解态。

    Raises:
        ValidationError: 输入类型或取值非法、``linear_solver`` 不可调用
            或其输出不是 ``StateOracle``。
    """
    require_instance(problem, PolynomialODE, "carleman.problem")
    finite_real(time, "carleman.time", minimum=0)
    if not callable(linear_solver):
        raise ValidationError("Carleman 需要可调用的线性求解器")
    generator, initial = (
        carleman_lift(problem, cutoff=cutoff),
        carleman_initial(problem, cutoff=cutoff),
    )
    state = linear_solver(generator, initial, time)
    require_instance(state, StateOracle, "carleman.linear_solver.output")
    selected = select_subspace(
        state, problem.width, 1 << ((cutoff - 1) * problem.width), label="carleman_level_one"
    )
    return StateOracle(
        tagged(
            selected.operation,
            "carleman_qode",
            cutoff=cutoff,
            truncation_assumption="Carleman tail pending",
        )
    )
