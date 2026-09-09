"""开放 input model 的 QODE/QPDE 生成器；截断与算法正确性统一待核验。"""

from __future__ import annotations

import cmath
import json
import math
from dataclasses import dataclass
from functools import partial

from ..builder import Builder
from ..combinators import adjoint_be, lcu, tensor
from ..ir import Bits, ValidationError
from ..library import BlockEncoding, _name, identity, product, scale
from ..oracles import (
    StateOracle,
    StatePreparation,
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from .solvers import apply_be_to_state, select_subspace


def tagged(operation, algorithm, **metadata):
    return annotate(
        operation,
        "unitary",
        algorithm=algorithm,
        correctness="pending",
        validation_stage="paradigm",
        **metadata,
    )


@dataclass(frozen=True)
class HermitianParts:
    """A=L+iH；Hermitian/半正定性质是输入模型声明，不由语言证明。"""

    hermitian: BlockEncoding
    h: BlockEncoding

    def __post_init__(self):
        if self.hermitian.width != self.h.width:
            raise ValidationError("Hermitian parts 宽度不匹配")

    @classmethod
    def from_operator(cls, a):
        adj = adjoint_be(a)
        return cls(lcu([(0.5, a), (0.5, adj)]), lcu([(-0.5j, a), (0.5j, adj)]))


@dataclass(frozen=True)
class LinearODE:
    parts: HermitianParts
    initial: StatePreparation
    label: str = "du_dt_equals_minus_A_u"

    def __post_init__(self):
        if self.parts.hermitian.width != self.initial.width:
            raise ValidationError("线性 ODE 初态和算子宽度不匹配")


@dataclass(frozen=True)
class QuadraturePlan:
    nodes: tuple[float, ...]
    weights: tuple[complex, ...]
    kernel: str = "user_supplied"

    def __post_init__(self):
        if not self.nodes or len(self.nodes) != len(self.weights):
            raise ValidationError("离散节点和权重长度不符")
        if any(not math.isfinite(x) for x in self.nodes):
            raise ValidationError("节点必须有限")

    @classmethod
    def cauchy(cls, cutoff=2, spacing=1.0):
        if cutoff < 0 or spacing <= 0:
            raise ValidationError("Cauchy 离散参数无效")
        nodes = tuple(k * spacing for k in range(-cutoff, cutoff + 1))
        return cls(nodes, tuple(spacing / (math.pi * (1 + k * k)) for k in nodes), "finite_cauchy")


@dataclass(frozen=True)
class ContourPlan:
    """QST Eq.12 主级数；辅助极点与有限截断余项在报告中明确保留。"""

    a: float = 1.0
    cutoff: int = 2
    poles: tuple[complex, ...] = (2j, -1 + 1j, 1j, 1 + 1j)

    def __post_init__(self):
        if not math.isfinite(self.a) or self.a <= 0 or self.cutoff < 0:
            raise ValidationError("CBMD a 必须为正且截断非负")
        if len(set(self.poles)) != len(self.poles) or any(
            p == -1j or p.imag == 0 for p in self.poles
        ):
            raise ValidationError("CBMD 当前要求非实互异简单辅助极点，且避开 -i")

    @property
    def nodes(self):
        return tuple(k / self.a for k in range(-self.cutoff, self.cutoff + 1))

    @property
    def weights(self):
        numerator = math.expm1(-2 * math.pi * self.a)
        return tuple(
            numerator
            / (
                self.a
                * 2
                * math.pi
                * 1j
                * (q + 1j)
                * math.prod((q - p) / (-1j - p) for p in self.poles)
            )
            for q in self.nodes
        )

    @property
    def auxiliary_coefficients(self):
        numerator = math.expm1(-2 * math.pi * self.a)
        return tuple(
            numerator
            / (
                (cmath.exp(-2 * math.pi * p * self.a * 1j) - 1)
                * math.prod((p - other) / (-1j - other) for other in self.poles if other != p)
            )
            for p in self.poles
        )

    def metadata(self):
        def pair(z):
            return [complex(z).real, complex(z).imag]

        return json.dumps(
            {
                "a": self.a,
                "cutoff": self.cutoff,
                "nodes": self.nodes,
                "poles": [pair(x) for x in self.poles],
                "weights": [pair(x) for x in self.weights],
                "auxiliary_coefficients": [pair(x) for x in self.auxiliary_coefficients],
                "omitted": ["auxiliary_nonhermitian_evolutions", "infinite_series_tail"],
                "assumption": "L>=0 and norm(integral L dt)<=2*pi*a",
                "source": "QST 11 035027 (2026), Eq.12-13; arXiv:2511.10267v3",
            },
            separators=(",", ":"),
        )


def taylor_hamiltonian(hamiltonian, time, *, degree=2):
    """可闭合的普通 Hamiltonian-function BE；可替换为 QSP/HamSim protocol。"""
    if degree < 0 or not math.isfinite(time):
        raise ValidationError("Taylor 阶数/演化时间无效")
    powers, current = [(1, identity(hamiltonian.width))], identity(hamiltonian.width)
    for k in range(1, degree + 1):
        current = product(hamiltonian, current)
        powers.append(((-1j * time) ** k / math.factorial(k), current))
    out = lcu(powers)
    return BlockEncoding(
        annotate(
            out.operation,
            "block_encoding",
            be_alpha=out.alpha,
            algorithm="truncated_taylor_hamiltonian_function",
            degree=degree,
            time=float(time),
            correctness="pending",
            success_condition="signal == 0",
        )
    )


def _lcu_dynamics(model, time, nodes, weights, hamiltonian_function, algorithm, **metadata):
    terms = []
    for node, weight in zip(nodes, weights, strict=True):
        hk = lcu([(1, model.parts.h), (node, model.parts.hermitian)])
        encoded = hamiltonian_function(hk, time)
        if not isinstance(encoded, BlockEncoding):
            raise ValidationError(
                "Hamiltonian-function protocol 必须返回 BlockEncoding，保留 alpha"
            )
        terms.append((weight, encoded))
    evolution = lcu(terms)
    state = apply_be_to_state(evolution, model.initial)
    return StateOracle(
        tagged(
            state.operation,
            algorithm,
            evolution_alpha=evolution.alpha,
            branch_count=len(nodes),
            input_assumption="L>=0; autonomous; homogeneous",
            **metadata,
        )
    )


def lchs_qode(model, time, *, plan=None, hamiltonian_function=taylor_hamiltonian):
    plan = plan or QuadraturePlan.cauchy()
    return _lcu_dynamics(
        model,
        time,
        plan.nodes,
        plan.weights,
        hamiltonian_function,
        "lchs_qode",
        quadrature_kernel=plan.kernel,
        quadrature_nodes=json.dumps(plan.nodes),
        remainder="finite quadrature pending",
    )


def cbmd_qode(model, time, *, plan=None, hamiltonian_function=taylor_hamiltonian):
    plan = plan or ContourPlan()
    return _lcu_dynamics(
        model,
        time,
        plan.nodes,
        plan.weights,
        hamiltonian_function,
        "cbmd_qode",
        contour_plan=plan.metadata(),
        remainder="auxiliary pole contribution and truncation pending",
    )


def cbmd_function(a, nodes, residue_weights, hermitian_function):
    """通用 f(A) 组装点：Hermitian function protocol 保持开放，不偷换为矩阵求逆。"""
    parts = HermitianParts.from_operator(a)
    terms = [
        (weight, hermitian_function(lcu([(node, parts.h), (1, parts.hermitian)])))
        for node, weight in zip(nodes, residue_weights, strict=True)
    ]
    result = lcu(terms)
    return BlockEncoding(
        annotate(
            result.operation,
            "block_encoding",
            be_alpha=result.alpha,
            algorithm="cbmd_matrix_function",
            correctness="pending",
            residue_sign_convention="caller supplies target-side weights from contour identity",
        )
    )


def qft(width):
    b = Builder("qft_" + str(width), {"target": Bits(width), "work": Bits(0)})
    for j in reversed(range(width)):
        b.h(b["target"][j])
        for k in reversed(range(j)):
            with b.control(b["target"][k]):
                b.gate("phase", b["target"][j], math.pi / (1 << (j - k)))
    for j in range(width // 2):
        b.swap(b["target"][j], b["target"][width - j - 1])
    return b.finish()


def fourier_momentum(width, period):
    """频率对角 BE，按位的投影求和；避免构造稠密矩阵。"""
    terms = []
    for bit in range(width):
        b = Builder(_name("momentum_bit", width, bit), {"target": Bits(width), "signal": Bits(1)})
        b.x(b["signal"])
        with b.control(b["target"][bit]):
            b.x(b["signal"])
        be = BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))
        coefficient = (1 << bit) * 2 * math.pi / period * (-1 if bit == width - 1 else 1)
        terms.append((coefficient, be))
    return lcu(terms)


@dataclass(frozen=True)
class SchrodingerPlan:
    auxiliary_width: int = 2
    period: float = 8.0
    selected_index: int = 1

    def __post_init__(self):
        if (
            self.auxiliary_width < 1
            or self.period <= 0
            or not 0 <= self.selected_index < (1 << self.auxiliary_width)
        ):
            raise ValidationError("Schrodingerization 辅助网格/通道无效")


def schrodinger_qode(
    generator, initial, time, *, plan=None, hamiltonian_function=taylor_hamiltonian
):
    """u'=Gu，G=H1+iH2；Fourier lift P⊗H1-I⊗H2。"""
    plan = plan or SchrodingerPlan()
    parts = HermitianParts.from_operator(generator)
    n, p = generator.width, plan.auxiliary_width
    momentum = fourier_momentum(p, plan.period)
    hamiltonian = lcu([(1, tensor(momentum, parts.hermitian)), (-1, tensor(identity(p), parts.h))])
    evolution = hamiltonian_function(hamiltonian, time)
    grid = [
        (j if j < (1 << (p - 1)) else j - (1 << p)) * plan.period / (1 << p) for j in range(1 << p)
    ]
    warp = gate_state_prep([math.exp(-abs(x)) for x in grid])
    transform = qft(p)
    prep = Builder(
        _name("schrod_warp_initial", initial.operation, plan),
        {"target": Bits(n + p), "work": Bits(initial.work_width)},
        resources_for(("initial", initial.operation)),
    )
    invoke(prep, initial.operation, "initial", target=prep["target"][:n], work=prep["work"])
    invoke(prep, warp.operation, target=prep["target"][n:], work=prep["work"][:0])
    invoke(prep, transform, target=prep["target"][n:], work=prep["work"][:0])
    state = apply_be_to_state(
        evolution, StatePreparation(annotate(prep.finish(), "state_prep_isometry", zero_input=True))
    )
    out = Builder(
        _name("schrod_inverse_fourier", state.operation),
        {"target": Bits(n + p), "signal": Bits(state.signal_qubits)},
        resources_for(("state", state.operation)),
    )
    invoke(out, state.operation, "state", target=out["target"], signal=out["signal"])
    with out.adjoint():
        invoke(out, transform, target=out["target"][n:], work=out["signal"][:0])
    selected = select_subspace(
        StateOracle(out.finish()),
        n,
        plan.selected_index,
        label="schrodingerization_physical_channel",
    )
    return StateOracle(
        tagged(
            selected.operation,
            "schrodingerization_qode",
            auxiliary_grid=json.dumps(grid),
            selected_p=grid[plan.selected_index],
            recovery_scale=math.exp(grid[plan.selected_index]),
            recovery_assumption="selected p in valid warped region; periodic truncation pending",
        )
    )


@dataclass(frozen=True)
class PolynomialODE:
    width: int
    coefficients: tuple[tuple[int, BlockEncoding], ...]
    initial: StatePreparation
    initial_norm: float = 1.0

    def __post_init__(self):
        if self.initial.width != self.width or self.initial_norm < 0:
            raise ValidationError("Carleman 初始数据布局无效")
        for order, coefficient in self.coefficients:
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
    generator, initial = (
        carleman_lift(problem, cutoff=cutoff),
        carleman_initial(problem, cutoff=cutoff),
    )
    state = linear_solver(generator, initial, time)
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


def linear_qode(method, *, hamiltonian_function=taylor_hamiltonian, **options):
    """通用 u'=Gu 接口，可直接注入既有 QHAM / make_qpde。"""
    if method == "schrodingerization":
        return partial(schrodinger_qode, hamiltonian_function=hamiltonian_function, **options)
    if method not in {"lchs", "cbmd"}:
        raise ValidationError("未知线性 QODE 方法")
    algorithm = lchs_qode if method == "lchs" else cbmd_qode

    def generate(generator, initial, time):
        model = LinearODE(HermitianParts.from_operator(scale(-1, generator)), initial)
        return algorithm(model, time, hamiltonian_function=hamiltonian_function, **options)

    return generate


@dataclass(frozen=True)
class PDEInput:
    """空间离散化后的输入模型，允许直接为 PolynomialODE；不要求稠密矩阵。"""

    model: object
    label: str = "open_spatial_discretization"


def qpde_solver(qode, spatial_discretizer=lambda problem: problem.model):
    def generate(problem, time):
        return qode(spatial_discretizer(problem), time)

    return generate
