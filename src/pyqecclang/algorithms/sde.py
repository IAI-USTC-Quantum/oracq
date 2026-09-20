"""Fokker–Planck/SDE 的线性 ODE 输入模型与纯 Python 经典见证。

一维 Ito 随机微分方程 ``dx = a(x) dt + sqrt(2 D(x)) dW`` 的概率密度满足
Fokker–Planck 方程 ``p' = -∂x(a p) + ∂xx(D p)``。本模块在均匀网格上做零通量
有限体积离散，得到列和为零的离散生成元 G；G 经显式 Pauli 展开成为
BlockEncoding 后即可组装 QODEProblem/LinearODE，交给现有线性 ODE 求解器
（LCHS 等）。生成元的耗散性是输入模型的声明，不由语言证明。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from pyqecclang.algorithms.block_encoding import matrix_pauli_encoding
from pyqecclang.algorithms.contracts import finite_real, positive_integer
from pyqecclang.algorithms.ode import QODEProblem
from pyqecclang.algorithms.ode_models import HermitianParts, LinearODE
from pyqecclang.algorithms.operators import scale
from pyqecclang.algorithms.oracles import (
    StatePreparation,
    gate_state_prep,
    qram_state_angles,
    qram_state_prep,
    uniform_state,
)
from pyqecclang.infrastructure.ir import ValidationError


def _coefficients(values, size, path, *, nonnegative=False):
    """把常量或逐点向量统一成有限实数元组。"""
    if type(values) in (int, float):
        values = (values,) * size
    result = tuple(values)
    if len(result) != size:
        raise ValidationError(path + " 长度必须等于网格点数")
    for value in result:
        finite_real(value, path, minimum=0 if nonnegative else None)
    return tuple(float(v) for v in result)


def _uniform_points(grid):
    """校验严格递增且均匀的网格坐标，返回 (points, spacing)。"""
    points = tuple(grid)
    if len(points) < 2:
        raise ValidationError("FokkerPlanckProblem.grid 至少需要两个网格点")
    for value in points:
        finite_real(value, "FokkerPlanckProblem.grid")
    spacing = (points[-1] - points[0]) / (len(points) - 1)
    if spacing <= 0:
        raise ValidationError("FokkerPlanckProblem.grid 必须严格递增")
    tolerance = 1e-9 * max(1.0, abs(spacing))
    for left, right in zip(points, points[1:], strict=False):
        if abs((right - left) - spacing) > tolerance:
            raise ValidationError("当前仅支持均匀网格，非均匀网格超出支持范围")
    return tuple(float(p) for p in points), float(spacing)


def _power_of_two_width(size, path):
    if size < 2 or size & (size - 1):
        raise ValidationError(path + " 需要二的幂个网格点以匹配量子寄存器宽度")
    return (size - 1).bit_length()


@dataclass(frozen=True)
class FokkerPlanckProblem:
    """守恒形式 Fokker–Planck 算子的零通量有限体积离散输入模型。

    Args:
        drift: 漂移系数 a(x)，可为常量或逐点向量。
        diffusion: 扩散系数 D(x) >= 0，可为常量或逐点向量。
        grid: 均匀递增的网格坐标序列。
    """

    drift: object
    diffusion: object
    grid: object

    def __post_init__(self):
        points, spacing = _uniform_points(self.grid)
        object.__setattr__(self, "points", points)
        object.__setattr__(self, "spacing", spacing)
        size = len(points)
        object.__setattr__(
            self, "drift", _coefficients(self.drift, size, "FokkerPlanckProblem.drift")
        )
        object.__setattr__(
            self,
            "diffusion",
            _coefficients(self.diffusion, size, "FokkerPlanckProblem.diffusion", nonnegative=True),
        )

    @property
    def size(self):
        """网格点数，也是离散生成元矩阵的维数。"""
        return len(self.points)

    @property
    def width(self):
        """网格对应的量子寄存器位数，等于 ``log2(size)``。

        网格点数不是二的幂时抛出 ``ValidationError``。"""
        return _power_of_two_width(self.size, "FokkerPlanckProblem")

    def generator_matrix(self):
        """返回离散生成元 G（列向量约定 ``p' = G p``，列和为零）。"""
        size, h = self.size, self.spacing
        g = [[0.0] * size for _ in range(size)]
        for face in range(size - 1):
            a_face = 0.5 * (self.drift[face] + self.drift[face + 1])
            d_face = 0.5 * (self.diffusion[face] + self.diffusion[face + 1])
            forward = (0.5 * a_face + d_face / h) / h
            backward = (0.5 * a_face - d_face / h) / h
            g[face][face] -= forward
            g[face + 1][face] += forward
            g[face][face + 1] -= backward
            g[face + 1][face + 1] += backward
        return tuple(tuple(row) for row in g)

    def generator_encoding(self):
        """小尺度显式 Pauli 展开的 BlockEncoding；不宣称量子加速。"""
        if self.width > 5:
            raise ValidationError("显式 Pauli 展开仅用于不超过 32 点的小网格，大实例需要访问 oracle")
        return matrix_pauli_encoding(self.generator_matrix())

    def qode_problem(self, initial=None):
        """组装 ``p' = G p`` 的 QODEProblem；dissipative 为调用方声明。"""
        initial = initial if initial is not None else uniform_state(self.width)
        if not isinstance(initial, StatePreparation):
            raise ValidationError("FokkerPlanckProblem.qode_problem 需要 StatePreparation 初态")
        if initial.width != self.width:
            raise ValidationError("初态宽度与 Fokker–Planck 网格宽度不匹配")
        return QODEProblem(
            self.generator_encoding(),
            initial,
            dissipative=True,
            evidence="fokker_planck_zero_flux_finite_volume; dissipative caller-declared",
        )

    def linear_ode(self, initial=None):
        """同一问题的 ``u' = -A u`` （A = -G）LinearODE 视图。"""
        initial = initial if initial is not None else uniform_state(self.width)
        return LinearODE(
            HermitianParts.from_operator(scale(-1, self.generator_encoding())),
            initial,
            label="fokker_planck_minus_A",
        )


def _probabilities(probabilities, path):
    values = tuple(probabilities)
    if len(values) < 2 or len(values) & (len(values) - 1):
        raise ValidationError(path + " 需要二的幂长度且至少两个点")
    for value in values:
        finite_real(value, path, minimum=0)
    total = math.fsum(values)
    if total <= 0:
        raise ValidationError(path + " 总和必须为正")
    return tuple(float(v) / total for v in values)


def sde_state_preparation(probabilities, *, implementation="gate", angle_width=8, work_width=0):
    """把离散初态分布编码为幅度等于 ``sqrt(p_i)`` 的 StatePreparation。

    Args:
        probabilities: 非负离散概率，长度为二的幂。
        implementation: ``"gate"`` 为普通复用旋转实现，``"qram"`` 为 QRAM 角表实现。
        angle_width: QRAM 实现的角表字宽。
        work_width: gate 实现附加的工作位宽。
    """
    values = _probabilities(probabilities, "sde_state_preparation.probabilities")
    amplitudes = [math.sqrt(v) for v in values]
    if implementation == "gate":
        return gate_state_prep(amplitudes, work_width=work_width)
    if implementation == "qram":
        positive_integer(angle_width, "sde_state_preparation.angle_width")
        return qram_state_prep((len(values) - 1).bit_length(), angle_width)
    raise ValidationError("未知态制备实现：" + repr(implementation))


def sde_state_angles(probabilities, *, angle_width=8):
    """QRAM 实现的角表绑定数据；键为旋转树节点地址。"""
    values = _probabilities(probabilities, "sde_state_angles.probabilities")
    positive_integer(angle_width, "sde_state_angles.angle_width")
    return qram_state_angles([math.sqrt(v) for v in values], angle_width)


def _square_matrix(matrix, path):
    result = tuple(tuple(row) for row in matrix)
    size = len(result)
    if size < 1 or any(len(row) != size for row in result):
        raise ValidationError(path + " 需要非空方阵")
    for row in result:
        for value in row:
            finite_real(value, path)
    return tuple(tuple(float(v) for v in row) for row in result)


def _matvec(matrix, vector):
    return [math.fsum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def _matmul(a, b):
    size = len(a)
    columns = [[b[i][j] for i in range(size)] for j in range(size)]
    return [[math.fsum(x * y for x, y in zip(row, col, strict=True)) for col in columns] for row in a]


def matrix_exponential(matrix, time=1.0):
    """缩放平方加 Taylor 的纯 Python 矩阵指数，仅用于小规模经典见证。"""
    finite_real(time, "matrix_exponential.time")
    a = [list(row) for row in _square_matrix(matrix, "matrix_exponential.matrix")]
    size = len(a)
    for i in range(size):
        for j in range(size):
            a[i][j] *= time
    norm = max((math.fsum(abs(v) for v in row) for row in a), default=0.0)
    halvings = max(0, math.ceil(math.log2(norm / 0.5))) if norm > 0.5 else 0
    if halvings:
        shrink = 2.0**halvings
        a = [[v / shrink for v in row] for row in a]
    result = [[float(i == j) for j in range(size)] for i in range(size)]
    term = [row[:] for row in result]
    for order in range(1, 200):
        term = [[v / order for v in row] for row in _matmul(term, a)]
        result = [[x + y for x, y in zip(rx, tx, strict=True)] for rx, tx in zip(result, term, strict=True)]
        if max(abs(v) for row in term for v in row) < 1e-17:
            break
    for _ in range(halvings):
        result = _matmul(result, result)
    return result


def evolve_distribution(matrix, initial, time, *, steps=None):
    """经典参考演化；steps 为 None 时用矩阵指数，否则显式 Euler。"""
    g = _square_matrix(matrix, "evolve_distribution.matrix")
    initial = tuple(initial)
    if len(initial) != len(g):
        raise ValidationError("evolve_distribution.initial 长度与矩阵不符")
    for value in initial:
        finite_real(value, "evolve_distribution.initial")
    finite_real(time, "evolve_distribution.time", minimum=0)
    vector = [float(v) for v in initial]
    if steps is None:
        return _matvec(matrix_exponential(g, time), vector)
    positive_integer(steps, "evolve_distribution.steps")
    dt = time / steps
    for _ in range(steps):
        delta = _matvec(g, vector)
        vector = [v + dt * dv for v, dv in zip(vector, delta, strict=True)]
    return vector


def distribution_moments(points, probabilities, orders=(1, 2)):
    """由网格坐标与概率向量计算各阶矩 ``<x^k>``，默认返回 ``(<x>, <x^2>)``。"""
    points = tuple(points)
    values = tuple(probabilities)
    if len(points) != len(values) or not points:
        raise ValidationError("distribution_moments 需要等长非空的坐标与概率")
    for value in (*points, *values):
        finite_real(value, "distribution_moments")
    if any(type(k) is not int or k < 1 for k in orders):
        raise ValidationError("distribution_moments.orders 需要正整数阶数")
    total = math.fsum(values)
    if total <= 0:
        raise ValidationError("distribution_moments 概率总和必须为正")
    return tuple(
        math.fsum(p * x**k for x, p in zip(points, values, strict=True)) / total for k in orders
    )


def _interface_ratio(problem, face):
    """零通量面的有效网格 Péclet 数 u = a h / (2 D)。"""
    if not isinstance(problem, FokkerPlanckProblem):
        raise ValidationError("需要 FokkerPlanckProblem")
    a_face = 0.5 * (problem.drift[face] + problem.drift[face + 1])
    d_face = 0.5 * (problem.diffusion[face] + problem.diffusion[face + 1])
    if d_face <= 0:
        raise ValidationError("稳态参考要求面上扩散系数严格为正")
    return a_face * problem.spacing / (2 * d_face)


def stationary_distribution(problem):
    """零通量离散的精确稳态；相邻点比值为 ``(1+u)/(1-u)``。"""
    weights = [1.0]
    for face in range(problem.size - 1):
        u = _interface_ratio(problem, face)
        if abs(u) >= 1:
            raise ValidationError("网格 Péclet 数过大，中心差分稳态不再为正")
        weights.append(weights[-1] * (1 + u) / (1 - u))
    total = math.fsum(weights)
    return tuple(w / total for w in weights)


def boltzmann_distribution(problem):
    """连续稳态参考 ``p ∝ exp(∫ a/D dx)``，与离散稳态相差 O(h^2)。"""
    weights = [1.0]
    for face in range(problem.size - 1):
        u = _interface_ratio(problem, face)
        weights.append(weights[-1] * math.exp(2 * u))
    total = math.fsum(weights)
    return tuple(w / total for w in weights)
