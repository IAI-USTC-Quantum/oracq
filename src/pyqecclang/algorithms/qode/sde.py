"""Fokker–Planck/SDE 的线性 ODE 输入模型与纯 Python 经典见证。

一维 Ito 随机微分方程 ``dx = a(x) dt + sqrt(2 D(x)) dW`` 的概率密度满足
Fokker–Planck 方程 ``p' = -∂x(a p) + ∂xx(D p)``。本模块在均匀网格上做零通量
有限体积离散，得到列和为零的离散生成元 G；G 经显式 Pauli 展开成为
BlockEncoding 后即可组装 QODEProblem/LinearODE，交给现有线性 ODE 求解器
（LCHS 等）。生成元的耗散性是输入模型的声明，不由语言证明。
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import cast

from pyqecclang.algorithms.input_model.block_encoding import matrix_pauli_encoding
from pyqecclang.algorithms.input_model.contracts import finite_real, positive_integer
from pyqecclang.algorithms.input_model.operators import BlockEncoding, scale
from pyqecclang.algorithms.input_model.oracles import (
    StatePreparation,
    gate_state_prep,
    qram_state_angles,
    qram_state_prep,
    uniform_state,
)
from pyqecclang.algorithms.qode.ode import QODEProblem
from pyqecclang.algorithms.qode.ode_models import HermitianParts, LinearODE
from pyqecclang.infrastructure.ir import ValidationError


def _coefficients(
    values: int | float | Iterable[float],
    size: int,
    path: str,
    *,
    nonnegative: bool = False,
) -> tuple[float, ...]:
    """把常量或逐点向量统一成有限实数元组。"""
    if type(values) in (int, float):
        values = cast("Iterable[float]", (values,) * size)
    result = tuple(cast("Iterable[float]", values))
    if len(result) != size:
        raise ValidationError(path + " 长度必须等于网格点数")
    for value in result:
        finite_real(value, path, minimum=0 if nonnegative else None)
    return tuple(float(v) for v in result)


def _uniform_points(grid: Iterable[float]) -> tuple[tuple[float, ...], float]:
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


def _power_of_two_width(size: int, path: str) -> int:
    """校验网格点数为不小于 2 的二的幂，并返回对应寄存器位数。"""
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

    drift: int | float | Iterable[float]
    diffusion: int | float | Iterable[float]
    grid: Iterable[float]

    def __post_init__(self) -> None:
        """校验均匀网格并把漂移、扩散系数规范化为逐点元组。"""
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
    def size(self) -> int:
        """网格点数，也是离散生成元矩阵的维数。"""
        return len(self.points)  # type: ignore[attr-defined]

    @property
    def width(self) -> int:
        """网格对应的量子寄存器位数，等于 ``log2(size)``。

        网格点数不是二的幂时抛出 ``ValidationError``。"""
        return _power_of_two_width(self.size, "FokkerPlanckProblem")

    def generator_matrix(self) -> tuple[tuple[float, ...], ...]:
        """返回离散生成元 G（列向量约定 ``p' = G p``，列和为零）。

        Returns:
            tuple[tuple[float, ...], ...]: 网格维数的方阵，按行嵌套的元组表示。
        """
        size, h = self.size, self.spacing  # type: ignore[attr-defined]
        g = [[0.0] * size for _ in range(size)]
        for face in range(size - 1):
            a_face = 0.5 * (
                cast("tuple[float, ...]", self.drift)[face]
                + cast("tuple[float, ...]", self.drift)[face + 1]
            )
            d_face = 0.5 * (
                cast("tuple[float, ...]", self.diffusion)[face]
                + cast("tuple[float, ...]", self.diffusion)[face + 1]
            )
            forward = (0.5 * a_face + d_face / h) / h
            backward = (0.5 * a_face - d_face / h) / h
            g[face][face] -= forward
            g[face + 1][face] += forward
            g[face][face + 1] -= backward
            g[face + 1][face + 1] += backward
        return tuple(tuple(row) for row in g)

    def generator_encoding(self) -> BlockEncoding:
        """小尺度显式 Pauli 展开的 BlockEncoding；不宣称量子加速。

        Returns:
            BlockEncoding: 离散生成元 G 的显式 Pauli 展开块编码。
        """
        if self.width > 5:
            raise ValidationError("显式 Pauli 展开仅用于不超过 32 点的小网格，大实例需要访问 oracle")
        return matrix_pauli_encoding(self.generator_matrix())

    def qode_problem(self, initial: StatePreparation | None = None) -> QODEProblem:
        """组装 ``p' = G p`` 的 QODEProblem；dissipative 为调用方声明。

        Args:
            initial: 初态制备句柄；缺省为网格宽度上的均匀分布态，宽度须与网格匹配。

        Returns:
            QODEProblem: 生成元 G 驱动的量子 ODE 问题，带零通量有限体积证据标注。
        """
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

    def linear_ode(self, initial: StatePreparation | None = None) -> LinearODE:
        """同一问题的 ``u' = -A u`` （A = -G）LinearODE 视图。

        Args:
            initial: 初态制备句柄；缺省为网格宽度上的均匀分布态。

        Returns:
            LinearODE: A 取 -G 的线性 ODE 视图，标签为 ``fokker_planck_minus_A``。
        """
        initial = initial if initial is not None else uniform_state(self.width)
        return LinearODE(
            HermitianParts.from_operator(scale(-1, self.generator_encoding())),
            initial,
            label="fokker_planck_minus_A",
        )


def _probabilities(probabilities: Iterable[float], path: str) -> tuple[float, ...]:
    """校验二的幂长度的非负概率序列并归一化。"""
    values = tuple(probabilities)
    if len(values) < 2 or len(values) & (len(values) - 1):
        raise ValidationError(path + " 需要二的幂长度且至少两个点")
    for value in values:
        finite_real(value, path, minimum=0)
    total = math.fsum(values)
    if total <= 0:
        raise ValidationError(path + " 总和必须为正")
    return tuple(float(v) / total for v in values)


def sde_state_preparation(
    probabilities: Iterable[float],
    *,
    implementation: str = "gate",
    angle_width: int = 8,
    work_width: int = 0,
) -> StatePreparation:
    """把离散初态分布编码为幅度等于 ``sqrt(p_i)`` 的 StatePreparation。

    Args:
        probabilities: 非负离散概率，长度为二的幂。
        implementation: ``"gate"`` 为普通复用旋转实现，``"qram"`` 为 QRAM 角表实现。
        angle_width: QRAM 实现的角表字宽。
        work_width: gate 实现附加的工作位宽。

    Returns:
        StatePreparation: 幅度为 ``sqrt(p_i)`` 的离散分布制备句柄。
    """
    values = _probabilities(probabilities, "sde_state_preparation.probabilities")
    amplitudes = [math.sqrt(v) for v in values]
    if implementation == "gate":
        return gate_state_prep(amplitudes, work_width=work_width)
    if implementation == "qram":
        positive_integer(angle_width, "sde_state_preparation.angle_width")
        return qram_state_prep((len(values) - 1).bit_length(), angle_width)
    raise ValidationError("未知态制备实现：" + repr(implementation))


def sde_state_angles(probabilities: Iterable[float], *, angle_width: int = 8) -> dict[int, int]:
    """QRAM 实现的角表绑定数据；键为旋转树节点地址。

    Args:
        probabilities: 非负离散概率，长度为二的幂，内部归一化。
        angle_width: 角表字的量化位宽，取正整数。

    Returns:
        dict[int, int]: 键为旋转树节点地址、值为量化后角度整数的绑定数据。
    """
    values = _probabilities(probabilities, "sde_state_angles.probabilities")
    positive_integer(angle_width, "sde_state_angles.angle_width")
    return qram_state_angles([math.sqrt(v) for v in values], angle_width)


def _square_matrix(matrix: Iterable[Iterable[float]], path: str) -> tuple[tuple[float, ...], ...]:
    """校验非空方阵且元素有限，并转成浮点元组表示。"""
    result = tuple(tuple(row) for row in matrix)
    size = len(result)
    if size < 1 or any(len(row) != size for row in result):
        raise ValidationError(path + " 需要非空方阵")
    for row in result:
        for value in row:
            finite_real(value, path)
    return tuple(tuple(float(v) for v in row) for row in result)


def _matvec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    """计算矩阵与向量的乘积，逐行用 ``math.fsum`` 累加。"""
    return [math.fsum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def _matmul(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> list[list[float]]:
    """计算两个同阶方阵的乘积，逐元素用 ``math.fsum`` 累加。"""
    size = len(a)
    columns = [[b[i][j] for i in range(size)] for j in range(size)]
    return [[math.fsum(x * y for x, y in zip(row, col, strict=True)) for col in columns] for row in a]


def matrix_exponential(matrix: Iterable[Iterable[float]], time: float = 1.0) -> list[list[float]]:
    """缩放平方加 Taylor 的纯 Python 矩阵指数，仅用于小规模经典见证。

    Args:
        matrix: 非空方阵，元素为有限实数。
        time: 作用时长缩放因子，取有限实数。

    Returns:
        list[list[float]]: exp(matrix·time) 的矩阵，按行嵌套的列表表示。
    """
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


def evolve_distribution(
    matrix: Iterable[Iterable[float]],
    initial: Iterable[float],
    time: float,
    *,
    steps: int | None = None,
) -> list[float]:
    """经典参考演化；steps 为 None 时用矩阵指数，否则显式 Euler。

    Args:
        matrix: 生成元方阵，约定 ``p' = matrix · p``。
        initial: 与矩阵同维的初值向量。
        time: 演化时长，取非负有限实数。
        steps: 显式 Euler 步数；缺省时一步调用矩阵指数求解。

    Returns:
        list[float]: 演化到给定时刻的分布向量。
    """
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


def distribution_moments(
    points: Iterable[float],
    probabilities: Iterable[float],
    orders: tuple[int, ...] = (1, 2),
) -> tuple[float, ...]:
    """由网格坐标与概率向量计算各阶矩 ``<x^k>``，默认返回 ``(<x>, <x^2>)``。

    Args:
        points: 网格坐标序列，与概率向量等长。
        probabilities: 各点概率权重，总和须为正，内部归一化。
        orders: 需要计算的矩阶数，取正整数元组。

    Returns:
        tuple[float, ...]: 与阶数一一对应的归一化矩 ``<x^k>``。
    """
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


def _interface_ratio(problem: FokkerPlanckProblem, face: int) -> float:
    """零通量面的有效网格 Péclet 数 u = a h / (2 D)。"""
    if not isinstance(problem, FokkerPlanckProblem):
        raise ValidationError("需要 FokkerPlanckProblem")
    a_face = 0.5 * (
        cast("tuple[float, ...]", problem.drift)[face]
        + cast("tuple[float, ...]", problem.drift)[face + 1]
    )
    d_face = 0.5 * (
        cast("tuple[float, ...]", problem.diffusion)[face]
        + cast("tuple[float, ...]", problem.diffusion)[face + 1]
    )
    if d_face <= 0:
        raise ValidationError("稳态参考要求面上扩散系数严格为正")
    return a_face * problem.spacing / (2 * d_face)  # type: ignore[attr-defined]


def stationary_distribution(problem: FokkerPlanckProblem) -> tuple[float, ...]:
    """零通量离散的精确稳态；相邻点比值为 ``(1+u)/(1-u)``。

    Args:
        problem: 已规范化的 Fokker–Planck 离散问题，各面扩散系数须为正。

    Returns:
        tuple[float, ...]: 与网格点一一对应、总和为一的离散稳态概率。
    """
    weights = [1.0]
    for face in range(problem.size - 1):
        u = _interface_ratio(problem, face)
        if abs(u) >= 1:
            raise ValidationError("网格 Péclet 数过大，中心差分稳态不再为正")
        weights.append(weights[-1] * (1 + u) / (1 - u))
    total = math.fsum(weights)
    return tuple(w / total for w in weights)


def boltzmann_distribution(problem: FokkerPlanckProblem) -> tuple[float, ...]:
    """连续稳态参考 ``p ∝ exp(∫ a/D dx)``，与离散稳态相差 O(h^2)。

    Args:
        problem: 已规范化的 Fokker–Planck 离散问题，各面扩散系数须为正。

    Returns:
        tuple[float, ...]: 与网格点一一对应、总和为一的连续稳态参考概率。
    """
    weights = [1.0]
    for face in range(problem.size - 1):
        u = _interface_ratio(problem, face)
        weights.append(weights[-1] * math.exp(2 * u))
    total = math.fsum(weights)
    return tuple(w / total for w in weights)
