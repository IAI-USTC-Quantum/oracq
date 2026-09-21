"""空间离散化结果与可替换 QODE 生成器的组装接口。"""

from __future__ import annotations

from dataclasses import dataclass

from pyqecclang.algorithms.input_model.operators import BlockEncoding
from pyqecclang.algorithms.input_model.oracles import (
    StatePreparation,
)


@dataclass(frozen=True)
class DiscretePDE:
    """空间离散化后的线性 PDE 输入对象。

    只保存离散生成元、初态制备和标签；网格、边界、量纲等应用元数据不由该类型携带。

    Attributes:
        generator: 离散生成元 G 的块编码，对应 ``u' = Gu``。
        initial: 初态制备。
        label: 问题标签。
    """

    generator: BlockEncoding
    initial: StatePreparation
    label: str = "linear_pde"


def make_qpde(qode, discretizer=lambda problem: problem):
    """把三参数线性 QODE 协议包装成 QPDE 生成函数。

    Args:
        qode: ``(generator, initial, time) -> StateOracle`` 形式的线性求解协议。
        discretizer: 把问题对象映射为具有 ``generator`` 与 ``initial`` 属性的对象（如 ``DiscretePDE``）的可调用，默认原样返回。

    Returns:
        callable: 形如 ``(problem, final_time) -> StateOracle`` 的生成函数。
    """
    def generate(problem, final_time):
        discrete = discretizer(problem)
        result = qode(discrete.generator, discrete.initial, final_time)
        return result

    return generate


@dataclass(frozen=True)
class PDEInput:
    """空间离散化后的输入模型，允许直接为 PolynomialODE；不要求稠密矩阵。"""

    model: object
    label: str = "open_spatial_discretization"


def qpde_solver(qode, spatial_discretizer=lambda problem: problem.model):
    """把 ``(model, time)`` 形式的求解协议包装成 QPDE 生成函数。

    与 ``make_qpde`` 不同，这里把离散化结果作为单个模型直接交给 ``qode``，适合 ``PolynomialODE`` 等模型级协议。

    Args:
        qode: ``(model, time) -> StateOracle`` 形式的求解协议。
        spatial_discretizer: 从问题对象提取模型的可调用，默认取 ``problem.model``。

    Returns:
        callable: 形如 ``(problem, time) -> StateOracle`` 的生成函数。
    """
    def generate(problem, time):
        return qode(spatial_discretizer(problem), time)

    return generate
