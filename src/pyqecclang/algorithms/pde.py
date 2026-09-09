"""空间离散化结果与可替换 QODE 生成器的组装接口。"""

from __future__ import annotations

from dataclasses import dataclass

from pyqecclang.algorithms.operators import BlockEncoding
from pyqecclang.algorithms.oracles import (
    StatePreparation,
)


@dataclass(frozen=True)
class DiscretePDE:
    generator: BlockEncoding
    initial: StatePreparation
    label: str = "linear_pde"


def make_qpde(qode, discretizer=lambda problem: problem):
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
    def generate(problem, time):
        return qode(spatial_discretizer(problem), time)

    return generate
