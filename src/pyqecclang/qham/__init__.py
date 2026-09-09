"""规则化 PDE → HAM → 量子适配线性化 → QODE 的生成接口。"""

from .linearization import Block, Coupling, HomotopyWeight, QHAMPlan
from .pde import Expr, Field, Known, PolynomialPDE
from .quantum import (
    PortBinding,
    QHAMBindings,
    QHAMInputModel,
    gate_bindings,
    open_qham_input,
    qham_input_model,
    taylor_qode,
)
from .reference import Discretization, Grid
from .stencils import structured_fd_bindings

__all__ = [
    "Block",
    "Coupling",
    "HomotopyWeight",
    "QHAMPlan",
    "Expr",
    "Field",
    "Known",
    "PolynomialPDE",
    "PortBinding",
    "QHAMBindings",
    "QHAMInputModel",
    "gate_bindings",
    "open_qham_input",
    "qham_input_model",
    "taylor_qode",
    "Discretization",
    "Grid",
]

__all__ += ["structured_fd_bindings"]
