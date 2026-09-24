"QHAM 的 PDE 表达式、推导和空间离散化支持。"

from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # 仅静态检查可见的真实类型；运行时仍走下方 _EXPORTS 惰性导入。
    # 冗余别名（X as X）标记为显式再导出，供 ruff F401 识别。
    from oracq.algorithms.input_model.qham import (
        PortBinding as PortBinding,
    )
    from oracq.algorithms.input_model.qham import (
        QHAMBindings as QHAMBindings,
    )
    from oracq.algorithms.input_model.qham import (
        QHAMInputModel as QHAMInputModel,
    )
    from oracq.algorithms.input_model.qham import (
        gate_bindings as gate_bindings,
    )
    from oracq.algorithms.input_model.qham import (
        open_qham_input as open_qham_input,
    )
    from oracq.algorithms.input_model.qham import (
        qham_input_model as qham_input_model,
    )
    from oracq.algorithms.input_model.qham import (
        taylor_qode as taylor_qode,
    )
    from oracq.applications.qham.linearization import (
        Block as Block,
    )
    from oracq.applications.qham.linearization import (
        Coupling as Coupling,
    )
    from oracq.applications.qham.linearization import (
        HomotopyWeight as HomotopyWeight,
    )
    from oracq.applications.qham.linearization import (
        QHAMPlan as QHAMPlan,
    )
    from oracq.applications.qham.pde import (
        Expr as Expr,
    )
    from oracq.applications.qham.pde import (
        Field as Field,
    )
    from oracq.applications.qham.pde import (
        Known as Known,
    )
    from oracq.applications.qham.pde import (
        PolynomialPDE as PolynomialPDE,
    )
    from oracq.applications.qham.reference import (
        Discretization as Discretization,
    )
    from oracq.applications.qham.reference import (
        Grid as Grid,
    )
    from oracq.applications.qham.stencils import (
        qram_coefficient_encoding as qram_coefficient_encoding,
    )
    from oracq.applications.qham.stencils import (
        qram_coefficient_memory as qram_coefficient_memory,
    )
    from oracq.applications.qham.stencils import (
        structured_fd_bindings as structured_fd_bindings,
    )

_EXPORTS = {
    "Block": ("oracq.applications.qham.linearization", "Block"),
    "Coupling": ("oracq.applications.qham.linearization", "Coupling"),
    "HomotopyWeight": ("oracq.applications.qham.linearization", "HomotopyWeight"),
    "QHAMPlan": ("oracq.applications.qham.linearization", "QHAMPlan"),
    "Expr": ("oracq.applications.qham.pde", "Expr"),
    "Field": ("oracq.applications.qham.pde", "Field"),
    "Known": ("oracq.applications.qham.pde", "Known"),
    "PolynomialPDE": ("oracq.applications.qham.pde", "PolynomialPDE"),
    "PortBinding": ("oracq.algorithms.input_model.qham", "PortBinding"),
    "QHAMBindings": ("oracq.algorithms.input_model.qham", "QHAMBindings"),
    "QHAMInputModel": ("oracq.algorithms.input_model.qham", "QHAMInputModel"),
    "gate_bindings": ("oracq.algorithms.input_model.qham", "gate_bindings"),
    "open_qham_input": ("oracq.algorithms.input_model.qham", "open_qham_input"),
    "qham_input_model": ("oracq.algorithms.input_model.qham", "qham_input_model"),
    "taylor_qode": ("oracq.algorithms.input_model.qham", "taylor_qode"),
    "Discretization": ("oracq.applications.qham.reference", "Discretization"),
    "Grid": ("oracq.applications.qham.reference", "Grid"),
    "structured_fd_bindings": ("oracq.applications.qham.stencils", "structured_fd_bindings"),
    "qram_coefficient_encoding": (
        "oracq.applications.qham.stencils",
        "qram_coefficient_encoding",
    ),
    "qram_coefficient_memory": (
        "oracq.applications.qham.stencils",
        "qram_coefficient_memory",
    ),
}
__all__ = list(_EXPORTS)


def __getattr__(name: str) -> object:
    """惰性导入 ``_EXPORTS`` 表中的成员并缓存到模块全局。"""
    if name not in _EXPORTS:
        raise AttributeError(name)
    module, member = _EXPORTS[name]
    value = getattr(import_module(module), member)
    globals()[name] = value
    return value
