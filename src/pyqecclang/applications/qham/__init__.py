"QHAM 的 PDE 表达式、推导和空间离散化支持。"

from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # 仅静态检查可见的真实类型；运行时仍走下方 _EXPORTS 惰性导入。
    # 冗余别名（X as X）标记为显式再导出，供 ruff F401 识别。
    from pyqecclang.algorithms.input_model.qham import (
        PortBinding as PortBinding,
    )
    from pyqecclang.algorithms.input_model.qham import (
        QHAMBindings as QHAMBindings,
    )
    from pyqecclang.algorithms.input_model.qham import (
        QHAMInputModel as QHAMInputModel,
    )
    from pyqecclang.algorithms.input_model.qham import (
        gate_bindings as gate_bindings,
    )
    from pyqecclang.algorithms.input_model.qham import (
        open_qham_input as open_qham_input,
    )
    from pyqecclang.algorithms.input_model.qham import (
        qham_input_model as qham_input_model,
    )
    from pyqecclang.algorithms.input_model.qham import (
        taylor_qode as taylor_qode,
    )
    from pyqecclang.applications.qham.linearization import (
        Block as Block,
    )
    from pyqecclang.applications.qham.linearization import (
        Coupling as Coupling,
    )
    from pyqecclang.applications.qham.linearization import (
        HomotopyWeight as HomotopyWeight,
    )
    from pyqecclang.applications.qham.linearization import (
        QHAMPlan as QHAMPlan,
    )
    from pyqecclang.applications.qham.pde import (
        Expr as Expr,
    )
    from pyqecclang.applications.qham.pde import (
        Field as Field,
    )
    from pyqecclang.applications.qham.pde import (
        Known as Known,
    )
    from pyqecclang.applications.qham.pde import (
        PolynomialPDE as PolynomialPDE,
    )
    from pyqecclang.applications.qham.reference import (
        Discretization as Discretization,
    )
    from pyqecclang.applications.qham.reference import (
        Grid as Grid,
    )
    from pyqecclang.applications.qham.stencils import (
        qram_coefficient_encoding as qram_coefficient_encoding,
    )
    from pyqecclang.applications.qham.stencils import (
        qram_coefficient_memory as qram_coefficient_memory,
    )
    from pyqecclang.applications.qham.stencils import (
        structured_fd_bindings as structured_fd_bindings,
    )

_EXPORTS = {
    "Block": ("pyqecclang.applications.qham.linearization", "Block"),
    "Coupling": ("pyqecclang.applications.qham.linearization", "Coupling"),
    "HomotopyWeight": ("pyqecclang.applications.qham.linearization", "HomotopyWeight"),
    "QHAMPlan": ("pyqecclang.applications.qham.linearization", "QHAMPlan"),
    "Expr": ("pyqecclang.applications.qham.pde", "Expr"),
    "Field": ("pyqecclang.applications.qham.pde", "Field"),
    "Known": ("pyqecclang.applications.qham.pde", "Known"),
    "PolynomialPDE": ("pyqecclang.applications.qham.pde", "PolynomialPDE"),
    "PortBinding": ("pyqecclang.algorithms.input_model.qham", "PortBinding"),
    "QHAMBindings": ("pyqecclang.algorithms.input_model.qham", "QHAMBindings"),
    "QHAMInputModel": ("pyqecclang.algorithms.input_model.qham", "QHAMInputModel"),
    "gate_bindings": ("pyqecclang.algorithms.input_model.qham", "gate_bindings"),
    "open_qham_input": ("pyqecclang.algorithms.input_model.qham", "open_qham_input"),
    "qham_input_model": ("pyqecclang.algorithms.input_model.qham", "qham_input_model"),
    "taylor_qode": ("pyqecclang.algorithms.input_model.qham", "taylor_qode"),
    "Discretization": ("pyqecclang.applications.qham.reference", "Discretization"),
    "Grid": ("pyqecclang.applications.qham.reference", "Grid"),
    "structured_fd_bindings": ("pyqecclang.applications.qham.stencils", "structured_fd_bindings"),
    "qram_coefficient_encoding": (
        "pyqecclang.applications.qham.stencils",
        "qram_coefficient_encoding",
    ),
    "qram_coefficient_memory": (
        "pyqecclang.applications.qham.stencils",
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
