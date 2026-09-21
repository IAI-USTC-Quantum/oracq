"QHAM 的 PDE 表达式、推导和空间离散化支持。"

from importlib import import_module

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


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(name)
    module, member = _EXPORTS[name]
    value = getattr(import_module(module), member)
    globals()[name] = value
    return value
