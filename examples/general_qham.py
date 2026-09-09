"""从普通 PDE 表达式到开放 QODE 输入的最小示例。"""

from functools import partial

from pyqecclang.algorithms.differential import linear_qode, taylor_hamiltonian
from pyqecclang.qham import (
    Discretization,
    Field,
    Grid,
    Known,
    PolynomialPDE,
    QHAMPlan,
    qham_input_model,
    structured_fd_bindings,
)

u = Field("u")
pde = PolynomialPDE.from_equations(
    {"u": 0.1 * u.d("x", 2) - u * u.d("x") + Known("forcing")},
    label="forced_burgers",
)
plan = QHAMPlan(pde, order=2)
grid = Grid(("x",), (4,), (1.0,), boundary="periodic")
space = Discretization(pde, grid, {"forcing": [0.05, 0, -0.05, 0]})
initial = space.encode_fields({"u": [0.1, 0.2, 0, -0.1]})
bindings = structured_fd_bindings(space, initial)
qode_input = qham_input_model(plan, bindings, eta=-0.4)

# Schrödingerization 接收一般生成元；恢复区间和时间近似仍需后续验证。
solver = linear_qode(
    "schrodingerization", hamiltonian_function=partial(taylor_hamiltonian, degree=1)
)
solution = qode_input.solve(solver, 0.01)

# 如改用 CBMD/LCHS，先显式施加使生成元耗散的整体移位。
shifted = qode_input.dissipative_shift()
cbmd_solution = shifted.solve(
    linear_qode("cbmd", hamiltonian_function=partial(taylor_hamiltonian, degree=1)), 0.01
)
