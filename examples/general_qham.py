"""A minimal example going from an ordinary PDE expression to an open QODE input."""

from __future__ import annotations

from functools import partial

from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.qode.ode import linear_qode
from oracq.applications.qham import (
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

# Schrödingerization accepts a general generator; the recovery interval and time approximation still need later validation.
solver = linear_qode(
    "schrodingerization", hamiltonian_function=partial(taylor_hamiltonian, degree=1)
)
solution = qode_input.solve(solver, 0.01)

# To switch to CBMD/LCHS instead, first apply an explicit global shift that makes the generator dissipative.
shifted = qode_input.dissipative_shift()
cbmd_solution = shifted.solve(
    linear_qode("cbmd", hamiltonian_function=partial(taylor_hamiltonian, degree=1)), 0.01
)
