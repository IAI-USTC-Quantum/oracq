# Generating QHAM input from a PDE expression

**English** · <a href="../zh/tutorials/qham.html">简体中文</a>

This example first defines the Burgers equation, then sets up a second-order HAM plan. Creating the plan expands no tensor block and builds no complete matrix. The unknown field is a {obj}`Field <oracq.applications.qham.pde.Field>`, the equation is frozen into a {obj}`PolynomialPDE <oracq.applications.qham.pde.PolynomialPDE>`, and the second-order plan is stored in a {obj}`QHAMPlan <oracq.applications.qham.linearization.QHAMPlan>`.

```{testcode}
from oracq.applications.qham import Field, PolynomialPDE, QHAMPlan

u = Field("u")
pde = PolynomialPDE.from_equations({"u": 0.1*u.d("x", 2) - u*u.d("x")})
plan = QHAMPlan(pde, order=2)
restored = QHAMPlan.loads(plan.dumps())
print([(port.name, port.arity) for port in pde.ports])
print(plan.block_count, plan.max_rank, plan.order)
assert restored == plan
assert any(port.arity == 2 for port in pde.ports)
```

```{testoutput}
[('L', 1), ('B_1', 2)]
8 3 2
```

The printed results show: the PDE decomposes into one linear port `L` and one second-order multilinear port `B_1`; the plan records only `8` tensor blocks, a maximum rank of `3`, and a truncation order of `2` — no block is expanded and no complete matrix is built. The nonlinear terms of the PDE are represented as multilinear ports. Once a grid and boundary are chosen, {obj}`structured_fd_bindings <oracq.applications.qham.stencils.structured_fd_bindings>` can produce BEs for these ports from shifts, coefficient multipliers, and same-point contractions; alternatively, declare the ports open first and wait for a data-access implementation.

The quantum assembly entry point {obj}`qham_input_model <oracq.algorithms.input_model.qham.qham_input_model>` produces the lifted operator and the initial state required by QODE. From there the linear solution method can be replaced, and a physical output channel representing the sum of the HAM orders can be chosen.

Run the complete example:

```bash
uv run python examples/general_qham.py
uv run python -m oracq.applications.qham --example burgers --order 2 --eta=-0.4
```

To check the derivation formulas, read the [QHAM mathematical derivation](../reference/qham-derivation.md). To bind concrete oracles, read the [QHAM implementation notes](../manual/qham.md). The algebraic representation of the truncated HAM and the convergence towards the original PDE must be validated separately.

The complete script and its step-by-step walkthrough are in [Scientific computing workflows](scientific-workflows.md), which also
contrasts the inputs, solver choices, and binding boundaries of Carleman, LCHS, and CBMD.

## Related pages

- Manual: [General QHAM automatic generation](../manual/qham.md)
- Specification: [QHAM mathematical derivation](../reference/qham-derivation.md)
- Algorithm pages: <a href="../zh/manual/algorithms/qham.html">QHAM</a>, <a href="../zh/manual/algorithms/qode-problem.html">QODE problem objects and protocol</a>
- API reference: [QHAM](../api/algorithms/input_model/qham.rst), [PDE models and adaptation](../api/applications/qham/pde.rst), [QHAM finite closure](../api/applications/qham/linearization.rst), [structured stencil ports](../api/applications/qham/stencils.rst)
