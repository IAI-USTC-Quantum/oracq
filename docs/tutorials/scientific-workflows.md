# Scientific computing workflows: input models, solvers, and reproduction

**English** · <a href="../../zh/tutorials/scientific-workflows.html">简体中文</a>

This chapter strings together the complete examples of QHAM, Carleman, LCHS, and CBMD, explaining which objects each step
produces, which interfaces remain open, and which conclusions the results can support. The paper keeps the mechanism and
representative experiments; full usage, input variants, and step-by-step walkthroughs continue to be maintained here.

## Start by distinguishing three layers

| Layer | What it decides | Typical objects |
|---|---|---|
| Mathematical input and approximation | PDE, discretization, truncation order, initial-state norm, physical output window | {obj}`PolynomialPDE <oracq.applications.qham.pde.PolynomialPDE>`, {obj}`QHAMPlan <oracq.applications.qham.linearization.QHAMPlan>`, {obj}`PolynomialODE <oracq.algorithms.qnlss.carleman.PolynomialODE>` |
| Quantum-algorithm generation | which solver and internal simulation kernel to choose, which module calls to generate | {obj}`QODESolver <oracq.algorithms.qode.ode.QODESolver>`, {obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>`, {obj}`carleman_qode <oracq.algorithms.qnlss.carleman.carleman_qode>` |
| Implementation binding and execution | implement the open slots with gate tables, QRAM, or arithmetic, and supply run-time data | the RIR {obj}`Program <oracq.infrastructure.ir.Program>`, {obj}`Binding <oracq.infrastructure.linking.Binding>`, memory snapshots |

Replacing the solver usually requires regenerating the program; swapping a compatible implementation into an existing
slot can use {obj}`bind <oracq.infrastructure.linking.bind>`.
Changing the public bit width or an alpha already used for assembly cannot be done by merely replacing a label.
Identical mathematical conditions and interface compatibility are also two different things: a generator being able
to provide a BE does not mean it satisfies the conditions of every QODE method.

## QHAM: from the equation to the lifted linear input

Read the [minimal QHAM tutorial](qham.md) first. What follows includes the actual example source directly, so imports,
equations, data, and solver configuration always match the repository scripts.

```{literalinclude} ../../examples/general_qham.py
:language: python
:start-at: from functools import partial
:caption: QHAM assembly of the forced Burgers equation and two solver choices.
```

Read this code through the lifecycle of its objects:

1. {obj}`Field <oracq.applications.qham.pde.Field>` and {obj}`Known <oracq.applications.qham.pde.Known>` build the unknown fields and the forcing terms referenced by name. The expression describes the field to be solved for;
   it does not apply an ordinary numeric function directly to quantum amplitudes.
2. {obj}`QHAMPlan(pde, order=2) <oracq.applications.qham.linearization.QHAMPlan>` fixes the finite HAM truncation and the quantum-adapted linearization rules.
   The plan can be saved, restored, and queried row by row; the full lifted matrix is never materialized up front.
3. {obj}`Grid <oracq.applications.qham.reference.Grid>` and {obj}`Discretization <oracq.applications.qham.reference.Discretization>` determine the periodic grid, the derivatives, and the known data. Boundary conditions belong to
   this layer; they are not guessed by the solver.
4. {obj}`structured_fd_bindings <oracq.applications.qham.stencils.structured_fd_bindings>` builds the port BEs from shifts, local coefficients, and multilinear contractions;
   {obj}`qham_input_model <oracq.algorithms.input_model.qham.qham_input_model>` then assembles the lifted generator and the initial state. The relative norms of the lifted
   initial state's tensor blocks must be preserved; the blocks must not be normalized independently.
5. {obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>` selects the linear-evolution method; {obj}`partial(taylor_hamiltonian, degree=1) <oracq.algorithms.common.hamiltonian.taylor_hamiltonian>`
   explicitly fixes this example's approximate simulation kernel. `solve` produces a state oracle carrying signal and
   recovery information; it does not directly return a certified-converged classical PDE solution.
6. `dissipative_shift` explicitly modifies the generator to meet another method's preconditions. The corresponding
   growth-recovery factor must not be dropped, and this step must not be treated as a free way to raise the success rate.

The original method is in [the QHAM paper by Xue et al.](https://doi.org/10.1007/s11433-024-2584-2).
The derivation of the finite closure, the tensor-block weights, and the physical output is in the [QHAM mathematical notes](../reference/qham-derivation.md);
general polynomials, forcing, and discretization restrictions are in the [full implementation documentation](../manual/qham.md).

## Carleman: open coefficient ports and a replaceable linear solver

The Carleman route starts from the coefficient ports of a polynomial ODE and builds the tensor lifting according to the
chosen cutoff.
Its truncation object differs from QHAM's; the two truncation orders must not be treated as one shared accuracy parameter.

```{literalinclude} ../../examples/ode_input_models.py
:language: python
:start-at: "    # 6. Ordinary PDE"
:end-before: "    # 7."
:dedent: 4
:caption: Carleman example fragment; imports, grids, and the linear-solver configuration are in the full script.
```

Here `concrete.coefficients` stores the coefficient BEs grouped by nonlinear degree. The loop establishes, for each
coefficient, an open declaration with the same width, signal layout, and alpha, while keeping the matching implementation
in `bindings`. `BurgersInitial` is likewise an open state-preparation slot; the original initial-state norm is carried
separately by the host.

{obj}`PolynomialODE <oracq.algorithms.qnlss.carleman.PolynomialODE>` therefore describes the mathematical input; {obj}`carleman_qode <oracq.algorithms.qnlss.carleman.carleman_qode>` performs the lifting; the passed-in
`linear_solver` handles the lifted linear problem. The script uses Schrödingerization
and an explicitly shifted LCHS respectively, and writes the recovery factors into the records. Running the generator
twice does not mean the two methods agree in success probability, truncation error, or circuit cost.

Mathematical conventions and interfaces are on the <a href="../../zh/manual/algorithms/carleman.html">Carleman page</a>.
Download the {download}`full input-model script <../../examples/ode_input_models.py>` with all the context.

## LCHS: gate-table and QRAM bindings on the same open graph

An LCHS configuration contains a quadrature plan and a simulation kernel for the Hermitian branch. The Cauchy nodes
and the first-order Taylor kernel in the example are concrete approximation choices, not a target accuracy that the
language guarantees by default.

```{literalinclude} ../../examples/ode_input_models.py
:language: python
:start-at: "    # 3. Same open graph"
:end-before: "    # 4."
:dedent: 4
:caption: Without changing the caller, only the angle-database and initial-state bindings differ.
```

`DiagonalAngles` stores rotation angle words, and {obj}`diagonal_block_encoding <oracq.algorithms.input_model.oracles.diagonal_block_encoding>` interprets them as
matrix coefficients. `Initial` is an independent initial-state preparation slot. The two `write` calls receive the same open
`state` but choose the gate-table and QRAM implementations respectively; the latter must additionally supply the two
tables in `memory`.

The script's `save_case` saves the open, the first partially bound, and the final closed RIR separately,
and produces modular backend text. Intermediate resource capture remains an explicit parameter; memory contents are
never stuffed into the RIR. For a smaller standalone example see [comparing implementations of an open program](algorithm-research.md);
method configuration and applicability conditions are on the <a href="../../zh/manual/algorithms/lchs.html">LCHS page</a>.

## CBMD: reusing QHAM input while handling the dissipative precondition explicitly

```{literalinclude} ../../examples/input_models.py
:language: python
:start-at: "    # ---- CBMD: QHAM"
:end-before: "    (args.output"
:dedent: 4
:caption: QHAM ports that were already built are handed to CBMD through an explicit adaptation.
```

CBMD's method originates in [contour matrix decomposition by Wang et al.](https://doi.org/10.1088/2058-9565/ae7b7e).
This code reuses the QHAM plan and the open coefficient input instead of reinterpreting the PDE inside CBMD.
The shift and recovery information produced by `model.dissipative_shift()` is a mathematical adaptation; only the final
`qram_dict` and `qram_memory` constitute the implementation and run-time data binding.

Nodes, residues, the finite truncation, and the omitted terms are on the <a href="../../zh/manual/algorithms/cbmd.html">CBMD page</a>.
Do not infer from `check().ok` that the contour truncation error or the physical recovery error has been certified.

## Input-model variants and reproduction

The {download}`input-model variant script <../../examples/input_models.py>` provides structured ports, spectral
constructions, and QRAM data paths for four families. It is a construction example, not a cross-method accuracy
leaderboard.

```bash
PYTHONPATH=src python examples/general_qham.py
PYTHONPATH=src python examples/ode_input_models.py
PYTHONPATH=src python examples/input_models.py
PYTHONPATH=src python examples/research_workflow.py
PYTHONPATH=src python tools/build_qlss_comparison.py
```

Real-backend checks must use an interpreter with both `pysparq` and `uniqc` installed:

```bash
PYTHONPATH=src /path/to/backend/python examples/research_workflow.py --native
PYTHONPATH=src /path/to/backend/python tools/run_verification.py \
  --group oracles --group blockencoding --group arithmetic \
  --group mathfunc --group qham_qfvm --group ode
```

Numerical reports should explain three things separately: the implementation error of the circuits against an
independent finite-precision reference, the method error of the finite precision or truncation relative to the original
problem, and the post-selected success rate together with physical norm recovery. "The program exports", "the backends
agree with each other", and "the method converges on the original PDE" are distinct conclusions. The full boundaries
are in [Applicability boundaries and validation status](../manual/limits.md).

## Format details, resource data, and extras

| Detail to look up | Where it is maintained |
|---|---|
| RIR nodes, the full serialization text, and the OriginIR-ext grammar | [RIR specification](../reference/rir.md) |
| Open declarations, batched binding, resource capture | [Open IR](../reference/open-ir.md) |
| Intermediate representation and serialization of mathematical functions | [Math IR](../reference/math-ir.md) |
| Toffoli/rotation/QRAM counts, generic scale sweeps | [Resource estimation](../manual/resource-estimation.md) |
| QFVM inputs, Roe arithmetic, and the recovery chain | [QFVM documentation](../manual/qfvm.md) |
| The wider algorithm catalog and validation status | <a href="../../zh/manual/algorithms/index.html">Algorithm catalog</a>, [validation coverage](../development/validation-coverage.md) |
| Classes, functions, parameter names, and returned objects | [API reference](../api/index.rst) |

These pages carry the full implementation and tutorial detail. Non-core algorithms remain maintained in the library,
but their count is not used as evidence of the framework's superiority in the paper. The cross-framework code-size
experiment scripts remain under `tools/expressiveness/`; they are not part of the evidence retained for the paper.
