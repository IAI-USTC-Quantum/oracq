# Composing PDE inputs and solve instances

**English** · <a href="../zh/tutorials/pde-instances.html">简体中文</a>

`PDEProblem` describes an equation, a grid, physical initial data, and the final
time. `qpde_solve` selects a method; `prepare()` checks its input requirements
and generates a reusable modular RIR circuit. Export, resource estimation, and
native execution all use that same prepared circuit. Neither problem creation
nor preparation imports a quantum backend.

## Prepare, inspect, and execute

```{testcode}
import oracq

problem = oracq.PDEProblem(
    type="heat",
    grid=oracq.UniformGrid1D(0, 2, 2),
    initial=oracq.ArrayInput([0.2, 0.1]),
    viscosity=0.1,
    final_time=0.01,
)
instance = oracq.qpde_solve(
    problem, method="QHAM",
    config=oracq.QHAMConfig(order=1, eta=-1, taylor_degree=1),
)
instance.prepare()
circuit = instance.circuit()
assert instance.prepared.report.ok
assert oracq.loads(circuit.dumps()) == circuit.program
assert circuit.resource_estimate(basic_gates="clifford+t+qram").complete
assert "DEF " in circuit.originir_ext().text
```

`prepare()` is cached. Inspection or execution before preparation raises an
error. For a different configuration, create another instance. An invalid
initial width, missing derivative table, unbound oracle, or missing memory
snapshot fails preparation before backend execution.

With real optional backends installed, execute either path:

```python
results = instance.run_pysparq()
results = instance.run_originir_ext()
print(results.success_probability)
print(results.physical_amplitudes())  # conditional normalized amplitudes
print(results.physical_values())      # magnitude of the generated approximation
```

Both return `PDEResult` with the public register state and field-major physical
readout. Padding addresses are excluded. A zero physical solution has zero
success probability, recoverable zero values, and no normalized conditional
state. `run_originir_ext` executes the exported artifact through UnifiedQuantum;
`originir_ext()` only exports text. Execution has explicit state, instruction,
or qubit budgets; preparing and estimating a circuit does not guarantee it fits
a simulator's budget.

The default QHAM configuration uses a finite Taylor polynomial of the lifted
linear generator. `physical_values()` recovers that polynomial's known LCU
normalization; it does not assert convergence to the continuous PDE. Spatial
discretization error, HAM truncation/convergence, angle quantization, and linear
solver error must be assessed independently.

## Replace the initial input independently

| Input | Access supplied | Physical norm |
| --- | --- | --- |
| `UniformInput(value)` | Constant amplitudes over physical nodes and fields | Computed |
| `ArrayInput(values, encoding="gates")` | Gate preparation of real or complex samples | Computed |
| `ArrayInput(values, encoding="qram")` | Quantized rotation-tree preparation of nonnegative real samples | Computed before angle quantization |
| `QRAMInput(bank, config)` | Existing rotation-tree angle words, with a base offset | Explicit |
| `OracleInput(preparation, norm, qrams=...)` | An existing reversible amplitude-preparation oracle | Explicit |

The same initial input can be combined with a uniform grid or a QRAM grid.
Arrays describe physical nodes only, with one contiguous array per field.
For a coupled PDE, use `ArrayInput({"u": u_samples, "v": v_samples})`;
register padding is inserted automatically. Custom signed or complex QRAM
loaders can be wrapped in `OracleInput` along with their memory banks.

Raw value access `|i,0> -> |i,data[i]>` does not itself prepare the amplitude
state `sum_i data[i]|i>`. `register_qram` therefore returns a host-side bank,
not an amplitude oracle. Its `database()` method supplies XOR value access.
An input's adapter explicitly chooses the preparation circuit and norm.

This example prepares uniform physical amplitudes from an externally stored
root angle, using a four-bit address bank for a one-bit target:

```{testcode}
import math

bank = oracq.register_qram({3: 16}, address_length=4, word_length=6,
                           name="initial_angles")
initial = oracq.QRAMInput(
    bank, oracq.QRAMInputConfig(state_width=1, norm=math.sqrt(0.08), offset=3),
)
prepared = initial.prepare_initial(oracq.InitialLayout(problem.grid, ("u",)))
state = oracq.simulate(prepared.preparation.operation.program(),
                      {bank.name: bank.snapshot()})
assert abs(state.amplitudes[(0, 0)] - 1 / math.sqrt(2)) < 1e-12
assert abs(state.amplitudes[(1, 0)] - 1 / math.sqrt(2)) < 1e-12
```

Rotation-tree addresses are `offset + 2**depth - 1 + prefix`; words encode Ry
angles in units of `2*pi/2**word_length`. The tree must represent the intended
physical state, including zero amplitudes at padded addresses. The declared
norm is a physical input promise, separate from normalized amplitudes.

## Choose the grid representation

`UniformGrid1D(0, 10, 10)` has ten physical samples and a four-bit address
register. Periodic coordinates cover `[0, 10)` with spacing `1`; zero Dirichlet
coordinates are interior nodes with spacing `(stop-start)/(points+1)`.

Power-of-two periodic grids use shift circuits. Other sizes and zero Dirichlet
boundaries currently use a spatial Pauli expansion of at most five qubits.
This bounded implementation supports ten samples but is not an efficient
large-grid implementation. It expands individual spatial derivatives rather
than a complete nonlinear port or lifted matrix. An application-defined
`PDEGrid` can provide larger-grid encodings directly.

For an unstructured mesh, store the **discrete derivative operators** after
choosing the geometric discretization and boundary treatment. Coordinates or
connectivity alone do not determine these operators. The provided layout uses
real matrices in an angle bank and supports nonsymmetric operators:

```{testcode}
d2 = (("x", 2),)
layout = oracq.UnstructuredGridConfig(
    size=2, derivative_offsets={d2: 0}, value_scale=3,
)
data = layout.encode_matrices({d2: [[-2, 2], [2, -2]]}, word_length=8)
mesh_bank = oracq.register_qram(data, address_length=2, word_length=8, name="mesh")
mesh = oracq.UnstructuredGrid(mesh_bank, layout)
mesh_problem = oracq.PDEProblem(type="heat", grid=mesh,
                                initial=oracq.UniformInput(0.2), final_time=0.01)
mesh_instance = oracq.qpde_solve(mesh_problem).prepare()
assert mesh_instance.circuit().memory
assert mesh_instance.circuit().resource_estimate().qram_total > 0
```

A derivative table uses `base + row*2**spatial_width + column`. Each word
represents `value_scale*cos(pi*word/2**word_length)`; in particular, a zero
matrix element uses a half-turn word, whereas a zero word represents positive
`value_scale`. `encode_matrices` inserts padding and saturates the negative
endpoint instead of wrapping its sign. Its per-entry error bound is
`value_scale*pi/2**word_length`. The coherent encoding uses two queries and
normalization `2**spatial_width*value_scale`; QRAM storage does not imply a
sparse-access speedup. Layouts using neighbor lists or other matrix-access
models should implement `PDEGrid` with their own verified block encoding.

Memory snapshots are immutable host data outside RIR. During preparation,
resource arguments are traced through the module graph to obtain all entry
bindings, including multiple aliases of a shared bank. This traversal does not
expand `Repeat`. Conflicting bank names or incompatible widths fail explicitly.

## Algorithms own their requirements

QHAM requires reversible, controlled block encodings for every multilinear
port and an initial preparation with matching width, adjoint/control support,
zero-input semantics, and clean work. The contract report is available as
`instance.prepared.report`; the lifted model is `instance.prepared.model`.

To use an existing PDE expression, pass a `PolynomialPDE` as `PDEProblem.type`.
The string presets are currently `burgers` and `heat`. To choose the linear
evolution separately, supply a three-argument QODE callable as
`QHAMConfig(linear_solver=...)`. Its Hermiticity or dissipation requirements
still apply; the wrapper does not silently shift a generator. Such a callable
does not automatically supply magnitude recovery to this instance, so
`physical_values()` raises while conditional amplitude readout remains usable.

Other PDE algorithms implement `PDEMethod.prepare(problem) -> PreparedPDE`.
They own their grid restrictions, oracle adaptations, input contract report,
and mathematical promises. Pass the method object to `qpde_solve`; it shares
the cached lifecycle, RIR export, resource estimation, backend execution, and
physical readout. QHAM is currently the only built-in named PDE method.

`clifford+t+qram` estimates include a seven-T Toffoli decomposition;
`toffoli+clifford+t+qram` retains Toffoli as a separate atom. Arbitrary rotations
remain pending synthesis and `t_total(epsilon)` estimates their T cost. QRAM
queries are reported separately from physical QRAM hardware costs.

Run the complete example, with native execution optional:

```bash
python examples/pde_instances.py
python examples/pde_instances.py --grid qram --initial qram --backend pysparq
python examples/pde_instances.py --backend originir_ext
```

See also [QHAM expression assembly](qham.md), [QHAM implementation](../manual/qham.md),
[PDE input API](../api/algorithms/input_model/pde.rst), and
[solve instance API](../api/algorithms/qpde/solver.rst).
