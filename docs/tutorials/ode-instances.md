# Composing linear ODE inputs and solve instances

**English** · <a href="../zh/tutorials/ode-instances.html">简体中文</a>

`ODEProblem` describes the autonomous homogeneous system `u' = G u`, physical
initial values, and a final time. `qode_solve` selects a method and returns an
unprepared instance. As with [PDE instances](pde-instances.md), `prepare()`
checks input contracts and generates one cached modular RIR circuit for
export, resource estimates, and native execution.

## Prepare, inspect, and execute

The matrix is the generator **G**, including its sign. For example, diagonal
entries `-1` and `-2` represent decay, following `u' = G u`:

```{testcode}
import oracq

problem = oracq.ODEProblem(
    generator=oracq.MatrixInput([[-1, 2], [0, -2]]),
    initial=oracq.ArrayInput([2, -1]),
    components=("u", "v"),
    final_time=0.1,
)
instance = oracq.qode_solve(
    problem, method="taylor", config=oracq.ODEConfig(taylor_degree=2),
).prepare()
circuit = instance.circuit()
assert instance.prepared.report.ok
assert oracq.loads(circuit.dumps()) == circuit.program
assert circuit.resource_estimate().complete
assert "DEF " in circuit.originir_ext().text
assert instance.prepare().circuit() is circuit
```

Preparation and export import neither PySparQ nor UnifiedQuantum. With those
optional backends installed, execute the same prepared circuit:

```python
result = instance.run_pysparq()
result = instance.run_originir_ext()
print(result.success_probability)
print(result.physical_amplitudes())
print(result.physical_values())  # approximately (1.64, -0.82) for this polynomial
```

`physical_amplitudes()` normalizes the physical success channel by default;
`normalize=False` returns the raw success amplitudes. `physical_values()` uses
the known LCU normalization and initial norm to recover the finite Taylor
polynomial's values. This is a quantum circuit evaluation of
`(I + t G + t² G² / 2) u(0)`, with independent Taylor truncation error; it does
not substitute a classical solver's output for backend amplitudes.

The vector contains physical components in the declared order. A scalar uses
one target bit; other dimensions use `ceil(log2(size))` bits, with a minimum of
one. Unused addresses are excluded from readout, and the linear method masks
both sides of G to prevent paths through padded components. A zero initial
vector produces zero physical values and zero success probability, with no
normalized conditional state.

## Replace the generator and initial inputs

| Generator input | Access model | Scope |
| --- | --- | --- |
| `MatrixInput(values)` | Pauli gate expansion | Finite real or complex matrices, at most five target bits |
| `MatrixInput(values, encoding="qram")` | Dense real matrix angle bank | Finite real entries bounded by `value_scale` |
| `QRAMMatrixInput(bank, config)` | Existing dense matrix angle bank | Explicit dimension, offset, scale, and immutable memory |
| `OracleGeneratorInput(encoding, size, qrams=...)` | Existing block-encoding provider | Caller supplies operator semantics and any leaf memory banks |

These choices are independent of `UniformInput`, `ArrayInput`, `QRAMInput`,
and `OracleInput`, the [shared initial inputs](../api/algorithms/input_model/initial.rst)
also used by PDEs. `ArrayInput` takes one flat vector of physical components;
named mappings take a single sample per component, such as
`{"u": [2], "v": [-1]}`. Signed and complex gate initial states are supported;
the built-in QRAM rotation tree accepts nonnegative real amplitudes.
`OracleInput` wraps application-specific state preparation with its physical
norm and memory banks.

An existing block encoding can replace the explicit matrix:

```{testcode}
from dataclasses import replace

encoded = problem.generator.prepare_generator(oracq.ODELayout(2))
oracle_problem = replace(
    problem,
    generator=oracq.OracleGeneratorInput(encoded.encoding, size=2),
)
assert oracq.qode_solve(oracle_problem).prepare().prepared.report.ok
```

A registered QRAM bank supplies raw word access, not amplitude preparation or
a block encoding by itself. The following adapter explicitly implements G:

```{testcode}
matrix_config = oracq.QRAMMatrixConfig(size=2, value_scale=2, offset=4)
words = matrix_config.encode_matrix([[-1, 1], [0.5, -2]], word_length=4)
bank = oracq.register_qram(words, address_length=4, word_length=4,
                           name="ode_generator")
qram_problem = replace(
    problem,
    generator=oracq.QRAMMatrixInput(bank, matrix_config),
    initial=oracq.ArrayInput([1, 1], encoding="qram", angle_width=3),
)
qram_instance = oracq.qode_solve(
    qram_problem, config=oracq.ODEConfig(taylor_degree=1),
).prepare()
assert qram_instance.circuit().memory
assert qram_instance.circuit().resource_estimate().qram_total > 0
```

Matrix words use `offset + row*2**width + column` and represent
`value_scale*cos(pi*word/2**word_length)`. A zero word represents positive
`value_scale`; zero matrix entries require half-turn words. `encode_matrix`
fills padded rows and columns and saturates the negative endpoint. The
per-entry quantization error is bounded by `value_scale*pi/2**word_length`.
The encoding uses two queries per invocation with normalization
`2**width*value_scale`; dense QRAM storage carries no sparse-access speedup.

Memory snapshots stay outside RIR. Preparation resolves resource aliases
through module calls without expanding `Repeat`. Missing snapshots, conflicting
bank names, unbound oracles, and incompatible widths fail before execution.

## Select an existing QODE algorithm

| Named method | Configuration | Preconditions and magnitude readout |
| --- | --- | --- |
| `taylor` (default) | `taylor_degree` | Finite polynomial; known normalization supports `physical_values()` |
| `lchs` | `QuadraturePlan` in `plan` | Requires `dissipative=True`; finite quadrature and Hamiltonian errors remain separate |
| `cbmd` | `ContourPlan` in `plan` | Requires `dissipative=True`; auxiliary-pole contributions and series tails remain omitted |
| `schrodingerization` | `SchrodingerPlan` in `plan` | Caller validates the auxiliary grid and physical recovery region |

For example, reuse LCHS with a caller-supplied quadrature:

```{testcode}
from oracq.algorithms.qode.lchs import QuadraturePlan

decay = oracq.ODEProblem(
    generator=oracq.MatrixInput([[-1, 0], [0, -2]]),
    initial=oracq.UniformInput(), final_time=0.01,
    dissipative=True, evidence="diagonal real nonpositive generator",
)
lchs_instance = oracq.qode_solve(
    decay, method="lchs",
    config=oracq.ODEConfig(plan=QuadraturePlan.cauchy(cutoff=1)),
).prepare()
assert lchs_instance.prepared.report.ok
```

`dissipative=True` records a mathematical promise; preparation checks the
declaration, not a proof of `Hermitian(G) <= 0`. No automatic generator shift
is applied. `ODEConfig.hamiltonian_function` replaces the Hamiltonian simulation
protocol used inside LCHS, CBMD, and Schrödingerization.

These three methods currently expose conditional amplitudes through the
instance, without a complete physical magnitude recovery contract;
`physical_values()` therefore raises. This remains true if an arbitrary
`QODESolver` is passed as `method`. An application-defined
`ODEMethod.prepare(problem) -> PreparedODE` can supply its own verified scale
and reuse the same lifecycle. Configure custom methods directly.

This entry covers autonomous homogeneous **linear** systems. Time-dependent
generators, source terms, and nonlinear ODEs require explicit adaptations.
The existing low-level `QODEProblem`, `QODESolver`, and `linear_qode` interfaces
remain available for circuit composition.

Run the complete example:

```bash
python examples/ode_instances.py
python examples/ode_instances.py --generator qram --initial qram --degree 1 --backend pysparq
python examples/ode_instances.py --backend originir_ext
python examples/ode_instances.py --method lchs
```

See the [ODE input API](../api/algorithms/input_model/ode.rst),
[ODE instance API](../api/algorithms/qode/solver.rst), and
[existing QODE assembly tutorial](differential-equations.md).
