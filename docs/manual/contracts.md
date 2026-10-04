# Conventions owned by algorithms: starting from one gate

**English** · <a href="../zh/manual/contracts.html">简体中文</a>

**Conventions belong to the algorithm library, not to the RIR language's type system.** An object can satisfy multiple Python protocols at once; higher-level algorithms check the methods, parameters, and invocation capabilities they need. A new algorithm can define new protocols in its own file without modifying the oracq syntax, the RIR, the serializers, or a global type catalog.

This chapter explains how algorithms declare, check, and adapt inputs. Complete runnable example: [examples/algorithm_contracts.py](../../examples/algorithm_contracts.py). For the longer mathematical assembly guide, see the [QPDE/QODE guide](differential-equations.md).

## 1. Protocols and algorithm generators

New code uses {obj}`AlgorithmContract <oracq.algorithms.input_model.contracts.AlgorithmContract>`, {obj}`QLSSSolver <oracq.algorithms.qlss.qlss.QLSSSolver>`, and {obj}`QODESolver <oracq.algorithms.qode.ode.QODESolver>`. The old names
{obj}`ProtocolContract <oracq.algorithms.input_model.contracts.ProtocolContract>`, {obj}`QLSSProtocol <oracq.algorithms.qlss.qlss.QLSSProtocol>`, and {obj}`QODEProtocol <oracq.algorithms.qode.ode.QODEProtocol>` are compatibility aliases of the same types.
A Python protocol describes the access interface; the solver is responsible for generation, the contract for checking, and open RIR modules
for holding implementations that have not been provided yet. These mechanisms act on the generation and linking boundaries, respectively.

Python's `typing.Protocol` states "which interface an object satisfies"; the QLSS/QODE protocols in the library mean "replaceable algorithm generators". The two are related but different:

```text
the input object satisfies some Python interface
    -> the algorithm generator checks the input and selects an implementation
    -> new oracles/operations are generated
    -> the generated result again satisfies some Python interfaces
```

These steps happen during the Python generation stage, not as type dispatch at quantum-circuit execution time. The RIR still stores registers, gates, {obj}`Call <oracq.infrastructure.ir.Call>`, {obj}`Repeat <oracq.infrastructure.ir.Repeat>`, and open declarations; it does not store Python callbacks.

## 2. One gate is simultaneously a unitary, a state preparation, and an LCU input

Here is the complete code snippet:

```python
from oracq import Bits, Builder, identity, requires
from oracq.algorithms.input_model.interfaces import UnitaryProtocol, StatePreparationProtocol, BlockEncodingProtocol
from oracq.algorithms.input_model.block_encoding import lcu

b = Builder("XGate", {"q": Bits(1)})
b.x(b["q"])
gate = b.finish()

requires(gate, UnitaryProtocol)
requires(gate, StatePreparationProtocol)
requires(gate, BlockEncodingProtocol)

initial = gate.state_preparation()   # U|0> = |1>
encoded_u = gate.block_encoding()  # U itself: alpha=1, zero signal bits
sum_encoding = lcu([(1, gate), (1, identity(1))])  # encodes U+I, alpha=2
```

`gate` is not given three string labels here. {obj}`Operation <oracq.infrastructure.builder.Operation>` provides three actual methods, so the structural protocols hold: `unitary()`, `state_preparation()`, and {obj}`block_encoding() <oracq.algorithms.input_model.operators.block_encoding>`.

**By default, all public registers of an Operation are treated as the target space of the full unitary.** Multiple registers are concatenated from the least significant to the most significant position in signature order. If an operation exposes a 2-bit target and 3 bits of work, treating the whole Operation as a unitary state preparation yields a 5-bit state; it is not silently assumed that those 3 bits can disappear.

If you know the work is clean after a zero-input preparation, you can explicitly narrow the target interpretation of the state:

```python
from oracq.algorithms.input_model.oracles import StatePreparation

# the operation's interface must be exactly q / tmp.
# prep = StatePreparation.from_unitary(operation, target="q", work="tmp", clean_work=True)
```

`clean_work=True` is your algorithmic promise about this implementation; the language does not prove it. The default all-public-registers adaptation needs no such promise; the single target of the current BE/state-prep convenience wrappers is still limited to 64 bits.

This also clarifies "whether an arbitrary unitary can serve as b": it always defines some `b=U|0>`. Whether that is the b of the problem you want to solve is defined by the application; whether it supports the inverse/controlled forms needed by QLSS is checked by the concrete QLSS.

## 3. QLSS accepts this gate directly as b

```python
from oracq import BlockSystem, LinearSystem, SpectralPromise, identity
from oracq.algorithms.qlss.qlss import CostaConfig, make_costa_qlss

problem = LinearSystem(
    block=BlockSystem(identity(1), gate, SpectralPromise(1.0, 1.0)),
    rhs_norm=1.0,
)
solver = make_costa_qlss(CostaConfig(steps=1))
report = solver.check(problem)
report.require()
result = solver(problem)
```

{obj}`BlockSystem <oracq.algorithms.qlss.qlss.BlockSystem>` obtains the gate as an initial-state interface through `state_preparation()`. Costa's requirements can be inspected directly:

```python
print(solver.contract.to_dict())
print(report.to_dict())
print([p.__name__ for p in solver.provides])  # ['StateOracleProtocol']
```

The current Costa composition entry requires A to provide a BE, or CKS sparse input that can be converted explicitly; b to provide a state preparation with zero input and a clean workspace; the inputs must have the same width and support the adjoint/controlled forms actually needed by the composition. Spectral bounds remain a mathematical declaration of the problem.

`result.state_oracle()` yields the output state oracle. The output state oracle carries a success signal; it does not automatically become "a clean initial-state preparation without post-selection". If the next step only needs its full physical unitary, you can explicitly take `result.operation`, in which case all of its signal registers also belong to the full space.

`check` does not run the QLSS kernel or a quantum simulation. It may call methods of the input objects to obtain/generate concrete access views, so custom adaptation methods should be deterministic and free of external side effects; when construction is expensive, the object itself may cache the result.

`contract.resolve(**inputs)` returns the views and reports obtained this time. Algorithms should use
`resolved.get(name, ViewType)` to obtain these already-checked views, avoiding calling the provider again after the check.
One solve entry reuses the views of that call; calling `check()` independently and then calling the solver counts as two operations
and does not share a hidden cache. Adaptation entries check callability and the parameterless signature that can be obtained; exceptions raised inside provider method bodies keep their original traceback.

## 4. Where A's parameters are read from

Existing BE wrapper classes keep their original API and add read-only aliases:

| Read | Meaning |
|---|---|
| `A.width` / `A.main_qubit` | bit width of the target register |
| `A.signal_qubits` / `A.anc_qubit` | number of public signal bits of the BE, excluding the resource peak of module-private locals |
| `A.alpha` | normalization constant of the current encoding |
| `A.capabilities.adjoint` / `.controlled` | effective invocation capabilities derived from the whole dependency graph |
| `A.spec` | JSON-serializable description snapshot covering registers, resources, and open/closed state |
| `A.type` | descriptive role name such as `block_encoding`; not a unique type or a dispatch basis |

If A is your own operator object, first obtain the concrete view this algorithm needs through `A.block_encoding()`, then read these parameters. Algorithms work against interface methods; there is no need to append type names to a central enumeration.

The CKS input is a collection of a position operation and an entry operation: `SparseAccess.sparse_access()` returns itself; its `spec.components` describe position/entry respectively, with sparsity and value_width recorded in the parameters. The collection as a whole has `anc_qubit=None`, avoiding the pretense that two different query operations are the total ancilla count of one and the same unitary.

## 5. Providing new inputs without inheriting from oracq

```python
from oracq import identity, requires
from oracq.algorithms.input_model.interfaces import BlockEncodingProtocol

class MyMatrix:
    def block_encoding(self):
        return identity(2)

A = MyMatrix()
requires(A, BlockEncodingProtocol)
```

The algorithm library performs structural checks with `@runtime_checkable typing.Protocol`; `MyMatrix` inherits from no base class and performs no registration.

{obj}`requires <oracq.algorithms.input_model.contracts.requires>` checks that the interface exists; it does not prove the method signature or the mathematical meaning of the return value. After calling `block_encoding()`, a concrete algorithm still checks the returned BE type, width, alpha, and required capabilities. For example, returning a string produces an `INPUT_ADAPTER` report; it is not accepted merely because "a method with the same name happens to exist".

New conventions can live entirely inside the application:

```python
from typing import Protocol, runtime_checkable
from oracq import requires

@runtime_checkable
class HasDiagonal(Protocol):
    def diagonal_values(self) -> tuple[float, ...]: ...

def my_diagonal_algorithm(operator):
    requires(operator, HasDiagonal, path="my_diagonal_algorithm.operator")
    values = operator.diagonal_values()
    # this algorithm continues to check length and value range, and generates its own operations.
    return values
```

`HasDiagonal` here is not a new language type; exporters never need to know about it.

## 6. What Hermitian, Trotterizable, and QSP each mean

Mathematical operators and the physical unitaries implementing them must be kept apart: a non-Hermitian A can still have a unitary block encoding `U_A`, but A must not be treated as a Hermitian Hamiltonian on that basis.

The current [hamiltonian.py](../api/algorithms/common/hamiltonian.rst) defines the interface owned by this family of algorithms:

| Protocol | What it provides | Who judges |
|---|---|---|
| {obj}`HermitianProtocol <oracq.algorithms.common.hamiltonian.HermitianProtocol>` | the `.hermitian` declaration | HamSim requires it to be True; mathematical truthfulness is the implementer's responsibility |
| {obj}`BlockEncodingProtocol <oracq.algorithms.input_model.interfaces.BlockEncodingProtocol>` | `.block_encoding()` | BE/QSP-style implementations obtain the access model and alpha |
| {obj}`TrotterizableProtocol <oracq.algorithms.common.hamiltonian.TrotterizableProtocol>` | `.trotter_list()` | Trotter implementations obtain the ordered, coefficient-carrying decomposition |
| {obj}`EvolvableProtocol <oracq.algorithms.common.hamiltonian.EvolvableProtocol>` | `.evolution(t)` | a single term provides the concrete unitary evolution of `exp(-it H_j)` |

One object can satisfy several of these protocols at once. The object below has only Trotter access and no BE:

```python
from oracq.algorithms.common.hamiltonian import PauliOperator, TrotterTerm, hamiltonian_simulation

class MyHamiltonian:
    hermitian = True

    def trotter_list(self):
        return (
            TrotterTerm(0.3, PauliOperator("I")),
            TrotterTerm(0.7, PauliOperator("X")),
        )

evolution = hamiltonian_simulation(MyHamiltonian(), 0.4, steps=3)
```

Here `H=Σ c_j H_j`; the generator calls `evolution(c_j*t/steps)` for each term in list order, then keeps the repeated structure with the RIR {obj}`Repeat(steps) <oracq.infrastructure.ir.Repeat>`. It returns a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` with alpha=1. In the generally non-commuting case this is a product-formula approximation; accuracy is the responsibility of the configuration and the algorithm analysis.

{obj}`TrotterTerm <oracq.algorithms.common.hamiltonian.TrotterTerm>` is more than a matrix name. The current implementation requires each single-term evolution to return a concrete `Operation` whose only non-zero public register is `target`, and all terms must have the same width. Exact HamSim that needs ancilla registers can extend this algorithm's contract later; a general Taylor BE with post-selection signals must not be passed off as an exact single-term unitary evolution.

{obj}`hamiltonian_simulation(method="auto") <oracq.algorithms.common.hamiltonian.hamiltonian_simulation>` currently prefers an available Trotter decomposition and otherwise tries the BE/QSP path. `method="qsp"` requires injecting an actual `qsp(BE,time)->BlockEncoding` generator; the library currently ships no general QSP-HamSim kernel and fails explicitly when it is missing. This is a local, replaceable selection strategy, not "the language automatically implements QSP for every Hermitian input". Additional premises of the concrete QSP are still checked by the injected implementation.

{obj}`EncodedOperator(encoding, hermitian=False) <oracq.algorithms.common.hamiltonian.EncodedOperator>` can represent a non-Hermitian operator and provide a BE; HamSim rejects it, while QODE can consume its BE and check the conditions with its own methods.

## 7. LCHS and Schrödingerization each own their conventions

```python
from oracq import QODEProblem, identity, scale
from oracq.algorithms.qode.ode import linear_qode

problem = QODEProblem(
    generator=scale(-1, identity(1)),
    initial=gate,
    dissipative=True,
    initial_norm=1.0,
)
lchs = linear_qode("lchs")
schrodinger = linear_qode("schrodingerization")
lchs.check(problem, time=0.1).require()
state = lchs.solve(problem, 0.1)
other = schrodinger.solve(problem, 0.1)
```

{obj}`QODEProblem <oracq.algorithms.qode.ode.QODEProblem>` and `QODESolver` are ordinary algorithm-library objects. The problem-level `solve` of LCHS/CBMD requires an explicit `dissipative=True`; Schrödingerization does not require this declaration, but the application must still guarantee a valid auxiliary grid and recovery region. `check().ok` means the structure and declarations satisfy the requirements — never that properties of the PDE have been proven.

The old `(G, initial, time)` call remains compatible, with the mathematical premises still borne by the caller as before; new applications are encouraged to use `.solve(QODEProblem(...),time)`, so that undeclared dissipation and explicitly non-dissipative inputs are not silently accepted. No "dissipative matrix" type was added to the RIR.

Carleman itself requires the layout of the multilinear coefficient ports `F_p` and the initial-state norm, then hands the generated linear system to the selected QODE protocol. The lifted system is not necessarily dissipative; shift explicitly when needed. These mathematical boundaries are still covered in the [QPDE/QODE guide](differential-equations.md).

## 8. requires, reports, and binding each do one thing

`requires(value, Interface)` is the minimal tool for writing an algorithm's own checks. When you need to aggregate errors and display contracts, use {obj}`InputRequirement <oracq.algorithms.input_model.contracts.InputRequirement>` and `AlgorithmContract`, passing in the Python protocol classes this algorithm accepts together with adaptation functions. They do not depend on a closed set of string types.

Errors in the report contain `code/path/expected/actual/message`. Common errors include `INPUT_PROTOCOL` (missing interface), `INPUT_ADAPTER` (failed to obtain a view), `INPUT_WIDTH`, `INPUT_CAPABILITY`, `INPUT_PROMISE`, and `OUTPUT_LAYOUT`. `report.require()` raises all discovered problems at once; the exception remains a subclass of {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` for compatibility with older code.

{obj}`bind <oracq.infrastructure.linking.bind>` continues to handle only the ABI, alpha, invocation capabilities, and resource wiring of existing RIR slots. It is not responsible for knowing `HermitianProtocol`, `TrotterizableProtocol`, or your new mathematical properties; the mathematical truthfulness of declarations remains with the algorithms and implementers. Same-signature/alpha allows late binding; otherwise rerun the host generator.

`A.spec` / `contract.to_dict()` / `report.to_dict()` are check reports that can be stored as JSON, but they are not another executable IR. Running implementations are still saved through `dumps/loads` on the RIR; entry operations can be restored with `Operation.from_program(program)`. Protocol objects, caches, and arbitrary host methods are never silently serialized.

## 9. Running and engineering boundaries

```bash
PYTHONPATH=src .venv/bin/python examples/algorithm_contracts.py
PYTHONPATH=src .venv/bin/python -m oracq inspect out/algorithm-contracts/qlss.closed.rir.yaml
```

The example output includes accept/reject reports, an open and a closed QLSS, an LCU of the same unitary, two QODE descriptions, and the modular OriginIR of the Trotter run.

Algorithm conventions extend at the Python layer while the RIR format stays at 0.3. For production-use boundaries and the full acceptance procedure see [engineering maturity](limits.md): the language and assembly checks pass engineering tests; the full numerical correctness of Costa/CKS, LCHS, Schrödingerization, Carleman, and QHAM remains experimental. The quantum solving stack as a whole cannot yet be called a certified production solver.
