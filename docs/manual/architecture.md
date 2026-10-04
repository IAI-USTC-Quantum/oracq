# Source code structure

**English** · <a href="../zh/manual/architecture.html">简体中文</a>

The source is divided into infrastructure, the algorithm library, and domain
applications. Directory names correspond to the responsibilities of the code,
and algorithm implementations live in separate files grouped by category.

```text
src/oracq/
├── infrastructure/       RIR, Builder, validation, serialization, execution, and backends
│   ├── backends/         OriginIR-ext, strict gate set, and PySparQ adapters
│   └── mathfunc/         pure math function frontend, MIR, and reversible lowering
├── algorithms/           quantum algorithms, oracle interfaces, and composition tools
│   ├── fourier.py        Fourier transform and Fourier addition
│   ├── search.py         Grover and amplitude amplification
│   ├── estimation.py     QPE, amplitude estimation, Hadamard/Swap tests
│   ├── hamiltonian.py    Hamiltonian evolution and swappable implementations
│   ├── qlss.py           linear-system contracts, Costa, and CKS
│   ├── lchs.py           LCHS
│   ├── schrodingerization.py
│   ├── carleman.py       Carleman linearization
│   └── ...
└── applications/         QFVM, Roe data, QHAM math support, and example directories
```

## Infrastructure

`infrastructure` provides the structural description and the execution
machinery; it does not choose circuits by algorithm name. The math function
frontend may call arithmetic generators from the algorithm library, and
convenience view methods on
{obj}`Operation <oracq.infrastructure.builder.Operation>` also call algorithm
adapters when needed. These calls happen during the Python generation stage
and add no new kinds of RIR instruction.

Optional backends are imported at the execution entry point (see
[Export and execution backends](backends.md)). Generating and serializing RIR
does not require PySparQ or UnifiedQuantum to be installed.

## Algorithm library

Each category in `algorithms` has its own file and input conventions. The
shared [`interfaces.py`](../api/algorithms/input_model/interfaces.rst)
provides common Python structural protocols;
[`contracts.py`](../api/algorithms/input_model/contracts.rst) provides check
reports. Applications may still define their own protocols.

[`operators.py`](../api/algorithms/input_model/operators.rst) holds
block-encoding wrapper classes and basic scaling operations, and
[`block_encoding.py`](../api/algorithms/input_model/block_encoding.rst)
provides LCU, tensor, and small-matrix constructions. They are the composition
tools of the algorithms. [`ode.py`](../api/algorithms/qode/ode.rst) provides
the swappable linear-solver interface; LCHS, Schrödingerization, CBMD, and
Carleman each maintain their own implementations.

## Domain applications

The physics parameters, Roe formulas, and classical data updates of
[QFVM](qfvm.md) live in `applications`. The quantum assembly entry point of
[QHAM](qham.md) is
[`algorithms/qham.py`](../api/algorithms/input_model/qham.rst); PDE
expressions, homotopy derivations, and spatial discretization support live in
`applications/qham/`. This way the quantum generation steps can be read
without unfolding the entire mathematical derivation code at the same time.

Early demonstration implementations are concentrated in
`applications/legacy.py` and `algorithms/legacy.py`. They exist for
compatibility with old cases; the current documentation does not treat them
as the default entry point for new applications.

## Compatibility and migration

Old import paths are forwarded centrally by `_compat.py` in the repository
root. The compatibility layer references the same function or class and does
not keep a second implementation. In-repository source, tests, and tutorials
use the canonical paths; external code can migrate step by step, see
[Import path migration](compatibility.md).
