# oracq

**English** · [简体中文](README.zh-CN.md)

[![CI](https://github.com/IAI-USTC-Quantum/oracq/actions/workflows/python-ci.yml/badge.svg)](https://github.com/IAI-USTC-Quantum/oracq/actions/workflows/python-ci.yml)
[![Docs](https://github.com/IAI-USTC-Quantum/oracq/actions/workflows/docs.yml/badge.svg)](https://github.com/IAI-USTC-Quantum/oracq/actions/workflows/docs.yml)
[![PyPI](https://img.shields.io/pypi/v/oracq.svg)](https://pypi.org/project/oracq/)
[![Python](https://img.shields.io/pypi/pyversions/oracq.svg)](https://pypi.org/project/oracq/)
[![License: MIT](https://img.shields.io/badge/License-MIT-informational.svg)](LICENSE)

oracq is a scientific-computing algorithm implementation framework for quantum
algorithm researchers. Composable algorithms are written in Python against
access models and compiled to a register-level intermediate representation
(RIR) that preserves module structure; oracle placeholders stay unimplemented
until you compare gate-network, QRAM, and reversible-arithmetic candidates,
bind them, and proceed to numerical validation and resource analysis.

The current package version is **0.1.0** and the RIR format is **0.1**
([full specification](docs/reference/rir.md)). Algorithms declare inputs and
outputs through plain Python protocols; no language-level type-system
extension is required.

## Core concepts at a glance

| Concept | One-liner | Manual / spec | API reference |
|---|---|---|---|
| RIR | Register-level IR; module calls and Repeat structures survive to the text form | [RIR 0.1 spec](docs/reference/rir.md) | [ir](docs/api/infrastructure/ir.rst) |
| Builder | `Builder`/`Operation` three-stage generation of modules and programs | [Core concepts](docs/manual/concepts.md#modules-and-oracle-placeholders) | [builder](docs/api/infrastructure/builder.rst) |
| Registers and views | bits/uint/qubit interpretations, slicing, reinterpretation; index 0 is the least significant bit | [Core concepts](docs/manual/concepts.md#registers-and-views) | [ir](docs/api/infrastructure/ir.rst) |
| Open oracles and binding | Declare first, bind later: capability conjunction, candidate comparison, QRAM capture | [Binding tutorial](docs/tutorials/oracle-binding.md) | [linking](docs/api/infrastructure/linking.rst) |
| Algorithm contracts | Checkable input protocols, capability specifications, and acceptance reports | [Contracts](docs/manual/contracts.md) | [contracts](docs/api/algorithms/input_model/contracts.rst) |
| Structural validation | Generation-time checks of bit widths, overlaps, and structural constraints | [Core concepts](docs/manual/concepts.md) | [validation](docs/api/infrastructure/validation.rst) |
| Serialization | RIR defaults to YAML with optional JSON; both text forms agree field by field | [RIR spec](docs/reference/rir.md) | [serialization](docs/api/infrastructure/serialization.rst) |
| Execution and readout | Register reference executor and host readout | [Backends](docs/manual/backends.md) | [execution](docs/api/infrastructure/execution.rst), [readout](docs/api/infrastructure/readout.rst) |
| QMem | Pointer-style reads/writes and multidimensional views over QRAM resources | [QMem](docs/manual/qmem.md) | [qmem](docs/api/infrastructure/qmem.rst) |
| QRAM data files | Loading and writing `*.qram.yaml` memory definitions | [QRAM memory](docs/reference/qram-memory.md) | [qram_schema](docs/api/infrastructure/qram_schema.rst) |
| Math-function frontend | Python callables compiled into reversible circuits | [Math functions](docs/manual/math-functions.md) | [mathfunc](docs/api/infrastructure/mathfunc.rst) |
| Resource estimation | Counting models and open analysis for T gates / rotations / QRAM accesses | [Resource estimation](docs/manual/resource-estimation.md) | [estimate](docs/api/infrastructure/estimate.rst) |
| Backend export | OriginIR-ext, strict netlist, PySparQ, quantikz | [Backends](docs/manual/backends.md) | [backends](docs/api/infrastructure/backends/originir.rst) |

## A minimal program

```python
from oracq import Bits, Builder, export_originir, simulate

b = Builder("bell_pair", {"pair": Bits(2)})
b.h(b["pair"][0])
b.xor(b["pair"][0], b["pair"][1])
program = b.finish().program()

print(simulate(program).amplitudes)
print(export_originir(program).text)
```

The result has equal amplitudes on `00` and `11`. Module calls and Repeat
structures are preserved in RIR and in the YAML/JSON text forms; they are not
expanded during generation.

APIs used in this example: [`Builder`](docs/api/infrastructure/builder.rst) ·
[`Bits`](docs/api/infrastructure/ir.rst) ·
[`simulate`](docs/api/infrastructure/execution.rst) ·
[`export_originir`](docs/api/infrastructure/backends/originir.rst).
For a step-by-step walkthrough see the
[first register program](docs/tutorials/first-program.md) tutorial.

## Typical workflow

| Step | Entry point | Further reading |
|---|---|---|
| 1. Declare open oracles | [`declare`](docs/api/algorithms/input_model/oracles.rst), capabilities and specs | [Algorithm contracts](docs/manual/contracts.md) |
| 2. Build the program | [`Builder`](docs/api/infrastructure/builder.rst) composing module calls | [Tutorial: first program](docs/tutorials/first-program.md) |
| 3. Validate structure | [`validate`](docs/api/infrastructure/validation.rst) | [Core concepts](docs/manual/concepts.md) |
| 4. Bind implementations | [`bind`](docs/api/infrastructure/linking.rst); `bind_with_report` compares candidates | [Tutorial: replacing oracles](docs/tutorials/oracle-binding.md), [algorithm-research tutorial](docs/tutorials/algorithm-research.md) |
| 5. Export / execute / estimate | [`export_originir`](docs/api/infrastructure/backends/originir.rst), [`simulate`](docs/api/infrastructure/execution.rst), [`estimate_resources`](docs/api/infrastructure/estimate.rst) | [Backends](docs/manual/backends.md), [resource estimation](docs/manual/resource-estimation.md) |

## Algorithm library

The algorithm library is organized into ten subpackages by purpose; the
complete catalog with selection guidance is the
[algorithm index](docs/manual/algorithms/index.md), with one manual page per
algorithm (interface, implementation notes, validation approach, and known
gaps). Representatives per category:

| Category | Representative algorithm pages |
|---|---|
| Query and search | [Grover search](docs/manual/algorithms/grover.md), [amplitude estimation](docs/manual/algorithms/qae.md), [quantum counting](docs/manual/algorithms/quantum-counting.md) |
| Basic query algorithms | [Deutsch–Jozsa](docs/manual/algorithms/deutsch-jozsa.md), [Simon](docs/manual/algorithms/simon.md), [Bernstein–Vazirani](docs/manual/algorithms/bernstein-vazirani.md) |
| Fourier and arithmetic | [QFT](docs/manual/algorithms/qft.md), [Fourier addition](docs/manual/algorithms/fourier-addition.md), [order finding](docs/manual/algorithms/order-finding.md) |
| Hamiltonian evolution | [Trotter](docs/manual/algorithms/trotter.md), [QSP phase synthesis](docs/manual/algorithms/qsp-phase-synthesis.md), [QSVT HamSim](docs/manual/algorithms/qsvt-hamiltonian-simulation.md) |
| Estimation and testing | [QPE](docs/manual/algorithms/qpe.md), [Hadamard test](docs/manual/algorithms/hadamard-test.md), [Swap test](docs/manual/algorithms/swap-test.md) |
| Quantum linear systems | [Costa walk](docs/manual/algorithms/costa-walk.md), [CKS](docs/manual/algorithms/cks.md), [VTAA-CKS](docs/manual/algorithms/vtaa-cks.md) |
| QODE / QPDE | [Schrödingerization](docs/manual/algorithms/schrodingerization.md), [LCHS](docs/manual/algorithms/lchs.md), [Carleman](docs/manual/algorithms/carleman.md) |
| Variational and optimization | [VQE](docs/manual/algorithms/vqe.md), [QAOA MaxCut](docs/manual/algorithms/qaoa-maxcut.md), [DQI](docs/manual/algorithms/dqi.md) |
| Quantum walks | [coined walk](docs/manual/algorithms/coined-cycle-walk.md), [Szegedy](docs/manual/algorithms/szegedy-walk.md), [MNRS](docs/manual/algorithms/mnrs-search.md) |
| Quantum machine learning | [QPCA](docs/manual/algorithms/qpca.md), [QCNN](docs/manual/algorithms/qcnn.md), KP recommendation ([API](docs/api/algorithms/qml/recommendation.rst); manual page pending) |
| Data loading and input models | [state preparation](docs/manual/algorithms/state-preparation.md), [XOR database](docs/manual/algorithms/xor-database.md), [Select-Swap QROM](docs/manual/algorithms/select-swap.md) |
| Quantum error correction | [repetition codes](docs/manual/algorithms/repetition-codes.md) |

## Input models and operators

Algorithms consume inputs through five composable access-model families:

- Operator wrappers and basic composition (`identity`/`product`/`scale`/LCU): [manual](docs/manual/operators.md), [API](docs/api/algorithms/input_model/operators.rst); block-encoding algebra in [block_encoding](docs/api/algorithms/input_model/block_encoding.rst).
- Oracle paradigms (XorDatabase, StatePreparation, StateOracle, SparseAccess): [API](docs/api/algorithms/input_model/oracles.rst).
- Quantum data structures QVector/QMatrix: [manual](docs/manual/qdata.md), [API](docs/api/algorithms/input_model/qdata.rst).
- Density matrices and Gibbs states, spectral and low-rank decompositions: [density](docs/api/algorithms/input_model/density.rst), [spectral](docs/api/algorithms/input_model/spectral.rst), [lowrank](docs/api/algorithms/input_model/lowrank.rst).

## Applications

- **QFVM** (quantum fluid solving): [manual](docs/manual/qfvm.md), [API](docs/api/applications/qfvm.rst), input-model review [spec](docs/reference/qfvm-input-models.md).
- **QHAM** (PDE → HAM → QHAM pipeline): [manual](docs/manual/qham.md), [derivation spec](docs/reference/qham-derivation.md), [API](docs/api/applications/qham/linearization.rst), [tutorial](docs/tutorials/qham.md).
- **Roe matrix elements**: [roe](docs/api/applications/roe.rst), [roe_formulas](docs/api/applications/roe_formulas.rst).
- **Case catalog**: 22 reference workloads ([catalog](docs/api/applications/catalog.rst)) and the gallery generator ([gallery](docs/api/applications/gallery.rst), [tutorial](docs/tutorials/gallery.md)).

## Documentation map

Pick a reading path by role:

- **Getting started**: [first program](docs/tutorials/first-program.md) → [core concepts](docs/manual/concepts.md) → [replacing oracles](docs/tutorials/oracle-binding.md).
- **Algorithm research**: [algorithm-research tutorial](docs/tutorials/algorithm-research.md) → [algorithm index](docs/manual/algorithms/index.md) → [validation coverage matrix](docs/development/validation-coverage.md); applicability boundaries in [validation and scope](docs/manual/limits.md).
- **Differential-equation applications**: [tutorial](docs/tutorials/differential-equations.md) → [manual](docs/manual/differential-equations.md) → [QHAM manual](docs/manual/qham.md).
- **Backend engineering**: [backend manual](docs/manual/backends.md) → [compatibility review](docs/reference/backend-compatibility.md) → [infrastructure API](docs/api/infrastructure/index.rst).
- **Language and specs**: [RIR spec](docs/reference/rir.md), [generation-layer boundary](docs/reference/language.md), [open IR](docs/reference/open-ir.md), [QRAM memory format](docs/reference/qram-memory.md), [math IR](docs/reference/math-ir.md).

## Install · build · check

```bash
uv sync --locked --extra dev --extra docs
uv run python tools/build_docs.py --lang all
```

This builds both language trees (warnings are errors): English to
`out/docs/en/html` and Chinese to `out/docs/zh/html`, plus doctest runs. Open
`out/docs/en/html/index.html` to browse the generated site. Algorithm
gallery and engineering checks:

```bash
uv run python examples/algorithm_gallery.py
uv run python tools/check_project.py --docs
```

The algorithm gallery generates RIR, OriginIR-ext, and readout notes for 22
small examples. The full native check requires a separate environment with
`pysparq` and `uniqc` installed:

```bash
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

The language core depends only on PyYAML for text serialization. Optional
backends are imported at their execution entry points; generated artifacts,
environments, and build files are never committed. See
[CONTRIBUTING.md](CONTRIBUTING.md) for the development workflow and the
[documentation style guide](docs/development/writing-docs.md); source-tree
classification and migration notes are in [architecture](docs/manual/architecture.md)
and [import paths](docs/manual/compatibility.md). Documentation is also
maintained in Chinese: see [README.zh-CN.md](README.zh-CN.md).

Advanced algorithms such as QLS/QODE/QHAM still have pending verification
items on numerical accuracy, success channels, or convergence (per-item status
in the [validation coverage matrix](docs/development/validation-coverage.md)).
A general-purpose QSP-HamSim kernel is still to be provided; the current
modular multiplication uses finite-size permutation synthesis, and classical
optimizers for VQE/QAOA are chosen by the application.
