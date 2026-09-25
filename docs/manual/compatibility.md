# Import path migration

**English** · [简体中文](../zh/manual/compatibility.html)

Release 0.8 moved the implementations into directories organized by
responsibility. The common exports of the root package remain available; old
module paths are forwarded through a centralized compatibility table and
reference the same implementations. New code should use the canonical paths.

| Old path | Canonical path |
|---|---|
| `oracq.ir` | [`oracq.infrastructure.ir`](../api/infrastructure/ir.rst) |
| `oracq.builder` | [`oracq.infrastructure.builder`](../api/infrastructure/builder.rst) |
| `oracq.validation` | [`oracq.infrastructure.validation`](../api/infrastructure/validation.rst) |
| `oracq.serialization` | [`oracq.infrastructure.serialization`](../api/infrastructure/serialization.rst) |
| `oracq.linking` | [`oracq.infrastructure.linking`](../api/infrastructure/linking.rst) |
| `oracq.execution` | [`oracq.infrastructure.execution`](../api/infrastructure/execution.rst) |
| `oracq.native` | [`oracq.infrastructure.native`](../api/infrastructure/native.rst) |
| `oracq.layout` | [`oracq.infrastructure.layout`](../api/infrastructure/layout.rst) |
| `oracq.readout` | [`oracq.infrastructure.readout`](../api/infrastructure/readout.rst) |
| `oracq.backends.originir` | [`oracq.infrastructure.backends.originir`](../api/infrastructure/backends/originir.rst) |
| `oracq.backends` | `oracq.infrastructure.backends` |
| `oracq.backends.pysparq` | [`oracq.infrastructure.backends.pysparq`](../api/infrastructure/backends/pysparq.rst) |
| `oracq.backends.basis` | [`oracq.infrastructure.backends.basis`](../api/infrastructure/backends/basis.rst) |
| `oracq.mathfunc.graph` | [`oracq.infrastructure.mathfunc.graph`](../api/infrastructure/mathfunc/graph.rst) |
| `oracq.mathfunc.frontend` | [`oracq.infrastructure.mathfunc.frontend`](../api/infrastructure/mathfunc/frontend.rst) |
| `oracq.mathfunc.numeric` | [`oracq.infrastructure.mathfunc.numeric`](../api/infrastructure/mathfunc/numeric.rst) |
| `oracq.mathfunc.lowering` | [`oracq.infrastructure.mathfunc.lowering`](../api/infrastructure/mathfunc/lowering.rst) |
| `oracq.mathfunc` | [`oracq.infrastructure.mathfunc`](../api/infrastructure/mathfunc.rst) |
| `oracq.mathfunc.roe_formulas` | [`oracq.applications.roe_formulas`](../api/applications/roe_formulas.rst) |
| `oracq.library` | [`oracq.algorithms.input_model.operators`](../api/algorithms/input_model/operators.rst) |
| `oracq.combinators` | [`oracq.algorithms.input_model.block_encoding`](../api/algorithms/input_model/block_encoding.rst) |
| `oracq.contracts` | [`oracq.algorithms.input_model.contracts`](../api/algorithms/input_model/contracts.rst) |
| `oracq.oracles` | [`oracq.algorithms.input_model.oracles`](../api/algorithms/input_model/oracles.rst) |
| `oracq.arithmetic` | [`oracq.algorithms.common.arithmetic`](../api/algorithms/common/arithmetic.rst) |
| `oracq.sparse_models` | [`oracq.algorithms.input_model.sparse`](../api/algorithms/input_model/sparse.rst) |
| `oracq.access` | [`oracq.algorithms.input_model.sparse`](../api/algorithms/input_model/sparse.rst) |
| `oracq.qlss` | [`oracq.algorithms.qlss.qlss`](../api/algorithms/qlss/qlss.rst) |
| `oracq.algorithms.costa` | [`oracq.algorithms.qlss.qlss`](../api/algorithms/qlss/qlss.rst) |
| `oracq.algorithms.cks` | [`oracq.algorithms.qlss.qlss`](../api/algorithms/qlss/qlss.rst) |
| `oracq.qode` | [`oracq.algorithms.qode.ode`](../api/algorithms/qode/ode.rst) |
| `oracq.qfvm` | [`oracq.applications.qfvm`](../api/applications/qfvm.rst) |
| `oracq.qfvm_sparse` | [`oracq.applications.qfvm`](../api/applications/qfvm.rst) |
| `oracq.flow_data` | [`oracq.applications.flow_data`](../api/applications/flow_data.rst) |
| `oracq.roe` | [`oracq.applications.roe`](../api/applications/roe.rst) |
| `oracq.applications` | `oracq.applications.legacy` |
| `oracq.workloads` | [`oracq.applications.catalog`](../api/applications/catalog.rst) |
| `oracq.qham.pde` | [`oracq.applications.qham.pde`](../api/applications/qham/pde.rst) |
| `oracq.qham.linearization` | [`oracq.applications.qham.linearization`](../api/applications/qham/linearization.rst) |
| `oracq.qham.reference` | [`oracq.applications.qham.reference`](../api/applications/qham/reference.rst) |
| `oracq.qham.quantum` | [`oracq.algorithms.input_model.qham`](../api/algorithms/input_model/qham.rst) |
| `oracq.qham` | `oracq.applications.qham` |
| `oracq.qham.examples` | [`oracq.applications.qham.examples`](../api/applications/qham/examples.rst) |
| `oracq.qham.report` | [`oracq.applications.qham.report`](../api/applications/qham/report.rst) |
| `oracq.qham.stencils` | [`oracq.applications.qham.stencils`](../api/applications/qham/stencils.rst) |

Package paths without a standalone API page (such as
`oracq.infrastructure.backends` and `oracq.applications.qham`) stay as-is;
open the API pages of their submodules to inspect them.

## Split-up legacy entry points

The query, search, QFT, QPE, and matrix-transform entry points of
`algorithms.elementary` moved into `oracle_algorithms`, `search`, `fourier`,
`estimation`, and `transforms`, respectively.

The concrete methods of `algorithms.differential` moved into `lchs`,
`schrodingerization`, `cbmd`, and `carleman`; the general entry point is
`ode`, and the PDE adapter is `pde`.

The state composition utilities of `algorithms.solvers` moved into
`state_preparation`, the Euler history into `ode`, and Trotter into
`hamiltonian`. The old factories remain available through the compatibility
entries.

## File formats

The RIR version is still 0.3, and the instruction set and execution semantics
are unchanged; the `binding_captures` field newly added to the linker holds
provenance information inside the existing scalar-attribute format, and the
schema and semantic checks now state its constraints as well. After the module
locations and generator compositions of the Python classes changed, module
hash names may change; do not use hash names as part of an application
protocol.

## Algorithm interface naming and the resource ledger

{obj}`AlgorithmContract <oracq.algorithms.input_model.contracts.AlgorithmContract>`,
{obj}`QLSSSolver <oracq.algorithms.qlss.qlss.QLSSSolver>`, and
{obj}`QODESolver <oracq.algorithms.qode.ode.QODESolver>` are the new canonical
names;
{obj}`ProtocolContract <oracq.algorithms.input_model.contracts.ProtocolContract>`,
{obj}`QLSSProtocol <oracq.algorithms.qlss.qlss.QLSSProtocol>`, and
{obj}`QODEProtocol <oracq.algorithms.qode.ode.QODEProtocol>` are kept as
same-type aliases. Problem builders accept provider protocols, and the
constructed fields still hold concrete role views.

The `rotations` of resource estimation changed from a list to a compact
read-only sequence. The previous length, indexing, and iteration reads keep
working; code that mutated the list directly should read `counts` instead or
use the report. In open analyses `qubits` is `None`; known lower bounds are in
`qubits_lower_bound`.
