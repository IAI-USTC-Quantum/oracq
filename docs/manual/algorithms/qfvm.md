# Quantum Finite Volume Method, QFVM

**English** · <a href="../../zh/manual/algorithms/qfvm.html">简体中文</a>

> Category C6 · Module [`oracq.applications.qfvm`](../../api/applications/qfvm.rst) · Stage V1

## Overview

QFVM constructs all the inputs needed for quantum linear-system solving (QLSS) from classical flow-field data: the matrix-element oracle, the sparse position oracle, and the right-hand-side residual state preparation, then hands the same problem to a replaceable QLSS (either the CKS base route declaring sparse input, or the Costa route declaring block-encoding input). The interface division follows the QFVM paper ([arXiv:2102.03557](https://arxiv.org/html/2102.03557v1); see the input-model review). The current implementation targets the Euler equations on a one-dimensional periodic grid: three conserved quantities (density, momentum, energy) plus the frozen-Roe Jacobian; it does not cover the full set of grids, boundaries, and physical models of the original paper.

The linear system takes the Hermitian dilation in physical coordinates `D = [[0, M], [M.T, 0]] + padding_value·I_pad` (M is the Roe block), and the physical-coordinate block is selected on the output side after solving; this dilation serves the linear solve and must not be taken as a QODE evolution equivalent to the original generator.

## Interface and input model

```python
roe_qfvm_inputs(*, cell_width=2, fmt=None, angle_width=8, prefix="RoeQfvm")
qfvm_sparse_access(inputs, *, padding_value=1.0, **entry_options)
roe_qfvm_block_encoding(inputs, *, amax=8.0, padding_value=1.0, **entry_options)
roe_qfvm_problem(inputs, *, spectrum, rhs_norm=None, amax=8.0, padding_value=1.0, **entry_options)
roe_qfvm_step(inputs, qlss, *, spectrum=None, rhs_norm=None, **options)
bind_qfvm(program, inputs)
qfvm_memories(inputs, flow, *, amax=8.0)
```

API entry points: {obj}`roe_qfvm_inputs <oracq.applications.qfvm.roe_qfvm_inputs>`, {obj}`qfvm_sparse_access <oracq.applications.qfvm.qfvm_sparse_access>`, {obj}`roe_qfvm_block_encoding <oracq.applications.qfvm.roe_qfvm_block_encoding>`, {obj}`roe_qfvm_problem <oracq.applications.qfvm.roe_qfvm_problem>`

- {obj}`roe_qfvm_inputs <oracq.applications.qfvm.roe_qfvm_inputs>` returns {obj}`RoeQfvmInputs <oracq.applications.qfvm.RoeQfvmInputs>`: six abstract databases (rho / momentum / energy / geometry / theta / residual values) plus an abstract right-hand-side state preparation. `cell_width` ranges 2..25 (at least four cells) and `angle_width` ranges 1..64; the matrix dimension is `2**width`, where `width = cell_width + 3` (3 conserved quantities + 1 padding component).
- Input-model layering: the flow-field, geometry, and angle data are QRAM (bound via {obj}`qram_database <oracq.algorithms.input_model.oracles.qram_database>`); the matrix-access artifact is SO ({obj}`qfvm_sparse_access <oracq.applications.qfvm.qfvm_sparse_access>` returns {obj}`SparseAccess <oracq.algorithms.input_model.oracles.SparseAccess>`: position + element oracles, sparsity 9); the right-hand side is SP (the signed residual tree of {obj}`rhs_qram_preparation <oracq.applications.qfvm.rhs_qram_preparation>`).
- {obj}`roe_qfvm_problem <oracq.applications.qfvm.roe_qfvm_problem>` assembles a {obj}`LinearSystem <oracq.algorithms.qlss.qlss.LinearSystem>`: `spectrum` must be an explicit {obj}`SpectralPromise <oracq.algorithms.qlss.qlss.SpectralPromise>` (merged with padding_value by taking the wider bound); `amax`, the spectral bounds, and the matrix properties are caller declarations that must hold for the actually quantized matrix; data assumptions are recorded in `data_assumptions` (coherent XOR QRAM, pinned snapshots, raw field quantities rather than matrix tables, etc.).
- {obj}`roe_qfvm_step <oracq.applications.qfvm.roe_qfvm_step>` requires a {obj}`QLSSProtocol <oracq.algorithms.qlss.qlss.QLSSProtocol>` that declares its input_model (bare callables are rejected) and returns `qlss.solve(problem)`; {obj}`bind_qfvm <oracq.applications.qfvm.bind_qfvm>` replaces the abstract slots with QRAM database bindings, and {obj}`qfvm_memories <oracq.applications.qfvm.qfvm_memories>` assembles the memory tables from {obj}`RoeFlowData <oracq.applications.flow_data.RoeFlowData>` snapshots.

`RoeQfvmInputs` has the attributes `width` (= cell_width + 3) and `geometry_width` (= width + 4 + cell_width + 7, packing the seven segments neighbor / reverse / source / rowvar / colvar / band / valid).

## Implementation notes

Matrix elements are not pre-stored as a table: {obj}`roe_entry <oracq.applications.qfvm.roe_entry>` queries the conserved-quantity databases of the source cell and its left/right neighbors, calls the reversibly compiled {obj}`roe_face <oracq.applications.roe.roe_face>` on the two interfaces (the pure function {obj}`frozen_roe_face <oracq.applications.roe_formulas.frozen_roe_face>` compiled to fixed-point arithmetic via mathfunc {obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>`, including sqrt / div / mul / select), selects the west / center / east flux difference by band, and the mass term is written only onto the component diagonal via the `component_equal` Boolean network; on arithmetic failure (status ≠ 0) the element totalizes to zero. The position oracle encodes only neighbors, slots, and component indices ({obj}`geometry_cells <oracq.applications.qfvm.geometry_cells>`; it contains no matrix values): the 9 structural slots per column (3 neighbors × 3 components) go through 9 coherent value transpositions ({obj}`value_transposition <oracq.algorithms.input_model.sparse.value_transposition>`) to produce the in-place index permutation CKS needs; the padding component (variable 3) is the diagonal `padding_value`, the remaining invalid slots return the zero element, and the boundaries are periodic. The right-hand side uses {obj}`ptheta_cells <oracq.applications.qfvm.ptheta_cells>` to quantize residual values into rotation angles `2·arccos(min(1, |v|/amax))`; `rhs_qram_preparation` prepares amplitudes from a QRAM angle table and writes the signs via Z kickback from the sign database.

The classical-side Riemann computation and local flow-field updates are managed by `RoeFlowData` (`applications/flow_data.py`): a single-point update recomputes only the adjacent interfaces and the affected tree nodes; modifying a bank on the native backend still requires rematerialization, which is not the constant-time physical write assumed by the paper. The module attribute `correctness="pending"`: what is currently witnessed is the input adaptation and the interchangeability of the two QLSS routes; numerical-solve correctness is not certified. Full assumptions and boundaries are in the [QFVM/QLSS input-model review](../../reference/qfvm-input-models.md).

## Validation approach

Category C6 (composition skeleton; acceptance criteria in `../../development/validation-plan.md` §2: component contracts satisfied + end-to-end small-instance semantic correctness + bind invariance). `validation-coverage.md` registers `qfvm.py`/`qham/` on the application-layer row (C6, V1, no gaps); QFVM's concrete witnesses sit in the evidence column of the ODE/PDE-solver row:

- Structure and numerical: `tests/core/test_differential.py:DifferentialStructureTests.test_qfvm_uses_raw_flow_and_modular_arithmetic` — with {obj}`FixedFormat(4, 1) <oracq.algorithms.common.arithmetic.FixedFormat>` and `amax=4` inputs, {obj}`roe_qfvm_block_encoding <oracq.applications.qfvm.roe_qfvm_block_encoding>` gives α = 36; the BE program contains exactly 4 unresolved slots, and after closing the resources are the four databases rho / momentum / energy / geometry; the program contains sqrt / div / mul / select arithmetic modules; the geometry table has size exactly `2**(width+4)`; the sparse encoding is annotated with sparsity 9 and the self-adjoint extension.
- Classical side: the same family's `test_local_riemann_updates_only_neighbor_faces_and_tree` — a single-point update of a 16-cell flow field recomputes only interfaces (4, 5) and cells (4, 5, 6), updates fewer rhs angle-tree entries than a full update, and materializes no matrix bank.
- Binding: inside `test_qfvm_uses_raw_flow_and_modular_arithmetic`, after `bind_qfvm` replaces all abstract slots there are no unresolved declarations.

## Known gaps and planned stages

Consistent with the application-layer row of `validation-coverage.md`: no registered gaps, stage V1. The `correctness="pending"` module attribute is a record of the implementation boundary (numerical-solve correctness certification is outside the current validation scope), not a gap on the coverage matrix.

## Numerical validation

Paper-grade numerical validation executed on 2026-09-16 by `tests/verification/verify_qham_qfvm.py` (real backends: the PySparQ native RIR interpreter and the UniQC full-amplitude state vector; no mocks, no skips). This machine has no C++ compiler, so {obj}`arithmetic_native_registry <oracq.algorithms.common.arithmetic.arithmetic_native_registry>` is unavailable and the Roe arithmetic executes through gate-level expansion on PySparQ RIR (about 1.6×10⁷ expanded gates per single-entry circuit); entry-level validation therefore uses superposition sampling instead of full-matrix enumeration.

### Experiment design

- **Fixed-point format choice**: `roe_face`'s compile_function requires the row/column inputs of {obj}`Index(2) <oracq.infrastructure.mathfunc.graph.Index>` to be representable, and the smallest usable format is `FixedFormat(5,2)`; `entropy_delta=0.5` makes 2δ, δ², and δ all exactly representable in that format. Note: if the fraction bits cannot represent 2δ (e.g. fmt=(4,1) with the default δ=0.125, where 2δ=0.25 truncates to raw 0), the entropy-correction branch divides by zero and the entry silently zeroes out per the documented totalize behavior — the integration test previously covered only the padding diagonal, not nonzero Roe entries.
- **Roe face flux**: the compiled `roe_face` is cross-checked **bit for bit** against an in-script independent fixed-point simulation (re-implemented per the documented toward_zero/modular_wrap semantics of {obj}`fixed_arithmetic <oracq.algorithms.common.arithmetic.fixed_arithmetic>`) over 16+8+8=32 superposition branches (all row/col pairs, a density sweep, a momentum sweep); the method error (against float64 `roe_formulas`) is reported separately.
- **Sparse entries/positions**: the entry oracle is checked bit by bit against the fixed-point simulation matrix under an 8-branch superposition on column 5 (three-band nonzeros + zero structure + padding zero); the padding diagonal is a separate case; the position oracle verifies, under a superposition of all 32 columns, the complete permutation and the 9-slot mapping (against independent geometry semantics).
- **RHS and classical identities**: the RHS residual-state preparation (QRAM angle tree + sign kickback) is checked against an independent residual computation; flow_data's F*(L,R)=left·U_L+right·U_R is checked against the matrix-inversion implementation {obj}`riemann_flux <oracq.applications.flow_data.riemann_flux>`; M·u−mass·u=−residual (implicit-FVM sign convention); the local update patch agrees with a full recomputation; the ptheta angle table is pointwise true. Uniform four-cell flow field ρ=(1,1,1.25,1.25), m=(0,0.25,0,−0.25), E=(1,1,1,1).

### Key metrics

| Case | Scale | Backend path | Value |
|---|---|---|---|
| roe_face bit-level cross-check | 32 branches, fmt=(5,2) | rir-pysparq | raw/status bit-identical (32/32); method error 0.427 |
| entry-matrix sampling | column 5, 8 branches | rir-pysparq | implementation error 0.0 (bit-identical); method error 0.481; zero structure correct |
| padding diagonal | coordinates (7,7)/(23,7) | rir-pysparq | diagonal raw=4 (1.0), 0 outside the dilation block |
| position oracle | 32 columns × 32 inputs | rir-pysparq | complete permutation per column; 288 slot mappings 0 mismatches; geometry table 0 pointwise mismatches; work cleaned |
| RHS preparation | 16 addresses, angle_width 8 | rir-pysparq | amplitude error 1.11e-16, signs 0 mismatches, total probability 1.0; angle tree/sign bank 0 pointwise mismatches |
| flux identity | 4 sampled states | classical independent | F* block form vs inversion implementation 3.33e-16; F*(U,U)=F(U) 2.22e-16 |
| single-step update | 4 cells | classical independent | M·u−mass·u=−residual 2.22e-16; residual bank 0 mismatches; local patch consistent with full recomputation |

The method error is the inherent quantization error of the 5-bit fixed-point pipeline (about 1–2 quanta), not an implementation defect: the implementation side (against the fixed-point-semantics simulation) is bit-exact.

### Reproduce

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_qham_qfvm.py
```

Artifact: `out/verification/qham_qfvm.json` (full metrics and criteria of the 15 cases).

## Related links

- Source: `src/oracq/applications/qfvm.py` (Roe arithmetic `applications/roe.py`, formulas `applications/roe_formulas.py`, classical data `applications/flow_data.py`)
- User guide: [QFVM input models and solver replacement](../qfvm.md)
- API reference: [QFVM application](../../api/applications/qfvm.rst)
- Same-group page: [General QHAM](qham.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)
