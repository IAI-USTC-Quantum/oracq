# QFVM input models and solver replacement

**English** · <a href="../../zh/manual/qfvm.html">简体中文</a>

The QFVM application builds a linear system from flow-field data and hands it to a replaceable QLSS. The current implementation targets the periodic one-dimensional Euler equations with three conserved quantities and the frozen-Roe Jacobian. It does not cover the full set of grids, boundaries, and physical models of the original paper. For the end-to-end workflow see the [tutorial: scientific computing workflows](../tutorials/scientific-workflows.md).

## Data and quantum access

The density, momentum, and energy of the flow field are kept in QRAM data. The geometry table records neighbors, slots, and component indices; no full Jacobian matrix is pre-stored.

The quantum matrix-element oracle queries the relevant cells, invokes reversible Roe arithmetic, and selects elements by row/column index. The position oracle gives the structural positions of each column through the in-place index permutation used by <a href="../../zh/manual/algorithms/cks.html">CKS</a>. The RHS preparation uses the angle table provided by the residual data structure.

Classical Riemann computation and local flow-field updates are managed by {obj}`RoeFlowData <oracq.applications.flow_data.RoeFlowData>`. A logical QRAM patch can update locally; the current native backend materialization may still rebuild modified banks, so the two must not be equated.

## Constructing the problem

```python
from oracq import SpectralPromise
from oracq.applications.qfvm import roe_qfvm_inputs, roe_qfvm_problem
from oracq.algorithms.qlss.qlss import CKSConfig, CostaConfig, make_cks_qlss, make_costa_qlss

inputs = roe_qfvm_inputs()
problem = roe_qfvm_problem(
    inputs,
    spectrum=SpectralPromise(norm_upper=10.0, sigma_min_lower=0.5),
    amax=4.0,
)
cks = make_cks_qlss(CKSConfig(order=2))
costa = make_costa_qlss(CostaConfig(steps=1))
cks.check(problem).require()
costa.check(problem).require()
```

`amax`, the spectral bounds, and the matrix properties are caller declarations and must apply to the actually quantized matrix. The values in the example illustrate the interface and do not replace a spectral analysis of a real flow field.

## Replacing the QLSS

The <a href="../../zh/manual/algorithms/cks.html">CKS</a> entry consumes sparse problems, and the <a href="../../zh/manual/algorithms/costa-walk.html">Costa</a> entry requests a BE. The current real-symmetric sparse adaptation constructs the corresponding BE via `T†ST` and records alpha. There is no general adaptation that recovers efficient sparse access from an arbitrary BE in reverse.

QFVM's non-symmetric physical matrix uses an explicit Hermitian dilation and selects the physical coordinate at output. This dilation serves the linear solve; it cannot be taken directly as a QODE evolution equivalent to the original generator.

The two QLSS routes can keep the same problem inputs and physical output interpretation, but their ancilla widths and internal modules differ. After switching implementations, the higher-level layout should be regenerated. Late binding into existing open slots is appropriate only when the ABI and alpha are unchanged.

## Output and norms

{obj}`SolveResult <oracq.algorithms.qlss.qlss.SolveResult>` returns the state oracle, matrix-norm probes, input alphas, spectral declarations, and adaptation records. It keeps the solve success probability separate from the conditional probability of the independent matrix probe, avoiding the misuse of filtering success rates as the solution-vector norm.

The solution-state direction, success channel, numerical error, and the outer CFD loop still need case-by-case validation. For the fuller mathematical assumptions and the comparison with the original paper see the [QFVM/QLSS input-model review](../reference/qfvm-input-models.md).

## Line-by-line implementation walkthrough (including the QRAM data structures)

This section takes the QFVM implementation completely apart: first the contracts of every QRAM bank (who writes, who reads, what is stored), then the geometry table and the residual tree field by field, then the three quantum circuits segment by segment (matrix-element oracle, sparse position oracle, RHS preparation), and finally a line-by-line checklist of the end-to-end assembly. The sources are `src/oracq/applications/qfvm.py` and `src/oracq/applications/flow_data.py`.

### QR1. QRAM bank overview

What {obj}`roe_qfvm_inputs <oracq.applications.qfvm.roe_qfvm_inputs>` declares are **abstract slots** (open modules); the data appears only at binding and execution time. The full contract of each bank after binding:

| bank | address width → data width | classical writer | quantum reader | contents |
|---|---|---|---|---|
| `rho` / `momentum` / `energy` | cell_width → fmt.width | `RoeFlowData.update` | {obj}`roe_entry <oracq.applications.qfvm.roe_entry>`: 3 neighbors × 3 fields per element = 9 queries | fixed-point encoding of the conserved quantities per cell |
| `geometry` | width+4 → geometry_width | {obj}`geometry_cells() <oracq.applications.qfvm.geometry_cells>` generated statically (independent of the flow values) | position oracle 9 times; entry oracle 9 times | the seven-segment packed pure-geometry word (QR2) |
| `rhs_values` | cell_width+2 → fmt.width | `RoeFlowData.update` | — (data plane of the classical norm tree) | quantized residual value per (cell, component) |
| `rhs_sign` | cell_width+2 → 1 | `RoeFlowData.update` | RHS preparation queries twice (write + un-query) | residual sign bit |
| `rhs_angles` | cell_width+2 → angle_width | `RoeFlowData.update` (tree cache) | {obj}`qram_state_prep <oracq.algorithms.input_model.oracles.qram_state_prep>` twice per level, 2·(cell_width+2) in total | rotation angle words of the internal nodes of the squared-norm tree |
| `theta` | fmt.width → angle_width | {obj}`ptheta_cells() <oracq.applications.qfvm.ptheta_cells>` generated statically | — (reserved data plane, see below) | angle words for θ=2·acos(min(1,|v|/amax)) |

Three key points:

1. **Matrix elements are never pre-stored**. The QRAM holds only raw field quantities and geometry; matrix elements are computed on the fly by reversible Roe arithmetic at query time (QR4).
2. **`theta` is reserved**. It is the conversion table "residual-value fixed-point word → rotation angle word", materialized together with {obj}`qfvm_memories <oracq.applications.qfvm.qfvm_memories>`; neither the CKS nor the Costa route has a circuit querying it today, and the documentation records this honestly instead of presenting it as a wired-in oracle.
3. **The memory snapshot must be pinned**. P, P†, controlled invocations, and norm probes within one coherent computation must read the same version of the data; this is written into the `data_assumptions` attribute of {obj}`roe_qfvm_problem <oracq.applications.qfvm.roe_qfvm_problem>`, and the language does not prove it on the hardware's behalf.

### QR2. The geometry table `geometry_cells`, field by field

First fix the coordinate layout: the full matrix coordinate is `half·2^width + 4·cell + var` — the lowest 2 bits are the component `var` (0..3, with 3 the padding component), the middle cell_width bits are the cell number, and the top bit `half` selects the upper/lower block of the Hermitian dilation `D=[[0,M],[Mᵀ,0]]`. The QRAM address = `row + slot·2^width`, where row packs (half, cell, var) and slot packs (band, other_var).

Reading `geometry_cells` line by line:

```python
cell, var, half = (row >> 2) % n, row % 4, row >> (cw + 2)
```

Unpacks the row coordinate: the low 2 bits are `var`, the middle segment `cell`, and the top bit `half` (upper block = M, lower block = Mᵀ).

```python
band, other_var = divmod(slot, 3)
valid = int(slot < 9 and var < 3)
```

Of the 16 slot values only the first 9 are structural slots (band∈{west,center,east} × other_var∈{0,1,2}); the padding component (var=3) and invalid slots are marked `valid=0`.

```python
neighbor_cell = (cell + band - 1) % n
neighbor = other_var + 4 * neighbor_cell + ((1 - half) << (cw + 2))
```

The partner coordinate: this row couples, through slot (band, other_var), to (neighbor_cell, other_var) in **the other half of the space**. Upper-block rows point at lower-block columns and lower-block rows point at upper-block columns — exactly the transpose coupling structure of D; `% n` implements the periodic boundary.

```python
reverse = (2 - band) * 3 + var if valid else 0
source = cell if half == 0 else neighbor_cell
rowvar, colvar, source_band = (
    (var, other_var, band) if half == 0 else (other_var, var, 2 - band)
)
```

`reverse` is the symmetric slot of the same matrix element; `source` is the source cell to use when computing this element (lower-block rows take the neighbor); the last line converts a lower-block entry into the same physical element of the upper block (row/column swapped, band replaced by 2−band).

```python
values = ((neighbor, w), (reverse, 4), (source, cw), (rowvar, 2),
          (colvar, 2), (source_band, 2), (valid, 1))
```

The seven segments are packed from low to high bits into one `geometry_width`-bit word: partner coordinate | symmetric slot | source cell | row component | column component | source band | valid bit. Both the position oracle and the entry oracle query this same word and slice it with the same segment widths from `_geometry_refs`.

### QR3. The QRAM signed residual tree for the RHS

**Classical side** (`RoeFlowData.update`) maintains a 1-based heap-array tree (`8n` entries) whose leaves store squared residual components:

```python
node = 4 * self.n + address          # leaves: heap 1..8n, leaves start at 4n+1 (address=4*cell+j)
self.tree[node] = value * value      # squared residual
while node > 1:                      # collect affected nodes up the ancestors
    node //= 2
    ancestors.add(node)
```

```python
self.tree[node] = self.tree[2*node] + self.tree[2*node + 1]      # recompute subtree sums bottom-up
angle = 2 * math.acos(math.sqrt(self.tree[2*node] / self.tree[node]))
changes["rhs_angles"][node - 1] = round(angle * (1 << aw) / (2 * math.pi)) % (1 << aw)
```

The angle formula comes from the RY half-angle convention: `RY(θ)|0> = cos(θ/2)|0> + sin(θ/2)|1>`; to make the left/right subtree amplitude ratio `sqrt(S_left/S_node)`, take `θ = 2·acos(sqrt(S_left/S_node))`. The angle is quantized into an `angle_width`-bit angle word and written to the `rhs_angles` bank (address = node−1). A single-point flow update touches only adjacent interfaces and O(log) tree nodes — this is where the "local update" claim comes from.

**Quantum-side amplitudes** (`qram_state_prep`, prepared level by level):

```python
for depth in range(width):
    bit = width - depth - 1                     # target bit decided at this level (from the top down)
    if depth:
        b.xor(b["target"][bit + 1:], addr[:depth])   # copy the settled prefix into the address
    offset = (1 << depth) - 1
    b.add_const(addr.reinterpret("uint"), offset)    # address = prefix + tree node offset
    b.qram("angles", addr, angle)                    # query this internal node's angle word
    for k in range(angle_width):
        with b.control(angle[k]):
            b.ry(b["target"][bit], 2 * math.pi * (1 << k) / (1 << angle_width))
    b.qram("angles", addr, angle)                    # un-query to clean angle
    b.add_const(addr.reinterpret("uint"), (-offset) % (1 << address_width))
```

Line-by-line highlights: the tree node addressed at depth `depth` is `(1<<depth)−1+prefix` (consistent with the table addressing of {obj}`qram_state_angles <oracq.algorithms.input_model.oracles.qram_state_angles>`); an `angle_width`-bit angle word need not be compiled into `2^angle_width` distinct rotation gates — bitwise weighting `2π·2^k/2^aw` synthesizes any angle with `angle_width` controlled RYs; one query and one un-query per level gives `2·width` total QRAM queries (versus O(2^n) queries for leaf-by-leaf preparation).

**Quantum-side signs** ({obj}`rhs_qram_preparation <oracq.applications.qfvm.rhs_qram_preparation>`): the amplitude tree only produces non-negative amplitudes; the signs are written through **phase kickback** using a separate 1-bit bank:

```python
invoke(b, prep.operation, "prep", target=b["target"][:n], work=b["work"][: n + aw])  # amplitudes
flag = b["work"][n + aw :]
invoke(b, sign.operation, "sign", address=b["target"][:n], data=flag)  # XOR the sign into flag
b.z(flag)                                                              # branches with flag=1 pick up a factor of −1
invoke(b, sign.operation, "sign", address=b["target"][:n], data=flag)  # un-query to clean flag
```

{obj}`Load <oracq.infrastructure.ir.Load>` is self-inverse, so the net effect of "write sign → Z → restore sign" is only a branch phase: negative residual components pick up a π phase, and flag and work are fully cleaned (`zero_input=True`).

### QR4. The matrix-element oracle `roe_entry`, segment by segment

Registers: `source` (cell), `row`/`col` (components), `band` (west/center/east), `value` (fixed-point output), `status` (arithmetic status). The resources are the three conserved-quantity banks.

```python
for offset in (-1, 0, 1):
    addr = g.local(cw)
    b.xor(b["source"], addr)                                   # addr = source
    b.add_const(addr.reinterpret("uint"), offset % (1 << cw))  # addr += offset (periodic wrap-around)
    addresses.append(addr)
    for name, op in fields:
        word = g.local()
        invoke(b, op, name, address=addr, data=word)           # query rho/momentum/energy
        state.append(word)
```

The three conserved-quantity banks are queried for source and its left/right neighbors — the inputs of one element computation are 3 cells × 3 fields, 9 QRAM queries in total.

```python
face = roe_face(fmt=fmt, gamma=gamma, entropy_delta=entropy_delta)
for i in range(2):
    b.call(face, rho_l=..., m_l=..., e_l=..., rho_r=..., m_r=..., e_r=...,
           row=b["row"], col=b["col"], left=left, right=right, status=flag)
```

{obj}`roe_face <oracq.applications.roe.roe_face>` is the reversible fixed-point circuit (with sqrt/div/mul/select) that {obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>` compiles from the pure function {obj}`frozen_roe_face <oracq.applications.roe_formulas.frozen_roe_face>`; two invocations yield the left/right flux characteristic components of the left and right interfaces. When `status≠0` (overflow/division by zero), the element totalizes to zero per the documented convention.

```python
west = -faces[0][0] / dx
center = (faces[1][0] - faces[0][1]) / dx + g.choose(equal, mass, 0)
east = faces[1][1] / dx
value = g.choose3(b["band"], [west, center, east])
```

The west/center/east candidates are assembled from the finite-volume flux differences; the mass term is added to the center-block diagonal only when `row == col` (the `component_equal` Boolean network); finally one of three is selected by `band`.

### QR5. The sparse position oracle (CKS in-place permutation), segment by segment

CKS requires "given a column and a sparse rank, obtain the row index in place". The locator of {obj}`qfvm_sparse_access <oracq.applications.qfvm.qfvm_sparse_access>` builds the full permutation with 9 coherent value transpositions:

```python
for rank in range(9):
    ...  # slot set to rank; query geometry(column, slot) for the seven-segment word; neighbor = data[:width]
    padded = locator.local(...); locator.xor(locator["column"], padded)
    locator.add_const(padded.reinterpret("uint"), rank)
    with locator.control(locator["column"][:2], 3):      # padding component (var=3)
        locator.xor(data[:n], neighbor)
        locator.xor(padded, neighbor)                    # partner coordinate rewritten as column+rank (padding diagonal)
    for bit in range(n):
        if (rank >> bit) & 1:
            locator.x(left[bit])                         # left starts at rank
    for previous in range(rank):
        locator.call(value_transposition(n), index=left, a=lefts[previous], b=neighbors[previous])
```

The padding component has no real neighbor structure: its partner coordinate is rewritten as `column+rank`, so the 9 ranks map exactly onto 9 distinct padding coordinates (the padding diagonal). `left` starts from rank and performs a value transposition against each earlier `(left_j, neighbor_j)` — if rank is already occupied it is swapped to a not-yet-used slot, guaranteeing that the 9 slot mappings form a **complete permutation**.

```python
setup = tuple(locator._frames[0])
for left, right in zip(lefts, neighbors, strict=True):
    locator.call(value_transposition(n), index=locator["index"], a=left, b=right)
locator.emit(Adjoint(setup))
```

Nine value transpositions are applied to `index` in order (a word in `index` equal to `left_j` is replaced by `neighbor_j`), then {obj}`Adjoint(setup) <oracq.infrastructure.ir.Adjoint>` cleans all work bits. The inverse mapping is the adjoint of the same circuit — exactly the in-place semantics the CKS interface needs; the cost is sparsity-many reversible comparisons/swaps, not a free identity gate.

The **entry oracle** uses the same geometry table: each of the 9 slots goes through the compute/uncompute pair "query → match `row` with {obj}`compare_words <oracq.algorithms.input_model.sparse.compare_words>` → XOR-accumulate the selected seven-segment word → un-compare → un-query"; the segments then call QR4's `roe_entry` to obtain `(value, status)`:

```python
with b.control(valid):
    with b.control(status, 0):
        b.xor(value, b["data"])                    # XOR-write elements that are valid and arithmetically successful
with b.control(fuse(same, b["column"][:2]), 7):    # diagonal with padding column component (3)
    ...                                            # write padding_value
```

An `Adjoint(forward)` layer wraps both ends to keep everything reversible. The final result is {obj}`SparseAccess(location, entry, width, fmt.width, 9) <oracq.algorithms.input_model.oracles.SparseAccess>`.

### QR6. End-to-end assembly, line by line

```python
from oracq import FixedFormat, SpectralPromise
from oracq.applications.flow_data import RoeFlowData
from oracq.applications.qfvm import (
    bind_qfvm, qfvm_memories, roe_qfvm_inputs, roe_qfvm_problem,
)
from oracq.algorithms.qlss.qlss import CKSConfig, CostaConfig, make_cks_qlss, make_costa_qlss

# Fixed-point format: numerical validation shows (5,2) is the smallest format in which
# all roe_face constants are exactly representable.
fmt = FixedFormat(5, 2)
# 1. Declare the input model: six abstract databases + one abstract RHS preparation
#    (all open slots, no data yet).
inputs = roe_qfvm_inputs(fmt=fmt, angle_width=6)
# 2. Classical-side flow: RoeFlowData writes the six banks rho/momentum/energy/rhs_*
#    and maintains the norm tree.
flow = RoeFlowData([(1, 0, 1), (1, 0.25, 1), (1.25, 0, 1), (1.25, -0.25, 1)],
                   fmt=fmt, angle_width=6)
# 3. Assemble the problem: Hermitian dilation + padding diagonal folded into the
#    spectral declaration; spectral bounds are the caller's duty.
problem = roe_qfvm_problem(
    inputs, spectrum=SpectralPromise(norm_upper=10.0, sigma_min_lower=0.5), amax=4.0,
)
# 4. Hand the same problem to both routes: CKS consumes the SparseSystem;
#    Costa triggers the sparse→BE adaptation.
cks = make_cks_qlss(CKSConfig(order=2))(problem)
costa = make_costa_qlss(CostaConfig(steps=1))(problem)
# 5. Close the open slots: five databases bind to concrete QRAM banks, the RHS slot
#    binds to the signed residual-tree preparation.
cks_rir = bind_qfvm(cks.operation.program(), inputs)
# 6. Materialize the runtime memory tables: the six flow-snapshot banks
#    + the geometry table + the theta table.
memories = qfvm_memories(inputs, flow, amax=4.0)
# 7. At execution hand memories to the backend by name (run_pysparq(cks_rir, memory=memories));
#    switching to Costa only changes the protocol in step 4; every other line is unchanged.
```

## Numerical validation

The paper-grade numerical validation of 2026-09-16 (`tests/verification/verify_qham_qfvm.py`, real backends with no mock substitutes) numerically confirmed the data and access constructions of this chapter:

- **Numerical flux**: the compiled Roe face circuit agrees bit for bit with an independent fixed-point simulation (32/32 branches), with a method error of 0.427 against the float64 Roe formula (inherent 5-bit fixed-point quantization); the sparse entry oracle agrees bit for bit on sampled columns, the padding diagonal yields exactly padding_value, and everything outside the structural domain is zero.
- **Single-step update**: F*(L,R)=left·U_L+right·U_R differs from the matrix-inversion implementation {obj}`riemann_flux <oracq.applications.flow_data.riemann_flux>` by 3.33e-16; the matrix identity M·u−mass·u=−residual (implicit FVM sign convention) has residual 2.22e-16; `RoeFlowData` local updates recompute only adjacent interfaces and affected tree nodes and agree bank by bank with a full recomputation.
- **Position access**: under a superposition of all 32 columns every column is a complete permutation, the 9 structural slot mappings show 0 mismatches against the independent geometry semantics, and the geometry QRAM table is point-wise true.
- **RHS preparation**: the residual-state amplitudes agree with independent residual/norm computations (error 1.11e-16), signs are written exactly via Z kickback, and the angle tree and sign bank are point-wise true.

The validation uses {obj}`FixedFormat(5,2) <oracq.algorithms.common.arithmetic.FixedFormat>` and `entropy_delta=0.5` (all constants exactly representable); under lower-precision formats 2δ may truncate to zero, causing division by zero in the entropy-correction branch and zeroing the entry per the documented behavior (see the <a href="../../zh/manual/algorithms/qfvm.html#数值验证">algorithm page's numerical validation</a>). The numerical solve accuracy of the two QLSS routes belongs to the solver side and is out of scope for this page. Artifacts: `out/verification/qham_qfvm.json`.

## QMem direct parallel path

`applications/qfvm_qmem.py` provides a parallel implementation of the same QFVM data access, with the entire data plane moved to [pointer-style QRAM access](qdata.md):

- The three conserved-quantity banks merge into one `(field, cell)` state table {obj}`QMem(b, "state", shape=(3, n)) <oracq.infrastructure.qmem.QMem>`; periodic neighbors are computed by modular addition on a cell_width-bit scratch (as in the old path) and then queried through the two-dimensional pointer `[field, addr]`.
- The geometry table is addressed two-dimensionally by `(slot, column)` (point-wise identical to the old path's {obj}`fuse(column, slot) <oracq.infrastructure.ir.fuse>` addressing).
- The residual state is prepared by the {obj}`QVector <oracq.algorithms.input_model.qdata.QVector>` squared-norm tree (replacing the `qram_state_prep` + sign-bank combination).
- Modules declare QRAM-form resources directly, bypassing abstract database slots and {obj}`bind_qfvm <oracq.applications.qfvm.bind_qfvm>`.

The reversible Roe arithmetic, the nine-slot geometry selection, the padding diagonal, and the in-place position permutation have circuit structure identical to this path. Equivalence evidence (`tests/core/test_qfvm_qmem.py`): superposition probes of geometry/state/periodic-neighbors agree word for word, the position oracle's full-column permutations agree amplitude by amplitude, residual-state amplitudes agree, and the physical layer containing compiled Roe arithmetic agrees bit for bit at single points. The slot-binding path described above remains the main path for QLSS integration; the QMem path is currently positioned as an equivalent rewrite of the data-access layer and a baseline for future evolution.
