# Algorithm catalog

**English** · <a href="../../../zh/manual/algorithms/index.html">简体中文</a>

The algorithm library is organized into ten subpackages by purpose: `input_model` (input models and data access), `common` (shared primitives), `qlss` (linear systems), `qnlss` (nonlinear systems), `qode` (ordinary differential equations), `qpde` (partial differential equations), `qml` (quantum machine learning), `optimization` (quantum optimization and variational methods), `basics` (basic example algorithms), and `qec` (quantum error correction). The table below lists the entry files, what is implemented, and the boundaries to keep in mind when using them. The API reference lists the full signatures; the per-algorithm pages are grouped by family — see the navigation below.

| Subpackage | File | Implementation and scope |
|---|---|---|
| input_model | [`contracts.py`](../../api/algorithms/input_model/contracts.rst) | Checkable input contracts: oracle capabilities/specifications, protocol contracts, and acceptance reports |
| input_model | [`operators.py`](../../api/algorithms/input_model/operators.rst) | BlockEncoding class and basic compositions such as identity/product/scale/LCU |
| input_model | [`oracles.py`](../../api/algorithms/input_model/oracles.rst) | The four oracle paradigms XorDatabase, StatePreparation, StateOracle, SparseAccess and the abstract/gate/qram factories |
| input_model | [`interfaces.py`](../../api/algorithms/input_model/interfaces.rst) | Shared protocols such as StatePreparation/Unitary/BlockEncoding and adaptation functions |
| input_model | [`block_encoding.py`](../../api/algorithms/input_model/block_encoding.rst) | Block-encoding algebraic composition: tensor, direct_sum, projector, lcu, pauli_word, and more |
| input_model | [`sparse.py`](../../api/algorithms/input_model/sparse.rst) | Adaptation from sparse access to block encodings, including the Chebyshev block and reversible table lookup |
| input_model | [`spectral.py`](../../api/algorithms/input_model/spectral.rst) | Spectral diagonal block encodings, sparse spectral block encodings, and spectral state preparation |
| input_model | [`lowrank.py`](../../api/algorithms/input_model/lowrank.rst) | Quantum-chemistry low-rank decomposition (DF/THC) Hamiltonian LCU/block encodings |
| input_model | [`data_loading.py`](../../api/algorithms/input_model/data_loading.rst) | Select-Swap QROM data loading and its cost model |
| input_model | [`qdata.py`](../../api/algorithms/input_model/qdata.rst) | QVector/QMatrix quantum data structures (squared-norm trees and sample-and-query) |
| input_model | [`density.py`](../../api/algorithms/input_model/density.rst) | Density-matrix purification access, Gibbs state preparation, and classical tools such as trace distance |
| input_model | [`qham.py`](../../api/algorithms/input_model/qham.rst) | Constructing QODE inputs from finite HAM closures and the physical output channel |
| input_model | [`graph_walks.py`](../../api/algorithms/input_model/graph_walks.rst) | Graph-oracle input models and the Szegedy/MNRS walk-search framework |
| common | [`prepare_select.py`](../../api/algorithms/common/prepare_select.rst) | The standard PREPARE–SELECT decomposition of an LCU and alias sampling |
| common | [`qsvt.py`](../../api/algorithms/common/qsvt.rst) | Standard QSVT transforms: QSP phases, matrix inversion, eigenstate filtering, fixed-point search |
| common | [`transforms.py`](../../api/algorithms/common/transforms.rst) | Assembly of qubitization, QSVT with explicit phase sequences, and oblivious amplification |
| common | [`state_preparation.py`](../../api/algorithms/common/state_preparation.rst) | Initial-state composition: extended initial states, physical subspace selection, applying a BE to a state |
| common | [`hamiltonian.py`](../../api/algorithms/common/hamiltonian.rst) | Pauli-term evolution, Trotter composition, Taylor BE, and the injectable QSP interface |
| common | [`fourier.py`](../../api/algorithms/common/fourier.rst) | Forward/inverse QFT, zero-width work adaptation, and Fourier addition modulo $2^n$ |
| common | [`arithmetic.py`](../../api/algorithms/common/arithmetic.rst) | Reversible fixed-point arithmetic: Boolean SSA networks, compute/uncompute, native registration |
| common | [`estimation.py`](../../api/algorithms/common/estimation.rst) | QPE, canonical amplitude estimation, Hadamard test, swap test |
| common | [`search.py`](../../api/algorithms/common/search.rst) | Grover, the Grover iterate, and amplitude amplification of a success subspace |
| common | [`walks.py`](../../api/algorithms/common/walks.rst) | Hadamard coined walks on periodic lattices |
| common | [`spectral_synthesis.py`](../../api/algorithms/common/spectral_synthesis.rst) | Operator-level synthesis optimization of spectral circuits (uniformly controlled preparation, fanout merging) |
| common | [`integration.py`](../../api/algorithms/common/integration.rst) | Heinrich quantum summation and numerical integration |
| qlss | [`qlss.py`](../../api/algorithms/qlss/qlss.rst) | Problem contracts, Costa walk/filter, and the CKS base Chebyshev/LCU routes |
| qlss | [`vtaa_cks.py`](../../api/algorithms/qlss/vtaa_cks.rst) | CKS §5 variable-time amplitude amplification: QSP verdict clock, banded inverse LCU, Ambainis nested amplification, and A' uncomputation |
| qnlss | [`newton.py`](../../api/algorithms/qnlss/newton.rst) | Quantum Newton method: M_F data structure, differenced Jacobian oracle |
| qnlss | [`carleman.py`](../../api/algorithms/qnlss/carleman.rst) | Finite-order Carleman tensor lifting for polynomial ODEs |
| qode | [`ode.py`](../../api/algorithms/qode/ode.rst) | QODE protocol and Euler history assembly; concrete methods live in separate files |
| qode | [`ode_models.py`](../../api/algorithms/qode/ode_models.rst) | Hermitian anti-Hermitian decomposition and the LinearODE shared input model |
| qode | [`cbmd.py`](../../api/algorithms/qode/cbmd.rst) | Contour-decomposition matrix functions and evolution assembly |
| qode | [`lchs.py`](../../api/algorithms/qode/lchs.rst) | Finite weighted sums of Hermitian branches for dissipative linear evolution |
| qode | [`schrodingerization.py`](../../api/algorithms/qode/schrodingerization.rst) | Quantum representation of non-unitary evolution via auxiliary coordinates + Fourier transform |
| qode | [`sde.py`](../../api/algorithms/qode/sde.rst) | Input models from Fokker–Planck/SDE to linear ODEs and classical witnesses |
| qode | `legacy.py` | Compatibility layer of the early evolution factories (make_lchs_qode and friends) |
| qpde | [`pde.py`](../../api/algorithms/qpde/pde.rst) | Thin QPDE input and solve wrapper on top of the QODE protocol |
| qml | [`recommendation.py`](../../api/algorithms/qml/recommendation.rst) | Kerenidis–Prakash quantum recommendation systems |
| qml | [`qpca.py`](../../api/algorithms/qml/qpca.rst) | Quantum principal component analysis via density-matrix exponentiation + phase estimation |
| qml | [`qsdp.py`](../../api/algorithms/qml/qsdp.rst) | Quantum semidefinite programming framework from Gibbs sampling + trace estimation + matrix multiplicative weights |
| qml | [`qcnn.py`](../../api/algorithms/qml/qcnn.rst) | Classical-side tensor/convolution/pooling logic of quantum convolutional neural networks |
| qml | [`qcnn_layer.py`](../../api/algorithms/qml/qcnn_layer.rst) | Quantum building blocks of quantum convolutional neural networks |
| optimization | [`dqi.py`](../../api/algorithms/optimization/dqi.rst) | DQI decoded quantum interference: GF(2) max-XORSAT |
| optimization | [`variational.py`](../../api/algorithms/optimization/variational.rst) | Parameterized ansätze, MaxCut QAOA, Pauli measurements, and VQE measurement-circuit collections |
| optimization | [`gradient.py`](../../api/algorithms/optimization/gradient.rst) | Jordan quantum gradient estimation |
| basics | [`oracle_algorithms.py`](../../api/algorithms/basics/oracle_algorithms.rst) | D-J, Bernstein–Vazirani, Simon sampling; Simon's GF(2) elimination runs on the classical side |
| basics | [`number_theory.py`](../../api/algorithms/basics/number_theory.rst) | Finite-size modular-multiplication permutations, QPE order finding, continued-fraction factor-candidate post-processing |
| qec | [`error_correction.py`](../../api/algorithms/qec/error_correction.rst) | Encoding and coherent recovery of the three-qubit bit/phase-flip repetition code |

## How to choose a starting point

If the input is a Boolean function, start with the query and search family. To
estimate a probability or an expectation value, start with `estimation`. Given
a decomposition into Hamiltonian terms, Trotter is available; with only an
operator block encoding, choose an implementation that consumes that access
model.

PDE algorithms first construct the spatial-discretization inputs, then choose
the QODE or the linearization route. Data living in QRAM does not mean the
matrix is already block encoded; the input-adaptation step still has to be
explicit.

## Running the gallery

```bash
uv run python examples/algorithm_gallery.py
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

The gallery contains 22 small examples, including different applications of
the same algorithm, for example amplitude estimation and quantum counting.
They illustrate interfaces and readout methods; they are not 22 unrelated
algorithms.

## Implementation basis

Fourier addition follows [Draper's QFT adder construction](https://arxiv.org/abs/quant-ph/0008033). Amplitude amplification and estimation follow [the framework of Brassard et al.](https://arxiv.org/abs/quant-ph/0005055). The QAOA cost/mixer layering follows [the original algorithm of Farhi et al.](https://arxiv.org/abs/1411.4028). Order finding and the classical factor post-processing follow [Shor's construction](https://arxiv.org/abs/quant-ph/9508027).

The current modular multiplication uses capped permutation synthesis, by
default at most 8 bits; it exercises the order-finding interface and circuit
and does not represent scalable Shor modular arithmetic. VQE and QAOA provide
the quantum circuits; the classical optimizer is chosen by the application.
QSVT receives caller-provided phase sequences and does not include a general
phase solver.

## Algorithm pages

One page per algorithm, describing the interface, input model, implementation
notes, and where to find the validation evidence. The algorithm pages are
organized into thirteen families, starting with data loading and input models
and ending with applications and infrastructure; the line under each page's
title records its category (C1–C6 or the application layer) and defining
module — the category definitions are in the
[validation plan](../../development/validation-plan.md).

```{toctree}
:maxdepth: 1
:caption: Algorithm families

groups/input-models
groups/basics
groups/search
groups/fourier-arithmetic
groups/estimation
groups/hamiltonian-evolution
groups/qlss
groups/differential-equations
groups/optimization
groups/walks
groups/qml
groups/qec
groups/applications
```

## Numerical validation

Paper-grade numerical experiments for quantum semidefinite programming
(`qsdp.py`, no dedicated page yet) and its inner primitives in `density.py`
(purification access, Gibbs state preparation) are in
`tests/verification/verify_misc_algorithms.py` (the misc_algorithms group),
all executed on real backends with fully independent classical oracles (numpy
eigendecomposition / `eigvalsh` / closed-form expressions).

**Experiment design**:
(a) purification witness — the purified state of a 2×2 complex density matrix
(off-diagonal entries 0.1±0.05j) is prepared by {obj}`gate_purification <oracq.algorithms.input_model.density.gate_purification>` and, after taking the partial trace over the environment, cross-checked against the original matrix element by element;
(b) Gibbs purification — two instances, diagonal $H = \mathrm{diag}(1, -0.5)$ ($\beta = 0.6$) and dense non-diagonal $H$ ($\beta = 0.8$), with `error = 0.05`; the normalized reduced state of the postselected signal == 0 branch is compared against the numpy Gibbs state $e^{-\beta H}/Z$, plus convergence sweeps at error = 0.4/0.2/0.1;
(c) trace estimation — probe readouts of {obj}`trace_estimate_circuit <oracq.algorithms.qml.qsdp.trace_estimate_circuit>` for the Z/X/Y Pauli observables on a complex density matrix against the numpy trace $\mathrm{Tr}(P\rho)$, and joint decoding of {obj}`trace_from_joint <oracq.algorithms.qml.qsdp.trace_from_joint>` along the Gibbs approximate-purification path;
(d) MMW driver — a four-constraint feasibility instance (equalities $r_z = 0.2$, $r_x = 0.1$ split into two-sided inequalities, $\varepsilon = 0.08$); after convergence, the returned average iterate $\bar\rho$ is independently audited for PSD/trace-1/violations;
(e) single quantum round — estimates produced on the reference executor by the Gibbs purification + per-constraint trace-estimation circuits generated by {obj}`iteration_circuits <oracq.algorithms.qml.qsdp.iteration_circuits>` are compared with a classical estimator constraint by constraint.

**Key metrics**:

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| purification-partial-trace | 2 qubits | four-path full-amplitude cross-check | max elementwise error | 1.1e-16 |
| gibbs-purification-diagonal | 9 qubits, qsp_degree 3 | reference + rir-pysparq + originir-ext | trace distance / postselection success rate | 3.0e-5 / 0.0965 |
| gibbs-purification-nondiagonal | 13 qubits, qsp_degree 3 | same | trace distance / postselection success rate | 1.9e-6 / 0.1130 |
| gibbs-error-scaling | error = 0.4/0.2/0.1 | reference | trace distance | 5.6e-4 / 8.7e-5 / 8.7e-5 (all ≤ error, monotone non-increasing) |
| trace-estimate-pauli-xyz | ≤5 qubits | three paths | Tr(Zρ)/Tr(Xρ)/Tr(Yρ) estimates | 0.4 / 0.2 / 0.1 (max error 5.3e-16) |
| trace-estimate-joint-gibbs | 10 qubits | reference + rir-pysparq | Tr(Zρ_gibbs) joint decoding | −0.42196 vs −0.42190 (error 6.0e-5) |
| qsdp-mmw-driver-feasibility | 4 constraints, 256 iterations | classical driver | convergence / max violation / min eigenvalue | ✓ / 0.0738 ≤ 0.08 / 0.4295 |
| qsdp-mmw-quantum-round | 2 constraints, single round | reference + rir-pysparq | max difference, quantum vs classical estimate | 2.3e-6 |

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

Artifacts: `out/verification/misc_algorithms.json` (24 cases, all passing; this section covers the eight `purification-*`, `gibbs-*`, `trace-estimate-*`, and `qsdp-*` cases). Numerical validation of the group's remaining algorithms is on their own pages (density-matrix exponentiation, QPCA, DQI, variational ansätze, VQE, QAOA, repetition codes).

## Numerical validation (application catalog and infrastructure modules)

All entries of `applications/catalog.py` and `applications/gallery.py` are
covered by the nt_qlss_sde group in
`tests/verification/verify_nt_qlss_sde.py`: every entry runs on real backends
(rir-pysparq, adapter-pysparq, some including the reference executor) and is
cross-checked amplitude by amplitude against the reference executor;
oracle-level numerical correctness is carried by the validation sections of
the underlying algorithms' own pages (number theory, CKS, and so on — see the
corresponding pages in this catalog).

**Experiment design**: the 22 gallery examples are cross-checked at full
amplitude on three backends (reference / rir-pysparq / adapter-pysparq); of the
33 catalog examples, 31 are cross-checked reference vs rir-pysparq
(`stateprep_qram` doubles as the regression sentinel for pysparq.rir work-bit
uncomputation and is also full-state cross-checked against the adapter),
while `qham_qode` and `qham_qpde` — with more than $10^6$ expansion steps and
roughly 60–90 seconds per backend — complete their reference vs adapter
cross-checks through spawned concurrent subprocesses. These entries are
**assembly examples**: their numerical content (modular-multiplication
permutations, order-finding distributions, CKS solution states, and so on) is
certified by the independent oracle validations on each algorithm's page; the
catalog level only pins down assembly determinism and does not repeat physical
oracles — which is why catalog/gallery entries are covered by "structural
tests + cross-backend comparison" rather than independent numerical oracles.

**Key metrics**:

| Case | Scale | Paths | Metric | Value |
|---|---|---|---|---|
| gallery-* (22 examples) | small instances | three backends | max cross-backend amplitude deviation | 0.0 (all) |
| catalog-* (29 examples) | small instances | reference + rir-pysparq | max cross-backend amplitude deviation | ≤ 1.4e-16 |
| catalog-stateprep_qram | QRAM angle table | three backends | full-state deviation / work-bit residue | 0.0 / 0.0 |
| catalog-qham_qode, qham_qpde | 39296 states | reference + adapter | max cross-backend amplitude deviation | 1.8e-9 |

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

Artifacts: `out/verification/nt_qlss_sde.json` (`gallery-*` 22 examples, `catalog-*` 33 examples).
