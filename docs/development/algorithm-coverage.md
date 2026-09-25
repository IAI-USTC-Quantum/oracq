# Quantum algorithm coverage workboard

**English** · <a href="../../zh/development/algorithm-coverage.html">简体中文</a>

This document tracks oracq's quantum-algorithm implementation coverage
relative to the open-source ecosystem, ordering implementation by the
principle of "input-model degrees of freedom first". Core principle:
algorithms are decoupled from input models — one algorithm definition serves
several input models (the abstract / gate / qram layers) through open
declarations and batched binding.

## Input model vocabulary

| Abbreviation | Input model | Repository location |
|---|---|---|
| SP | State preparation / amplitude encoding | `algorithms/oracles.py` StatePreparation (gate/qram) |
| QRAM | QRAM/QROM data lookup | `infrastructure/ir.py` Load, `oracles.py` XorDatabase/banked |
| SO | Sparse-matrix oracle (position + element separated) | `oracles.py` SparseAccess |
| BE | Block encoding / LCU | `algorithms/operators.py`, `block_encoding.py` |
| UO | Unitary-operator oracle (controlled U) | `interfaces.py` UnitaryProtocol |
| HAM | Hamiltonian representation (Pauli/sparse/low-rank) | `algorithms/hamiltonian.py` |
| FO | Function oracle (phase/amplitude/Boolean) | `infrastructure/mathfunc/`, `oracles.py` |
| CP | Direct parameterization by classical data (QUBO, variational parameters) | `algorithms/variational.py` etc. |
| ODE | Differential-equation input (access model for A + initial state + b) | `algorithms/ode_models.py`, `ode.py`, `qham.py` |
| DM | Density matrix / Gibbs state | Not implemented |
| EST | Sampling/measurement access (shadows, tomography) | Not implemented |

## Workboard

Status: ✅ implemented; 🚧 in progress; 📄 candidate paper (no reliable
open-source implementation, or worth rewriting with this repository's
abstractions).

### Class A: frameworks explicitly decoupling algorithms from input models (highest priority)

| Entry | Algorithm | Input model | Status | Notes |
|---|---|---|---|---|
| A1 | QSVT standard-transform library: phase synthesis (parity decomposition), fixed-point search, matrix inversion/sign function, eigenstate filtering, Hamiltonian-simulation phases | BE | ✅ | Gilyén et al. 2019; implemented in `algorithms/qsvt.py`, assembled on the `transforms.py` sequence skeleton |
| A2 | PREPARE–SELECT standard decomposition: independent prepare/select oracles, coefficients bound via gate/QRAM/alias, wired through to the qubitization walk | BE/HAM | ✅ | Low–Chuang 2019, Babbush et al. 2018; implemented in `algorithms/prepare_select.py`, including alias sampling |
| A3 | Quantum-walk framework: graph-oracle input model, Szegedy/MNRS marked-vertex search frameworks | FO (adjacency oracle) | ✅ | MNRS 2011; implemented in `algorithms/graph_walks.py`, including the three-layer adjacency oracle and hitting-time tools |
| — | Element distinctness / collision finding | FO | 📄 | Ambainis; depends on the A3 framework + Johnson-graph walks (needs nonuniform setup and extra data registers) |
| — | Non-Abelian HSP | FO | 📄 | High mathematical barrier, deferred |

### Class B: new input models opening new algorithm families

| Entry | Algorithm | Input model | Status | Notes |
|---|---|---|---|---|
| B1 | Gibbs/thermal-state preparation, quantum Metropolis | DM + HAM | ✅ | Implemented in `algorithms/density.py`: the three-layer purification-access paradigm + Gibbs preparation via the QSVT purification route (the DM input model is established) |
| B2 | Quantum SDP / convex optimization (Brandão–Svore, van Apeldoorn–Gilyén) | DM/BE | ✅ | Implemented in `algorithms/qsdp.py`: trace-estimation probe circuits + MMW driving + per-round quantum-subroutine generation |
| B3 | Heinrich quantum summation/integration | FO + QRAM | ✅ | Implemented in `algorithms/integration.py`: comparator linear readout + QAE, cross-checked across the three binding layers |
| B4 | Jordan gradient estimation | FO (probabilistic oracle) | ✅ | Implemented in `algorithms/gradient.py`: three-layer phase oracle + single-query d-component readout |
| B5 | Alias-sampling state preparation, Select-Swap QROM | CP→SP / QRAM | ✅ | Alias in `prepare_select.py`; Select-Swap implemented in `algorithms/data_loading.py` (with the {obj}`qrom_cost <oracq.algorithms.input_model.data_loading.qrom_cost>` resource cross-check) |

### Class C: algorithms that this repository's abstraction level can rewrite more correctly

| Entry | Algorithm | Input model | Status | Notes |
|---|---|---|---|---|
| C1 | Quantum-chemistry low-rank decomposition (DF/THC/SF → LCU → BE) | HAM (low-rank tensors) → BE | ✅ | Implemented in `algorithms/lowrank.py`: DF (two-level synthetic rotations + controlled-rotation diagonal encoding) and THC; the products can feed qubitization_walk |
| C2 | Recent nonlinear/stochastic ODE advances (Fokker–Planck, SDE) | ODE | ✅ | Implemented in `algorithms/sde.py`: conservative-form FP discretization feeding the QODEProblem contract, with analytic OU-moment cross-checks |
| C3 | DQI (decoding quantum interferometry optimization) | CP (encoding constraints) | ✅ | Implemented in `algorithms/dqi.py`: Dicke states + open decoder declarations, with pointwise distribution cross-checks against Krawtchouk closed forms |
| C4 | QRAM-based QML (QPCA, recommendation systems) | QRAM + SP | ✅ | QPCA implemented in `algorithms/qpca.py`: LMR density-matrix exponentiation + QPE, with a first-order convergence-rate witness; recommendation systems deferred |

### Implemented baseline (cross-check catalog in `applications/catalog.py`)

| Area | Algorithms |
|---|---|
| Oracle/search/estimation | Deutsch–Jozsa, Bernstein–Vazirani, Simon, Grover, amplitude amplification, QAE/quantum counting, QPE, Hadamard/SWAP test |
| Fourier/number theory | QFT, QFT addition, modular multiplication, order finding |
| Simulation | Trotter, Taylorization, qubitization walk, QSVT sequence skeleton, OAA |
| Linear systems | CKS (Chebyshev), Costa (discrete adiabatic walk), the CKS §5 VTAA variable-time solver |
| Differential equations | LCHS, CBMD, Schrödingerization, Carleman, Euler history states, structured FD, QHAM, QFVM |
| Other | Variational (QAOA-MaxCut/VQE), repetition-code error correction, fixed-point arithmetic, mathfunc front end |

## Maintenance rules

- New algorithms go into the matching category file; do not write legacy-import
  compatibility layers; register names, bit widths, and views are preserved
  down to backend lowering.
- When editing this workboard, update the status column in sync; once an
  algorithm is implemented, note the implementation file on its row.
- Criteria for picking new papers: the paper defines its input oracles
  explicitly; it can be transcribed mechanically into the protocol-contract
  system; its correctness can be endorsed by the contracts checks.
- Per-algorithm implementation and validation details live on the <a href="../../zh/manual/algorithms/index.html">algorithm
  pages</a>; update the corresponding pages
  when editing the workboard.
