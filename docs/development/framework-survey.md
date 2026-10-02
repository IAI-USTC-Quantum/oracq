# Framework survey notes (stage 1)

**English** · <a href="../../zh/development/framework-survey.html">简体中文</a>

Working notes collected for the review-driven revision. Each section ends in
one-sentence conclusions so the paper edits never need to re-verify facts.
Everything below is documentation- and paper-level evidence gathered on
2026-10-02; no benchmarking was performed. "Not found in docs" means absent
from the public documentation (and, for Qualtran, the repository layout) as
of that date — it is not a claim about unreleased code.

## Qrisp BlockEncoding paper (arXiv:2604.18276)

Bibliographic record (from the arXiv abstract page): Matic Petrič and René
Zander, "Block-encodings as programming abstractions: The Eclipse Qrisp
BlockEncoding Interface", 2026, arXiv:2604.18276 [quant-ph], v1 of
2026-04-20, DOI 10.48550/arXiv.2604.18276. Companion documentation:
[Qrisp BlockEncoding reference](https://qrisp.eu/reference/Block%20Encodings/BlockEncoding.html).

Content, point by point against oracq:

- Construction interface: `BlockEncoding(alpha, ancillas, unitary, num_ops,
  is_hermitian)` with factories `from_array`, `from_operator`
  (QubitOperator/FermionicOperator), `from_lcu` (PREP/SELECT with a
  balanced-binary-tree `q_switch`), `from_projector`, `from_eye`. Ancillas
  are declared as templates and the operator shape is dynamic at application
  time. oracq's closest constructors are the explicit Pauli expansion
  (`matrix_pauli_encoding`, ≤ 5 qubits), `projector`, and the sparse/low-rank
  input models — oracq has no dense-matrix or physics-operator constructor.
- Algebra: operator overloading `+`/`−` (LCU), scalar `*` (α → |c|·α), `@`
  (product, α multiplicative), `kron`, `dagger`, plus `qubitization()` and
  `chebyshev(k)`. This is a direct counterpart of oracq's block-encoding
  algebra (`linear_combination`, `product`, `scale`, `tensor`, `adjoint_be`,
  `qubitization_walk`, `chebyshev_block`), including the tracked α.
- Matrix inversion: `inv(eps, kappa, method={"QET","QSVT","GQSVT"})` with the
  stated complexity 𝒪(κ² log(κ/ε)) = degree-𝒪(κ log(κ/ε)) polynomial × 𝒪(κ)
  post-selection repetitions; phase synthesis is delegated to internal
  solvers and **no degree cap or precision limit is documented**. oracq's
  `qsvt_matrix_inversion` documents its degree guard (40) and validates every
  synthesis by a round-trip self-check. Qrisp additionally ships turnkey
  solvers (`CKS`, `dalzell_inversion`, `lanczos_alg`), `pseudo_inv`, `svt`,
  `poly`, and `sim(t, N)`.
- Resource estimation: `resources(*operands)` returns gate counts, depth,
  and qubits for one execution of the encoding unitary, built on Jasp
  tracing; whole-program counting via `@count_ops`. oracq's
  `estimate_resources` returns a fault-tolerant ledger (Toffoli/Clifford+T
  atoms, rotations pending synthesis, QRAM queries/writes, per-oracle call
  counts, an incompleteness flag for open slots) without executing anything.
- Execution: Qrisp BEs are directly runnable (`apply`, `apply_rus`,
  `expectation_value`) and JAX-integrated; oracq separates generation from
  execution and reads out through independent executor paths.

Differences that remain oracq's own: persisted open programs with
declared-but-unbound oracles and checked late binding; named runtime QRAM
resources; static contract checking with structured issues; the
Python-formula-to-reversible-circuit compiler. One-sentence conclusion: the
Qrisp BlockEncoding interface is a direct comparator for the component and
analysis layers — same algebra with tracked α, inversion, and resource
queries — while the differences are the persisted open-program workflow and
the scientific-convention/QRAM capture around it.

## Qualtran and Catalyst evidence check

For the softened related-work claims. Sources: the
[Qualtran Bloq reference](https://qualtran.readthedocs.io/en/latest/reference/qualtran/Bloq.html),
[Register](https://qualtran.readthedocs.io/en/latest/reference/qualtran/Register.html)
and
[Signature](https://qualtran.readthedocs.io/en/latest/reference/qualtran/Signature.html)
pages, the
[resource counting module](https://qualtran.readthedocs.io/en/latest/reference/qualtran/resource_counting.html),
the serialization
[module page](https://qualtran.readthedocs.io/en/latest/reference/qualtran/serialization.html)
plus `qualtran/serialization/bloq.py` and `qualtran/protos/bloq.proto` in the
repository, the Qualtran paper (arXiv:2409.04643), and the Catalyst
[index](https://docs.pennylane.ai/projects/catalyst/en/stable/index.html),
[architecture](https://docs.pennylane.ai/projects/catalyst/en/stable/dev/architecture.html),
and
[dialects](https://docs.pennylane.ai/projects/catalyst/en/stable/dev/dialects.html)
pages.

- Qualtran **has**: register signatures (name, dtype, shape, side);
  hierarchical decomposition (`decompose_bloq`, `CompositeBloq`,
  `BloqBuilder` with register-mismatch errors at construction); symbolic
  call counts without a complete decomposition (`build_call_graph`,
  SymPy-typed counts); protobuf serialization. Undecomposed bloqs survive
  serialization only as shallow name+signature+attributes entries whose
  deserialization resolves against a registry of concrete Python classes.
- Qualtran: **no documented** late binding (replacing a sub-bloq's
  realization after the program object exists, under an interface re-check;
  CompositeBloq is documented as immutable) and **no documented** named
  runtime QRAM resources (QROM/QROAM are circuit bloqs).
- Catalyst **has** a high-level structured hybrid IR (MLIR dialects:
  Quantum, Gradient, Catalyst, …), but it is a transient compiler artifact:
  the pipeline lowers to LLVM/QIR and produces a binary, and intermediate
  stages are exposed only through debug utilities. **No documented**
  declared-but-unimplemented quantum oracle slots resolved later under an
  interface check (the docs index contains no "oracle" entries;
  `quantum.custom` is a fully specified named gate).

One-sentence conclusion: the mechanisms exist but stop at different points —
Qualtran carries signatures, hierarchy, symbolic counts, and shallow
serialization of abstract bloqs, and Catalyst carries a high-level hybrid IR
inside the compiler; neither documents the complete workflow of persisting an
open program, checking scientific conventions at binding, capturing named
runtime QRAM resources, and repeating partial binding.

## QSP phase-generation interfaces

Sources: the [pyqsp README](https://github.com/ichuang/pyqsp/blob/master/README.md)
and source docstrings (`pyqsp/angle_sequence.py`, `pyqsp/sym_qsp_opt.py`,
`pyqsp/response.py`), and the
[pyLIQTR Features.md](https://github.com/isi-usc-edu/pyLIQTR/blob/main/docs/Features.md)
plus its `phase_factors` package.

- pyqsp v0.2.0 (MIT; deps numpy/scipy/matplotlib; `pip install pyqsp`):
  inputs are **Chebyshev-basis** coefficients; the main entry
  `QuantumSignalProcessingPhases(poly, method=…, chebyshev_basis=True)`
  returns d+1 phases for a degree-d target. Phase convention: the Wx signal
  operator, U = e^{iφ₀Z} ∏ W·e^{iφ_k Z}, φ₀ the last-applied rotation, no
  π/4 offsets; for `method="sym_qsp"` the target appears as Im⟨0|U|0⟩.
  Methods: `laurent` (root-finding completion, machine precision at moderate
  degree, unstable at high degree) and `sym_qsp` (Newton iteration on reduced
  phases with an FFT Jacobian, after arXiv:2002.11649 and arXiv:2307.12468;
  documented as stable "well into the thousands of phases"). Response
  self-evaluation: `pyqsp.response.ComputeQSPResponse` / `PlotQSPResponse`.
- pyLIQTR does not wrap pyqsp: its `Angler_opt` is an L-BFGS symmetric-QSP
  generator "pulled very directly" from QSPPACK, targeting **Re⟨0|U|0⟩** with
  explicit π/4 end offsets — its phases are not interchangeable with pyqsp's
  without conversion. Its older scipy/mpsolve root-finding path is being
  deprecated.
- PennyLane correction: the "externally supplied phase angles" reading is
  outdated — `qml.poly_to_angles` (root-finding to ~degree 1000, iterative
  L-BFGS-B beyond) now exists, documented to validate circuits against P(x)
  at ~1e-10.

oracq's contract choice, implemented in this stage: a PhaseSynthesizer
receives ascending real coefficients of the real target f (monomial basis)
plus an optional pinned imaginary completion h and returns time-ordered
phases in the reflection convention; pinned mode requires p = f + i·h
exactly, free mode pins only Re p = f. The pyqsp adapter (free mode only)
converts the target to the Chebyshev basis and maps Wx phases into the
reflection convention (reverse the list, subtract π/4 from the end phases
and π/2 from the middle ones, add a d mod 4 dependent offset to the
first-applied phase, negate the target when d ≡ 3 mod 4); the mapping is
independently verified by `qsp_response` grid checks on every call. pyLIQTR's
Re-target π/4-offset convention was not adopted. One-sentence conclusion: the
phase-convention mapping is exact and machine-checked, and the adapter lifts
only the synthesis-degree guard — the monomial coefficient pipeline has its
own reliability limit around degree 40, so high-degree targets must be
constructed and validated in the Chebyshev basis (the κ=8, ε=1e-2 inversion
target is reached this way at degree 585 with round-trip residual ≈ 1e-13).

## Maturity grading material

Four levels: **L1** declared interface · **L2** circuit generation · **L3**
closed execution (simulate/run end-to-end) · **L4** specified numerical
precision with validation. Documentation-level grading, sources linked per
framework; "—" means not found in the reviewed documentation.

| Framework | L1 interface | L2 circuits | L3 execution | L4 precision-validated | QLS | ODE/PDE |
|---|---|---|---|---|---|---|
| Qualtran | yes ([BE bloq](https://qualtran.readthedocs.io/en/latest/bloqs/block_encoding/block_encoding.html)) | yes ([GeneralizedQSP](https://qualtran.readthedocs.io/en/latest/bloqs/qsp/generalized_qsp.html)) | partial (classical bloq simulation; purpose is resource estimation) | — | — | — |
| pyLIQTR | yes ([Features.md](https://github.com/isi-usc-edu/pyLIQTR/blob/main/docs/Features.md)) | yes (QSVT circuit classes, OpenQASM export) | partial (cirq in internals; deliverable is Clifford+T estimates) | — | — | classical only (`clam`) |
| PennyLane | yes ([qml.QSVT](https://docs.pennylane.ai/en/stable/code/api/pennylane.QSVT.html), `poly_to_angles`, FABLE, Qubitization) | yes (BlockEncode is simulator-only) | yes (default.qubit doc examples) | yes (`poly_to_angles` ~1e-10 validation; QSVT demos vs classical) | QSVT inversion demos | — |
| Classiq | yes ([qsvt_inversion library](https://docs.classiq.io/latest/explore/algorithms/quantum_linear_solvers/qsvt_matrix_inversion/qsvt_matrix_inversion), `qsp_approximate`) | yes (closed synthesis engine) | yes (simulators and cloud backends) | yes (QSVT/HHL examples vs classical) | HHL, QSVT, VQLS, adiabatic | Poisson, LCHS, time-marching, CFD examples |
| UnitaryLab | yes ([algorithms manual](https://docs.unitarylab.com/en/docs/unitarylab-algorithms-user-manual/)) | yes (`.run()` builds circuits) | yes (numpy/torch backends) | partial — Schrödingerization `method='block'` is documented to fall back to classical computation ([manual page](https://docs.unitarylab.com/en/docs/unitarylab-algorithms-user-manual/schrodingerization/)) | HHL, QSVT-QLSA, VQLS, AQC | Schrödingerization heat/advection (trotter + classical) |
| Qrisp | yes ([BlockEncoding](https://qrisp.eu/reference/Block%20Encodings/BlockEncoding.html)) | yes | yes (`apply_rus`, sim + hardware) | yes ([inv() QSLP example](https://www.qrisp.eu/reference/Block%20Encodings/methods/inv.html) vs classical) | CKS, `inv()` | — |

oracq's own row for calibration: L1 yes (typed oracle interfaces, contracts),
L2 yes, L3 yes (reference executor plus backend paths), L4 yes for the
witnessed fixtures (self-verifying phase synthesis, κ=8/ε=1e-2 inversion
reached through the replaceable synthesizer at the phase and 2×2
block-encoding level), QLS yes (CKS, VTAA-CKS), ODE/PDE yes (LCHS,
Schrödingerization, CBMD, Carleman, QHAM/QFVM families).

One-sentence conclusion: templates, demonstrations, and validated solvers
are cleanly separable by the four levels; the reviewed documentation shows
full L1–L4 chains for QSVT/QLS at PennyLane (simulator), Classiq (closed
engine), Qrisp, and oracq, partial chains at Qualtran and pyLIQTR
(resource-estimation oriented), and a documented classical fallback at
UnitaryLab's Schrödingerization `block` option; ODE/PDE solver families
remain oracq's distinguishing coverage.
