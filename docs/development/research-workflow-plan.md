# Algorithm research workflow improvements

**English** · <a href="../zh/development/research-workflow-plan.html">简体中文</a>

The core audience is researchers who implement, validate, and compare quantum
scientific-computing algorithms from papers. Implementation advances in
dependency order; each stage synchronizes specifications, examples, and
validation results.

| Stage | Deliverables | Status |
|---|---|---|
| 1 | Single-input adaptation, structured diagnostics, naming-compatible migration of solvers and contracts | Implemented, regressions pass |
| 2 | Partially bound contracts, structured binding reports, resource and capability regressions | Implemented, regressions pass |
| 3 | Compact rotation counting, QRAM actual-argument accounting, an open-oracle call ledger | Implemented, regressions pass |
| 4 | Implementation-selection cases, Roe/QFVM research reports, traceable validation | Four-way cross-checks of 9 new cases and 193 domain cases pass |
| 5 | Algorithm-author tutorials, paper mainline and evidence updates, full acceptance | Done |

Reproduced issues: a single generation repeatedly called provider methods;
open estimation still reported a closed error; rotation counts materialized
with Repeat; QRAM queries merged by internal formal parameters instead of the
entry actual arguments.

RIR keeps the concrete bit widths and assembly constants; what stays open is
the module implementation. Changing the public layout or alpha requires
explicit adaptation or regeneration. Python providers do not enter the
serialization result; analysis reports and RIR are stored separately. Unknown
oracle costs and private workspaces must not be booked as zero cost.

Core acceptance uses unittest, the schema, ruff, and mypy; semantic execution
uses real uniqc/pysparq. Sphinx HTML and doctest treat warnings as errors,
and the paper is built with LaTeX to check citations and figures. Generated
artifacts are written only into out; neighboring backend repositories are not
modified, and nothing is committed or pushed.

## Native validation environment

This round uses `out/native-validation/bin/python`, created and installed via
uv. PySparQ 0.1.1 on PyPI has no native RIR interface yet, so a snapshot was
extracted from QRAM-Simulator's committed revision
`4e4c9f16fe3e912a828915231ecb4a6663048bd2` into `out/native-sources/` in this
repository and built with uv as `0.0.0+g4e4c9f16`. The corresponding UniQC
distribution is `unified-quantum==0.1.1`, with the C++ simulator
`uniqc-cppsimulator==1.0.1`. Uncommitted changes in neighboring repositories
are not loaded.

The new workflow cases compare gate tables, QRAM tables, and arithmetic at
width=2,3,4 with repetitions=3: the maximum amplitude error of the four actual
execution paths against an independent reference is `2.5e-16`. The detailed
report lives at `out/research-workflow-native/report.json`.

QMatrix native tests are layered by execution budget: the public registers of
the 10-bit angle-word row preparation take 14 bits and the private workspace
for address computation takes 16 bits, 30 bits in total; two real PySparQ
paths are kept for cross-checking. A new 4-bit angle-word case, 24 bits in
total, completes a four-way cross-check with the reference executor, two
PySparQ paths, and OriginIR. This layering was confirmed by the user; no mock
substitutes or skips were used.

## Acceptance results

- `tools/check_project.py --docs`: all checks pass; core and schema total 391
  tests and 278 subtests pass.
- Real backend `tests/integration`: 41 pass; the original 10-bit QMatrix case
  is retained.
- `tools/run_verification.py`: six groups — oracles, blockencoding,
  arithmetic, mathfunc, qham_qfvm, ode — total 193 cases pass. This result
  does not replace the version coverage of the historical thirteen-group
  reports.
- New implementation-comparison cases: 9 configurations, all four execution
  paths pass, maximum amplitude error `2.5e-16`.
- ruff, mypy, Sphinx HTML, and 55 doctests pass; the paper is generated at
  `out/paper/main.pdf`.

Workflow checks are recorded in `out/checks/research-workflow/report.json`;
the native execution logs are `out/research-native-integration.log` and
`out/research-native-verification.log`. The six JSON reports produced by the
latter share the same source fingerprint, verification-script fingerprint,
and backend versions. The paper mainline is `paper/main.tex`; the old
Overleaf copy was not synchronized and is not a delivery entry for this
round.
