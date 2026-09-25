# Repository reorganization and the algorithm expansion panel

**English** · [简体中文](../zh/development/roadmap.html)

This round first reorganizes code and documentation, then expands a batch of
algorithms that can actually generate circuits. All work happens inside the
oracq repository, preserving existing uncommitted changes.

| Stage | Scope | Status |
|---|---|---|
| R1 | Layering into infrastructure / algorithms / applications; split up elementary, solvers, differential while keeping legacy import compatibility | Done |
| R2 | Sphinx + MyST + API documentation; manuals separated from tutorials, historical records archived | Done |
| R3 | Oracle queries, Fourier arithmetic, search and amplitude amplification, estimation, variational algorithms, quantum walks, number theory, simple error correction | Done |
| R4 | Mathematical witnesses for new algorithms, regression of old cases, real-backend validation, strict documentation builds, packaging and migration acceptance | Done |

Algorithm files are named by purpose: oracle_algorithms.py, fourier.py,
search.py, estimation.py, variational.py, walks.py, number_theory.py,
error_correction.py, plus the existing hamiltonian.py, qlss.py, lchs.py,
schrodingerization.py, cbmd.py, carleman.py, qham.py.

The first batch of new implementations covers Bernstein–Vazirani; Simon
sampling with GF(2) post-processing; QFT addition; general amplitude
amplification; Hadamard/Swap tests; standard amplitude estimation; MaxCut
QAOA; VQE measurement circuits; parameterized ansatz; the periodic coined
walk; modular multiplication for order finding and factor post-processing;
and three-bit bit/phase-flip encoding and decoding. Gate-level output and
each algorithm's scope of applicability are recorded separately — the small
modular-multiplication-table implementation is not treated as scalable Shor
arithmetic.

The documentation states the problem, inputs, steps, and results in complete
Chinese sentences. Tutorials are organized around user tasks; specifications
give explicit contracts; the API is generated from source. Builds must be
free of Sphinx warnings, and tutorials must run. Existing stage reports moved
into archive, keeping their references without mixing them into the main
reading path.
