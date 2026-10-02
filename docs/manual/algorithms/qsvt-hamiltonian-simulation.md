# QSVT Hamiltonian Simulation

**English** · <a href="../../../zh/manual/algorithms/qsvt-hamiltonian-simulation.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.common.qsvt`](../../api/algorithms/common/qsvt.rst) · Stage V2

## Overview

Construct a block encoding of $e^{itA/\alpha}$ from a block encoding of a Hermitian matrix $A$ ($\alpha$ is the input BE normalization). The implementation follows the QSVT framework (Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838), cited in the module docstring) and the Jacobi–Anger expansion of $e^{itx}$, $e^{itx} = \cos(tx) + i\sin(tx)$: the even branch approximates $\cos(tx)$ and the odd branch approximates $\sin(tx)$; each branch goes through QSP phase synthesis (see [QSP phase synthesis](qsp-phase-synthesis.md)) and real-part extraction via $(U_\Phi + U_{-\Phi})/2$, and the two are finally combined in an LCU with coefficients 1 and $i$. The zero-signal block of the returned BE is approximately $e^{itA/\alpha}/\text{sim\_scale}$.

## Interface and input model

```python
qsvt_hamiltonian_simulation(a, t, *, error=0.01)
```

API entry: {obj}`qsvt_hamiltonian_simulation <oracq.algorithms.common.qsvt.qsvt_hamiltonian_simulation>`

- `a`: a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`, the block encoding of the simulated Hamiltonian (input model BE).
- `t`: the evolution time, which must be a nonzero finite real number.
- `error`: the polynomial approximation error, which must lie in $(0, 1)$.

Returns a `BlockEncoding`. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"qsvt_hamiltonian_simulation"` |
| `be_alpha` | normalization of the output BE (an LCU of the two branches with coefficients 1 and $i$; 2.0) |
| `time` / `error` | echo of the call parameters |
| `qsp_degree` | the Bessel truncation degree $K$ (the degrees actually synthesized for the two branches are the even / odd degrees not exceeding $K$) |
| `sim_scale` | the output scaling $2s$ (the zero-signal block is approximately $e^{itA/\alpha}/\text{sim\_scale}$) |

## Implementation notes

The truncation degree is adaptive: $J_k(t)$ is computed by power series (pure Python), $K$ is the smallest degree for which the sum of the absolute values of the Bessel tail terms does not exceed `error`/4, capped at an upper bound of 40. The two branches share the uniform scaling $s = 1.5\max(\lVert f_c\rVert_\infty, \lVert f_s\rVert_\infty, 10^{-3})$, leaving headroom for the imaginary-part completion, and $\text{sim\_scale} = 2s$. A parity-check violation (the cos branch must be even, the sin branch odd) raises {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`; the candidate families of imaginary-part completions (a constant plus a single even power for the cos branch, a single odd power for the sin branch) are tried one by one, and the first synthesizable candidate is taken.

Applicability boundary: $t = 0$ is rejected (a trivial case); for very large $|t|$, the required truncation degree reaches the cap of 40 before the error condition does, in which case the `error` bound may not be met without an error being raised, and the caller must verify this on its own. The output is scaled by `sim_scale` and must be divided back out after readout.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../../development/validation-plan.md` §2): the acting operator must approximate the target continuous function within tolerance. Three layers of evidence:

- Structure: `tests/core/test_qsvt.py:TransformWitnessTests.test_transform_input_validation` (negative cases such as t = 0 and out-of-range error) and `PhaseSynthesisTests.test_input_validation` (the input validation of the underlying synthesis).
- Numerical: `tests/core/test_qsvt.py:TransformWitnessTests.test_hamiltonian_simulation_block` — a 2×2 diagonal matrix (eigenvalues 0.5 and $-0.25$, $\alpha = 0.5$, spectral variables 1.0 and $-0.5$), t = 0.7, error = 0.01: the zero-signal block amplitudes are read out column by column and cross-checked against $e^{itx}/\text{sim\_scale}$, delta = 5e-4.
- Binding: this algorithm has no independent binding witness (the input already requires a concrete BE).

## Known gaps and planned stages

Consistent with the `qsvt.py` row of the validation coverage matrix: convergence scans are missing (batch curves of the error versus t / degree are not automated, pending the stage V2 convergence-scan framework); the ill-conditioned negative cases (degree ≳ 16 or too-low precision) are partially covered by the parametrized input_validation, and a separate batch parametrization of "ill-conditioned inputs must raise" is still to be added.

## Related links

- Source: `src/oracq/algorithms/common/qsvt.py`
- API reference: [QSVT standard transforms](../../api/algorithms/common/qsvt.rst)
- Same-family pages: [QSP phase synthesis](qsp-phase-synthesis.md), [QSVT matrix inversion](qsvt-matrix-inversion.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_hamiltonian.py` (real-backend execution, no mock substitutes). Experiment design: the input BE is a non-diagonal Hermitian matrix $A = \begin{pmatrix} 0.5 & 0.2 \\ 0.2 & -0.3 \end{pmatrix}$ encoded by {obj}`matrix_pauli_encoding <oracq.algorithms.input_model.block_encoding.matrix_pauli_encoding>` ($\alpha = 0.7$), with $t \in \{0.7, 2.0\}$ and `error = 0.01`; on the three paths reference / rir-pysparq / originir-ext, the full $2\times2$ zero-signal block is read out column by column and cross-checked against $e^{itA/\alpha}/\text{sim\_scale}$ computed by `scipy.linalg.expm`. The implementation error (circuit block vs the $e^{itx}$ reference) and the method error (the Jacobi–Anger truncation-tail bound $2\sum_{j>K}|J_j(t)|$, evaluated independently by `scipy.special.jv`) are reported separately.

| Case | Scale | Backend paths | Metric | Value |
|---|---|---|---|---|
| `qsvt-hamsim-t0.7` | $K = 4$, sim_scale ≈ 3.0 | reference, rir-pysparq, originir-ext | impl_error | 2.72e-5 |
| same as above | — | — | method tail bound / total error bound | 9.10e-5 / 1.18e-4 |
| `qsvt-hamsim-t2.0` | $K = 6$ | reference, rir-pysparq, originir-ext | impl_error | 5.42e-5 |
| same as above | — | — | method tail bound / total error bound | 4.00e-4 / 4.54e-4 |

The total error bounds at the two time points (1.2e-4 and 4.5e-4) are both far below the requested `error = 0.01`, and the implementation error is significantly below the method error, indicating that the truncation-degree choice is the dominant contributor to the error while phase synthesis and circuit assembly sit at the noise level.

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_hamiltonian.py
```

Artifacts: `out/verification/hamiltonian.json`.
