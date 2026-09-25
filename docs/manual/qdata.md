# Quantum data structures and the KP recommendation system

**English** · [简体中文](../zh/manual/qdata.html)

{obj}`QVector <oracq.algorithms.input_model.qdata.QVector>` and {obj}`QMatrix <oracq.algorithms.input_model.qdata.QMatrix>` package the quantum data structures shared by [QFVM](qfvm.md) (arXiv:2102.03557 Eqs. 15–22 residual sum-of-squares tree) and the Kerenidis–Prakash recommendation system (arXiv:1603.08675 Thm 5.1 + Appendix A.1), with the entire data plane accessed through [QRAM pointers](qmem.md). The API is documented in [quantum data structures (qsample and sample-and-query)](../api/algorithms/input_model/qdata.rst).

## QVector: the qsample vector

A squared-norm binary tree: leaves store squared components and internal nodes cache RY rotation angle words (θ = 2·acos(√(S_left/S_node))). The quantum side prepares the normalized state level by level with two QRAM queries per level; signs are written via phase kickback on a 1-bit bank. On the classical side, a single-point update recomputes only the leaf-to-root path.

```{doctest}
>>> from oracq import QVector, simulate
>>> from oracq.algorithms.common.arithmetic import FixedFormat
>>> vector = QVector((3.0, 4.0), fmt=FixedFormat(8, 4), angle_width=12)
>>> state = simulate(vector.preparation().operation.program(), vector.snapshot())
>>> [round(a.real, 3) for _, a in sorted(state.amplitudes.items())]
[0.6, 0.8]
```

## QMatrix: the sample-and-query matrix

The entry bank supports queries in arbitrary superposition (addresses use a two-dimensional (row, column) pointer); one QVector tree per row gives Ũ: |i⟩|0⟩→|i⟩|Ā_i⟩; the row-norm root tree gives Ṽ: |0⟩|j⟩→|Ã⟩|j⟩. Entries are restricted to non-negative fixed-point values.

```{doctest}
>>> from oracq import QMatrix, simulate
>>> from oracq.algorithms.common.arithmetic import FixedFormat
>>> matrix = QMatrix([[0.5, 0.25], [0.5, 0.25]], fmt=FixedFormat(8, 4), angle_width=12)
>>> state = simulate(matrix.query().operation.program(),
...                  {"entries": matrix.snapshot()["entries"]}, initial={"address": 0b10})
>>> sorted(state.amplitudes.items())
[((2, 8), (1+0j))]
```

## The KP quantum recommendation system

{obj}`kp_recommendation <oracq.algorithms.qml.recommendation.kp_recommendation>` builds W = Ũ R₁ Ũ⁻¹ · Ṽ R₀ Ṽ⁻¹ per Lemma 5.3 of the paper (R₀/R₁ are reflections about the zero basis vectors) and runs phase estimation on it: cos(θᵢ/2) = σᵢ/‖A‖_F. σ̂ is estimated phase word by phase word and the flag flips against the threshold (the deterministic projection of Alg 2); after inverse phase estimation, (flag, item) is measured; the entry distribution of the flag=1 branch is exactly the recommendation sampling distribution. The mirror phase 2^p−t and t correspond to the same singular value (the eigenvalues of W come in pairs e^{±iθ}).

```{doctest}
>>> from oracq import KPRecommendationConfig, QMatrix, simulate, kp_recommendation
>>> from oracq.algorithms.common.arithmetic import FixedFormat
>>> matrix = QMatrix([[0.5, 0.25], [0.5, 0.25]], fmt=FixedFormat(8, 4), angle_width=12)
>>> result = kp_recommendation(matrix, 0, KPRecommendationConfig(precision=4, sigma=0.5))
>>> success, distribution = result.readout(simulate(result.operation.program(), result.memories()))
>>> round(success, 3), {k: round(v, 3) for k, v in sorted(distribution.items())}
(1.0, {0: 0.8, 1: 0.2})
```

For a rank-one matrix the recommendation distribution lands exactly on the principal singular vector v₁²=(0.8, 0.2) (the residual is the angle-quantization error). On matrices with small perturbations, the threshold σ filters out small-singular-value components, the success probability drops below 1, and the conditional distribution approaches v₁² — numerical validation is in `tests/core/test_recommendation.py`. Declared assumptions: non-negative fixed-point entries, angle-word quantization, and no amplitude amplification (the variable-time amplification of KP §6 is noted as future work).

## The QMem direct path of QFVM

`applications/qfvm_qmem.py` rewrites QFVM's data access onto the same syntactic sugar: the three conserved quantities merge into a (field, cell) state table {obj}`QMem(b, "state", shape=(3, n)) <oracq.infrastructure.qmem.QMem>`, neighbors cell±1 implement the periodic boundary through pointer arithmetic modulo 2^cell_width; the geometry table is addressed two-dimensionally by (slot, column); the residual state is prepared by a QVector tree. Modules declare QRAM-form resources directly, no longer going through abstract slots and {obj}`bind <oracq.infrastructure.linking.bind>`. The circuit semantics (reversible Roe arithmetic, nine-slot geometry, padding diagonal, in-place position permutation) match the existing `applications/qfvm.py` path; equivalence is cross-checked amplitude by amplitude on {obj}`simulate <oracq.infrastructure.execution.simulate>` by `tests/core/test_qfvm_qmem.py`. For the line-by-line walkthrough of the slot-binding main path see [QFVM input models and solver replacement](qfvm.md#qmem-direct-parallel-path).
