# Oblivious Amplitude Amplification

**English** · <a href="../../../zh/manual/algorithms/oblivious-amplification.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.common.transforms`](../../api/algorithms/common/transforms.rst) · Stage V1

## Overview

Applies oblivious amplitude amplification (OAA) to a block encoding $U$: an iterative boost of the zero-signal block amplitude that does not rely on knowledge of the initial state. With $\Pi = |0\rangle\langle 0|_{\text{signal}} \otimes I_{\text{target}}$ and $R = I - 2\Pi$, the operator with $m$ iterations is

$$
W_m \;=\; U\,\big(R\,U^\dagger R\,U\big)^{m}
$$

($m=1$ is the standard three-query form $U\,R\,U^\dagger R\,U$ from the literature; the two global minus signs of $R$ relative to $2\Pi - I$ cancel). For an input whose zero-signal block is $V/2$ ($V$ a partial isometry), the zero-signal block of $W_m$ is $(-1)^m \sin\!\big((2m+1)\theta\big)\, V$ ($\sin\theta = 1/2$) — the magnitude evolves along the standard OAA narrative $\sin\theta \to \sin 3\theta \to \cdots$, and a single iteration restores the amplitude to $O(1)$. General block encodings satisfy the Chebyshev identity $\Pi W_1 \Pi = B\,(4B^\dagger B - 3I)$. This construction is the standard companion component of the qubitization framework (framework reference: Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838); for walk-operator background see [qubitization walk](qubitization-walk.md)).

> Historical defect record: versions before 2026-09 assembled $[R\,U^\dagger R\,U]^m$ (missing the trailing $U$), degenerating the zero-signal block to $2B^\dagger B - I$ and performing no amplification; it was found by a validation round and fixed to the formula above, then pinned by the "library operator = standard sequence" regression case below — see [Numerical validation](#numerical-validation).

## Interface and input model

```python
oblivious_amplification(a, iterations=1)
```

API entry point: {obj}`oblivious_amplification <oracq.algorithms.common.transforms.oblivious_amplification>`

- `a`: a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`, the block encoding being amplified (input model BE).
- `iterations`: the iteration count, default 1.

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` (a bare operation, not a `BlockEncoding`), with registers `target` (width `a.width`) and `signal` (width `a.signal_qubits`). Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"oaa"` |

BE normalization is not propagated through attributes; the semantics of the amplified amplitude are interpreted by the caller.

## Implementation notes

The program's time order is: {obj}`invoke <oracq.algorithms.input_model.oracles.invoke>` first mounts the input BE (the initial $U$), then each iteration makes four structural calls: {obj}`reflect_zero(signal) <oracq.algorithms.input_model.block_encoding.reflect_zero>` (default `positive=False`: the signal-zero branch takes $-1$ and the other branches take $+1$, i.e. $R = I - 2\Pi$) → adjoint `invoke` → `reflect_zero(signal)` again → the trailing `invoke`. The iteration count is expressed through the IR's {obj}`Repeat <oracq.infrastructure.ir.Repeat>` structure and is not unconditionally expanded during generation or text serialization.

Applicability boundary: the input must be a block encoding. The amplification guarantee relies on the structural assumption that "the zero-signal block is close to a $1/2$-scaled partial-isometry operator", which this function does not check; iterating a block encoding whose normalization is far from $1/2$ does not converge in the OAA sense and must first be adjusted by scaling or LCU composition.

## Validation approach

Category C2 (approximate continuous semantics; acceptance criterion in `../../development/validation-plan.md` §2). Current evidence (matching the `transforms.py` row of the validation coverage matrix):

- Structure: `tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` and `test_trotter_only_input_does_not_need_block_encoding` cover the protocols and register contracts of the same module's assembly chain (including same-kind assertions that the Repeat structure is preserved).
- Numerical: `tests/core/test_qsvt.py:PhaseSynthesisTests.test_convention_matches_qsvt_sequence` pins the signal/reflection convention shared by this module (delta = 1e-10); OAA and the [qubitization walk](qubitization-walk.md) share the `invoke` + `reflect_zero` components and are thereby covered indirectly.
- Binding: no independent binding witness (the input already requires a concrete BE).

## Known gaps and planned stages

No known gaps (matching the `transforms.py` row of the validation coverage matrix), stage V1. The amplification semantics have an independent numerical witness: the $V/2$ fixture goes from magnitude $0.5 \to 1.0$ in one iteration ($\sin 3\theta$); the Chebyshev identity for general blocks and the "library operator = standard sequence" regression pin are in [Numerical validation](#numerical-validation).

## Related links

- Source: `src/oracq/algorithms/common/transforms.py`
- API reference: [Matrix transform sequences](../../api/algorithms/common/transforms.rst)
- Same-family pages: [qubitization walk](qubitization-walk.md), [QSVT phase sequences](qsvt-sequence.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_fourier.py` (real-backend execution, no mock substitutes). The input BE is a small gate-level bound instance. Assembly correctness: the full unitary for iterations $it = 1, 2, 3$ is extracted via OriginIR-ext + UniQC `to_matrix` and compared element by element with $U\,[R\,U^\dagger R\,U]^{it}$ ($R = I - 2\Pi$). Behavioral measurement: a minimal gate-level instance whose zero-signal block is exactly $X/2$ is built (the "$V/2$" fixture with $\sin\theta = 1/2$); after one iteration the zero-signal block magnitude recovers to $1.0$ ($-\!X$, i.e. $(-1)^1\sin 3\theta \cdot X$), and $it = 2, 3$ give $0.5$ and $0.5$ ($\sin 5\theta$, $\sin 7\theta$) — the magnitudes agree step by step with the standard three-query OAA narrative, and the global sign $(-1)^{it}$ is recorded as an informational metric (observable at the unitary level, indistinguishable at the measurement level).

| Case | Scale | Backend path | Metric | Value |
|---|---|---|---|---|
| `oaa-unitary-it1` | 8-dimensional, 1 iteration | originir-ext + UniQC to_matrix | max_error | 2.3e-16 |
| `oaa-unitary-it2` | 8-dimensional, 2 iterations (Repeat structure) | originir-ext + UniQC to_matrix | max_error | 4.4e-16 |
| `oaa-unitary-it3` | 8-dimensional, 3 iterations (Repeat structure) | originir-ext + UniQC to_matrix | max_error | 4.6e-16 |
| `oaa-half-block-it1` | zero block = $X/2$ fixture | originir-ext + UniQC to_matrix | zero block vs $-\sin 3\theta\,X$ | 6.1e-17 |
| `oaa-half-block-it2/3` | same as above (Repeat structure) | originir-ext + UniQC to_matrix | zero block vs $\pm\sin 5\theta\,X$, $\sin 7\theta\,X$ | 8.9e-16 / 1.1e-15 |

Reproduction command:

```bash
PYTHONPATH=src <interpreter with pysparq+uniqc> tests/verification/verify_fourier.py
```

Artifacts: `out/verification/fourier.json`.

### Supplementary validation: Chebyshev identity for general blocks and the standard-sequence regression pin (verify_hamiltonian.py)

`tests/verification/verify_hamiltonian.py` (27 cases, all PASS) provides two complementary pieces of evidence that corroborate the measurements of the fourier group above: first, for a **general non-isometric block** ($2\times2$ matrix BE, $B = \Pi U \Pi = A/\alpha$) it verifies the Chebyshev amplification identity of the library operator $W_1 = U\,R\,U^\dagger R\,U$, $\Pi W_1 \Pi = B\,(4B^\dagger B - 3I)$, on the three paths reference / rir-pysparq / originir-ext (max_error 3.6e-16); second, it independently assembles the literature's standard sequence from the library's public components (`invoke` + `reflect_zero` + `adjoint`) and exactly recovers $-V$ on the zero-signal block of the $V/2$ fixture ($V = R_z(0.4)R_y(0.9)$) (max_error 1.1e-16), with the target amplitude amplified from 0.4502 to 0.9004 (exactly ×2.0), and **the library operator and the script-assembled standard sequence agree amplitude by amplitude on all three paths** (library_vs_script_error = 0.0) — this is the regression pin for "library = standard three-query sequence", preventing a regression to the historical defect that missed the trailing $U$ (that defect was found and fixed by a validation round; see the historical record in the overview above).

| Case | Scale | Backend paths | Metric | Value |
|---|---|---|---|---|
| `oaa-iterate-block-identity` | general block, 1 iteration | reference, rir-pysparq, originir-ext | Chebyshev identity max_error | 3.6e-16 |
| `oaa-standard-sequence-amplification` | $V/2$ fixture, $URU^\dagger RU$ | reference, rir-pysparq, originir-ext | max_error vs $-V$ | 1.1e-16 |
| same as above | — | — | per-amplitude deviation, library vs script assembly | 0.0 |
| same as above | — | — | amplitude amplification (×2.0) | 0.4502 → 0.9004 |

Reproduction command:

```bash
PYTHONPATH=src <interpreter with pysparq+uniqc> tests/verification/verify_hamiltonian.py
```

Artifacts: `out/verification/hamiltonian.json`.
