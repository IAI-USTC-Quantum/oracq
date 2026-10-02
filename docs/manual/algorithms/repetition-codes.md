# Repetition Codes

**English** · <a href="../../../zh/manual/algorithms/repetition-codes.html">简体中文</a>

> Category C1 · Module [`oracq.algorithms.qec.error_correction`](../../api/algorithms/qec/error_correction.rst) · Stage V1

## Overview

Encodes the one-bit logical state $\alpha|0\rangle + \beta|1\rangle$ into the three-physical-bit repetition code $\alpha|000\rangle + \beta|111\rangle$ and coherently recovers from any single error of the specified type. `error="bit"` corrects a single X error; `error="phase"` is the H-conjugated version of the same code and corrects a single Z error. Encoding and recovery are both purely unitary circuits, with no measurement or reset, so the logical amplitudes (including relative phases) are preserved pointwise across the encode–error–recover cycle.

## Interface and input model

```python
repetition_encode(*, error="bit")
repetition_recover(*, error="bit")
```

API entry points: {obj}`repetition_encode <oracq.algorithms.qec.error_correction.repetition_encode>`, {obj}`repetition_recover <oracq.algorithms.qec.error_correction.repetition_recover>`

- `error`: `"bit"` or `"phase"`; the encoder and the recoverer must take the same value; other values raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time.
- The input model is CP: there is no oracle input; the error type and the code structure are given as classical parameters.
- The two entry points share registers: `target: Bits(1)` (the logical bit) and `syndrome: Bits(2)`; encoding requires the `syndrome` input to be zero. The three physical bits are laid out as `target`, `syndrome[0]`, `syndrome[1]`.

Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"repetition_encode"` / `"repetition_recover"` |
| `error_kind` | `"bit"` or `"phase"` |
| `syndrome_policy` | recoverer only: `"retained; host reset required before reuse"` |

## Implementation notes

The encoder applies two CNOTs (`target` → `syndrome[0]`, `target` → `syndrome[1]`) to obtain the three-bit repetition; the `"phase"` variant then applies H to all three physical bits, moving the codewords into the phase-flip basis. The recoverer first XORs `target` into the two syndrome bits to obtain the syndrome: when the error is on the `target` bit, both bits are 1, and a controlled X (controlled on `syndrome == 3`) flips `target` back to the logical value; when the error is on one of the syndrome bits, the syndrome is not all 1s, and `target` was never affected. The `"phase"` variant first applies H back to the bit-flip basis before recovery, reusing the same logic.

Design constraint: after recovery the syndrome deterministically retains the error location (0 when there is no error), but it is **not uncomputed back to clean** — the host must reset before reusing the registers, which is exactly what the `syndrome_policy` attribute states. Applicability boundary: only the correction of at most one error consistent with `error` is guaranteed; two or more same-type errors, or mixed errors, are beyond the three-bit code's ability.

## Validation approach

Category C1 (exact discrete semantics; acceptance criterion in `../../development/validation-plan.md` §2): the composite action of encode–error–recover must equal the identity operator pointwise on arbitrary logical amplitudes. Evidence:

- Structure: construction and attribute assertions in `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`; `test_bad_inputs_fail_at_generation` covers the generation-time violation of {obj}`repetition_encode(error="unknown") <oracq.algorithms.qec.error_correction.repetition_encode>`.
- Numerical: `AlgorithmExpansionTests.test_repetition_codes_preserve_arbitrary_logical_amplitudes` — for both `error` kinds × three physical locations, an arbitrary logical state is prepared with Ry(0.73)/Rz(0.29); after injecting a single error and recovering, the `target` amplitudes are exactly restored to $e^{-0.145j}\cos 0.365$ and $e^{0.145j}\sin 0.365$ (places = 10), and the syndrome converges to a single deterministic value.
- Binding: this algorithm has no independent binding witness (no open declaration entry point).

## Known gaps and planned stages

No known gaps (the gap column of the validation matrix reads —), stage V1.

## Related links

- Source: `src/oracq/algorithms/qec/error_correction.py`
- API reference: [Repetition codes and error recovery](../../api/algorithms/qec/error_correction.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_misc_algorithms.py` (the misc_algorithms group), all executed on a real backend.

**Experiment design**: the arbitrary logical state $\mathrm{Ry}(0.73)\mathrm{Rz}(0.29)\lvert 0\rangle$ (numpy independently provides the exact amplitudes) goes through "encode–inject–recover" for both `error` kinds × 4 injection locations (no error / target / syndrome[0] / syndrome[1]): pointwise restoration of the target amplitudes and the deterministic syndrome value (location mapping 0/3/1/2) are checked; in addition, for the error-free encode–recover composite, the syndrome = 0 input block and the leakage are extracted via `originir-ext + UniQC Circuit.to_matrix`. Backend paths: `reference`, `rir-pysparq`, `adapter-pysparq`, `originir-ext` (3-qubit full-amplitude cross-check).

**Key metrics**:

| Case | Scale | Paths | Amplitude-restoration error | syndrome | unitary block error / leakage |
|---|---|---|---|---|---|
| repetition-bit-flip-injection | 3 bits × 4 locations | four paths | 0 | all deterministic and equal to the error location ✓ | — |
| repetition-phase-flip-injection | 3 bits × 4 locations | four paths | ≤ 4.6e-16 | all deterministic and equal to the error location ✓ | — |
| repetition-bit-encode-recover-unitary | 3 bits | to_matrix | — | — | 0 / 0 |
| repetition-phase-encode-recover-unitary | 3 bits | to_matrix | — | — | 4.4e-16 / 2.4e-17 |

**Reproduction**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

Artifacts: `out/verification/misc_algorithms.json` (24 cases all passing; this page corresponds to the four `repetition-*` cases).
