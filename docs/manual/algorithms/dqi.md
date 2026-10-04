# Decoded Quantum Interferometry

**English** · <a href="../../zh/manual/algorithms/dqi.html">简体中文</a>

> Category C4 · Module [`oracq.algorithms.optimization.dqi`](../../api/algorithms/optimization/dqi.rst) · Stage V1

## Overview

DQI solves GF(2) max-XORSAT: given a constraint system $Bx=v$ ($B$ as the sparse row representation of an $m\times n$ 0/1 matrix), find the assignment that satisfies the most constraints. The implementation follows the circuit skeleton of Jordan et al. 2024 ([arXiv:2408.08292](https://arxiv.org/abs/2408.08292), the single-weight version of figure 4): prepare a Dicke state of weight $l$ in the $m$-qubit error register, apply the right-hand-side phase $(-1)^{v\cdot y}$, compute $B^T y$ reversibly into the $n$-qubit syndrome register, use a reversible classical decoder to uncompute the error register back to $\lvert 0\rangle$, and finally apply a Hadamard transform to the syndrome and measure. After postselecting the branch with zero error, the probability of measuring an assignment $x$ is proportional to $K_l(u(x))^2$ — $u(x)$ is the number of unsatisfied constraints and $K_l$ is the Krawtchouk polynomial, so the sampling is biased toward assignments satisfying more constraints.

## Interface and input model

```python
dqi(instance, decoder, *, weight)
XorSatInstance(rows, rhs, num_variables)
abstract_decoder(name, syndrome_width, error_width)
table_decoder(syndrome_width, error_width, table, *, name=None)
bruteforce_decoder(instance, *, max_weight=None, name=None)
dicke_state(m, weight)
```

API entry points: {obj}`dqi <oracq.algorithms.optimization.dqi.dqi>`, {obj}`XorSatInstance <oracq.algorithms.optimization.dqi.XorSatInstance>`, {obj}`abstract_decoder <oracq.algorithms.optimization.dqi.abstract_decoder>`, {obj}`table_decoder <oracq.algorithms.optimization.dqi.table_decoder>`

- `instance`: an {obj}`XorSatInstance <oracq.algorithms.optimization.dqi.XorSatInstance>`, the sparse row representation of the constraints (each row is the indices of the variables participating in that constraint, with no duplicates within a row; `rhs` takes 0 or 1). `num_variables` is the syndrome bit width $n$, and the constraint count `num_constraints` is the error bit width $m$. The constraint instance is classical data parameterized directly, so the input model is CP (it encodes the constraints).
- `decoder`: a {obj}`DecoderOracle <oracq.algorithms.optimization.dqi.DecoderOracle>` with semantics `|syndrome, error> → |syndrome, error XOR D(syndrome)>`, an FO-class open input (a reversible classical function, the `reversible_function` paradigm). {obj}`abstract_decoder <oracq.algorithms.optimization.dqi.abstract_decoder>` declares the slot and binds it in batches via {obj}`bind <oracq.infrastructure.linking.bind>`; {obj}`table_decoder <oracq.algorithms.optimization.dqi.table_decoder>` uses an explicit lookup table as the gate witness (unlisted syndromes map to the zero error); {obj}`bruteforce_decoder <oracq.algorithms.optimization.dqi.bruteforce_decoder>` enumerates all $2^m$ error patterns and returns the lightest error of weight at most `max_weight`; it only accepts $m \le 16$.
- `weight`: the Dicke-state weight $l$, range 0..m; the decoding radius must cover this weight, and the decoder bit widths must match the instance's n/m.

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with registers `error`(m) and `syndrome`(n). Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` / `field` | `"dqi"` / `"GF(2)"` |
| `num_constraints` / `num_variables` | m / n |
| `dicke_weight` | l |
| `readout_register` | `"syndrome"` |
| `postselection` | `"error_zero"` |

{obj}`dicke_state(m, weight) <oracq.algorithms.optimization.dqi.dicke_state>` is exported separately and prepares $\lvert D_l^m\rangle$ (explicit-amplitude construction, $m \le 16$); `XorSatInstance.satisfied_count(assignment)` provides the satisfied-count statistic for classical-side cross-checks.

## Implementation notes

Generation chain: Dicke preparation ({obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>` with explicit amplitudes, using a zero-width view of the syndrome as work) → right-hand-side phase (a Z on each error bit with `rhs[i] = 1`) → syndrome computation (an `xor(error_i → syndrome_j)` for each variable index of each constraint, the total being $B^T y$) → decoder call → Hadamard on the syndrome. On the branches where decoding succeeds, error returns to $\lvert 0\rangle$; postselection (discarding branches with nonzero error) is performed by the caller according to the `postselection` attribute.

Design decision: the decoder is part of the input model rather than an internal algorithmic detail — syndromes that cannot be decoded may be mapped to any error pattern (the corresponding branches are eliminated in postselection), and efficient classical decoding (e.g. a reversible implementation of belief propagation) plugs in via `abstract_decoder` plus `bind`, consistent with the repository's three-layer open-declaration paradigm. Applicability boundary: only GF(2) is currently supported; the GF(q) case needs a q-ary discrete Fourier transform and generalized Dicke states, with registers organized into $\log_2 q$-bit sub-units, which is left as an extension; the explicit-amplitude Dicke state and brute-force decoding only serve small-instance witnesses with $m \le 16$, and larger instances should switch to dedicated Dicke circuits (e.g. the O(l·m) construction of Bartschi–Eidenbenz, plugged in via an open declaration through `state_prep_isometry`) and efficient decoders.

## Validation approach

Category C4 (acceptance criteria in `../../development/validation-plan.md` §2: strictly better than the random baseline, and reaching the known optimum / theoretical fraction on small instances). Three layers of evidence in `tests/core/test_dqi.py:DqiTests`:

- Structure: `test_invalid_inputs_fail_at_generation` covers 12 classes of generation-time violations (out-of-range/duplicate indices, rhs other than 0/1, empty constraints, out-of-range Dicke weight, out-of-range decoding table, bit-width mismatch, wrong decoder type, negative weight, etc.).
- Numerical: `test_planted_instance_beats_random_guessing` — a planted instance with 7 constraints and 3 variables (the right-hand side planted from the assignment x\* = 0b101); at weight 1 the syndrome distribution is cross-checked pointwise against the Krawtchouk closed form $K_l(u(x))^2$ (places = 10), and the expected satisfied count 6.5 is strictly better than the random baseline 3.5; `test_identity_instance_with_weight_two` — an m = n = 5 instance with B the identity at weight 2, where after the distribution cross-check the expected satisfied count returns to the random baseline m/2 (the assertion uses `assertGreaterEqual` with a 1e-9 tolerance); `test_dicke_state_weight_and_uniformity` — $\lvert D_2^5\rangle$ has exactly C(5,2) = 10 equal-amplitude nonzero components (1/√10, places = 12) and the work returns clean.
- Binding: `test_abstract_decoder_binds_to_witness` — the abstract decoder is the program's only unresolved slot; after `bind` binds the brute-force implementation there are no unresolved declarations, and the syndrome distribution agrees pointwise with the gate witness (places = 10).

The shared `check_distribution` also pins down that the error register is deterministically restored to zero when decoding succeeds (P(error = 0) = 1, places = 12).

## Known gaps and planned stages

Consistent with `validation-coverage.md`: no outstanding gaps, stage V1. An identity boundary issue that previously arose — when the decoder provides no advantage the expectation exactly equals the random baseline, so a strict-inequality assertion fails at the boundary — has been fixed per the determinism policy of validation-plan §4 into ≥/≤ with a numerical tolerance (`test_identity_instance_with_weight_two` is the witness of that fix).

## Related links

- Source: `src/oracq/algorithms/optimization/dqi.py`
- API reference: [DQI decoded quantum interferometry optimization](../../api/algorithms/optimization/dqi.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical experiments are in `tests/verification/verify_misc_algorithms.py` (misc_algorithms group), all executed on real backends.

**Experimental design**: (a) a planted instance (7 constraints, 3 variables, right-hand side planted from $x^* = \mathtt{0b101}$, Dicke weight 1, brute-force decoder); the syndrome distribution against the paper's Krawtchouk closed form $K_l(u(x))^2$ (independently computed with `math.comb`), and the optimization quality against all 8 assignments of a classical brute-force enumeration; (b) the abstract decoder, after `bind` binds the brute-force witness, cross-checked amplitude by amplitude against the direct witness; (c) the amplitude uniformity of the Dicke state $\lvert D_2^5\rangle$. Backend paths: `reference`, `rir-pysparq`, `adapter-pysparq`, `originir-ext` (a 10-qubit instance, distribution cross-check across four paths).

**Key metrics**:

| Case | Scale | Path | Metric | Value |
|---|---|---|---|---|
| dqi-planted-krawtchouk | m = 7, n = 3, 10 qubits | four paths | distribution TVD / pointwise error | 4.0e-16 / 6.7e-16 |
| same | — | — | P(error = 0) (decoding restored to zero) | 1.0000 |
| same | — | — | expected satisfied count (random baseline 3.5) | 6.5000 |
| same | — | — | probability-peak assignment = brute-force optimum (satisfying 7/7) | $x^*$ = 5 ✓ |
| dqi-abstract-decoder-bind | same | rir-pysparq | max amplitude deviation of bind vs the direct witness | 0 |
| dicke-state-uniformity | m = 5, l = 2 | reference + originir-ext | support / amplitude error | 10 basis states / 5.6e-17 |

**Reproduce**:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

Artifact: `out/verification/misc_algorithms.json` (all 24 cases pass; this page corresponds to the three `dqi-*` and `dicke-state-*` cases).
