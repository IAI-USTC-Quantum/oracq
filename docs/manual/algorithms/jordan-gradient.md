# Jordan Quantum Gradient Estimation

**English** · <a href="../../zh/manual/algorithms/jordan-gradient.html">简体中文</a>

> Category C3 · Module [`oracq.algorithms.optimization.gradient`](../../api/algorithms/optimization/gradient.rst) · Stage V1

## Overview

Estimates the gradient $\nabla f$ of a $d$-dimensional real function $f$ at a given point. Implementation basis: Jordan 2005 (PRL 95, 050501), with the scaling convention matching the generalization of Gilyén–Arunachalam–Wiebe 2019: the phase oracle implements

$$
O\,|x\rangle = e^{2\pi i \cdot N \cdot f(x)}\,|x\rangle, \qquad N = 2^{m},
$$

where each coordinate takes an $m$-bit fixed-point grid (the integer value $k$ corresponds to the grid point $x = k/N$). When $f$ is approximately linear on the grid, a single oracle call writes $N \cdot \partial f/\partial x_i$ into the Fourier basis of the $i$-th coordinate register; after a per-coordinate inverse QFT all $d$ components are read out at once, while a classical deterministic evaluation of the same gradient requires $O(d)$ function queries.

## Interface and input model

```python
gradient_estimation(oracle, *, dimension, grid_bits)
```

API entry point: {obj}`gradient_estimation <oracq.algorithms.optimization.gradient.gradient_estimation>`

- `oracle`: a {obj}`PhaseOracle <oracq.algorithms.optimization.gradient.PhaseOracle>` or a phase-oracle operation with a `target` signature, input model FO (phase oracle); its `phase_scale` attribute must equal $2^{\text{grid\_bits}}$, otherwise it raises as an INPUT_PROMISE violation.
- `dimension`: the grid dimension $d$, range 1..16.
- `grid_bits`: the bits per coordinate $m$, range 1..32; the total width $d \cdot m$ must not exceed 64.

The oracle's three construction entry points: {obj}`abstract_phase_oracle(name, width, *, phase_scale) <oracq.algorithms.optimization.gradient.abstract_phase_oracle>` (an open declaration, with the implementation left for {obj}`bind <oracq.infrastructure.linking.bind>`), {obj}`gate_phase_oracle(width, angles, *, phase_scale, name=None) <oracq.algorithms.optimization.gradient.gate_phase_oracle>` (an explicit phase table, whose length must be $2^{\text{width}}$), and {obj}`function_phase_oracle(source, *, dimension, grid_bits, fmt=None, scale=None, ...) <oracq.algorithms.optimization.gradient.function_phase_oracle>` (compiles $f$ via mathfunc arithmetic and then applies phase kickback).

Returns an {obj}`Operation <oracq.infrastructure.builder.Operation>` with register `target` (width $d \cdot m$). After readout, decode with {obj}`gradient_from_readout(value, *, dimension, grid_bits) <oracq.algorithms.optimization.gradient.gradient_from_readout>`: the $i$-th component interprets the bit segment $[i \cdot m, (i+1) \cdot m)$ as two's complement and divides by $2^m$. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"jordan_gradient"` |
| `readout_register` / `decoder` | `"target"` / `"gradient_from_readout"` |
| `dimension` / `grid_bits` | the dimension $d$ and the bits per coordinate $m$ |
| `oracle_queries` / `classical_queries` | `1` / `O(dimension)` |

## Implementation notes

The generation chain is: prepare a uniform superposition over all coordinate registers → a single call to the phase oracle → apply {obj}`inverse_qft(grid_bits) <oracq.algorithms.common.fourier.inverse_qft>` to each coordinate bit segment. The register layout is the single `target`, with the $i$-th coordinate occupying the segment $[i \cdot m, (i+1) \cdot m)$; there is no work register, and the arithmetic work bits inside `function_phase_oracle` are restored after the phase kickback by the inverse call.

Design decision: the phase-scaling convention (`phase_scale == 2**grid_bits`) is checked at generation time rather than implicitly assumed, and the oracle width must also equal $d \cdot m$; both mismatches raise {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>` at generation time. `function_phase_oracle` defaults to the fixed-point format {obj}`FixedFormat(grid_bits+12, grid_bits+8) <oracq.algorithms.common.arithmetic.FixedFormat>` and requires a signed format with fraction ≥ grid_bits so that the grid coordinates are represented exactly.

Applicability boundary: the readout components are distinguishable only mod 1 (the phase $e^{2\pi i N a j}$ is the same for $a$ and $a \pm 1$), so the witnesses use components with $|\partial f/\partial x_i| < 1/2$; larger gradients require the caller to pre-scale $f$. When $f$ deviates from linearity, the readout peak broadens around the true value, and the failure probability decays at an approximately quadratic rate as the grid is refined (see the validation approach).

## Validation approach

Category C3 (probability-distribution semantics; acceptance criteria in `../development/validation-plan.md` §2): the output distribution must equal the closed-form expectation. Three layers of evidence:

- Structure: `tests/core/test_gradient.py:GradientTests.test_invalid_inputs_fail_at_generation` covers all entry-point violations — out-of-range dimension / bit counts, oracle-width mismatch, phase-table length mismatch, non-positive phase_scale, out-of-range readout values, etc.
- Numerical: `test_linear_function_exact_two_dimensions` — for a linear function (components $3/16$ and $-2/16$, $m = 4$) the readout distribution lands deterministically on the encoded values (probability 1.0, places = 12) and the decoding equals the gradient exactly; `test_function_phase_oracle_from_mathfunc` — the mathfunc path $f(x) = 0.25x$ recovers the component $0.25$ exactly; `test_perturbed_linear_concentrates_with_grid_bits` — for $f(x) = ax + x^2/N^2$ ($a = 3/8$) at $m = 3, 4, 5$ the peak position stays at the true value, the success probability rises monotonically, and the failure-probability decay rate is asserted as $q_{m+1} \le 0.34 \cdot q_m$ (under the phase perturbation $o(1/N)$ the out-of-peak leakage converges approximately quadratically; measured ratios 0.292 / 0.268).
- Binding: `test_abstract_oracle_binds_to_gate_implementation` — after the abstract declaration is bound via `bind` to the gate phase table, the readout decoding agrees with the direct construction, covering declaration/implementation two-layer consistency.

## Known gaps and planned stages

No known gaps; the stage-V1 witnesses are complete (exact linear readout + the mathfunc path + the failure-probability decay-rate assertion $q_{m+1} \le 0.34 \cdot q_m$, measured 0.292 / 0.268 + binding consistency).

## Related links

- Source: `src/oracq/algorithms/optimization/gradient.py`
- API reference: [quantum gradient estimation](../../api/algorithms/optimization/gradient.rst)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

`tests/verification/verify_estimation.py` performs readout-level numerical validation of this interface on real backends (11 cases in total, all passing). Experiment design:

- Exact linear readout: `gate_phase_oracle` phase tables, $d=1$ ($m=3,4,5$, gradient components $3/8$, $-5/16$, $9/32$), $d=2$ ($m=4$), $d=3$ ($m=3$); the readout distribution must land deterministically on the encoded values.
- mathfunc phase oracle: `function_phase_oracle` compiles $f(x)=0.25x$ ($d=1,m=4$) and $f(x_0,x_1)=0.25x_0-0.125x_1$ ($d=2,m=3$), also compared against central finite differences (for linear functions they agree with the truth exactly, gap $0$). The workspace of this path is budgeted at $1768$ / $1638$ qubits by {obj}`workspace_table <oracq.infrastructure.layout.workspace_table>`, exceeding the 24-qubit OriginIR budget, so it runs only on the three register-level paths reference / rir-pysparq / adapter-pysparq.
- Perturbed-linear convergence: $f(x)=ax+x^2/N^2$ ($a=3/8$), $m=3,4,5$; the peak position stays at the true value, and the failure-probability decay rate is compared against the approximately-quadratic-convergence criterion; the informational metric gives the central finite difference at $x=1/2$ (which differs from the linear coefficient by $O(1/N^2)$ when $f$ is nonlinear).
- Backend paths: reference, rir-pysparq, adapter-pysparq, originir-ext (gate phase-table cases, total width $\le 9$ qubits).

| Case | Scale | Path | Value |
|---|---|---|---|
| linear readout (gate table) | $d=1$, $m=3,4,5$ | all four paths | decode error $=0$, peak probability $=1$ |
| linear readout (gate table) | $d=2$, $m=4$; $d=3$, $m=3$ | all four paths | decode error $=0$, peak probability $=1$ |
| mathfunc oracle | $d=1,m=4$; $d=2,m=3$ | ref/rir/adp | decode error $=0$, peak probability $=1$, consistent with finite differences |
| perturbed convergence | $d=1$, $m=3,4,5$ | all four paths | success probabilities $0.9588/0.9880/0.9968$, decay rates $0.2922/0.2678$ ($<0.34$); finite differences $0.3906/0.3789/0.3760\to 0.375$ |

Reproduction command:

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_estimation.py
```

Artifact path: `out/verification/estimation.json` (case names prefixed `jordan-`).
