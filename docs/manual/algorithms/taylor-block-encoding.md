# Truncated Taylor Block Encoding

**English** · <a href="../../../zh/manual/algorithms/taylor-block-encoding.html">简体中文</a>

> Category C2 · Module [`oracq.algorithms.common.hamiltonian`](../../api/algorithms/common/hamiltonian.rst) · Stage V2

## Overview

Assemble the truncated Taylor series of a matrix exponential into a block encoding: given a block encoding of the generator and an order $d$,

$$
e^{-iHt} \;\approx\; \sum_{k=0}^{d} \frac{(-it)^k}{k!}\, H^k,
$$

form a linear combination of the $k$-th powers $H^k$ (the $k$-fold product of the input BE) with analytic coefficients. The module docstring positions it as "a closable ordinary Hamiltonian-function BE", for QODE solvers to use as an injectable `hamiltonian_function`; it does not require Hermiticity and does not claim to be an efficiently optimal HamSim algorithm, and can be replaced by the QSP/HamSim protocol (see [Hamiltonian simulation protocol](hamiltonian-simulation.md)).

## Interface and input model

```python
taylor_hamiltonian(hamiltonian, time, *, degree=2)
```

API entry: {obj}`taylor_hamiltonian <oracq.algorithms.common.hamiltonian.taylor_hamiltonian>`

- `hamiltonian`: a {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`, the block encoding of the evolution generator (input model BE); the mathematical semantics (such as Hermiticity) are the caller's responsibility.
- `time`: the evolution time, a finite real number.
- `degree`: the truncation order $d$, a non-negative integer, default 2.

Returns a `BlockEncoding` whose zero-signal block approximately equals $\sum_k (-it)^k H^k / k!$. Module attributes:

| Attribute | Meaning |
|---|---|
| `algorithm` | `"truncated_taylor_hamiltonian_function"` |
| `degree` / `time` | echo of the call parameters |
| `be_alpha` | the LCU normalization $\sum_{k=0}^{d} (\alpha t)^k / k!$ ($\alpha$ is the input BE normalization) |
| `correctness` / `success_condition` | `"pending"` / `"signal == 0"` |

## Implementation notes

The power sequence is generated iteratively: starting from $(1, I)$, {obj}`product(hamiltonian, current) <oracq.algorithms.input_model.operators.product>` right-multiplies step by step to obtain $H^k$, the coefficients are taken analytically as $(-it)^k / k!$, and the whole is then combined by {obj}`lcu(powers) <oracq.algorithms.input_model.block_encoding.lcu>`. The register layout is decided by the LCU: `target = Bits(width)`, and `signal` is the selector ($\lceil\log_2(d+1)\rceil$ bits) concatenated with each term's signal space; the phases of complex coefficients are realized by global phase gates.

Applicability boundary: the gate count grows linearly with $d$ (a deep power chain); this function performs no truncation-error-bound analysis, and the error is to be assessed by the caller from $\lVert H\rVert t$ and $d$. Generation-time validation checks that `degree` is a non-negative integer and `time` a finite real number.

## Validation approach

Category C2 (approximately continuous semantics; the acceptance criteria are in `../development/validation-plan.md` §2). The statement consistent with the validation coverage matrix is: the `hamiltonian.py` row lists no separate numerical witness for this entry point, and this function's evidence comes from being consumed by the ODE solvers —

- Structure: `tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` injects it with `degree=1` into the four routes lchs / cbmd / schrodingerization / carleman, checking that the input oracle slots are kept and that the IR round-trips through JSON (registered in the ode row of the validation coverage matrix).
- Numerical: `tests/core/test_sde.py:SolverContractTests.test_qode_problem_accepted_by_lchs` likewise injects the LCHS contract with `degree=1` and solves (registered in the sde row). A direct cross-check of the truncation accuracy itself, "error vs order", is missing.
- Binding: no independent binding witness; the binding semantics of the LCU combinator are covered by `tests/core/test_language.py:BlockEncodingTests.test_alpha_survives_ir_serialization` (see [block encoding algebra](block-encoding-algebra.md)).

## Known gaps and planned stages

Consistent with the `hamiltonian.py` row of the validation coverage matrix: the module's registered gap is the Trotter order-error-rate sweep (stage V2 convergence-scan framework); a direct cross-check of this entry point's own Taylor truncation error is the same kind of convergence evidence, to be added once that framework lands.

## Related links

- Source: `src/oracq/algorithms/common/hamiltonian.py`
- Same-family pages: [Hamiltonian simulation protocol](hamiltonian-simulation.md), [Trotter product-formula simulation](trotter.md), [block encoding algebra](block-encoding-algebra.md)
- API reference: [Hamiltonian evolution](../../api/algorithms/common/hamiltonian.rst)
- Concepts: [Oracles and operator representations](../operators.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)

## Numerical validation

Paper-grade numerical validation script: `tests/verification/verify_blockencoding.py` (real-backend execution, no mocks, no skips; all 73 cases passed as of 2026-09-16), with artifacts in `out/verification/blockencoding.json`. This page corresponds to the 4 cases `taylor-block-encoding-d{1,2,3,4}`, closing the registered gap of "a direct error-vs-order cross-check of the truncation accuracy itself is missing".

Experiment design: the input generator is the 2×2 Hermitian matrix $H=\begin{pmatrix}1.0&0.4\\0.4&-0.6\end{pmatrix}$ (encoded by {obj}`matrix_pauli_encoding <oracq.algorithms.input_model.block_encoding.matrix_pauli_encoding>`, $\alpha_{\rm in}=1.4$), with evolution time $t=0.7$ and orders $d=1..4$. The block semantics are $\big(\sum_{k\le d}(-itH)^k/k!\big)/\alpha$, with $\alpha=\sum_{k\le d}(\alpha_{\rm in}t)^k/k!$. Two classes of error are reported separately: **implementation error** — the zero-signal block extracted column by column on reference against the same-order truncated series (computed independently with numpy); **method error** — the block against $e^{-iHt}/\alpha$ (computed independently with a numpy eigendecomposition); additionally, `be_alpha` is checked against the closed form of the series, and rir-pysparq / adapter-pysparq cross-check column by column.

| Case | Degree | Implementation error | Method error (informational) | $\alpha$ deviation | Backend cross-check |
|---|---|---|---|---|---|
| `taylor-block-encoding-d1` | 1 | 8.3e-17 | 1.41e-1 | 0 | 0 |
| `taylor-block-encoding-d2` | 2 | 1.8e-16 | 2.81e-2 | 0 | 0 |
| `taylor-block-encoding-d3` | 3 | 1.7e-16 | 5.20e-3 | 0 | 0 |
| `taylor-block-encoding-d4` | 4 | 1.9e-16 | 7.75e-4 | 0 | 0 |

The method error decreases monotonically with the order (on the order of $(\alpha_{\rm in}t)^{d+1}/(d+1)!$, with $\alpha_{\rm in}t=0.98$), while the implementation error stays at machine precision — the truncation error is indeed controlled by the caller through $d$, and the assembly itself introduces no extra error.

Reproduction command:

```bash
PYTHONPATH=src /path/to/backend/python tests/verification/verify_blockencoding.py
```

Artifacts: `out/verification/blockencoding.json`.
