# QCNN

**English** · <a href="../../../zh/manual/algorithms/qcnn.html">简体中文</a>

> Category C4 · Module [`oracq.algorithms.qml.qcnn`](../../api/algorithms/qml/qcnn.rst) and [`oracq.algorithms.qml.qcnn_layer`](../../api/algorithms/qml/qcnn_layer.rst) · Paper arXiv:1911.01117 (ICLR 2020)

## Overview

QCNN moves one layer of a classical CNN onto quantum circuits: convolution becomes a matrix product A^ℓ·F^ℓ=Y^{ℓ+1} via im2col (Eq. 5–7); the quantum side loads row/column vectors from QRAM and estimates inner products with Hadamard-type circuits (Eq. 16–19); sampling and l∞ tomography draw high-valued pixels weighted by squared pixel values (Eq. 28–38); pooling is completed in online QRAM updates (§5.2.2). This implementation covers all four forward operation classes (convolution, activation, pooling, and sampled readout), with small-scale semantics verified elementwise against a classical CNN.

The implementation boundary recorded by `correctness="pending"`: amplitude-estimation register encoding (Eq. 22–25) and amplitude amplification (Eq. 30–33) are modeled by their exact distribution (sampling probability f(Y)²/Σf², Eq. 34); the mathfunc `asin` primitive needed by the angle arithmetic of the conditional rotations (θ=2·asin(√(f/C))) is available, but the controlled-power + inverse QFT encoding of AE has not been built into circuits; backpropagation (Thm 6.1) is not implemented.

## Interface

```python
ConvSpec(input_shape, kernel_shape, cap, pool, pool_kind)  # layer spec
im2col(x, spec)                 # input tensor -> A^ℓ (Eq. 6)
kernel_columns(kernel, spec)    # rank-4 tensor kernel -> F^ℓ columns (Eq. 6)
cap_relu(value, cap)            # capReLU (§5.1.6)
convolution_forward(x, kernel, spec)   # classical mirror: convolution+capReLU+pooling
QCNNQRAM(rows)                  # row-tree QRAM (§4.3/§5.2, with pooling overwrite)

vector_angle_tables(rows, angle_width)  # row/column -> angle tree + sign bank (runtime table)
qcnn_vector_prep(count, width, angle_width)   # |p>|0>->|p>|A_p/||A_p||> (Eq. 15)
qcnn_inner_product(rows, cols, width, aw)     # Hadamard inner-product circuit (Eq. 16–19)
quantized_prepared_state(vector, aw)          # angle-quantized mirror (classical reproduction of the circuit semantics)
qcnn_sampled_layer(x, kernel, spec, *, samples, eta, seed)  # sampled driver (Eq. 34–39)
```

API entries: {obj}`ConvSpec <oracq.algorithms.qml.qcnn.ConvSpec>`, {obj}`im2col <oracq.algorithms.qml.qcnn.im2col>`, {obj}`kernel_columns <oracq.algorithms.qml.qcnn.kernel_columns>`, {obj}`cap_relu <oracq.algorithms.qml.qcnn.cap_relu>`

## Implementation notes

**Row/column preparation ({obj}`qcnn_vector_prep <oracq.algorithms.qml.qcnn_layer.qcnn_vector_prep>`)**: each node of the multiplexed rotation tree is driven by a QRAM angle bank — the address is {obj}`fuse(node, index) <oracq.infrastructure.ir.fuse>` (bank key `(index<<width)|node`), in three steps: query, bitwise-weighted composition of controlled RY (angle 2π·2^k/2^aw), unquery; negative components are written into the phase via a Z kickback from the sign bank. Swapping the data swaps only the memory table, not the circuit (QRAM semantics). A tree of depth d has 2^d nodes and 2·width total queries.

**Hadamard inner product ({obj}`qcnn_inner_product <oracq.algorithms.qml.qcnn_layer.qcnn_inner_product>`)**: the three registers p, q, flag are put in uniform superposition; the flag=0 branch loads the row vector and the flag=1 branch loads the column vector (two sets of controlled preparations), and a Hadamard is finally applied to flag. The probability of measuring (p,q,flag=0) is (1+⟨A_p|F_q⟩)/(2·rows·cols) — an exact reproduction of Eq. 20; the inner product may be positive or negative, while the probability is always positive.

**Sampled driver ({obj}`qcnn_sampled_layer <oracq.algorithms.qml.qcnn_layer.qcnn_sampled_layer>`)**: the quantized mirror first recovers all Y_pq=(2P_pq−1)‖A_p‖‖F_q‖ and capReLU values, then samples according to f²/Σf² (the distribution of Eq. 34); each sample yields a triple (p,q,f(Y)), pixels with f<η count as unsampled and are set to zero (Eq. 36); after mapping into pooling regions per Eq. 39, values are aggregated under the QRAM overwrite rules (max keeps the high value, average spreads it).

## Validation

`tests/core/test_qcnn.py` (9 cases, all green):

| Layer | Assertions | Result |
|---|---|---|
| Classical basis | im2col inner product vs naive triple loop; forward+pooling vs naive implementation; QRAM rows/norms/overwrites | exact elementwise (1e-12) |
| Row preparation | circuit amplitudes vs rows/‖·‖ (including negative components); circuit vs quantized mirror | worst 3e-4 (aw=10 angle-quantization lower bound) |
| Inner-product circuit | measurement probabilities vs Eq. 20 | worst 4e-5 |
| Sampled driver | max pooling recovers the region maximum from 40000 samples; η→∞ all zero; average converges to the f²-weighted expectation | 5e-3 / exact / 2e-2 |

## Known boundaries

- The register encoding of amplitude estimation (AE + median purification) has not been built into circuits — Y recovery currently goes through the quantized mirror (consistent with the circuit semantics), while conditional rotations and amplitude amplification are modeled by distribution. Once built, the mirror part of `qcnn_sampled_layer` should be replaced by register readout.
- Backpropagation (paper §6) is not implemented.
- Row/column counts must be powers of two (a small-scale verification constraint; general shapes are reachable via padding).

## Related links

- Source: `src/oracq/algorithms/qml/qcnn.py` (overall flow) and `src/oracq/algorithms/qml/qcnn_layer.py` (quantum building blocks)
- API reference: [Quantum convolutional neural network](../../api/algorithms/qml/qcnn.rst), [Quantum building blocks of the QCNN](../../api/algorithms/qml/qcnn_layer.rst)
- User manual: [QCNN line-by-line implementation walkthrough](../qcnn-walkthrough.md)
- Validation matrix: [Validation coverage matrix](../../development/validation-coverage.md)
