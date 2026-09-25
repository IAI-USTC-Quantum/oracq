# QCNN line-by-line implementation walkthrough

**English** · [简体中文](../zh/manual/qcnn-walkthrough.html)

This page explains the QCNN (arXiv:1911.01117, ICLR 2020) implementation line by line: first an end-to-end minimal example, then a segment-by-segment reading of the sources of the classical basis (`qcnn.py`) and the quantum building blocks (`qcnn_layer.py`). Every code block is runnable (`tests/core/test_qcnn.py` is built from these snippets). For the algorithm overview see the [quantum convolutional neural network (QCNN)](../zh/manual/algorithms/qcnn.html) page; for the API reference see [qcnn](../api/algorithms/qml/qcnn.rst) and [qcnn_layer](../api/algorithms/qml/qcnn_layer.rst).

## 1. End-to-end minimal example

```python
import random
from oracq.algorithms.qml.qcnn import ConvSpec, convolution_forward
from oracq.algorithms.qml.qcnn_layer import qcnn_sampled_layer

# 5×5 single-channel input, 2×2 kernel ×2 output channels, capReLU cap 3.0, 2×2 max pooling.
spec = ConvSpec(input_shape=(5, 5, 1), kernel_shape=(2, 2, 1, 2), cap=3.0, pool=2)
random.seed(7)
x = [random.uniform(-1, 1) for _ in range(25)]            # input image (H×W×D flattened)
kernel = [random.uniform(-1, 1) for _ in range(2 * 2 * 1 * 2)]  # 4-tensor kernel flattened
classical = convolution_forward(x, kernel, spec)          # classical mirror (reference for comparison)
pooled, stats = qcnn_sampled_layer(x, kernel, spec, samples=40000, eta=0.0, seed=1)
```

## 2. The classical basis line by line (`qcnn.py`)

### ConvSpec — the layer specification (dimension rules of paper §3.3–3.4)

```python
@dataclass(frozen=True)
class ConvSpec:
    input_shape: tuple      # (H^l, W^l, D^l): the three dimensions of the input tensor
    kernel_shape: tuple     # (h, w, D^l, D^{l+1}): the 4-tensor kernel
    cap: float = 5.0        # capReLU cap C (§5.1.6: provides a known max f for the conditional rotation)
    pool: int = 1           # pooling window P (Eq. 39)
    pool_kind: str = "max"  # "max" or "average" (the two overwrite rules of §5.2.2)
```

`output_shape` follows Eq. (4) of the paper: `H^{l+1} = H^l − h + 1` (no padding/stride); `pooled_shape` divides the output spatial dimensions by the pooling window (`__post_init__` enforces divisibility — a 5×5 input with a 2×2 kernel gives a 4×4 output, exactly divisible by 2).

### im2col — the matrix form of Eq. (6)

```python
def im2col(x, spec):
    h, w, _ = spec.input_shape
    kh, kw, d, _ = spec.kernel_shape
    rows = []
    for i in range(h - kh + 1):        # each output row corresponds to one receptive field
        for j in range(w - kw + 1):
            row = []
            for dd in range(d):        # stack channel by channel (column-major vectorization, paper §3.4)
                for ki in range(kh):
                    for kj in range(kw):
                        row.append(tensor_get(x, spec.input_shape, i + ki, j + kj, dd))
            rows.append(row)
    return rows                        # A^l: row-major, one receptive field per row
```

Each row `A^l_p` is the vectorization of the receptive field at position (i,j); {obj}`kernel_columns <oracq.algorithms.qml.qcnn.kernel_columns>` expands each kernel into one column of `F^l` in the same (d, ki, kj) order — the two orders agree, so the inner product is exactly the convolution output `Y^{l+1}_{p,q}`.

### convolution_forward — the classical mirror

```python
a = im2col(x, spec)
f = kernel_columns(kernel, spec)
y = [[cap_relu(sum(ar * fc for ar, fc in zip(row, col, strict=True)), spec.cap)
      for col in f] for row in a]     # convolution + capReLU (§5.1.6)
return pool_tensor(y, spec)           # the pooling of Eq. 39
```

`cap_relu = min(max(x, 0), cap)` — the paper uses the capped ReLU so that `max_pq f(Y_pq)` is known a priori (the cap constant), allowing the conditional rotation (Eq. 28) to use a fixed normalization.

### QCNNQRAM — the row-tree data structure (§4.3 + §5.2.2)

```python
self.leaf_start = 1
while leaf_start < self.width:      # complete binary tree, leaves start at leaf_start
    leaf_start *= 2
self.tree_size = 2 * leaf_start     # 1-based heap: root is 1, internal nodes store squared partial sums
```

`_write(p, r, value)` writes one leaf and refreshes the squared partial sums up the ancestors — a single-point update costs logarithmically (the O(log² n) update of Theorem 4.3 in the paper). `update_with_pooling` implements the online overwrite: `max` writes only when the new value is higher (paper §5.2.2, "update the leaf if the new sampled value is higher"), and `average` amortizes with a state count ("replace the actual value by the new averaged value").

## 3. The quantum building blocks line by line (`qcnn_layer.py`)

### vector_angle_tables — from data to runtime banks

```python
angles[(index << width) | node] = round(angle * scale) % (1 << angle_width)
```

Line by line, node by node: `angle = 2·atan2(√S_right, √S_left)` (the RY half-angle convention, amplitude ratio √(S_right/S_node)); quantized into an `angle_width`-bit angle word. Signs go into a separate `signs` bank (key `(index<<width)|leaf`). **Changing to another layer's data only recomputes these two tables — the circuit object is completely unchanged** (QRAM semantics).

### qcnn_vector_prep — the row-loading circuit of Eq. (15)

```python
for depth in range(vector_width):          # tree level: bit from high to low
    for prefix in range(1 << depth):
        node = (1 << depth) - 1 + prefix   # tree node number (addressing consistent with qram_state_prep)
        node_register = _constant_register(b, f"node_{depth}_{prefix}", vector_width, node)
        if depth:
            with b.control(b["target"][bit + 1:], prefix):   # prefix control on the settled high bits
                _prep_node_rotation(b, angles, node_register, bit, angle_width)
        else:
            _prep_node_rotation(b, angles, node_register, bit, angle_width)
        # symmetrically restore node_register (X flips it back)
```

The three steps of `_prep_node_rotation`: query `angles` (address {obj}`fuse(node_register, index) <oracq.infrastructure.ir.fuse>` — **index in the high bits, consistent with the bank key**); synthesize controlled RYs with bitwise weighting (each set angle-word bit k applies `ry(2π·2^k/2^aw)`); un-query to clean the work. The sign kickback at the end: query `signs` addressed by the already-prepared `target`, apply a Z on the flag, then un-query — negative components are left with only a −1 phase.

### qcnn_inner_product — Eqs. (16)–(19)

```python
for register in (b["p"], b["q"], b["flag"]):
    b.h(register)                                  # uniform superposition Σ|p>|q> (1/√2)(|0>+|1>)
with b.control(b["flag"], 0):
    invoke(b, row_prep, "row", index=b["p"], target=b["vec"], work=b["work"])
with b.control(b["flag"], 1):
    invoke(b, column_prep, "col", index=b["q"], target=b["vec"], work=b["work"])
b.h(b["flag"])                                     # amplitude √P_pq lands on the flag=0 branch
```

This is exactly the three-line derivation of Eq. (16) in the paper: branch loading of |A_p⟩/|F_q⟩ → Hadamard interference → measurement probability P_pq=(1+⟨A_p|F_q⟩)/2 on flag=0. Validation against Eq. (20): the joint probability of (p,q,flag=0) after summing out (p,q,vec,work) = (1+⟨A_p|F_q⟩)/(2·rows·columns), with a measured worst-case deviation of 4×10⁻⁵.

### qcnn_sampled_layer — the sampling driver of Eqs. (34)–(39)

```python
values[(p, q)] = cap_relu(inner, spec.cap)          # Y_pq=(2P_pq-1)||A_p||||F_q|| after activation
weight = {key: value * value for key, value in values.items()}
# sampling probability ∝ f(Y)^2 (Eq. 34); each sample yields the triple (p, q, f(Y)):
if value < eta:      # Eq. 36: pixels below η are treated as unsampled
    continue
i, j = p // ow, p % ow                       # Eq. 7: recover the row number (i^{l+1}, j^{l+1})
index = (i // spec.pool) * pw + (j // spec.pool)   # Eq. 39: pooling-region mapping
pooled[index][q] = max(pooled[index][q], value)    # QRAM online overwrite (§5.2.2)
```

{obj}`quantized_prepared_state <oracq.algorithms.qml.qcnn_layer.quantized_prepared_state>` is the classical mirror of the circuit semantics (level-by-level quantized angle rotations + signs); the driver uses it to recover all of Y — the circuit and the mirror are compared amplitude by amplitude in tests (1e-9), guaranteeing that the distribution modeling introduces no extra approximation.

## 4. Validation matrix (all passing)

| Experiment | Scale | Result |
|---|---|---|
| im2col vs naive convolution | 5×5×1, kernel 2×2×1×2 | element-wise 1e-12 |
| forward + pooling vs naive | same + 2×2 max | element-wise 1e-12 |
| QRAM rows/norms/overwrite | 16 rows | exact; max keeps the higher value, average gives (2+4)/2=3 |
| row-preparation circuit vs normalized row | with negative components | worst 3e-4 (= the lower bound of aw=10 angle quantization) |
| circuit vs quantized mirror | same | 1e-9 |
| inner-product probability vs Eq. 20 | 2×2 rows/columns | 4e-5 |
| max-pooling sampled recovery | 40000 samples | differs from the classical maximum by ≤5e-3 |
| η→∞ | 1000 samples | all zero (Eq. 36) |
| average convergence | 60000 samples | to the f²-weighted expectation within ≤2e-2 |
