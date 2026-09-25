"""Implementation of the quantum convolutional neural network (QCNN, arXiv:1911.01117, ICLR 2020).

The paper writes the CNN convolution as an im2col matrix product
A^l F^l = Y^{l+1} (Eq. 5-7); on the quantum side, QRAM loads the row/column
vectors, a Hadamard-type circuit estimates the inner product (Eq. 16-19),
amplitude estimation encodes P_pq into a register and arithmetic recovers
Y_pq=(2P-1)||A_p|| ||F_q|| (Eq. 22-25), Boolean circuits implement capReLu
(a ReLU with cap C), conditional rotation moves pixel values into amplitudes
(Eq. 28-29), and after amplitude amplification and l_inf tomography sampling
the values are written back into the next layer's QRAM under the online
overwrite rule while pooling completes at the same time (§5.2.2).

This file provides the classical substrate first: the layer specification,
im2col, the classical forward mirror, and the row-tree QRAM (including the
pooling overwrite rule); the quantum building blocks are implemented
equation by equation in qcnn_layer.py per §5.1 of the paper. The semantic
validation of both files lives in tests/core/test_qcnn.py.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import cast


def tensor_get(x: Sequence[float], shape: tuple[int, int, int], i: int, j: int, d: int) -> float:
    """Read the element at (i, j, d) of a row-major flat tensor with channels last (H×W×D).

    Args:
        x: Tensor data stored flat in row-major order.
        shape: Shape tuple of the form ``(H, W, D)``.
        i: Row coordinate.
        j: Column coordinate.
        d: Channel coordinate.

    Returns:
        The element value, i.e. ``x[(i*W+j)*D+d]``.
    """
    h, w, _ = shape
    return x[(i * w + j) * shape[2] + d]


@dataclass(frozen=True)
class ConvSpec:
    """Single convolution layer specification: input H×W×D, kernel h×w×D×D', cap C, and P×P pooling."""

    input_shape: tuple[int, int, int]
    kernel_shape: tuple[int, int, int, int]
    cap: float = 5.0
    pool: int = 1
    pool_kind: str = "max"

    def __post_init__(self) -> None:
        """Validate kernel/input channel consistency and the value constraints of cap, pool, and pool_kind."""
        _h, _w, d, _dp = self.kernel_shape
        if self.input_shape[2] != d:
            raise ValueError("the kernel channel count does not match the input channel count")
        if self.cap <= 0 or self.pool < 1:
            raise ValueError("cap must be positive and pool must be at least 1")
        if self.pool_kind not in {"max", "average"}:
            raise ValueError("pool_kind must be max or average")

    @property
    def output_shape(self) -> tuple[int, int, int]:
        """Output shape of the valid convolution, i.e. ``(H-kh+1, W-kw+1, D')``."""
        h, w = self.input_shape[:2]
        kh, kw = self.kernel_shape[:2]
        return (h - kh + 1, w - kw + 1, self.kernel_shape[3])

    @property
    def pooled_shape(self) -> tuple[int, int, int]:
        """Output shape after pooling; identical to ``output_shape`` when ``pool`` is 1.

        Raises:
            ValueError: The output spatial dimensions are not divisible by the
                pooling window.
        """
        if self.pool == 1:
            return self.output_shape
        h, w, d = self.output_shape
        if h % self.pool or w % self.pool:
            raise ValueError("the output spatial dimensions must be divisible by the pooling window")
        return (h // self.pool, w // self.pool, d)


def im2col(x: Sequence[float], spec: ConvSpec) -> list[list[float]]:
    """Unroll the input tensor into A^l (Eq. 6): each row is one receptive field stacked by channel.

    Args:
        x: Input tensor data stored flat in row-major order (H×W×D).
        spec: Convolution layer specification, determining the receptive field
            size and the channel counts.

    Returns:
        list[list[float]]: The unrolled matrix, one receptive field per row,
        stacked in (d, ki, kj) order.
    """
    h, w, _ = spec.input_shape
    kh, kw, d, _ = spec.kernel_shape
    rows = []
    for i in range(h - kh + 1):
        for j in range(w - kw + 1):
            row = []
            for dd in range(d):
                for ki in range(kh):
                    for kj in range(kw):
                        row.append(tensor_get(x, spec.input_shape, i + ki, j + kj, dd))
            rows.append(row)
    return rows


def _kernel_value(
    kernel: Sequence[float], spec: ConvSpec, ki: int, kj: int, dd: int, q: int
) -> float:
    """Look up a kernel coefficient by the kernel index (ki, kj, dd, q)."""
    kh, kw, d, dp = spec.kernel_shape
    return kernel[((ki * kw + kj) * d + dd) * dp + q]


def kernel_columns(kernel: Sequence[float], spec: ConvSpec) -> list[list[float]]:
    """Columns of F^l (Eq. 6): each column is one kernel vectorized in (d, ki, kj) order.

    Args:
        kernel: Flat sequence of convolution kernel coefficients, matching the
            kernel shape of the specification.
        spec: Convolution layer specification, determining the kernel count
            and the vectorization order.

    Returns:
        list[list[float]]: Coefficient vectors, one column per output-channel
        kernel.
    """
    kh, kw, d, dp = spec.kernel_shape
    columns = []
    for q in range(dp):
        columns.append(
            [
                _kernel_value(kernel, spec, ki, kj, dd, q)
                for dd in range(d)
                for ki in range(kh)
                for kj in range(kw)
            ]
        )
    return columns


def cap_relu(value: float, cap: float) -> float:
    """capReLU (paper §5.1.6): min(max(x, 0), cap).

    Args:
        value: Scalar to activate.
        cap: Activation cap, must be positive.

    Returns:
        float: The activation value clamped to [0, cap].
    """
    return min(max(value, 0.0), cap)


def convolution_forward(
    x: Sequence[float], kernel: Sequence[float], spec: ConvSpec
) -> list[list[float]]:
    """Classical mirror: A F = Y followed by capReLU, then pooling per Eq. 39.

    Args:
        x: Input tensor data stored flat in row-major order.
        kernel: Flat sequence of convolution kernel coefficients.
        spec: Convolution layer specification, determining the receptive
            field, the activation cap, and the pooling kind.

    Returns:
        list[list[float]]: The pooled output matrix, rows for spatial
        positions and columns for output channels.
    """
    a = im2col(x, spec)
    f = kernel_columns(kernel, spec)
    _rows, _cols = spec.output_shape[0] * spec.output_shape[1], spec.kernel_shape[3]
    y = [
        [cap_relu(sum(ar * fc for ar, fc in zip(row, col, strict=True)), spec.cap)
         for col in f]
        for row in a
    ]
    return pool_tensor(y, spec)


def pool_tensor(y_rows: Sequence[Sequence[float]], spec: ConvSpec) -> list[list[float]]:
    """Pool the (H^{l+1}W^{l+1})×D^{l+1} row-major matrix per Eq. 39.

    Args:
        y_rows: Convolution output matrix, one spatial position per row.
        spec: Convolution layer specification; returned as-is when pool is 1.

    Returns:
        list[list[float]]: The pooled matrix; max takes the window high value
        and average the window mean.
    """
    if spec.pool == 1:
        return cast("list[list[float]]", y_rows)
    oh, ow, d = spec.output_shape
    ph, pw = spec.pooled_shape[:2]
    pooled = [[0.0] * d for _ in range(ph * pw)]
    counts = [[0] * d for _ in range(ph * pw)]
    for p, row in enumerate(y_rows):
        i, j = p // ow, p % ow
        pi, pj = i // spec.pool, j // spec.pool
        for q in range(d):
            index = pi * pw + pj
            if spec.pool_kind == "max":
                pooled[index][q] = max(pooled[index][q], row[q])
            else:
                pooled[index][q] += row[q]
            counts[index][q] += 1
    if spec.pool_kind == "average":
        for index in range(ph * pw):
            for q in range(d):
                pooled[index][q] /= max(counts[index][q], 1)
    return pooled


class QCNNQRAM:
    """The QRAM data structure of paper §5.2 (classical side): one partial-norm tree per row.

    The tree is a 1-based complete binary heap: leaves store row elements and
    internal nodes store squared partial sums; it supports the row-value and
    row-norm queries needed by the ``|p>|0> → |p>|A_p>`` semantics, as well as
    the online pooling overwrite of §5.2.2 (max keeps the high value, average
    spreads it).
    """

    def __init__(self, rows: Sequence[Sequence[float]]) -> None:
        """Build the QRAM from a sequence of equal-length rows, writing one squared-partial-sum tree per row."""
        if not rows:
            raise ValueError("the QRAM requires at least one row")
        self.width: int = len(rows[0])
        if any(len(r) != self.width for r in rows):
            raise ValueError("all rows must have the same length")
        self.count: int = len(rows)
        leaf_start = 1
        while leaf_start < self.width:
            leaf_start *= 2
        self.leaf_start: int = leaf_start
        self.tree_size: int = 2 * leaf_start
        self.trees: list[list[float]] = [
            [0.0] * self.tree_size for _ in range(self.count)
        ]
        self.leaves: list[list[float]] = [[0.0] * self.width for _ in range(self.count)]
        for p, row in enumerate(rows):
            for r, value in enumerate(row):
                self._write(p, r, value)

    def _write(self, p: int, r: int, value: float) -> None:
        """Write the value at ``(p, r)`` and update that row's squared-partial-sum tree bottom-up."""
        self.leaves[p][r] = value
        node = self.leaf_start + r
        self.trees[p][node] = value * value
        node //= 2
        while node:
            self.trees[p][node] = self.trees[p][2 * node] + self.trees[p][2 * node + 1]
            node //= 2

    def row(self, p: int) -> list[float]:
        """Return the list of values of row ``p`` (a reference to internal storage, not a copy).

        Args:
            p: Row index, in 0..number of rows − 1.

        Returns:
            list[float]: Current values of each column in that row.
        """
        return self.leaves[p]

    def norm(self, p: int) -> float:
        """Euclidean norm of row ``p``, the square root of the root of the squared-partial-sum tree.

        Args:
            p: Row index, in 0..number of rows − 1.

        Returns:
            float: Euclidean norm of that row vector.
        """
        return math.sqrt(self.trees[p][1])

    def update_with_pooling(
        self, p: int, r: int, value: float, kind: str, state: int | None = None
    ) -> int | None:
        """Online overwrite (§5.2.2): max writes only when the new value is higher, average spreads it.

        state is the accumulated count at (p, r) (needed by average pooling);
        the new state is returned.

        Args:
            p: Row index, in 0..number of rows − 1.
            r: Column index, in 0..row width − 1.
            value: Candidate value to overwrite with.
            kind: Pooling kind, ``max`` or ``average``.
            state: Accumulated count at this position; used by average
                pooling only.

        Returns:
            int | None: The updated accumulated count; returned unchanged
            under max pooling.
        """
        if kind == "max":
            if value > self.leaves[p][r]:
                self._write(p, r, value)
            return state
        count = (state or 0) + 1
        current = self.leaves[p][r] * (state or 0)
        self._write(p, r, (current + value) / count)
        return count
