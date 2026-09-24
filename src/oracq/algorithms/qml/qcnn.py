"""量子卷积神经网络（QCNN，arXiv:1911.01117，ICLR 2020）的实现。

论文把 CNN 的卷积写为 im2col 矩阵乘 A^l F^l = Y^{l+1}（Eq. 5-7），
量子侧用 QRAM 加载行/列向量、Hadamard 型电路估计内积（Eq. 16-19），
幅度估计把 P_pq 编码进寄存器后经算术恢复 Y_pq=(2P-1)||A_p|| ||F_q||
（Eq. 22-25），布尔电路实现 capReLu（带上限 C 的 ReLU），条件旋转把
像素值转入幅度（Eq. 28-29），幅度放大与 l_inf 层析采样后按在线覆写规则
写回下一层 QRAM 并同时完成池化（§5.2.2）。

本文件先给出经典基底：层规格、im2col、经典前向镜像与行树 QRAM（含
池化覆写规则）；量子构件在 qcnn_layer.py 中按论文 §5.1 逐式实现。
两个文件的语义验证见 tests/core/test_qcnn.py。
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import cast


def tensor_get(x: Sequence[float], shape: tuple[int, int, int], i: int, j: int, d: int) -> float:
    """读取通道在后（H×W×D）行主序扁平张量在 (i, j, d) 处的元素。

    Args:
        x: 按行主序扁平存储的张量数据。
        shape: 形如 ``(H, W, D)`` 的形状元组。
        i: 行坐标。
        j: 列坐标。
        d: 通道坐标。

    Returns:
        对应元素值，即 ``x[(i*W+j)*D+d]``。
    """
    h, w, _ = shape
    return x[(i * w + j) * shape[2] + d]


@dataclass(frozen=True)
class ConvSpec:
    """单层卷积规格：输入 H×W×D，核 h×w×D×D'，cap C 与 P×P 池化。"""

    input_shape: tuple[int, int, int]
    kernel_shape: tuple[int, int, int, int]
    cap: float = 5.0
    pool: int = 1
    pool_kind: str = "max"

    def __post_init__(self) -> None:
        """校验核/输入通道一致性、cap、pool 与 pool_kind 的取值约束。"""
        _h, _w, d, _dp = self.kernel_shape
        if self.input_shape[2] != d:
            raise ValueError("核通道数与输入通道数不匹配")
        if self.cap <= 0 or self.pool < 1:
            raise ValueError("cap 必须为正，pool 至少为 1")
        if self.pool_kind not in {"max", "average"}:
            raise ValueError("pool_kind 必须是 max 或 average")

    @property
    def output_shape(self) -> tuple[int, int, int]:
        """有效卷积的输出形状，即 ``(H-kh+1, W-kw+1, D')``。"""
        h, w = self.input_shape[:2]
        kh, kw = self.kernel_shape[:2]
        return (h - kh + 1, w - kw + 1, self.kernel_shape[3])

    @property
    def pooled_shape(self) -> tuple[int, int, int]:
        """池化后的输出形状；``pool`` 为 1 时与 ``output_shape`` 相同。

        Raises:
            ValueError: 输出空间维度不能被池化窗口整除。
        """
        if self.pool == 1:
            return self.output_shape
        h, w, d = self.output_shape
        if h % self.pool or w % self.pool:
            raise ValueError("输出空间维度必须能被池化窗口整除")
        return (h // self.pool, w // self.pool, d)


def im2col(x: Sequence[float], spec: ConvSpec) -> list[list[float]]:
    """把输入张量展开为 A^l（Eq. 6）：每行是一个感受野的按通道堆叠。

    Args:
        x: 按行主序扁平存储的输入张量数据（H×W×D）。
        spec: 卷积层规格，决定感受野尺寸与通道数。

    Returns:
        list[list[float]]: 每行一个感受野、按 (d, ki, kj) 顺序堆叠的展开矩阵。
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
    """按卷积核下标 (ki, kj, dd, q) 查询核系数。"""
    kh, kw, d, dp = spec.kernel_shape
    return kernel[((ki * kw + kj) * d + dd) * dp + q]


def kernel_columns(kernel: Sequence[float], spec: ConvSpec) -> list[list[float]]:
    """F^l 的列（Eq. 6）：每列一个核的按 (d, ki, kj) 顺序向量化。

    Args:
        kernel: 卷积核系数的扁平序列，与规格的核形状一致。
        spec: 卷积层规格，决定核数量与向量化顺序。

    Returns:
        list[list[float]]: 每列对应一个输出通道核的系数向量。
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
    """capReLU（论文 §5.1.6）：min(max(x, 0), cap)。

    Args:
        value: 待激活的标量。
        cap: 激活上限，必须为正。

    Returns:
        float: 截断到 [0, cap] 的激活值。
    """
    return min(max(value, 0.0), cap)


def convolution_forward(
    x: Sequence[float], kernel: Sequence[float], spec: ConvSpec
) -> list[list[float]]:
    """经典镜像：A F = Y 后 capReLU，再按 Eq. 39 池化。

    Args:
        x: 按行主序扁平存储的输入张量数据。
        kernel: 卷积核系数的扁平序列。
        spec: 卷积层规格，决定感受野、激活上限与池化方式。

    Returns:
        list[list[float]]: 池化后的输出矩阵，行对应空间位置、列对应输出通道。
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
    """把 (H^{l+1}W^{l+1})×D^{l+1} 的行主序矩阵按 Eq. 39 池化。

    Args:
        y_rows: 卷积输出矩阵，每行一个空间位置。
        spec: 卷积层规格，pool 为 1 时原样返回。

    Returns:
        list[list[float]]: 池化后的矩阵，max 取窗口高值、average 取窗口均值。
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
    """论文 §5.2 的 QRAM 数据结构（经典侧）：每行一棵部分范数树。

    树为 1-based 完全二叉堆：叶存行元素，内部节点存平方部分和；支持
    ``|p>|0> → |p>|A_p>`` 语义所需的行值与行范数查询，以及 §5.2.2 的在线
    池化覆写（max 保高值、average 均摊）。
    """

    def __init__(self, rows: Sequence[Sequence[float]]) -> None:
        """由等长行序列构建 QRAM，为每行写入一棵平方部分和树。"""
        if not rows:
            raise ValueError("QRAM 至少需要一行")
        self.width: int = len(rows[0])
        if any(len(r) != self.width for r in rows):
            raise ValueError("所有行必须等长")
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
        """写入 ``(p, r)`` 处的值并自底向上更新该行的平方部分和树。"""
        self.leaves[p][r] = value
        node = self.leaf_start + r
        self.trees[p][node] = value * value
        node //= 2
        while node:
            self.trees[p][node] = self.trees[p][2 * node] + self.trees[p][2 * node + 1]
            node //= 2

    def row(self, p: int) -> list[float]:
        """返回第 ``p`` 行的取值列表（内部存储引用，非拷贝）。

        Args:
            p: 行编号，取 0..行数−1。

        Returns:
            list[float]: 该行各列的当前取值。
        """
        return self.leaves[p]

    def norm(self, p: int) -> float:
        """第 ``p`` 行的欧几里得范数，取平方部分和树根节点的平方根。

        Args:
            p: 行编号，取 0..行数−1。

        Returns:
            float: 该行向量的欧几里得范数。
        """
        return math.sqrt(self.trees[p][1])

    def update_with_pooling(
        self, p: int, r: int, value: float, kind: str, state: int | None = None
    ) -> int | None:
        """在线覆写（§5.2.2）：max 仅当新值更高时写，average 均摊。

        state 为该 (p, r) 的累计计数（average 池化需要）；返回新状态。

        Args:
            p: 行编号，取 0..行数−1。
            r: 列编号，取 0..行宽−1。
            value: 覆写候选值。
            kind: 池化方式，取 ``max`` 或 ``average``。
            state: 该位置的累计计数；仅 average 池化使用。

        Returns:
            int | None: 更新后的累计计数；max 池化时原样返回。
        """
        if kind == "max":
            if value > self.leaves[p][r]:
                self._write(p, r, value)
            return state
        count = (state or 0) + 1
        current = self.leaves[p][r] * (state or 0)
        self._write(p, r, (current + value) / count)
        return count
