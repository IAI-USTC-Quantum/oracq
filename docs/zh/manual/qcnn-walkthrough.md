# QCNN 逐行实现讲解

<a href="../../en/manual/qcnn-walkthrough.html">English</a> · **简体中文**

本页把 QCNN（arXiv:1911.01117，ICLR 2020）的实现逐行讲清：先给端到端最小示例，再逐段读经典基底（`qcnn.py`）与量子构件（`qcnn_layer.py`）的源码。所有代码块均可运行（`tests/core/test_qcnn.py` 即由这些片段构成）。算法概述见[量子卷积神经网络（QCNN）](algorithms/qcnn.md)；API 参考见 [qcnn](../api/algorithms/qml/qcnn.rst) 与 [qcnn_layer](../api/algorithms/qml/qcnn_layer.rst)。

## 1. 端到端最小示例

```python
import random
from oracq.algorithms.qml.qcnn import ConvSpec, convolution_forward
from oracq.algorithms.qml.qcnn_layer import qcnn_sampled_layer

# 5×5 单通道输入，2×2 核 ×2 个输出通道，capReLU 上限 3.0，2×2 max 池化。
spec = ConvSpec(input_shape=(5, 5, 1), kernel_shape=(2, 2, 1, 2), cap=3.0, pool=2)
random.seed(7)
x = [random.uniform(-1, 1) for _ in range(25)]            # 输入图像（H×W×D 展平）
kernel = [random.uniform(-1, 1) for _ in range(2 * 2 * 1 * 2)]  # 4-张量核展平
classical = convolution_forward(x, kernel, spec)          # 经典镜像（对照参考）
pooled, stats = qcnn_sampled_layer(x, kernel, spec, samples=40000, eta=0.0, seed=1)
```

## 2. 经典基底逐行（`qcnn.py`）

### ConvSpec——层规格（论文 §3.3–3.4 的维度规则）

```python
@dataclass(frozen=True)
class ConvSpec:
    input_shape: tuple      # (H^l, W^l, D^l)：输入张量三维度
    kernel_shape: tuple     # (h, w, D^l, D^{l+1})：4-张量核
    cap: float = 5.0        # capReLU 上限 C（§5.1.6：为条件旋转提供已知的 max f）
    pool: int = 1           # 池化窗口 P（Eq. 39）
    pool_kind: str = "max"  # "max" 或 "average"（§5.2.2 的两种覆写规则）
```

`output_shape` 按论文 Eq. (4)：`H^{l+1} = H^l − h + 1`（无 padding/stride）；`pooled_shape` 按池化窗口整除输出空间维度（`__post_init__` 强制整除，5×5 输入配 2×2 核得 4×4 输出，恰好被 2 整除）。

### im2col——Eq. (6) 的矩阵化

```python
def im2col(x, spec):
    h, w, _ = spec.input_shape
    kh, kw, d, _ = spec.kernel_shape
    rows = []
    for i in range(h - kh + 1):        # 每个输出行对应一个感受野
        for j in range(w - kw + 1):
            row = []
            for dd in range(d):        # 逐通道堆叠（列优先向量化，论文 §3.4）
                for ki in range(kh):
                    for kj in range(kw):
                        row.append(tensor_get(x, spec.input_shape, i + ki, j + kj, dd))
            rows.append(row)
    return rows                        # A^l：行主序，每行一个感受野
```

每一行 `A^l_p` 是位置 (i,j) 感受野的向量化；{obj}`kernel_columns <oracq.algorithms.qml.qcnn.kernel_columns>` 用同一 (d, ki, kj) 顺序把每个核展开为 `F^l` 的一列——两者顺序一致，内积即为卷积输出 `Y^{l+1}_{p,q}`。

### convolution_forward——经典镜像

```python
a = im2col(x, spec)
f = kernel_columns(kernel, spec)
y = [[cap_relu(sum(ar * fc for ar, fc in zip(row, col, strict=True)), spec.cap)
      for col in f] for row in a]     # 卷积 + capReLU（§5.1.6）
return pool_tensor(y, spec)           # Eq. 39 的池化
```

`cap_relu = min(max(x, 0), cap)`——论文用带 cap 的 ReLU 保证 `max_pq f(Y_pq)` 先验已知（cap 常数），条件旋转（Eq. 28）因此可用固定归一化。

### QCNNQRAM——行树数据结构（§4.3 + §5.2.2）

```python
self.leaf_start = 1
while leaf_start < self.width:      # 完全二叉树，叶从 leaf_start 起
    leaf_start *= 2
self.tree_size = 2 * leaf_start     # 1-based 堆：根为 1，内部节点存平方部分和
```

`_write(p, r, value)` 写一个叶并沿祖先向上刷新平方部分和——单点更新的代价是对数级（论文 Theorem 4.3 的 O(log² n) 更新）。`update_with_pooling` 实现在线覆写：`max` 只在新值更高时写（论文 §5.2.2 "update the leaf if the new sampled value is higher"），`average` 以状态计数做均摊（"replace the actual value by the new averaged value"）。

## 3. 量子构件逐行（`qcnn_layer.py`）

### vector_angle_tables——数据到运行时 bank

```python
angles[(index << width) | node] = round(angle * scale) % (1 << angle_width)
```

逐行逐节点：`angle = 2·atan2(√S_right, √S_left)`（RY 半角约定，幅度比为 √(S_right/S_node)）；量化为 `angle_width` 位角字。符号单独入 `signs` bank（键 `(index<<width)|leaf`）。**换一层的数据只需重算这两个表——线路对象完全不变**（QRAM 语义）。

### qcnn_vector_prep——Eq. (15) 的行装载电路

```python
for depth in range(vector_width):          # 树层：bit 从高位到低位
    for prefix in range(1 << depth):
        node = (1 << depth) - 1 + prefix   # 树节点编号（与 qram_state_prep 编址一致）
        node_register = _constant_register(b, f"node_{depth}_{prefix}", vector_width, node)
        if depth:
            with b.control(b["target"][bit + 1:], prefix):   # 已确定的高位做前缀控制
                _prep_node_rotation(b, angles, node_register, bit, angle_width)
        else:
            _prep_node_rotation(b, angles, node_register, bit, angle_width)
        # 对称恢复 node_register（X 反转回去）
```

`_prep_node_rotation` 的三步：查询 `angles`（地址 {obj}`fuse(node_register, index) <oracq.infrastructure.ir.fuse>`——**index 在高位，与 bank 键一致**）；按位加权合成受控 RY（每个置位角字位 k 施加 `ry(2π·2^k/2^aw)`）；反查询复净 work。末尾的符号反冲：以已制备的 `target` 为地址查 `signs`，Z 旗标后反查——负分量只留下 −1 相位。

### qcnn_inner_product——Eq. (16)–(19)

```python
for register in (b["p"], b["q"], b["flag"]):
    b.h(register)                                  # 均匀叠加 Σ|p>|q> (1/√2)(|0>+|1>)
with b.control(b["flag"], 0):
    invoke(b, row_prep, "row", index=b["p"], target=b["vec"], work=b["work"])
with b.control(b["flag"], 1):
    invoke(b, column_prep, "col", index=b["q"], target=b["vec"], work=b["work"])
b.h(b["flag"])                                     # 幅度 √P_pq 落在 flag=0 分支
```

这正是论文 Eq. (16) 的三行推导：分支装载 |A_p⟩/|F_q⟩ → Hadamard 干涉 → flag=0 的测量概率 P_pq=(1+⟨A_p|F_q⟩)/2。验证对照 Eq. (20)：对 (p,q,vec,work) 求和后的 (p,q,flag=0) 联合概率 = (1+⟨A_p|F_q⟩)/(2·行数·列数)，实测最差偏差 4×10⁻⁵。

### qcnn_sampled_layer——Eq. (34)–(39) 的采样驱动

```python
values[(p, q)] = cap_relu(inner, spec.cap)          # Y_pq=(2P_pq-1)||A_p||||F_q|| 后激活
weight = {key: value * value for key, value in values.items()}
# 采样概率 ∝ f(Y)^2（Eq. 34），每次采样得三元组 (p, q, f(Y))：
if value < eta:      # Eq. 36：低于 η 的像素视为未采样
    continue
i, j = p // ow, p % ow                       # Eq. 7：行号还原 (i^{l+1}, j^{l+1})
index = (i // spec.pool) * pw + (j // spec.pool)   # Eq. 39：池化区域映射
pooled[index][q] = max(pooled[index][q], value)    # QRAM 在线覆写（§5.2.2）
```

{obj}`quantized_prepared_state <oracq.algorithms.qml.qcnn_layer.quantized_prepared_state>` 是电路语义的经典镜像（逐层量化角旋转 + 符号），驱动用它恢复全部 Y——电路与镜像在测试中逐振幅对照（1e-9），保证分布建模不引入额外近似。

## 4. 验证矩阵（全部通过）

| 实验 | 规模 | 结论 |
|---|---|---|
| im2col vs 朴素卷积 | 5×5×1, 核 2×2×1×2 | 逐元素 1e-12 |
| 前向+池化 vs 朴素 | 同上 + 2×2 max | 逐元素 1e-12 |
| QRAM 行/范数/覆写 | 16 行 | 精确；max 保高、average (2+4)/2=3 |
| 行制备电路 vs 归一化行 | 含负分量 | 最差 3e-4（=aw=10 角量化下界） |
| 电路 vs 量化镜像 | 同上 | 1e-9 |
| 内积概率 vs Eq. 20 | 2×2 行列 | 4e-5 |
| max 池化采样恢复 | 40000 样本 | 与经典最大值差 ≤5e-3 |
| η→∞ | 1000 样本 | 全零（Eq. 36） |
| average 收敛 | 60000 样本 | 到 f² 加权期望 ≤2e-2 |
