# 量子卷积神经网络（QCNN）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C4 · 模块 [`oracq.algorithms.qml.qcnn`](../../api/algorithms/qml/qcnn.rst) 与 [`oracq.algorithms.qml.qcnn_layer`](../../api/algorithms/qml/qcnn_layer.rst) · 论文 arXiv:1911.01117（ICLR 2020）

## 概述

QCNN 把经典 CNN 的一层搬到量子线路：卷积经 im2col 化为矩阵乘 A^ℓ·F^ℓ=Y^{ℓ+1}（Eq. 5–7），量子侧用 QRAM 装载行/列向量、Hadamard 型电路估计内积（Eq. 16–19），采样与 l∞ 层析按像素值平方加权抽取高值像素（Eq. 28–38），池化在 QRAM 在线更新中完成（§5.2.2）。本实现覆盖前向的全部四类操作（卷积、激活、池化、采样读出），小规模语义对照经典 CNN 逐元素验证。

`correctness="pending"` 记录的实现边界：幅度估计寄存器编码（Eq. 22–25）与幅度放大（Eq. 30–33）按其精确分布建模（采样概率 f(Y)²/Σf²，Eq. 34），条件旋转的角度算术（θ=2·asin(√(f/C))）依赖的 mathfunc `asin` 原语已具备，但 AE 的受控幂+逆 QFT 编码未落成线路；反向传播（Thm 6.1）未实现。

## 接口

```python
ConvSpec(input_shape, kernel_shape, cap, pool, pool_kind)  # 层规格
im2col(x, spec)                 # 输入张量 -> A^ℓ（Eq. 6）
kernel_columns(kernel, spec)    # 4-张量核 -> F^ℓ 列（Eq. 6）
cap_relu(value, cap)            # capReLU（§5.1.6）
convolution_forward(x, kernel, spec)   # 经典镜像：卷积+capReLU+池化
QCNNQRAM(rows)                  # 行树 QRAM（§4.3/§5.2，含池化覆写）

vector_angle_tables(rows, angle_width)  # 行/列 -> 角度树+符号 bank（运行时表）
qcnn_vector_prep(count, width, angle_width)   # |p>|0>->|p>|A_p/||A_p||>（Eq. 15）
qcnn_inner_product(rows, cols, width, aw)     # Hadamard 内积电路（Eq. 16–19）
quantized_prepared_state(vector, aw)          # 角度量化镜像（电路语义的经典复现）
qcnn_sampled_layer(x, kernel, spec, *, samples, eta, seed)  # 采样驱动（Eq. 34–39）
```

API 入口：{obj}`ConvSpec <oracq.algorithms.qml.qcnn.ConvSpec>`、{obj}`im2col <oracq.algorithms.qml.qcnn.im2col>`、{obj}`kernel_columns <oracq.algorithms.qml.qcnn.kernel_columns>`、{obj}`cap_relu <oracq.algorithms.qml.qcnn.cap_relu>`

## 实现要点

**行/列制备（{obj}`qcnn_vector_prep <oracq.algorithms.qml.qcnn_layer.qcnn_vector_prep>`）**：多路旋转树的每个节点由 QRAM 角度 bank 驱动——地址为 {obj}`fuse(node, index) <oracq.infrastructure.ir.fuse>`（bank 键 `(index<<width)|node`），查询、按位加权合成受控 RY（角度 2π·2^k/2^aw）、反查询三步；负分量经符号 bank 的 Z 反冲写入相位。换数据只换内存表，线路不变（QRAM 语义）。深度 d 的节点数为 2^d，总查询 2·width 次。

**Hadamard 内积（{obj}`qcnn_inner_product <oracq.algorithms.qml.qcnn_layer.qcnn_inner_product>`）**：p、q、flag 三寄存器均匀叠加，flag=0 分支装载行向量、flag=1 分支装载列向量（两套受控制备），最后对 flag 施加 Hadamard。测量 (p,q,flag=0) 的概率为 (1+⟨A_p|F_q⟩)/(2·行数·列数)——Eq. 20 的精确复现；内积可正可负，概率恒正。

**采样驱动（{obj}`qcnn_sampled_layer <oracq.algorithms.qml.qcnn_layer.qcnn_sampled_layer>`）**：先经量化镜像恢复全部 Y_pq=(2P_pq−1)‖A_p‖‖F_q‖ 与 capReLU 值，再按 f²/Σf² 采样（Eq. 34 的分布）；每次采样得三元组 (p,q,f(Y))，f<η 的像素视为未采样置零（Eq. 36）；按 Eq. 39 映射进池化区域后以 QRAM 覆写规则（max 保高、average 均摊）聚合。

## 验证

`tests/core/test_qcnn.py`（9 例，全绿）：

| 层 | 断言 | 结果 |
|---|---|---|
| 经典基底 | im2col 内积 vs 朴素三重循环；前向+池化 vs 朴素实现；QRAM 行/范数/覆写 | 逐元素精确（1e-12） |
| 行制备 | 电路幅度 vs 行/‖·‖（含负分量）；电路 vs 量化镜像 | 最差 3e-4（aw=10 角量化下界） |
| 内积电路 | 测量概率 vs Eq. 20 | 最差 4e-5 |
| 采样驱动 | max 池化 40000 样本恢复区域最大值；η→∞ 全零；average 收敛到 f² 加权期望 | 5e-3 / 精确 / 2e-2 |

## 已知边界

- 幅度估计的寄存器编码（AE + 中值提纯）未落成线路——当前 Y 恢复走量化镜像（与电路语义一致），条件旋转与幅度放大按分布建模。落成后 `qcnn_sampled_layer` 的镜像部分应替换为寄存器读出。
- 反向传播（论文 §6）未实现。
- 行数/列数要求 2 的幂（小规模验证约束；一般形状经填充可达）。

## 相关链接

- 源码：`src/oracq/algorithms/qml/qcnn.py`（整体流程）与 `src/oracq/algorithms/qml/qcnn_layer.py`（量子构件）
- API 参考：[量子卷积神经网络](../../api/algorithms/qml/qcnn.rst)、[量子卷积神经网络的量子构件](../../api/algorithms/qml/qcnn_layer.rst)
- 使用手册：[QCNN 逐步指南](../qcnn-walkthrough.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
