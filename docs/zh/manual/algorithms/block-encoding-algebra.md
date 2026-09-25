# 块编码组合代数（Block Encoding Algebra）

<a href="../../../en/index.html">English</a> · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.input_model.block_encoding`](../../api/algorithms/input_model/block_encoding.rst) · 阶段 V1

## 概述

块编码（block encoding）把矩阵 $A$ 藏进更大酉的角块：$U$ 的零信号块满足 $(\langle 0| \otimes I)\, U\, (|0\rangle \otimes I) = A/\alpha$，归一化 $\alpha$ 作为 RIR 模块属性 `be_alpha` 随操作保存，精度不属于语言核心。本模块提供 BE 的组合算术：线性组合（LCU）、张量、伴随、信号空间扩张、直和、Kronecker 和、布尔嵌入与小矩阵显式 Pauli 展开。上层算法（QSVT、qubitization、QLSS、低秩分解、ODE 求解器）都在这套代数上组装输入；最低层的 {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` 类型与二元基元 {obj}`identity <oracq.algorithms.input_model.operators.identity>` / {obj}`pauli_x <oracq.algorithms.input_model.operators.pauli_x>` / {obj}`zero <oracq.algorithms.input_model.operators.zero>` / {obj}`scale <oracq.algorithms.input_model.operators.scale>` / {obj}`linear_combination <oracq.algorithms.input_model.operators.linear_combination>` / {obj}`product <oracq.algorithms.input_model.operators.product>` 位于 `operators.py`。

## 接口与输入模型

Input model 统一为 BE：所有入口接收并返回 `BlockEncoding`（{obj}`reflect_zero <oracq.algorithms.input_model.block_encoding.reflect_zero>` 是 {obj}`Builder <oracq.infrastructure.builder.Builder>` 辅助）。

```python
lcu(terms)
tensor(a, b)
adjoint_be(a)
pad_signal(a, width)
kronecker_sum(a, b=None)
direct_sum(a, b)
projector(width, accepted)
truncated_shift(width, last)
pauli_word(word)
matrix_pauli_encoding(matrix, *, drop_tolerance=1e-12)
reflect_zero(builder, register, *, positive=False)
```

API 入口：{obj}`lcu <oracq.algorithms.input_model.block_encoding.lcu>`、{obj}`tensor <oracq.algorithms.input_model.block_encoding.tensor>`、{obj}`adjoint_be <oracq.algorithms.input_model.block_encoding.adjoint_be>`、{obj}`pad_signal <oracq.algorithms.input_model.block_encoding.pad_signal>`

主要入口的语义与输出归一化：

| 入口 | 语义 | 输出 alpha |
|---|---|---|
| {obj}`lcu(terms) <oracq.algorithms.input_model.block_encoding.lcu>` | $\sum_j c_j A_j$ 的 PREPARE/SELECT 组合（复相位经全局相位；零系数剔除，单项退化为 `scale`） | $\sum_j \lvert c_j\rvert\, \alpha_j$ |
| {obj}`tensor(a, b) <oracq.algorithms.input_model.block_encoding.tensor>` | $A \otimes B$ | $\alpha_A \alpha_B$ |
| {obj}`adjoint_be(a) <oracq.algorithms.input_model.block_encoding.adjoint_be>` | $A^\dagger$ | $\alpha_A$ |
| {obj}`kronecker_sum(a, b) <oracq.algorithms.input_model.block_encoding.kronecker_sum>` | $A \otimes I + I \otimes B$ | $\alpha_A + \alpha_B$ |
| {obj}`direct_sum(a, b) <oracq.algorithms.input_model.block_encoding.direct_sum>` | 同宽矩阵直和（单选择位） | $\alpha_A + \alpha_B$ |
| {obj}`pauli_word(word) <oracq.algorithms.input_model.block_encoding.pauli_word>` / {obj}`matrix_pauli_encoding(matrix) <oracq.algorithms.input_model.block_encoding.matrix_pauli_encoding>` | Pauli 词 BE / 小矩阵显式 Pauli LCU | 1 / LCU 归一化 |

`matrix_pauli_encoding` 按 $\operatorname{Tr}(P^\dagger M)/2^n$ 逐列相位累积解析计算展开系数，丢弃幅值低于 `drop_tolerance` 的项，全零时退化为 `zero`；仅接受 2 的幂维、至多 5 位的方阵，定位为小型应用的显式门实现，不宣称矩阵输入或经典展开具有量子加速。

## 实现要点

`lcu` 是核心组合子：selector 位宽 $\lceil\log_2 L\rceil$，权重 $\sqrt{\lvert c_j\rvert \alpha_j / \alpha}$ 经 {obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>` 制备，各项在自己的信号分区上受控调用，末尾逆制备；`operators.py` 的二元 `linear_combination` 是单 selector 的同构特例，两者共同覆盖组合语义。{obj}`pad_signal <oracq.algorithms.input_model.block_encoding.pad_signal>` 只允许放大信号空间（拒绝缩小），供同签名晚绑定。{obj}`projector <oracq.algorithms.input_model.block_encoding.projector>` / {obj}`truncated_shift <oracq.algorithms.input_model.block_encoding.truncated_shift>` 构造布尔嵌入的 $\alpha = 1$ BE，供 `direct_sum` 等结构化组装；`reflect_zero` 实现信号零态的正/负反射，是行走类算子的共享构件。

适用边界：所有组合只保证角块语义与 alpha 演算，不检查输入 BE 的数学声明；`be_alpha` 必须是有限正数，非法值（0、负数、`True`、inf）在构造期拒绝。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。与验证覆盖矩阵 `sparse.py` / `block_encoding.py` 行一致，三层证据都在 `tests/core/test_language.py:BlockEncodingTests`：

- 结构：类内构造与属性断言；`test_invalid_alpha` 覆盖 alpha 有效性（0、−1、`True`、inf 均抛 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`）。
- 数值：`test_unequal_normalization_and_signed_sum`（$\mathrm{lc}(-0.5,\ 2I,\ 2,\ 3X)$ 归一化 7、矩阵 $[[-1, 6], [6, -1]]$）、`test_complex_coefficients`（复系数相位进入角块）、`test_product_order_and_independent_signal_spaces`（乘积顺序与信号空间分区，$k = 2\cos(0.4)/\sqrt{2}$ 逐元素对拍）、`test_zero_and_zero_coefficient`（零系数与零矩阵退化）。见证技术为逐列模拟初始 `target` 基态、读零信号块幅度并乘回 alpha（与 `tests/core/witness.py` 的 `assert_block_equals` 同口径）。
- 绑定：`test_alpha_survives_ir_serialization`——`scale(4, identity(2))` 经 JSON {obj}`dumps <oracq.infrastructure.serialization.dumps>`/{obj}`loads <oracq.infrastructure.serialization.loads>` 往返后 `be_alpha` 保持 4。

## 已知缺口与计划阶段

无已知缺口，阶段 V1（验证覆盖矩阵 `sparse.py` / `block_encoding.py` 行缺口列为空）。

## 相关链接

- 源码：`src/oracq/algorithms/input_model/block_encoding.py`（组合子）与 `src/oracq/algorithms/input_model/operators.py`（`BlockEncoding` 类型与二元基元）
- 同族页面：[稀疏矩阵块编码](sparse-block-encoding.md)、[QSVT 矩阵求逆](qsvt-matrix-inversion.md)、[截断 Taylor 块编码](taylor-block-encoding.md)
- API 参考：[Block encoding 组合](../../api/algorithms/input_model/block_encoding.rst)
- 概念：[Oracle 与算子表示](../operators.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值验证脚本：`tests/verification/verify_blockencoding.py`（真实后端执行，无 mock、无 skip；2026-09-16 共 73 个案例全部通过），产物 `out/verification/blockencoding.json`。本页对应 `matrix-pauli-encoding-*`、`be-algebra-*`、`pauli-word-YZX`、`projector-w2-03`、`truncated-shift-w2-l3` 共 12 个案例。

实验设计：固定种子的随机 Hermitian / 非 Hermitian 小矩阵经 `matrix_pauli_encoding` 编码，七种组合子（tensor / adjoint / pad_signal / lcu 复系数 / {obj}`kronecker_sum <oracq.algorithms.input_model.block_encoding.kronecker_sum>` / {obj}`direct_sum <oracq.algorithms.input_model.block_encoding.direct_sum>` / product）各自组装；程序全幺正经 OriginIR-ext 导出 + UniQC `Circuit.to_matrix` 取出，`harness.effective_block` 提取零信号有效块并报告失败分支（信号 ≠ 0）泄漏上界；期望矩阵由 numpy 独立组装（Pauli 字约定 `word[0]` 作用于最低位）。`pauli_word`（无信号位，泄漏恰为 0）与布尔嵌入 `projector` / `truncated_shift` 同法与显式经典矩阵逐元对拍。各案例另以 reference / rir-pysparq / adapter-pysparq 逐列交叉（偏差均恰为 0），α 与独立计算的 Pauli l1 上界对拍。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `matrix-pauli-encoding-2x2-general` | 2×2 非 Hermitian | originir-ext + to_matrix，三后端交叉 | max_error | 1.3e-16 |
| `matrix-pauli-encoding-4x4-hermitian` | 4×4 Hermitian | 同上 | max_error | 4.6e-16 |
| 同上 | — | — | α − 独立 Pauli l1 | 0（恰等于上界） |
| `be-algebra-tensor` | 2×2 ⊗ 2×2 | originir-ext + to_matrix | max_error | 1.3e-16 |
| `be-algebra-adjoint` | 2×2 非 Hermitian | 同上 | max_error | 2.3e-16 |
| `be-algebra-pad-signal` | signal +2 位 | 同上 | max_error | 1.2e-16 |
| `be-algebra-lcu-complex` | 复系数 0.6+0.3j / −0.8 | 同上 | max_error | 2.4e-16 |
| `be-algebra-kronecker-sum` | A⊗I+I⊗B | 同上 | max_error | 2.5e-16 |
| `be-algebra-direct-sum` | 同宽直和 | 同上 | max_error | 1.1e-16 |
| `be-algebra-product` | A·C（非 Hermitian） | 同上 | max_error | 2.2e-16 |
| `pauli-word-YZX` | 3 位，α=1 | originir-ext + to_matrix | max_error / leakage | 0 / 0 |
| `projector-w2-03` / `truncated-shift-w2-l3` | 2 位 | 同上 | max_error | 0 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_blockencoding.py
```

产物：`out/verification/blockencoding.json`。
