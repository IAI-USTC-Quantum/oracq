# 量子化行走（Qubitization Walk）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.common.transforms`](../../api/algorithms/common/transforms.rst) · 阶段 V1

## 概述

给定厄米矩阵 $A$ 的块编码（block encoding）$U$，构造量子化行走算子

$$
W = (2\Pi - I)\,U, \qquad \Pi = |0\rangle\langle 0|_{\text{signal}} \otimes I_{\text{target}},
$$

即先作用 $U$，再对信号寄存器的 $|0\rangle$ 态做正反射。$W$ 在 $A$ 每个谱值 $x$ 对应的二维不变子空间上作用为旋转，是把块编码转化为谱函数变换的基本迭代单元（Low–Chuang 2019 的 qubitization 构造，见 `../../development/algorithm-coverage.md` A2；QSVT 框架见 Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)）。本页只覆盖单步行走算子的组装；相位序列的交替调用见 [QSVT 相位序列](qsvt-sequence.md)。

## 接口与输入模型

```python
qubitization_walk(a)
```

API 入口：{obj}`qubitization_walk <oracq.algorithms.common.transforms.qubitization_walk>`

- `a`：{obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`，被行走矩阵的块编码（input model 为 BE）。

返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`，寄存器为 `target`（宽度 `a.width`）与 `signal`（宽度 `a.signal_qubits`）。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"qubitization_walk"` |

输入 BE 的归一化 `alpha` 不经属性传播：返回值是裸 `Operation` 而非 `BlockEncoding`，谱变量 $x$ 相对输入归一化 $\alpha$ 定义，即 $x = \lambda/\alpha$。

## 实现要点

生成策略为两次结构调用：{obj}`invoke <oracq.algorithms.input_model.oracles.invoke>` 把输入 BE 挂到 `target | signal` 上，随后 {obj}`reflect_zero(b, signal, positive=True) <oracq.algorithms.input_model.block_encoding.reflect_zero>` 实现 $2\Pi - I$（`positive=True` 使 $|0\rangle_{\text{signal}}$ 分支取 $+1$、其余分支取 $-1$）。重复调用由调用方（如 {obj}`qsvt_sequence <oracq.algorithms.common.transforms.qsvt_sequence>` 或低秩分解流水线）编排，本函数不内置 {obj}`Repeat <oracq.infrastructure.ir.Repeat>` 结构。

适用边界：输入必须是块编码；稀疏 oracle 或低秩分解产物需先经 `sparse.py` / `lowrank.py` 等适配为 BE。行走算子本身不施加任何多项式变换，其谱性质只在配合相位序列或投影测量时显现。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../../development/validation-plan.md` §2）。当前证据（与验证覆盖矩阵 `transforms.py` 行一致）：

- 结构：`tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` 与 `test_trotter_only_input_does_not_need_block_encoding` 覆盖同模块组装链的协议与寄存器契约。
- 数值：`tests/core/test_qsvt.py:PhaseSynthesisTests.test_convention_matches_qsvt_sequence` 在随机相位下钉死同一信号/反射约定的电路零信号块与 {obj}`qsp_response <oracq.algorithms.common.qsvt.qsp_response>` 逐点一致（delta = 1e-10），行走步所嵌入的序列骨架由此间接见证。
- 绑定：无独立绑定见证（输入已要求具体 BE）。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（序列约定逐点一致 1e-10）。单步行走算子的独立旋转性质（二维不变子空间上的转角）暂无直接见证，由序列级约定测试覆盖。

## 相关链接

- 源码：`src/oracq/algorithms/common/transforms.py`
- API 参考：[矩阵变换序列](../../api/algorithms/common/transforms.rst)
- 同族页面：[QSVT 相位序列](qsvt-sequence.md)、[Oblivious 振幅放大](oblivious-amplification.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_fourier.py`（真实后端执行，无模拟替身）。输入 BE 用 gate 级绑定小实例：{obj}`matrix_pauli_encoding <oracq.algorithms.input_model.block_encoding.matrix_pauli_encoding>` 编码 1 比特厄米矩阵 $A = \begin{pmatrix} 0.5 & 0.3 \\ 0.3 & -0.1 \end{pmatrix}$（独立显式 Pauli 展开得 $\alpha = 0.8$，谱值 $\lambda/\alpha \approx 0.780, -0.280$）。行走算子全幺正经 OriginIR-ext + UniQC `to_matrix` 取出，与 $(2\Pi - I)U$ 逐元素对比，其中 $U$ 为同一 BE 程序的全幺正（并与 reference 路径逐基态列组装的幺正交叉对拍）；谱性质独立检验：行走幺正本征角须落入 $\{\pm\arccos(\lambda/\alpha)\} \cup \{0, \pi\}$，零信号块须等于 $A/\alpha$。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `qubitization-walk-unitary-spectrum` | target 1 + signal 2（8 维） | originir-ext + UniQC to_matrix，reference 交叉 | unitary_max_error | 1.7e-16 |
| 同上 | — | — | be_path_crosscheck（两路径 BE 幺正偏差） | 1.1e-16 |
| 同上 | — | — | zero_block_max_error（$\Pi W\Pi$ vs $A/\alpha$） | 1.8e-16 |
| 同上 | — | — | spectrum_max_deviation（本征角 vs $\pm\arccos(\lambda/\alpha)$） | 4.4e-16 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_fourier.py
```

产物：`out/verification/fourier.json`。

## 数值验证

实验设计：构造零信号块为 $x=\cos(\pi/6)$ 的厄米（Householder 型）块编码 $U=R_y(2\theta)Z\otimes I_{\rm target}$（$\theta=\pi/6$，target/signal 各 1 位）。三层对照：

1. 旋转性质（补上"二维不变子空间转角无直接见证"的缺口）：行走重复 $k=1\ldots4$ 次，零信号概率对照 Chebyshev 闭式 $T_k(x)^2=\cos^2(k\arccos x)$，即序列 $0.75,\ 0.25,\ 0,\ 0.25$（$k=3$ 为精确零点），四条后端路径（reference、rir-pysparq、adapter-pysparq、originir-ext）；
2. 幺正分解：UniQC `Circuit.to_matrix` 核对 $W=(2\Pi-I)U$；
3. 谱：$W$ 的本征值应为 $\mathrm e^{\pm i\pi/6}$（两个 target 扇区各贡献一对，二重简并）。

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| qubitization-walk-chebyshev | 2 qubits, k=1..4 | 4 路径 + to_matrix | 零信号概率最大误差 | 5.6e-16 |
| 〃 | 〃 | 〃 | $\|W-(2\Pi-I)U\|_{\max}$ | 2.1e-16 |
| 〃 | 〃 | 〃 | 本征相位误差（vs $\pm\pi/6$） | 3.3e-16 |

前提说明：$T_k$ 旋转关系要求 BE 幺正厄米；对非厄米 $U$，正确的 qubitized 迭代需交替 $U$ 与 $U^\dagger$（QSVT 序列即如此处理），简单幂 $W^k$ 不服从 Chebyshev——本验证因此采用厄米构造。

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

产物：`out/verification/search_walks.json`（案例 `qubitization-walk-chebyshev`）。
