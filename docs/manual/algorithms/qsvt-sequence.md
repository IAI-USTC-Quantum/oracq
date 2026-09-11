# QSVT 相位序列（QSVT Phase Sequence）

> 类别 C2 · 模块 `pyqecclang.algorithms.transforms` · 阶段 V1

## 概述

给定矩阵 $A$ 的块编码（block encoding）$U$ 与一列实数相位 $\Phi = (\varphi_0, \dots, \varphi_d)$，显式组装 QSVT 序列：在 $A$ 的每个奇异值 $x$ 对应的二维不变子空间上，零信号块实现多项式

$$
p(x) = \bigl[S(\varphi_0)\, W(x)\, S(\varphi_1)\, W(x) \cdots W(x)\, S(\varphi_d)\bigr]_{00},
\qquad
W(x) = \begin{pmatrix} x & s \\ s & -x \end{pmatrix},\ s = \sqrt{1 - x^2},\
S(\varphi) = \operatorname{diag}(e^{i\varphi}, e^{-i\varphi}),
$$

相位按时间顺序排列（$\varphi_0$ 最先作用），共 $d$ 次 BE 调用、$d + 1$ 个相位。可实现的 $(P, Q)$ 满足充要条件：$\deg P \le d$、$\deg Q \le d - 1$、$P$ 的奇偶性为 $d \bmod 2$、$Q$ 的奇偶性为 $(d-1) \bmod 2$，且多项式恒等式 $P\bar{P} + (1-x^2)Q\bar{Q} \equiv 1$ 成立（特别地必有 $|P(\pm 1)| = 1$）。序列约定取自 QSVT 框架（Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)，见 `qsvt.py` 模块 docstring）；本函数是框架的电路组装层，合法相位由 [QSP 相位合成](qsp-phase-synthesis.md) 计算并自检。

## 接口与输入模型

```python
qsvt_sequence(a, phases)
```

- `a`：`BlockEncoding`，被变换矩阵的块编码（input model 为 BE）。
- `phases`：可迭代的实数相位，时间正序（$\varphi_0$ 最先作用）；内部转换为浮点元组，本身不做正确性校验。

返回 `Operation`（裸操作而非 `BlockEncoding`），寄存器为 `target`（宽度 `a.width`）与 `signal`（宽度 `a.signal_qubits`）。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"qsvt_sequence"` |
| `validation_stage` | `"paradigm"`（见证级别标注） |

BE 归一化不经属性传播；谱变量 $x$ 相对输入归一化 $\alpha$ 定义（$x = \lambda/\alpha$），幅度语义由调用方补全。

## 实现要点

每个相位实现为一对全局相位门：`global_phase(-φ)` 加上 `signal` 零分支条件下的 `global_phase(2φ)`，净效果是信号零分支取 $e^{+i\varphi}$、其余分支取 $e^{-i\varphi}$，即 $S(\varphi)$。相位之间交替调用输入 BE：下标为偶数的相位之后正向 `invoke`，下标为奇数的相位之后伴随（`adjoint`）调用，共 $d$ 次。`a.signal_qubits == 0` 的退化情形有专门分支（相位退化为无条件全局相位）。

适用边界：输入必须是块编码，稀疏 oracle 或 QRAM 需先经适配层构造 BE；本函数接受任意实相位，序列的数值正确性由调用方负责（合法相位的计算与往返自检见 [QSP 相位合成](qsp-phase-synthesis.md)）。序列直接以复多项式 $P$ 为零信号块、不含实部提取；实目标的实部提取由 `qsvt.py` 标准变换族经 $(U_\Phi + U_{-\Phi})/2$ 的 LCU 组合完成。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../../development/validation-plan.md` §2）。当前证据（与验证覆盖矩阵 `transforms.py` 行一致）：

- 结构：`tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` 与 `test_trotter_only_input_does_not_need_block_encoding` 覆盖同模块组装链的协议与寄存器契约。
- 数值：`tests/core/test_qsvt.py:PhaseSynthesisTests.test_convention_matches_qsvt_sequence`——固定种子随机相位（6 个）下 2×2 对角 BE 的电路零信号块与 `qsp_response` 逐点一致（delta = 1e-10），把电路组装层与数学约定互相钉死。
- 绑定：无独立绑定见证（输入已要求具体 BE）。

## 已知缺口与计划阶段

无已知缺口（与验证覆盖矩阵 `transforms.py` 行一致），阶段 V1 见证已齐（序列约定逐点一致 1e-10）。

## 相关链接

- 源码：`src/pyqecclang/algorithms/transforms.py`
- API 参考：[矩阵变换序列](../../api/algorithms/transforms.rst)
- 同族页面：[量子化行走](qubitization-walk.md)、[Oblivious 振幅放大](oblivious-amplification.md)、[QSP 相位合成](qsp-phase-synthesis.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
