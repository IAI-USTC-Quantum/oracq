# 量子化行走（Qubitization Walk）

> 类别 C2 · 模块 `pyqecclang.algorithms.transforms` · 阶段 V1

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

- `a`：`BlockEncoding`，被行走矩阵的块编码（input model 为 BE）。

返回 `Operation`，寄存器为 `target`（宽度 `a.width`）与 `signal`（宽度 `a.signal_qubits`）。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"qubitization_walk"` |

输入 BE 的归一化 `alpha` 不经属性传播：返回值是裸 `Operation` 而非 `BlockEncoding`，谱变量 $x$ 相对输入归一化 $\alpha$ 定义，即 $x = \lambda/\alpha$。

## 实现要点

生成策略为两次结构调用：`invoke` 把输入 BE 挂到 `target | signal` 上，随后 `reflect_zero(b, signal, positive=True)` 实现 $2\Pi - I$（`positive=True` 使 $|0\rangle_{\text{signal}}$ 分支取 $+1$、其余分支取 $-1$）。重复调用由调用方（如 `qsvt_sequence` 或低秩分解流水线）编排，本函数不内置 Repeat 结构。

适用边界：输入必须是块编码；稀疏 oracle 或低秩分解产物需先经 `sparse.py` / `lowrank.py` 等适配为 BE。行走算子本身不施加任何多项式变换，其谱性质只在配合相位序列或投影测量时显现。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../../development/validation-plan.md` §2）。当前证据（与验证覆盖矩阵 `transforms.py` 行一致）：

- 结构：`tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` 与 `test_trotter_only_input_does_not_need_block_encoding` 覆盖同模块组装链的协议与寄存器契约。
- 数值：`tests/core/test_qsvt.py:PhaseSynthesisTests.test_convention_matches_qsvt_sequence` 在随机相位下钉死同一信号/反射约定的电路零信号块与 `qsp_response` 逐点一致（delta = 1e-10），行走步所嵌入的序列骨架由此间接见证。
- 绑定：无独立绑定见证（输入已要求具体 BE）。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（序列约定逐点一致 1e-10）。单步行走算子的独立旋转性质（二维不变子空间上的转角）暂无直接见证，由序列级约定测试覆盖。

## 相关链接

- 源码：`src/pyqecclang/algorithms/transforms.py`
- API 参考：[矩阵变换序列](../../api/algorithms/transforms.rst)
- 同族页面：[QSVT 相位序列](qsvt-sequence.md)、[Oblivious 振幅放大](oblivious-amplification.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
