# 振幅估计（Amplitude Estimation）

> 类别 C3 · 模块 `pyqecclang.algorithms.estimation` · 阶段 V2

## 概述

估计态制备 $A$ 的好子空间概率：给定标记集合 $M\subset\{0,1\}^n$，$a=\sum_{i\in M}|\langle i|A|0\rangle|^2$。实现依据 Brassard 等人的振幅放大/估计框架（[arXiv:quant-ph/0005055](https://arxiv.org/abs/quant-ph/0005055)，另见算法目录"实现依据"节）：对 Grover 迭代算子

$$
Q = A\,(2|0\rangle\langle 0|-I)\,A^\dagger\, S_{\mathrm{good}}
$$

做相位估计。$Q$ 的本征值为 $e^{\pm 2i\theta}$，$\sin^2\theta = a$；读出相位栅格值 $y$ 后按 $\hat a=\sin^2(\pi y/2^p)$ 解码。查询复杂度 $O(1/\varepsilon)$，相对经典 Monte Carlo 的 $O(1/\varepsilon^2)$ 呈二次改进。

## 接口与输入模型

```python
amplitude_estimation(preparation, marked, *, precision=3)
amplitude_from_phase(value, precision)
```

- `preparation`：零输入、复净工作区的态制备（SP），且必须支持伴随与受控调用——生成期经 `checked_state_preparation(adjoint=True, controlled=True)` 检查契约。
- `marked`：目标空间中好状态的整数编号集合。
- `precision`：相位寄存器位数，范围 1..63。

返回 `Operation`，寄存器为 `target`、`work`、`phase`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"amplitude_estimation"` |
| `readout_register` | `"phase"` |
| `decoder` | `"sin(pi*phase/2**precision)**2"` |

读出 `phase` 后调用 `amplitude_from_phase(value, precision)` 解码：返回 $\sin^2(\pi\cdot\text{value}/2^p)\in[0,1]$，`value` 取 $0..2^p-1$。

## 实现要点

生成链为制备 → `phase_estimation(grover_iterate(preparation, marked))`：Grover 迭代由 `search.py` 的 `grover_iterate` 组装（标记来自 `phase_marks` 生成的相位 oracle，见 [Grover 搜索](grover.md)），相位估计见[量子相位估计](qpe.md)。

镜像对称：本征相位 $+2\theta$ 与 $-2\theta$ 对应栅格点 $y$ 与 $2^p-y$，解码值相同——相位读出本身不区分两者。有限精度读出可能对应多个近似概率，采样与统计处理（含置信区间）由宿主负责，这是 docstring 的明确口径。

积分模块（`integration.py` 的 Heinrich 求和/积分）以本入口为读出核：比较器构造使好状态概率恰为 $E[v]/2^w$，再做标准振幅估计。

## 验证方案

类别 C3（判定准则见 `docs/development/validation-plan.md` §2）。三层证据与 `docs/development/validation-coverage.md` 的 `estimation.py` 行一致：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言；制备的零输入/复净/伴随/受控契约违例在生成期抛错。
- 数值：`AlgorithmExpansionTests.test_amplitude_estimation_half_probability`——`uniform_state(1)` 上标记 `(1,)`（$a=1/2$ 恰在栅格上）：相位分布全部质量落在 2 与 6 两点（$p(2)+p(6)=1$，places=10），且 `amplitude_from_phase(2, 3)` 与 `amplitude_from_phase(6, 3)` 都精确等于 0.5。
- 绑定：模块行登记为"—"。

见证技术为精确态矢量模拟取相位边际分布，逐点对拍闭式（无采样断言）。

## 已知缺口与计划阶段

与 `docs/development/validation-coverage.md` 一致：QAE 置信区间声明未见证——现有见证只验证了点估计（$a=1/2$ 的栅格精确情形）。validation-plan §2 将"估计值落在理论置信区间"列为 QAE 类的标准见证技术，有限精度相位读出的误差分布与置信区间对拍尚未落地。`estimation.py` 模块因此处于阶段 V2；补齐后须同步覆盖矩阵。

## 相关链接

- 源码：`src/pyqecclang/algorithms/estimation.py`
- API 参考：[相位、振幅与重叠估计](../../api/algorithms/estimation.rst)
- 同组页面：[量子相位估计](qpe.md)、[量子计数](quantum-counting.md)、[Grover 搜索](grover.md)、[振幅放大](amplitude-amplification.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
