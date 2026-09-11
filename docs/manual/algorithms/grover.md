# Grover 搜索（Grover Search）

> 类别 C3 · 模块 `pyqecclang.algorithms.search` · 阶段 V1

## 概述

在 $2^{\text{width}}$ 维基态空间中搜索被相位 oracle 标记的目标：每次 Grover 迭代把态向标记子空间旋转，若干次后命中概率达到峰值。实现依据 Brassard 等人的振幅放大/估计框架（[arXiv:quant-ph/0005055](https://arxiv.org/abs/quant-ph/0005055)，见算法目录"实现依据"节）；迭代次数由调用方给定，本入口不自动选次。

## 接口与输入模型

```python
grover(phase_oracle, width, *, iterations=1, preparation=None)
```

- `phase_oracle`：标记目标基态的相位 `Operation`（FO，相位函数 oracle），接口为 `target`，可另有 `work`。
- `width`：搜索空间的目标位宽。
- `iterations`：非负迭代次数（默认 1）。
- `preparation`：初态制备（SP，须支持伴随调用）；省略时用 `uniform_state(width)`。宽度与 `width` 不一致抛 `ValidationError`。

返回 `StateOracle`：`target` 为搜索结果，`signal` 保留相位查询和制备的工作空间。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"grover"` |
| `iterations` | 迭代次数回显 |
| `validation_stage` | `"paradigm"` |

同模块的相关入口：

- `phase_from_database(database)`：把单输出位的 `XorDatabase`（QRAM 查找）适配为相位 oracle——加载、work 上 Z、卸载，产物标注 `phase_oracle`。
- `grover_iterate(preparation, marked)`：返回迭代算子 $Q=A(2|0\rangle\langle 0|-I)A^\dagger S_{\mathrm{good}}$（标记来自 `phase_marks`），供 [振幅估计](qae.md) 组装 QPE。

## 实现要点

单次迭代为相位 oracle → $A^\dagger$ → 关于全零态的正反射（`reflect_zero`，`positive=True`）→ $A$；初态反射同时包括 `target` 与制备工作区（docstring 口径）。寄存器布局：`signal` 前段是相位 oracle 的 work（`phase_work` 位），后段是制备的 work；相位 oracle 无 work 时前段宽度为 0。迭代次数以 `repeat(iterations)` 保存，不在生成或 JSON 序列化阶段展开。

产物包装为 `StateOracle`，可继续作为态 oracle 组合（例如 `amplify_success` 的输入）。适用边界：相位 oracle 的 work 须在其调用后自身复原（反射只覆盖 target 与制备工作区）；迭代次数过多会越过峰值降低命中概率，与 [振幅放大](amplitude-amplification.md) 的选次约束相同。

## 验证方案

类别 C3（判定准则见 `docs/development/validation-plan.md` §2：输出分布等于闭式期望）。三层证据与 `docs/development/validation-coverage.md` 的 `search.py` 行一致：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言（迭代次数非负、初态宽度一致等在生成期校验）。
- 数值：模块行登记 `AlgorithmExpansionTests.test_amplification_of_known_quarter_probability`（同一迭代原理的直接见证，见 [振幅放大](amplitude-amplification.md)）；`grover()` 完整电路未单列数值条目——其迭代算子 `grover_iterate` 经 QAE 的 `test_amplitude_estimation_half_probability` 间接见证（`estimation.py` 行），均匀初态特例由展示目录 `grover` 条目演示（`grover(phase_marks(2, [3]), 2)`，读出 target=3、signal=0）。
- 绑定：模块行登记为"—"。

## 已知缺口与计划阶段

无已知缺口（`search.py` 行缺口列为"—"），阶段 V1。

## 相关链接

- 源码：`src/pyqecclang/algorithms/search.py`
- API 参考：[搜索与振幅放大](../../api/algorithms/search.rst)
- 同组页面：[振幅放大](amplitude-amplification.md)、[振幅估计](qae.md)、[量子计数](quantum-counting.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
