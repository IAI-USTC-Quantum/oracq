# 振幅放大（Amplitude Amplification）

<a href="../../../en/index.html">English</a> · **简体中文**

> 类别 C3 · 模块 [`oracq.algorithms.common.search`](../../api/algorithms/common/search.rst) · 阶段 V1

## 概述

对任意态 oracle 的成功子空间做相干放大：设 {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>` 的制备 $A$ 以概率 $a$ 产出成功态（成功条件 signal==0），迭代算子

$$
Q = A\,(2|0\rangle\langle 0|-I)\,A^\dagger\, S_{\mathrm{good}}
$$

每次把成功幅度旋转约 $2\theta$（$\sin^2\theta=a$）。实现依据 Brassard 等人的振幅放大框架（[arXiv:quant-ph/0005055](https://arxiv.org/abs/quant-ph/0005055)，见算法目录"实现依据"节）。它是 Grover 搜索的推广：初态不必均匀，标记信息由 signal 寄存器而非相位 oracle 承担。

## 接口与输入模型

```python
amplify_success(state, *, iterations=1)
```

API 入口：{obj}`amplify_success <oracq.algorithms.common.search.amplify_success>`

- `state`：支持伴随调用的 `StateOracle`（`target`/`signal` 接口的态 oracle；signal 指示成功子空间，属 FO 幅度型访问）。非 `StateOracle` 输入在生成期抛 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`。
- `iterations`：非负放大次数（默认 1）。次数需要结合输入成功概率选择；过多迭代可能降低成功概率（docstring 口径），本入口不自动选次。

返回 `StateOracle`：公开接口与输入相同，成功条件仍为 signal==0。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"amplitude_amplification"` |
| `success_condition` | `"signal == 0"` |

## 实现要点

先正向调用一次 $A$ 制备初态；随后每次迭代为 {obj}`reflect_zero(signal) <oracq.algorithms.input_model.block_encoding.reflect_zero>`（对 signal==0 的好子空间相位翻转，充当 $S_{\mathrm{good}}$）→ $A^\dagger$ → 关于 target+signal 全零态的正反射 → $A$。与 {obj}`grover_iterate <oracq.algorithms.common.search.grover_iterate>` 的算子形状一致，差别只在标记来源：这里 $S_{\mathrm{good}}$ 由 signal 寄存器本身承担，因此适用于任意来源的成功指示（如 QRAM 查询结果的判定位、算法工作区的检验位）。迭代以 `repeat(iterations)` 保存，不在生成或文本序列化阶段展开。

## 验证方案

类别 C3（判定准则见 `docs/development/validation-plan.md` §2：输出分布等于闭式期望）。三层证据与 `docs/development/validation-coverage.md` 的 `search.py` 行一致：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言（{obj}`require_instance <oracq.algorithms.input_model.contracts.require_instance>` 校验 StateOracle、迭代次数非负）。
- 数值：`AlgorithmExpansionTests.test_amplification_of_known_quarter_probability`——构造 signal 上 `ry(2π/3)` 的态 oracle，成功概率 $P(\text{signal}=0)=\cos^2(\pi/3)=1/4$；默认一次迭代后 $P(\text{signal}=0)=1$（$\theta=\pi/6$，一次迭代转到 $\sin^2(3\theta)=1$），places=10。
- 绑定：模块行登记为"—"。

见证技术为精确态矢量模拟取 signal 边际分布后对拍闭式，无采样断言。

## 已知缺口与计划阶段

无已知缺口（`search.py` 行缺口列为"—"），阶段 V1。

## 相关链接

- 源码：`src/oracq/algorithms/common/search.py`
- API 参考：[搜索与振幅放大](../../api/algorithms/common/search.rst)
- 同组页面：[Grover 搜索](grover.md)（均匀初态 + 相位 oracle 的特例）、[振幅估计](qae.md)（同一迭代算子的对偶读出）
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

实验设计：构造成功子空间为 `signal==0` 的态 oracle（`target` 两位均匀叠加，信号位 $R_y(3\pi/4)$，初始成功概率 $a=\cos^2(3\pi/8)=\sin^2(\pi/8)\approx0.1464$，故 $\theta_a=\pi/8$），对 {obj}`amplify_success <oracq.algorithms.common.search.amplify_success>` 在 $k=0\ldots3$ 上测零信号概率，对照 $\sin^2((2k+1)\theta_a)$，四条后端路径（reference、rir-pysparq、adapter-pysparq、originir-ext）。序列覆盖放大（$k=1,2$ 达 $0.8536$）与过冲回落（$k=3$ 回 $0.1464$）。

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| amplify-success-curve | target 2 + signal 1，k=0..3 | 4 路径 | 成功概率最大误差 | 1.9e-15 |

各迭代数理论值：0.146447、0.853553、0.853553、0.146447，逐点复现。

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

产物：`out/verification/search_walks.json`（案例 `amplify-success-curve`）。
