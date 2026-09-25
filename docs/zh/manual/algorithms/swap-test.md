# SWAP 检验（Swap Test）

[English](../../../index.html) · **简体中文**

> 类别 C3 · 模块 [`oracq.algorithms.common.estimation`](../../api/algorithms/common/estimation.rst) · 阶段 V2

## 概述

测量两个纯态的重叠模方 $F=|\langle a|c\rangle|^2$：受控交换把两态的对称分量干涉到探针上，

$$
P(\text{probe}=0) = \frac{1+F}{2}.
$$

该操作不执行测量；采样与统计由宿主完成。

## 接口与输入模型

```python
swap_test(first, second)
```

API 入口：{obj}`swap_test <oracq.algorithms.common.estimation.swap_test>`

- `first` / `second`：两个同宽态的制备（SP），要求零输入与干净工作区（{obj}`checked_state_preparation <oracq.algorithms.input_model.interfaces.checked_state_preparation>` 契约，生成期检查）。宽度不一致抛 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`。

返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`，寄存器为 `left`、`right`、`left_work`、`right_work`、`probe`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"swap_test"` |
| `readout_register` | `"probe"` |

## 实现要点

两个制备只正向调用；条件交换（`swap` 于 `control(probe)` 之下）由算法生成。`left`/`right` 各带独立的 work 寄存器，复净由 SP 契约保证——标注 `clean_work=False` 的脏制备在生成期被拒绝。读出按上式从 probe 零概率解码 $F$。

## 验证方案

类别 C3（判定准则见 `docs/development/validation-plan.md` §2：输出分布等于闭式期望）。三层证据与 `docs/development/validation-coverage.md` 的 `estimation.py` 行一致：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`，其中 `test_measurement_algorithms_require_clean_preparation` 覆盖脏制备拒绝路径（`clean_work=False` 抛 `ValidationError`）。
- 数值：`AlgorithmExpansionTests.test_swap_test_overlap`——{obj}`basis_state(1, 0) <oracq.algorithms.input_model.oracles.basis_state>` 与 `basis_state(1, 1)` 两两组合：同一基态（$F=1$）时 $P(\text{probe}=0)=1$，正交基态（$F=0$）时为 $1/2$，places=10。
- 绑定：模块行登记为"—"。

见证技术为精确态矢量模拟取 probe 边际分布后对拍闭式，无采样断言。

## 已知缺口与计划阶段

本算法无单独缺口；`estimation.py` 登记的唯一缺口属于 QAE 的置信区间声明（见[振幅估计](qae.md)），模块整体处于阶段 V2。

## 相关链接

- 源码：`src/oracq/algorithms/common/estimation.py`
- API 参考：[相位、振幅与重叠估计](../../api/algorithms/common/estimation.rst)
- 同组页面：[Hadamard 检验](hadamard-test.md)（单酉期望的对应读出）
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

`tests/verification/verify_estimation.py` 在真实后端上对本接口做概率级数值验证（共 4 个案例，全部通过）。实验设计：

- 实例：{obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>` 制备的纯态对，覆盖 $F=1$（同态）、$F\approx 0$（正交）、$F=0.36$（部分重叠）与双比特态 $F=0.125$；重叠 $F=|\langle a|c\rangle|^2$ 由 numpy 对输入幅度向量独立计算。
- 后端路径：reference、rir-pysparq、adapter-pysparq、originir-ext 四条；$P(\text{probe}=0)$ 对照 $(1+F)/2$。

| 案例 | 规模 | 路径 | 指标值 |
|---|---|---|---|
| 同态 $(0.6,0.8)$ | 各 1 位 | 全部四条 | $F=1$，$P(0)=1$，误差 $2.2\times10^{-16}$ |
| 正交 $(0.6,0.8)$ vs $(0.8,-0.6)$ | 各 1 位 | 全部四条 | $F\approx 0$，$P(0)=0.5$，误差 $1.1\times10^{-16}$ |
| 部分重叠 $|0\rangle$ vs $(0.6,0.8)$ | 各 1 位 | 全部四条 | $F=0.36$，$P(0)=0.68$，误差 $2.2\times10^{-16}$ |
| 双比特态 | 各 2 位 | 全部四条 | $F=0.125$，$P(0)=0.5625$，误差 $1.1\times10^{-16}$ |

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_estimation.py
```

产物路径：`out/verification/estimation.json`（案例名前缀 `swap-test-`）。
