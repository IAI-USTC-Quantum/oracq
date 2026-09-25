# MaxCut QAOA（QAOA for MaxCut）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C4 · 模块 [`oracq.algorithms.optimization.variational`](../../api/algorithms/optimization/variational.rst) · 阶段 V3

## 概述

给定无向加权图 $G=(V,E)$，MaxCut 要求把顶点划分成两部分，最大化两端分属不同子集的边权之和。QAOA 从均匀叠加出发，交替施加 cost 演化与 mixer 演化，层数与角度 $(\gamma,\beta)$ 由经典外层选择；cost/mixer 分层依据 Farhi 等人的原始算法（[arXiv:1411.4028](https://arxiv.org/abs/1411.4028)，本手册算法目录"实现依据"已引用）。cost 算子按位串 $z$ 取值

$$
C(z) = \sum_{(u,v,w)\in E} \frac{w}{2}\,(1 - z_u z_v),
$$

即被割开的边权之和；mixer 取逐位 $X$。本函数生成固定角度的量子电路，读取 `target` 得到一个割的候选位串；不运行经典优化器。

## 接口与输入模型

```python
qaoa_maxcut(width, edges, gammas, betas)
```

API 入口：{obj}`qaoa_maxcut <oracq.algorithms.optimization.variational.qaoa_maxcut>`

- `width`：图的顶点数，也是 target 位宽，范围 1..64。
- `edges`：`(u, v, weight)` 三元组；`u`、`v` 必须落在 `[0, width)`，权重为非负有限实数，不接受自环。
- `gammas` / `betas`：各层 cost / mixer 演化角，两个列表非空且等长，角度须为有限实数。

图的边表与变分角度都是经典数据直接参数化，input model 为 CP。返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`，寄存器只有 `target`（width 位），无工作位。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"qaoa_maxcut"` |
| `layers` | QAOA 层数，即 `len(gammas)` |

## 实现要点

初始化为整根 `target` 上的 Hadamard 均匀叠加。每个 QAOA 层先做 cost 演化、再做 mixer 演化：

- cost 层逐边施加 XOR–`Rz(-gamma·weight)`–XOR 共轭，把相位写到 $Z_u Z_v$ 分量上，配合全局相位 $-\gamma w/2$，单边门恰为 $e^{-i\gamma w(1-Z_uZ_v)/2}$；各边的 ZZ 门相互可交换，整层实现 $e^{-i\gamma C}$。
- mixer 层对每个比特施加 `rx(2·beta)`。

寄存器布局：顶点 i 对应 `target` 第 i 位。参数在生成期全部校验（宽度上限、自环、负权重、gamma/beta 长度不一致或为空均抛 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`）。适用边界：只生成给定角度的电路，角度初始化、期望割值评估与外层优化循环由调用方组织；本模块不含经典优化器。

## 验证方案

类别 C4（启发式/优化语义，判定准则见 `../../development/validation-plan.md` §2：严格优于随机基线，且小实例达到已知最优/理论分数）。三层证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言；同类 `test_bad_inputs_fail_at_generation` 覆盖自环边 `(0, 0, 1)` 等生成期违例。
- 数值：`AlgorithmExpansionTests.test_qaoa_single_edge_optimal_layer`——2 顶点单边图（权重 1.0）取 $\gamma=\pi/2$、$\beta=\pi/8$ 的单层 QAOA，两个最优割态 `|01>` 与 `|10>` 的概率之和等于 1（places = 10），即单层即达最优割。
- 绑定：本模块无独立绑定见证（输入为经典参数，无抽象 oracle 槽位）。

## 已知缺口与计划阶段

与 `validation-coverage.md` 的 `variational.py` 行一致："达到小实例最优割"的强见证缺失——当前只有单边图的解析最优层，多顶点图上与穷举最优割值的分布级对拍未做。计划在 V3（catalog 登记新算法 + tests/integration 真实后端对拍扩展）接入。

## 相关链接

- 源码：`src/oracq/algorithms/optimization/variational.py`
- 同模块页面：[硬件高效拟设](variational-ansatz.md)、[VQE 测量电路](vqe.md)
- API 参考：[变分算法电路](../../api/algorithms/optimization/variational.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_misc_algorithms.py`（misc_algorithms 组），全部在真实后端上执行。

**实验设计**：(a) 单边图（2 顶点、单位权）单层 QAOA 的解析锚点 $\gamma = \pi/2$、$\beta = \pi/8$；(b) 4 顶点环 C4（单位权）单层 QAOA——角度由 numpy 精确模拟的网格搜索选定（$\gamma = 0.85$、$\beta = 0.45$，经典外层按模块契约自理），量子分布与 numpy 精确模拟（cost 对角相位 + mixer 逐位 $e^{-i\beta X}$）逐点对拍，最优割为 $\lvert 0101\rangle$ 与 $\lvert 1010\rangle$（随机基线 2/16 = 0.125）。后端路径：`reference`、`originir-ext` 全振幅对拍。

**关键指标**：

| 案例 | 规模 | 路径 | 分布 TVD | 最优割概率 | 相对随机基线 |
|---|---|---|---|---|---|
| qaoa-single-edge-optimal | 2 比特、1 层 | reference | — | 1.0000（解析最优） | — |
| qaoa-c4-distribution | 4 比特、1 层 | reference + originir-ext | 1.8e-16 | 0.5379 | 4.30× |

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

产物：`out/verification/misc_algorithms.json`（24 个案例全过，本页对应 `qaoa-*` 两个案例）。
