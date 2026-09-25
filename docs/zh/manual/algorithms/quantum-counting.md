# 量子计数（Quantum Counting）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C3 · 模块 [`oracq.algorithms.common.estimation`](../../api/algorithms/common/estimation.rst) · 阶段 V2

## 概述

估计大小为 $N=2^n$ 的搜索空间中被标记元素的个数 $t$。量子计数不是独立电路，而是振幅估计的直接应用（同一 Brassard 框架，[arXiv:quant-ph/0005055](https://arxiv.org/abs/quant-ph/0005055)）：均匀叠加上好状态概率恰为 $a=t/N$，做标准振幅估计后乘回 $N$：

$$
\hat t = N\cdot\hat a.
$$

本仓库将其登记为展示目录条目（`applications/gallery.py` 的 `quantum_counting`），复用 {obj}`amplitude_estimation <oracq.algorithms.common.estimation.amplitude_estimation>` 入口。

## 接口与输入模型

没有单独的入口函数；配方是：

```python
operation = amplitude_estimation(uniform_state(n), marked, precision=p)
```

- {obj}`uniform_state(n) <oracq.algorithms.input_model.oracles.uniform_state>`：均匀态制备（SP），使 $a=t/N$。
- `marked`：被标记基态的整数编号集合，$t=|M|$。
- `precision`：相位寄存器位数，范围 1..63。

返回对象的寄存器与属性同[振幅估计](qae.md)（`target`/`work`/`phase`，`decoder` 为 `"sin(pi*phase/2**precision)**2"`）：读出 `phase` 后经 {obj}`amplitude_from_phase <oracq.algorithms.common.estimation.amplitude_from_phase>` 解码、再乘 $2^n$ 得 $\hat t$。展示实例取 $n=2$、`marked=(3,)`、`precision=4`（$t=1$，$a=1/4$）。

## 实现要点

电路生成与 QAE 完全相同，差别全在经典侧解码（乘 $N$）。$\hat t$ 不必是整数，四舍五入到合法范围由宿主决定；估计精度由 QAE 栅格密度（`precision`）决定，$t/N$ 恰落栅格时点估计精确。

## 验证方案

类别 C3（判定准则见 `docs/development/validation-plan.md` §2）。量子计数在 `docs/development/validation-coverage.md` 中随 `estimation.py` 行登记（算法列含"量子计数"），无单独见证条目：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`。
- 数值：电路见证即 QAE 的 `AlgorithmExpansionTests.test_amplitude_estimation_half_probability`（见[振幅估计](qae.md)）；展示实例由 `applications/gallery.py` 的 `quantum_counting` 条目演示读出（"将振幅估计乘以 4，得到标记个数的估计"）。
- 绑定：模块行登记为"—"。

## 已知缺口与计划阶段

随 QAE 一致：置信区间声明未见证（只验证了点估计），点估计乘 $N$ 的误差传播同样未对拍；`estimation.py` 模块阶段 V2。

## 相关链接

- 源码：`src/oracq/algorithms/common/estimation.py`（电路）、`src/oracq/applications/gallery.py`（展示条目）
- API 参考：[相位、振幅与重叠估计](../../api/algorithms/common/estimation.rst)
- 同组页面：[振幅估计](qae.md)、[Grover 搜索](grover.md)、[量子相位估计](qpe.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

实验设计：$n=3$（$N=8$）、marked $=\{1,5,7\}$（$t=3$，$a=3/8$）、precision $=4$（栅格 $M=16$）。按上文配方生成 `amplitude_estimation(uniform_state(3), marked, precision=4)`，在四条后端路径（reference、rir-pysparq、adapter-pysparq、originir-ext）上读出 phase 寄存器的完整分布，对照独立 Dirichlet 核参考（本征相位 $\pm\theta/\pi$、$\theta=\arcsin\sqrt{3/8}$ 上的 QPE 栅格卷积），并按 $\hat t=N\sin^2(\pi y/M)$ 解码计数。

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| quantum-counting-n3-t3 | N=8, t=3, M=16 | 4 路径 | 相位分布 TVD | 5.7e-15 |
| 〃 | 〃 | 〃 | 峰值估计 $\hat t$（误差） | 2.4693（0.531，≤1 栅格步） |
| 〃 | 〃 | 〃 | 中心栅格质量（QAE 下界 $8/\pi^2$） | 0.8529（≥ 0.8106） |

真值 $\theta M/\pi\approx3.357$ 落在栅格点 3 与 4 之间，分布双峰（$y=3,13$ 各约 0.325）与理论位置一致；点估计精度由栅格密度决定，与"$\hat t$ 不必为整数"的口径吻合。

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

产物：`out/verification/search_walks.json`（案例 `quantum-counting-n3-t3`）。
