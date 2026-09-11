# 量子计数（Quantum Counting）

> 类别 C3 · 模块 `pyqecclang.algorithms.estimation` · 阶段 V2

## 概述

估计大小为 $N=2^n$ 的搜索空间中被标记元素的个数 $t$。量子计数不是独立电路，而是振幅估计的直接应用（同一 Brassard 框架，[arXiv:quant-ph/0005055](https://arxiv.org/abs/quant-ph/0005055)）：均匀叠加上好状态概率恰为 $a=t/N$，做标准振幅估计后乘回 $N$：

$$
\hat t = N\cdot\hat a.
$$

本仓库将其登记为展示目录条目（`applications/gallery.py` 的 `quantum_counting`），复用 `amplitude_estimation` 入口。

## 接口与输入模型

没有单独的入口函数；配方是：

```python
operation = amplitude_estimation(uniform_state(n), marked, precision=p)
```

- `uniform_state(n)`：均匀态制备（SP），使 $a=t/N$。
- `marked`：被标记基态的整数编号集合，$t=|M|$。
- `precision`：相位寄存器位数，范围 1..63。

返回对象的寄存器与属性同[振幅估计](qae.md)（`target`/`work`/`phase`，`decoder` 为 `"sin(pi*phase/2**precision)**2"`）：读出 `phase` 后经 `amplitude_from_phase` 解码、再乘 $2^n$ 得 $\hat t$。展示实例取 $n=2$、`marked=(3,)`、`precision=4`（$t=1$，$a=1/4$）。

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

- 源码：`src/pyqecclang/algorithms/estimation.py`（电路）、`src/pyqecclang/applications/gallery.py`（展示条目）
- API 参考：[相位、振幅与重叠估计](../../api/algorithms/estimation.rst)
- 同组页面：[振幅估计](qae.md)、[Grover 搜索](grover.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
