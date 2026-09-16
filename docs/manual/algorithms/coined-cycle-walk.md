# 周期格点硬币行走（Coined Cycle Walk）

> 类别 C1 · 模块 `pyqecclang.algorithms.walks` · 阶段 V1

## 概述

离散时间硬币量子行走（coined quantum walk）的最小实例：在周期为 $2^n$ 的环形格点 $\mathbb{Z}_{2^n}$ 上，一位 Hadamard 硬币决定每步移动 $+1$ 或 $-1$。硬币行走是量子行走的标准模型之一，本模块给出其一维周期版本；任意图上以 oracle 输入驱动的行走见 [Szegedy 量子行走](szegedy-walk.md)。

## 接口与输入模型

```python
cycle_walk(width, *, steps=1)
```

- `width`：`position` 位宽 $n$，周期长度为 $2^n$，上限 64。
- `steps`：非负步数。
- 输入模型为 CP：位宽与步数经典参数化，无 oracle 输入；行走空间初态由调用方制备，零输入对应位置 0、硬币 0。

返回 `Operation`，寄存器为 `position: Bits(n)` 与 `coin: Bits(1)`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"coined_cycle_walk"` |
| `steps` | 调用参数回显 |

## 实现要点

每步先更新硬币再条件移动：对 `coin` 作用 H，随后 `coin == 0` 分支把 `position` 加 1、`coin == 1` 分支加 $2^n - 1$（即 $-1 \bmod 2^n$）；`position` 以无符号解释运算，越界自然回绕。步数循环用 IR 的 Repeat 结构生成，不在生成期展开。适用边界：周期固定为 2 的幂、硬币固定为 Hadamard、移动固定为 $\pm 1$；需要一般图或可替换 oracle 时应改用 `graph_walks.py` 的框架。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../../development/validation-plan.md` §2）：作用酉须与行走定义逐点相等。证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言。
- 数值：`AlgorithmExpansionTests.test_ansatz_zero_angles_and_walk_one_step` 的行走部分——`cycle_walk(2)` 从零输入走一步后，态为 $(|1,0\rangle + |3,1\rangle)/\sqrt{2}$（基态按 `(position, coin)` 记），两个基态幅度均为 $1/\sqrt{2}$。
- 绑定：本算法无独立绑定见证（无开放声明入口）。

## 已知缺口与计划阶段

无已知缺口（验证矩阵缺口列为 —），阶段 V1。

## 相关链接

- 图上的行走：[图邻接 oracle](adjacency-oracle.md)、[Szegedy 量子行走](szegedy-walk.md)、[MNRS 量子行走搜索](mnrs-search.md)
- 源码：`src/pyqecclang/algorithms/walks.py`
- API 参考：[量子行走](../../api/algorithms/walks.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

实验设计：周期 $N=4$（width=2）与 $N=8$（width=3）的环上 Hadamard coined walk，从 $|0\rangle_{\rm position}|0\rangle_{\rm coin}$ 出发扫描步数 $s=0\ldots8$，四条后端路径（reference、rir-pysparq、adapter-pysparq、originir-ext）逐振幅对照独立 numpy 参考（coin H 后按 0/1 条件 $\pm1$ 移位的逐步模拟），位置边缘分布另算 TVD。注：$s=0$ 时程序不含门，UniQC 对零门线路无 qubit mapping，故零步情形不含 originir-ext 路径。

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| cycle-walk-w2 | N=4, s=0..8 | 4 路径 | 振幅最大误差 / 位置分布最大 TVD | 5.6e-17 / 1.1e-16 |
| cycle-walk-w3 | N=8, s=0..8 | 4 路径 | 振幅最大误差 / 位置分布最大 TVD | 5.6e-17 / 1.1e-16 |

信息性指标（$s=8$ 的平均环距离，量子弹道输运 vs 经典扩散）：$N=8$ 时量子 3.0 vs 经典对称随机游走 1.875；$N=4$ 时量子 0.0 vs 经典 1.0（小环上 8 步恰好回到原点，呈 revival）。两者定性差异符合离散时间量子行走的已知行为。

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

产物：`out/verification/search_walks.json`（案例 `cycle-walk-w2`、`cycle-walk-w3`）。
