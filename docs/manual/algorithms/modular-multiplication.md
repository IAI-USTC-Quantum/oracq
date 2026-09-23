# 模乘置换（Modular Multiplication）

> 类别 C1 · 模块 `pyqecclang.algorithms.basics.number_theory` · 阶段 —

## 概述

给定与模数 $m$ 互素的乘数 $a$，在 $w$ 位寄存器上实现可逆映射

$$
U_{a,m}\,|x\rangle =
\begin{cases}
|a \cdot x \bmod m\rangle, & x < m,\\
|x\rangle, & x \ge m,
\end{cases}
$$

即 Shor 求阶所需的模乘酉的有限规模版本。它是[量子求阶](order-finding.md)的底层构件；经典因子后处理见该页。

## 接口与输入模型

```python
modular_multiply(multiplier, modulus, *, width=None, max_width=8)
```

- `multiplier`：与 `modulus` 互素的正整数（input model 为 CP，乘数与模数直接以经典参数给出，无 oracle 输入）。
- `modulus`：模数，不小于 2。
- `width`：目标位宽；省略时取能容纳 `modulus - 1` 的最小位宽。
- `max_width`：置换合成预算，默认 8，最多允许 12。

返回 `Operation`，寄存器为 `target: Bits(width)`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"modular_multiply"` |
| `implementation_scope` | `"bounded permutation synthesis"` |
| `modulus` | 调用参数回显 |

## 实现要点

生成策略是枚举 $2^w$ 个基态上的完整置换，分解为不相交循环后逐对转置（`_transposition`）合成，因此任意基态上的作用与定义逐点一致。乘数在生成期先对模数约化；`modulus > 2^width` 或 $\gcd(a, m) \neq 1$ 时直接抛出 `ValidationError`。

适用边界：置换枚举个数随位宽指数增长，`max_width` 默认 8、上限 12 是有意的合成预算。**当前模乘使用有上限的置换合成，默认最多 8 位；它用于检查求阶接口与线路，不代表已实现可扩展的 Shor 模算术**（沿用算法目录页的既有口径）。需要大规模模算术的场景应换用可扩展构造，本模块不提供。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../../development/validation-plan.md` §2）：作用酉须与真值表逐点相等。证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言；`test_bad_inputs_fail_at_generation` 覆盖乘数与模数不互素（`modular_multiply(2, 4)`）等生成期违例。
- 数值：`AlgorithmExpansionTests.test_modular_multiplication_total_permutation_and_order`——对 `modular_multiply(2, 5)`（位宽 3）穷举全部 8 个基态，逐点验证 $x<5$ 映到 $2x\bmod 5$、$x\ge 5$ 保持不变，目标基态幅度精确为 1。
- 绑定：本算法无独立绑定见证（无开放声明入口）。

## 已知缺口与计划阶段

缺口沿用验证矩阵口径：小规模见证**已注明不可外推为完整 Shor**（保留原有口径）；可扩展模算术构造未实现。阶段列为 —（不属 V1–V4 的任何一档）。

## 相关链接

- 同模块：[量子求阶与因子后处理](order-finding.md)
- 源码：`src/pyqecclang/algorithms/basics/number_theory.py`
- API 参考：[模乘与求阶](../../api/algorithms/basics/number_theory.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_nt_qlss_sde.py`（nt_qlss_sde 组），全部在真实后端执行；经典预言机为独立的数论置换真值（$x \mapsto ax \bmod m$ 当 $x<m$，其余基态不变），不经过被测实现的任何辅助函数。

**实验设计**：(a) 对 5 个实例 $(a, m) \in \{(2,5), (3,5), (2,7), (2,9), (7,100)\}$（位宽 3–7），用叠加态（target 加 H 后调用模乘）在**一次运行内穷举全部 $2^w$ 个基态输入**，输出与经典置换重标记的均匀叠加逐振幅对拍；路径为 reference、rir-pysparq、adapter-pysparq 与 OriginIR-ext（UniQC 全振幅态向量，工作区 0）。(b) 对 $(2,5)$ 做 OriginIR-ext `Circuit.to_matrix` 全幺正检查：工作区复净后的有效块须等于经典置换矩阵，且向工作区的泄漏为零。

**关键指标**：

| 案例 | 规模 | 路径 | max_error |
|---|---|---|---|
| modmul-superposition-2-5 / 3-5 / 2-7 | $2^3$ 基态穷举 | 四路径 | 5.6e-17 |
| modmul-superposition-2-9 | $2^4$ 基态穷举 | 四路径 | 8.3e-17 |
| modmul-superposition-7-100 | $2^7$ 基态穷举 | 四路径 | 4.2e-17 |
| modmul-unitary-2-5 | $8 \times 8$ 全幺正 | originir-ext + to_matrix | 0.0（泄漏 0.0） |

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

产物：`out/verification/nt_qlss_sde.json`（`modmul-*` 共 6 个案例）。
