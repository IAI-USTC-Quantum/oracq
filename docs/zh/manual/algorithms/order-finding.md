# 量子求阶与因子后处理（Order Finding）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C1 · 模块 [`oracq.algorithms.basics.number_theory`](../../api/algorithms/basics/number_theory.rst) · 阶段 —

## 概述

对与模数 $m$ 互素的乘数 $a$，求满足 $a^r \equiv 1 \pmod m$ 的最小正整数 $r$（阶）。实现依据 Shor 的构造（[arXiv:quant-ph/9508027](https://arxiv.org/abs/quant-ph/9508027)）：把[模乘置换](modular-multiplication.md) $U_{a,m}$ 作用于本征态 $|1\rangle$ 的叠加分量，对相位寄存器做量子相位估计，读出 $s/r$ 的近似分数，再由经典连分数还原候选阶并尝试分解 $m$。阶为偶数且 $a^{r/2} \not\equiv -1 \pmod m$ 时，$\gcd(a^{r/2} \pm 1, m)$ 给出非平凡因子。

## 接口与输入模型

```python
order_finding(multiplier, modulus, *, precision=3, max_width=8)
factors_from_phase(value, precision, multiplier, modulus)
```

API 入口：{obj}`order_finding <oracq.algorithms.basics.number_theory.order_finding>`、{obj}`factors_from_phase <oracq.algorithms.basics.number_theory.factors_from_phase>`

- `multiplier` / `modulus`：与模乘相同的 CP 参数（经典整数直接参数化，无 oracle 输入）；底层模乘位宽预算由 `max_width` 控制。
- `precision`：QPE 相位寄存器位宽。
- {obj}`factors_from_phase <oracq.algorithms.basics.number_theory.factors_from_phase>` 为纯经典函数：`value` 是相位寄存器的整数读出（$0 \le \text{value} < 2^{\text{precision}}$），返回已验证的因子对（升序元组），或 `None` 表示该样本未能给出因子。

{obj}`order_finding <oracq.algorithms.basics.number_theory.order_finding>` 返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`，寄存器为 `target: Bits(n)`（$n$ 为模乘位宽）与 `phase: Bits(precision)`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"order_finding"` |
| `modulus` / `multiplier` | 调用参数回显 |
| `readout_register` | `"phase"`（读出后需经典后处理） |

## 实现要点

生成链为 {obj}`modular_multiply <oracq.algorithms.basics.number_theory.modular_multiply>` → {obj}`phase_estimation <oracq.algorithms.common.estimation.phase_estimation>`：先在 `target` 最低位加 X 制备 $|1\rangle$，再调用 QPE。$|1\rangle$ 在各阶循环本征态上均匀展开，故相位读出近似 $s/r$（$s$ 均匀），单个样本不保证给出阶——`factors_from_phase` 用 `Fraction(value, 2**precision).limit_denominator(modulus)` 取候选阶，逐一校验偶数性和 $a^r \equiv 1$，再对 $\gcd(a^{r/2}\pm1, m)$ 做验证，全部失败返回 `None`。乘数与模数不互素且公因子非平凡（$1 < \gcd(a,m) < m$）时跳过量子路径，直接返回由该公因子拆出的升序因子对。

适用边界：底层模乘是默认最多 8 位的有界置换合成（见[模乘置换](modular-multiplication.md)），因此求阶同样只面向小实例的接口与线路检查，**不代表可扩展的 Shor 模算术**。`precision` 不足时连分数可能还原不出正确的阶，调用方需按 $r < m$ 选择相位位宽并接受多次采样。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../../development/validation-plan.md` §2）：输出分布与经典后处理结果须与精确构造逐点相等。证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言（寄存器布局 `target/phase`、属性回显）。
- 数值：`AlgorithmExpansionTests.test_modular_multiplication_and_order` 同名见证 `test_modular_multiplication_total_permutation_and_order`——`order_finding(2, 3, precision=2)` 的相位分布逐点对拍：$p(0) = p(2) = 0.5$（阶 $r=2$，$s \in \{0,1\}$ 各占一半）；后处理路径由 `factors_from_phase(1, 2, 2, 15) == (3, 5)` 与 `factors_from_phase(0, 2, 2, 15) is None` 覆盖。
- 绑定：本算法无独立绑定见证（无开放声明入口）。

## 已知缺口与计划阶段

缺口沿用验证矩阵口径：小规模见证**已注明不可外推为完整 Shor**（保留原有口径）；多采样成功率、相位位宽与还原成功率的定量扫描未做。阶段列为 —（不属 V1–V4 的任何一档）。

## 相关链接

- 同模块：[模乘置换](modular-multiplication.md)
- 源码：`src/oracq/algorithms/basics/number_theory.py`
- API 参考：[模乘与求阶](../../api/algorithms/basics/number_theory.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_nt_qlss_sde.py`（nt_qlss_sde 组），全部在真实后端执行。

**实验设计**：(a) 对 7 个 $(a, m, p)$ 实例（$m \in \{3,5,7,15\}$，阶 $r \in \{2,3,4\}$，precision 3–5），在 reference 与 rir-pysparq 上运行完整求阶程序，phase 边缘分布对照独立理论：$|1\rangle$ 在 $r$ 个循环本征态上均布，每个本征相位 $s/r$ 经 QPE 产生 Dirichlet 核峰，总分布为 $\frac{1}{r}\sum_{s=0}^{r-1} D^2(y - 2^p s/r)$（numpy 独立计算）；$(2,3)$ 例另加 adapter-pysparq 与 OriginIR-ext 全振幅路径。同时报告连分数还原出完整阶的概率。(b) `factors_from_phase` 穷举：$m \in \{15,21,33,35\}$ × 全部 $a \in [2, m-2]$ × 全部 $2^8$ 个相位读出（precision=8，共 23552 次求值），任何返回因子对必须满足乘积为 $m$ 且两因子非平凡；成功率按 QPE 理论分布加权，并与独立教科书预言（$s$ 在 $[1,r)$ 均匀、候选阶 $d = r/\gcd(s,r)$、要求 $d$ 偶且 $a^d \equiv 1$ 且 $\gcd(a^{d/2}\pm1, m)$ 非平凡）逐模数对比——实测低于预言的差额是 precision=8 对较大的阶的连分数还原分辨率损失（真实物理，非实现缺陷）。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| order-finding-2-3-p3 | $r=2$ | 四路径 | TVD vs Dirichlet 理论 / 阶还原概率 | 4.4e-16 / 0.500 |
| order-finding-2-5-p4、3-5-p4 | $r=4$ | 双路径 | TVD / 阶还原概率 | 6.7e-16 / 0.500 |
| order-finding-2-7-p4、4-7-p4 | $r=3$ | 双路径 | TVD / 阶还原概率 | 8.6e-16 / 0.459 |
| order-finding-2-15-p5、7-15-p5 | $r=4$ | 双路径 | TVD / 阶还原概率 | 7.8e-16 / 0.500 |
| factors-from-phase-exhaustive | 4 模数 × 56 单位群 | 经典 oracle | 无效因子对 / 加权成功率均值 / 与教科书预言最小比值 | 0 / 0.2925 / 0.791 |

逐模数成功率（实测 vs 教科书预言）：$m=15$：0.500 vs 0.500；$m=21$：0.214 vs 0.233；$m=33$：0.185 vs 0.233；$m=35$：0.271 vs 0.318。奇阶或 $a^{r/2}\equiv-1$ 的固有失败单位占 16/56，与 Shor 构造的理论预期一致。

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

产物：`out/verification/nt_qlss_sde.json`（`order-finding-*`、`factors-from-phase-exhaustive` 共 8 个案例）。
