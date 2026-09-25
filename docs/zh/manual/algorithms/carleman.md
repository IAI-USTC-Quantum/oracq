# Carleman 线性化（Carleman Linearization）

<a href="../../../en/index.html">English</a> · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.qnlss.carleman`](../../api/algorithms/qnlss/carleman.rst) · 阶段 V2

## 概述

把多项式非线性 ODE 提升为有限维线性系统。写

$$
u'=\sum_{p=0}^{D}F_p\,u^{\otimes p},\qquad u^{\otimes 0}=1,
$$

其中 $F_0$ 是常量强迫向量、$F_1$ 是线性算子、$F_p$（$p\ge 2$）是对整个张量输入线性的收缩。令 $y_k=u^{\otimes k}$，乘积法则给出

$$
y_k'=\sum_{p=0}^{D}\sum_{j=0}^{k-1}\bigl(I^{\otimes j}\otimes F_p\otimes I^{\otimes(k-j-1)}\bigr)\,y_{k+p-1}.
$$

截断阶 $K$ 保留 $y_0,\dots,y_K$，丢弃来源阶数大于 $K$ 的项，得到线性生成元 $G_K$；有限矩阵构造是确定的，逼近原非线性方程的程度取决于截断与适用条件（量子收敛/效率依据见 [Carleman 论文 arXiv:2011.03185v4](https://arxiv.org/abs/2011.03185v4)，另见[微分方程章 §6](../differential-equations.md)）。

## 接口与输入模型

```python
PolynomialODE(width, coefficients, initial, initial_norm=1.0)
carleman_lift(problem, *, cutoff=2)
carleman_initial(problem, *, cutoff=2)
carleman_qode(problem, time, linear_solver, *, cutoff=2)
```

API 入口：{obj}`PolynomialODE <oracq.algorithms.qnlss.carleman.PolynomialODE>`、{obj}`carleman_lift <oracq.algorithms.qnlss.carleman.carleman_lift>`、{obj}`carleman_initial <oracq.algorithms.qnlss.carleman.carleman_initial>`、{obj}`carleman_qode <oracq.algorithms.qnlss.carleman.carleman_qode>`

- `width`：物理寄存器位数 $n$（维数 $d=2^n$，上限 64）。
- `coefficients`：`((p, F_p 的 BE), ...)`，次数非负且不重复；$F_p$ 的 BE 宽度必须为 $\max(1,p)\cdot n$（补齐到方形的零信号块契约：前 $d$ 行有效、其余行声明为零）。
- `initial` / `initial_norm`：物理初态的 {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>`（SP）与原始向量范数；范数决定提升初态各阶的相对幅值，不是装饰属性。
- `linear_solver`：三参数协议 `(generator, initial, time) -> StateOracle`，通常注入 {obj}`linear_qode(...) <oracq.algorithms.qode.ode.linear_qode>`。

三个入口分别返回：提升生成元 $G_K$ 的 {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`（属性 `algorithm="carleman_lift"`、`cutoff`、`coefficient_assumption`）、提升初态的 `StatePreparation`、以及最终 {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>`（属性 `algorithm="carleman_qode"`、`cutoff`、`truncation_assumption="Carleman tail pending"`）。

## 实现要点

寄存器布局为补齐布局：`target = 数据区（K·n 位，低）+ level 寄存器（K.bit_length() 位，高）`；阶 $k$ 的有效数据只占前 $k\cdot n$ 位，其余为零——不是紧凑的 $1+d+\dots+d^K$ 地址，资源差别保留在报告中。

{obj}`carleman_lift <oracq.algorithms.qnlss.carleman.carleman_lift>` 对每个 $(k,p,\text{position})$（$1\le k\le K$、来源阶 $k+p-1\le K$、位置 $0..k-1$）生成一个放置 BE：level 寄存器从来源阶路由到输出阶（XOR 翻转差位），$F_p$ 作用在 position 起的连续 $p$ 个 $n$ 位组上，$p>1$ 时输出保留在组内最低 $n$ 位、零行移到数据区高端；所有项 LCU 求和得 $G_K$。

{obj}`carleman_initial <oracq.algorithms.qnlss.carleman.carleman_initial>` 制备归一化提升向量

$$
\frac{(1,\,u_0,\,u_0^{\otimes 2},\,\dots,\,u_0^{\otimes K})}{\sqrt{\sum_{k=0}^{K}\lVert u_0\rVert^{2k}}},
$$

按 level 受控地多次调用初态 oracle（work 预留 $K\times$ 初态 work 宽度）。{obj}`carleman_qode <oracq.algorithms.qnlss.carleman.carleman_qode>` 把 $G_K$ 与提升初态交给注入的线性求解器，再以 {obj}`select_subspace <oracq.algorithms.common.state_preparation.select_subspace>` 选中 level 1 通道（选择值 `1 << ((cutoff-1)*width)` 恰落在 level 寄存器最低位，同时要求高位数据位为零），选择条件并入 signal。

适用边界：$G_K$ 不保证耗散，即使原 PDE 或 $F_1$ 耗散。接 Schrödingerization 可直接传 $G_K$（另判恢复区）；接 LCHS 必须确认 $\mathrm{Hermitian}(-G_K)\succeq 0$，或做作用在整个提升系统（含 level 0）的整体移位 $\mu\ge\alpha_{G_K}$，移位代价应在应用配置中可见。张量尾项的截断误差未量化（`truncation_assumption` 标注 pending）。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。三层证据：

- 结构：`tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` 的 Carleman 分支——{obj}`PolynomialODE(1, ((1, a), (2, f2)), initial) <oracq.algorithms.qnlss.carleman.PolynomialODE>` 以抽象 BE 为系数、CBMD 为线性求解器时，开放槽 `input_A`/`input_b`/`input_F2` 全部保留，且存在 `carleman_row_level == 1`、`carleman_column_level == 2` 的放置模块，序列化往返成立。
- 数值：`tests/core/test_differential.py` 定位为结构测试；截断收敛的数值对拍归入统一的 ODE 收敛基准（见已知缺口）。
- 绑定：同类注册与契约检查——系数 BE 宽度、次数唯一性与初态布局在 `PolynomialODE.__post_init__` 生成期校验，违例抛 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`。

## 已知缺口与计划阶段

统一的"解析可解 ODE 族"收敛基准缺失（同一问题过全部求解器，含 Carleman 截断阶扫描），归入阶段 V2 的收敛性扫描框架；与 `validation-coverage.md` 的 carleman.py 行一致。

## 数值验证

论文级数值实验见 `tests/verification/verify_ode.py`（ode 组，本文件覆盖 `carleman.py` 的部分）。经典参考全部独立：提升生成元 $G_K$ 按乘积法则在补齐布局上用 numpy 逐基组装、提升初态按张量幂解析构造、截断线性系统用 `scipy.linalg.expm`、非线性精确解用 `scipy.integrate.solve_ivp`（Riccati，rtol=1e-12）；量子程序在真实后端 reference、rir-pysparq、adapter-pysparq 与 OriginIR-ext（UniQC 全振幅/`to_matrix`）上执行。

**实验设计**：(a) 提升块——Riccati 问题 $u'=-u+u^2$（$F_1=-I$，$F_2$ 为收缩矩阵 $C e_{i_0+2i_1}=\delta_{i_0,i_1}e_{i_0}$，$d=2$），截断 $K=2$，提取 $G_K$ 的零信号块（幺正 `to_matrix` + `effective_block`，并与 reference 基态扫描交叉）对照独立组装；(b) 提升初态——$(1,u_0,u_0^{\otimes2})/Z$（$u_0=(0.6,0.8)$、$r=0.5$）四路径全振幅穷举；(c) 端到端——`carleman_qode` 注入三参数协议求解器（验证脚本内的最小截断 Taylor $e^{Gt}$ 块编码组装，真实量子程序，degree 3），$t=0.2$，物理通道对全堆叠 numpy 仿真；(d) 截断趋势——$K=1,2$ 量子运行的截断误差对 solve_ivp 精确 Riccati 解，$K=3$ 经典补点。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| carleman-lift-block | $K$=2、10 量子位 | originir+to_matrix、reference | $G_K$ 块逐元素最大偏差 | 3.4e-16 |
| carleman-initial-state | 4 量子位 | 四路径 | 全振幅最大偏差 | 0.0 |
| carleman-riccati-endtoend | 18 量子位 | reference、rir | 实现误差（对独立仿真） | 2.8e-17 |
| 同上 | | | Taylor 余项（degree 3，信息性） | 1.4e-4 |
| 同上 | | | Carleman 截断误差（对精确非线性解） | 1.1e-3 |
| 同上 | | | 量子方向 vs 截断线性精确解方向 | 4.8e-5 |
| carleman-cutoff-trend | $t=0.2$ | 量子（K=1,2）+ 经典（K=3） | 截断误差 K=1 / K=2 / K=3 | 9.3e-3 / 1.1e-3 / 1.05e-4 |

截断误差每升一阶约降一个量级，是 Carleman 截断收敛的直接数值证据。成功子空间概率（信号全零）同时与经典值逐点对拍一致。

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

产物：`out/verification/ode.json`（`carleman-*` 共 5 个案例）。

## 相关链接

- 源码：`src/oracq/algorithms/qnlss/carleman.py`
- 教程：[为同一个线性问题替换 QODE 方法](../../tutorials/differential-equations.md)
- API 参考：[Carleman 线性化](../../api/algorithms/qnlss/carleman.rst)
- 相关页：[QODE 问题对象与协议](qode-problem.md) · [LCHS](lchs.md) · [Schrödingerization](schrodingerization.md)（典型组合 Carleman → 线性求解器） · [CBMD](cbmd.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
