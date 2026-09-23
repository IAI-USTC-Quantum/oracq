# Fokker–Planck 输入模型（Fokker–Planck Input Model）

> 类别 C2 · 模块 `pyqecclang.algorithms.qode.sde` · 阶段 V3

## 概述

一维 Itô 随机微分方程 $dx=a(x)\,dt+\sqrt{2D(x)}\,dW$ 的概率密度满足 Fokker–Planck 方程

$$
p'=-\partial_x(a\,p)+\partial_{xx}(D\,p).
$$

本模块在均匀网格上做零通量有限体积离散，得到列和为零的离散生成元 $G$（约定 $p'=Gp$）；$G$ 经显式 Pauli 展开成为 `BlockEncoding` 后组装 `QODEProblem` / `LinearODE`，交给现有线性 ODE 求解器（LCHS 等）。生成元的耗散性是输入模型的声明，不由语言证明。模块同时提供纯 Python 的经典见证工具（矩阵指数、Euler 演化、矩、离散与连续稳态参考），用于数值对拍。

## 接口与输入模型

```python
FokkerPlanckProblem(drift, diffusion, grid)
problem.generator_matrix() / generator_encoding() / qode_problem(initial=None)
problem.linear_ode(initial=None)
sde_state_preparation(probabilities, *, implementation="gate", angle_width=8, work_width=0)
sde_state_angles(probabilities, *, angle_width=8)
matrix_exponential(matrix, time=1.0)
evolve_distribution(matrix, initial, time, *, steps=None)
distribution_moments(points, probabilities, orders=(1, 2))
stationary_distribution(problem)
boltzmann_distribution(problem)
```

- `FokkerPlanckProblem`：input model 为 CP（漂移、扩散与网格坐标直接参数化），产出 BE + SP。`drift` / `diffusion` 可为常量或逐点向量（扩散要求 $D\ge 0$），`grid` 须是严格递增的均匀网格；派生属性 `points`、`spacing`、`size`、`width`（网格点数须为二的幂）。
- `generator_encoding()`：显式 Pauli 展开的小网格 BE，仅支持 `width <= 5`（不超过 32 点），大实例需要访问 oracle，不宣称量子加速。
- `qode_problem(initial=None)`：返回 `QODEProblem`（`dissipative=True`，`evidence="fokker_planck_zero_flux_finite_volume; dissipative caller-declared"`），默认初态为均匀叠加。
- `linear_ode(initial=None)`：同一问题 $u'=-Au$（$A=-G$）的 `LinearODE` 视图，label 为 `"fokker_planck_minus_A"`。
- `sde_state_preparation`：把离散分布编码为幅度 $\sqrt{p_i}$ 的 `StatePreparation`，`"gate"` 为复用旋转实现、`"qram"` 为 QRAM 角表实现（配 `sde_state_angles` 生成绑定数据）。
- 经典见证工具：`matrix_exponential`（缩放平方 + Taylor）、`evolve_distribution`（矩阵指数精确或显式 Euler）、`distribution_moments`、`stationary_distribution` / `boltzmann_distribution`。

## 实现要点

离散化在相邻网格点之间的面上进行，面上系数取两点平均；前向 / 后向速率分别为 $\bigl(\tfrac{a}{2}+\tfrac{D}{h}\bigr)/h$ 与 $\bigl(\tfrac{a}{2}-\tfrac{D}{h}\bigr)/h$，按守恒形式写入 $G$，使每一列的和恰好为零（零通量边界）。离散稳态有精确闭式：相邻点比值为 $(1+u)/(1-u)$，其中 $u=ah/(2D)$ 为面上的有效网格 Péclet 数；$|u|\ge 1$ 时中心差分稳态不再为正，直接抛错。连续稳态参考 $p\propto\exp(\int a/D\,dx)$ 的离散相邻比值为 $e^{2u}$，与离散稳态相差 $O(h^2)$。

耗散声明与初值范数经 `QODEProblem` 的 `evidence` / `dissipative` 字段进入求解协议（见 [QODE 问题对象与协议](qode-problem.md)）；QRAM 态制备只接受非负实幅度，带符号或复相位需另行实现。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）：与经典参考（解析矩 / 闭式稳态）在容差内比较。三层证据位于 `tests/core/test_sde.py`：

- 结构：`GeneratorDiscretizationTests` / `StatePreparationTests` / `SolverContractTests` 的构造断言；`SolverContractTests.test_invalid_inputs_rejected` 与 `StatePreparationTests.test_invalid_probabilities_rejected` 覆盖负扩散、非均匀网格、点数非二的幂、非法初态宽度与非法概率等违例。
- 数值：`GeneratorDiscretizationTests.test_ou_moments_match_closed_form`——16 点 OU 过程（$\theta=1$、$D=0.5$、$t=0.4$）演化后均值对闭式 $\mu_0 e^{-\theta t}$ 容差 2e-3、方差对 $\sigma_0^2 e^{-2\theta t}+\frac{D}{\theta}(1-e^{-2\theta t})$ 容差 2e-2；`test_stationary_approaches_boltzmann`——均匀初态弛豫 30 个时间单位后与离散稳态逐点差 1e-6，离散稳态与 Boltzmann 参考逐点差 2e-2；`test_column_sums_vanish`（列和为零，places=12）。辅助见证：`test_evolution_conserves_probability`（精确与 Euler 演化保概率 places=9、逐点 delta 1e-4）与 `test_matrix_exponential_matches_euler_limit`（Euler 一阶收敛率 $e_{\text{fine}}<0.6\,e_{\text{coarse}}$）。
- 绑定：`SolverContractTests.test_qode_problem_accepted_by_lchs`——4 点 OU 网格的 `QODEProblem` 以 Boltzmann 分布为初态，经 `linear_qode("lchs")`（Cauchy cutoff=0 计划、Taylor 1 阶）的 `check().require()` 与 `solve`，输出宽度 2、`qode_dissipative_promise` 透传正确、程序序列化往返；`StatePreparationTests.test_gate_preparation_amplitudes_are_root_probabilities` 与 `test_qram_preparation_declares_angle_table` 覆盖 gate / QRAM 两类态制备实现。

## 已知缺口与计划阶段

SDE 与 LCHS 的端到端 catalog 案例缺失（离散化 → QODEProblem → 求解 → 读出的目录级登记），归入阶段 V3 的 catalog 登记计划；与 `validation-coverage.md` 的 sde.py 行一致。

## 相关链接

- 源码：`src/pyqecclang/algorithms/qode/sde.py`
- API 参考：[SDE/Fokker–Planck 输入模型](../../api/algorithms/qode/sde.rst)
- 相关页：[QODE 问题对象与协议](qode-problem.md) · [LCHS](lchs.md)（本输入模型的消费协议）
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_nt_qlss_sde.py`（nt_qlss_sde 组）。经典预言机全部独立：`scipy.linalg.expm` 传播子、OU 闭式矩、固定种子 Euler–Maruyama Monte Carlo（numpy，seed=20260916）、$G$ 的零空间向量（numpy 特征分解）；量子部分在真实后端（reference、rir-pysparq、adapter-pysparq、OriginIR-ext）上验证。

**实验设计**：(a) 输入模型——16 点 OU 网格（$\theta=1$、$D=0.5$）生成元的列和、演化保概率、谱横坐标，库矩阵指数对照 scipy.expm；(b) 显式 Euler（`evolve_distribution(steps=k)`）对照 scipy 精确传播子的时间收敛阶；(c) OU 矩三方对拍——Fokker–Planck 网格矩 vs 闭式 $\mu_0 e^{-\theta t}$ / $\sigma_0^2 e^{-2\theta t}+\frac{D}{\theta}(1-e^{-2\theta t})$ vs $2^{18}$ 粒子 Monte Carlo（$dt=10^{-3}$）；(d) 稳态——离散稳态对照零空间向量、长时间松弛收敛、与 Boltzmann 连续参考的 $O(h^2)$ 差距；(e) 量子通道——`sde_state_preparation` 的 gate 实现振幅恰为 $\sqrt{p_i}$（四路径全振幅），QRAM 角表实现对照独立量化角旋转树（实现误差与方法误差分开报告），4 点 OU 生成元的 Pauli 块编码在基态驱动下验证 $(\text{signal}=0)$ 块振幅乘 $\alpha$ 恰为 $G$（四路径）。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| sde-generator-conservation | 16 点 | scipy oracle | 列和 / 保概率 / 矩阵指数偏差 | 0.0 / 4.4e-16 / 6.9e-16 |
| sde-euler-convergence-order | $k$=2000→8000 | scipy oracle | 误差 / 实测阶 | 3.12e-6→7.80e-7 / 1.0001、1.0001 |
| sde-ou-moments-vs-monte-carlo | $t=0.4$ | numpy MC | 均值：FPE vs 闭式 / MC vs 闭式 | 7.0e-9 / 4.2e-4 |
| 同上 | | | 方差：FPE vs 闭式 / MC vs 闭式 | 1.94e-2（网格截断）/ 1.1e-3（统计涨落） |
| sde-stationary-boltzmann | 16 点、$T=30$ | numpy/scipy | TVD：稳态 vs 零空间 / 弛豫 / Boltzmann | 2.8e-15 / 4.2e-14 / 1.2e-3 |
| sde-state-preparation-gate | 2 量子位 | 四路径 | 振幅与 $\sqrt{p}$ 最大偏差 | 5.6e-17 |
| sde-state-preparation-qram | angle_width 8 | 双路径 | 实现误差 / 量化方法误差 | 2.2e-16 / 2.9e-3 |
| sde-generator-encoding | 4 点、$\alpha=4$ | 四路径 | 块振幅 × $\alpha$ vs $G$ | 3.3e-16 |

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

产物：`out/verification/nt_qlss_sde.json`（`sde-*` 共 7 个案例）。

## 数值验证（ode 组）

本节为 ode 组对 `QODEProblem` 消费链路的补充验证（与上一节 nt_qlss_sde 组的输入模型验证并列），见 `tests/verification/verify_ode.py`。经典参考全部独立：`scipy.linalg.expm` 与库内纯 Python `sde.matrix_exponential` 互证、numpy 逐分支 LCHS 仿真；量子程序在真实后端 reference 与 OriginIR-ext 上执行。

**实验设计**：4 点 OU 网格（$\theta=1$、$D=0.5$、网格 $(-1.5,-0.5,0.5,1.5)$）零通量离散生成元经 `generator_encoding()`（显式 Pauli 展开 BE）与 `qode_problem(initial)`（Boltzmann 根振幅初态、`dissipative=True`）组装为 `QODEProblem`，由 `linear_qode("lchs")`（Cauchy cutoff=2、Taylor degree 2）的 `solve` 求解 $t=0.3$；后选择块对照 $L+iH$ 分解的独立仿真与精确传播子。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| lchs-fokker-planck-ou | 16 量子位 | reference、originir | 实现误差（对独立仿真） | 1.1e-16 |
| 同上 | | | 方法误差（求积+Taylor 余项，信息性） | 7.9e-2 |
| 同上 | | | 两个经典参考互证（scipy vs sde 矩阵指数） | 1.1e-16 |
| 同上 | | | `qode_dissipative_promise` 透传 | True |

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

产物：`out/verification/ode.json`（`lchs-fokker-planck-ou` 案例）。
