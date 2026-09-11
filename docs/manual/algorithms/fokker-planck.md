# Fokker–Planck 输入模型（Fokker–Planck Input Model）

> 类别 C2 · 模块 `pyqecclang.algorithms.sde` · 阶段 V3

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

- 源码：`src/pyqecclang/algorithms/sde.py`
- API 参考：[SDE/Fokker–Planck 输入模型](../../api/algorithms/sde.rst)
- 相关页：[QODE 问题对象与协议](qode-problem.md) · [LCHS](lchs.md)（本输入模型的消费协议）
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
