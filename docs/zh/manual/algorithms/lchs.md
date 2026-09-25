# LCHS 线性组合哈密顿模拟（Linear Combination of Hamiltonian Simulations）

<a href="../../../en/index.html">English</a> · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.qode.lchs`](../../api/algorithms/qode/lchs.rst) · 阶段 V2

## 概述

对自治齐次方程 $u'=-Au$、$A=L+iH$ 且 $L\succeq 0$（不要求 $[L,H]=0$），把非酉演化写成 Hermitian 演化的加权积分（LCHS 原文 Theorem 1，[arXiv:2303.01029](https://arxiv.org/html/2303.01029v2)，另见[微分方程章 §4](../differential-equations.md)）：

$$
e^{-At}=\int_{\mathbb R}\frac{e^{-it(H+kL)}}{\pi(1+k^2)}\,dk .
$$

本实现用有限节点 $\{k_j\}$ 与已含积分核的复权重 $\{w_j\}$ 离散该积分，得到有限和 $\widetilde V=\sum_j w_j\,e^{-it(H+k_jL)}$，再以 LCU 组合各分支的哈密顿模拟。与通用接口 $u'=Gu$ 的对齐方式是 $A=-G$（由 {obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>` 的路由完成，见 [QODE 问题对象与协议](qode-problem.md)）。

## 接口与输入模型

```python
lchs_qode(model, time, *, plan=None, hamiltonian_function=taylor_hamiltonian)
QuadraturePlan(nodes, weights, kernel="user_supplied")
QuadraturePlan.cauchy(cutoff=2, spacing=1.0)
```

API 入口：{obj}`lchs_qode <oracq.algorithms.qode.lchs.lchs_qode>`、{obj}`QuadraturePlan <oracq.algorithms.qode.lchs.QuadraturePlan>`

- `model`：{obj}`LinearODE(HermitianParts(L, H), initial) <oracq.algorithms.qode.ode_models.LinearODE>`，input model 为 ODE；`parts.hermitian` 存 $L$、`parts.h` 存 $H$，初态为 SP。
- `plan`：节点为有限实数、权重为有限复数且不全为零、两者等长；`cauchy` 工厂生成节点 $k=-J\cdot h,\dots,J\cdot h$ 与权重 $w=h/[\pi(1+k^2)]$，核标记 `"finite_cauchy"`。
- `hamiltonian_function`：`(K, t) -> BlockEncoding` 的可替换协议，默认 {obj}`taylor_hamiltonian <oracq.algorithms.common.hamiltonian.taylor_hamiltonian>`（截断 Taylor 的非酉 BE）。
- `time`：非负实数。

返回 {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>`（`target`/`signal`），模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"lchs_qode"` |
| `evolution_alpha` | 有限和 LCU 的归一化常数 $\sum_j|w_j|\,\alpha_{E_j}$ |
| `branch_count` | 求积节点数 |
| `quadrature_kernel` / `quadrature_nodes` | 计划核名与 JSON 节点表 |
| `input_assumption` | `"L>=0; autonomous; homogeneous"` |
| `remainder` | `"finite quadrature pending"` |

## 实现要点

逐节点构造 $K_j = H + k_j L$（{obj}`lcu([(1, parts.h), (k_j, parts.hermitian)]) <oracq.algorithms.input_model.block_encoding.lcu>`，零系数项由 `lcu` 自动剔除，$\alpha_{K_j}=\alpha_H+|k_j|\,\alpha_L$），交给 `hamiltonian_function` 得到近似演化 $E_j$ 的 BE（必须保留 alpha），再把 $\{(w_j, E_j)\}$ 用 `lcu` 组合成 $\widetilde V$ 并经 {obj}`apply_be_to_state <oracq.algorithms.common.state_preparation.apply_be_to_state>` 作用到初态。寄存器布局继承输入 oracle：target 宽度即 $L/H$ 宽度，signal 为 LCU 选择位与各分支信号位的拼接。成功子空间为 signal 全零，读出的是归一化方向 $\widetilde V|u_0\rangle/\alpha_V$；物理幅值还含初值范数，协议不自动恢复。

适用边界：$L\succeq 0$、自治、齐次均为输入声明；`remainder` 明示有限求和没有尾积分保证，不要把权重强行归一化后仍称为同一个近似算子。原文的时间依赖与非齐次结果不在当前接口内。求解前由 {obj}`operator_state_contract <oracq.algorithms.input_model.interfaces.operator_state_contract>` 检查 $L$、$H$ 两个 BE 与初态的能力契约。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）：作用算子需在容差内逼近目标连续函数。`tests/core/test_differential.py` 定位为结构与开放输入差异测试（见其 docstring），数值对拍目前由输入模型侧承担。三层证据：

- 结构：`tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` 的 lchs 子测试——抽象 BE/SP 输入下开放槽（`input_A`/`input_b`）保留、`algorithm` 属性正确、RIR 序列化往返。
- 数值：`tests/core/test_sde.py:SolverContractTests.test_qode_problem_accepted_by_lchs`——4 点 OU 网格的 {obj}`QODEProblem <oracq.algorithms.qode.ode.QODEProblem>`（`QuadraturePlan.cauchy(cutoff=0)`、Taylor 1 阶）经 `linear_qode("lchs")` 的 `check`/`solve` 全链路，输出宽度与耗散声明透传正确、序列化往返。
- 绑定：同类注册与契约检查；`FokkerPlanckProblem.qode_problem` 即以本协议为消费方（见 [Fokker–Planck 输入模型](fokker-planck.md)）。

## 已知缺口与计划阶段

统一的"解析可解 ODE 族"收敛基准缺失：同一解析可解问题过全部 ODE 求解器的批量对拍未建立，归入阶段 V2 的收敛性扫描框架；与 `validation-coverage.md` 的 lchs.py 行一致。

## 数值验证

论文级数值实验见 `tests/verification/verify_ode.py`（ode 组，本文件覆盖 `lchs.py` 与 `_dynamics.py` 组装骨架的部分）。经典参考全部独立：numpy 逐分支 $K_j=H+k_jL$ 截断 Taylor 求和、`scipy.linalg.expm` 精确解、解析解（$e^{-t}$、热方程 Fourier 特征值）、`sde.matrix_exponential` 纯 Python 矩阵指数互证；量子程序在真实后端 reference、rir-pysparq、adapter-pysparq、OriginIR-ext 上执行。实现误差（量子后选择块对独立仿真）与方法误差（有限求积+Taylor 余项对精确解，库中标注 pending）分开报告。

**实验设计**：$u'=-Au$、$A=L+iH$ 五类输入模型的端到端对拍（plan 均 `QuadraturePlan.cauchy(cutoff=2, spacing=1.0)`）：(1) 整体 $G=-I$ 的 BE（解析解）；(2) 直接 Hermitian parts（$[L,H]\neq0$ 的 $2\times2$ 矩阵，Pauli 展开 BE）；(3) 对角谱角数据库 $A=\mathrm{diag}(1,\cos\pi/4)$，同一开放 RIR 分别绑定 gate 表与 QRAM 表；(4) Fokker–Planck OU 4 点零通量离散生成元（Pauli 展开）经 `QODEProblem.solve`（Boltzmann 根振幅初态）；(5) 周期 4 点热方程的结构化移位差分 BE。另有求积收敛扫描（标量 $L=I$，cutoff 2/8 量子 + 2..32 经典核趋势线）。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| lchs-given-be-scalar-decay | 7 量子位 | 四路径 | 实现误差 / 方法误差 | 1.2e-16 / 1.2e-2 |
| lchs-given-parts-noncommuting | 15 量子位 | 三路径 | 实现误差 / 方法误差 | 1.1e-16 / 7.3e-2 |
| lchs-diagonal-spectral-gate / -qram | 2×程序 | reference | 实现误差（两者） / gate-QRAM 振幅差 | 5.8e-17 / 0.0 |
| lchs-fokker-planck-ou | 16 量子位 | reference、originir | 实现误差 / 方法误差 / 经典参考互证 | 1.1e-16 / 7.9e-2 / 1.1e-16 |
| lchs-heat-structured-stencil | 15 量子位 | reference、originir | 实现误差 / 方法误差 / Fourier-expm 互证 | 1.7e-17 / 0.19 / 1.4e-17 |
| lchs-quadrature-scan | cutoff 2 / 8 | reference | 实现误差 | 5.6e-17 / 1.1e-16 |
| lchs-quadrature-kernel-trend | cutoff 2→32 | 经典核求积 | 核求积误差 | 7.1e-2→3.7e-3（含振荡尾部） |

成功概率（信号全零子空间）与经典值一致（偏差 < 1e-16）。方法误差为 Cauchy 求积余项：截断 $K_{\max}$ 与步长 $h$ 的尾部分别呈 $O(1/K_{\max})$ 振荡下降与混叠地板 $2e^{-(2\pi/h-\lambda t)}$；深嵌套 LCU 程序上 pysparq 两实现对 junk 分支有约 1e-7 数值地板（物理块不受影响，见组报告）。

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

产物：`out/verification/ode.json`（`lchs-*` 共 12 个案例）。

## 相关链接

- 源码：`src/oracq/algorithms/qode/lchs.py`
- 教程：[为同一个线性问题替换 QODE 方法](../../tutorials/differential-equations.md)
- API 参考：[LCHS](../../api/algorithms/qode/lchs.rst)
- 相关页：[QODE 问题对象与协议](qode-problem.md) · [Carleman 线性化](carleman.md)（提升后接 LCHS）· [Fokker–Planck 输入模型](fokker-planck.md) · [CBMD](cbmd.md)、[Schrödingerization](schrodingerization.md)（其他线性路线）
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
