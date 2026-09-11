# LCHS 线性组合哈密顿模拟（Linear Combination of Hamiltonian Simulations）

> 类别 C2 · 模块 `pyqecclang.algorithms.lchs` · 阶段 V2

## 概述

对自治齐次方程 $u'=-Au$、$A=L+iH$ 且 $L\succeq 0$（不要求 $[L,H]=0$），把非酉演化写成 Hermitian 演化的加权积分（LCHS 原文 Theorem 1，[arXiv:2303.01029](https://arxiv.org/html/2303.01029v2)，另见[微分方程章 §4](../differential-equations.md)）：

$$
e^{-At}=\int_{\mathbb R}\frac{e^{-it(H+kL)}}{\pi(1+k^2)}\,dk .
$$

本实现用有限节点 $\{k_j\}$ 与已含积分核的复权重 $\{w_j\}$ 离散该积分，得到有限和 $\widetilde V=\sum_j w_j\,e^{-it(H+k_jL)}$，再以 LCU 组合各分支的哈密顿模拟。与通用接口 $u'=Gu$ 的对齐方式是 $A=-G$（由 `linear_qode` 的路由完成，见 [QODE 问题对象与协议](qode-problem.md)）。

## 接口与输入模型

```python
lchs_qode(model, time, *, plan=None, hamiltonian_function=taylor_hamiltonian)
QuadraturePlan(nodes, weights, kernel="user_supplied")
QuadraturePlan.cauchy(cutoff=2, spacing=1.0)
```

- `model`：`LinearODE(HermitianParts(L, H), initial)`，input model 为 ODE；`parts.hermitian` 存 $L$、`parts.h` 存 $H$，初态为 SP。
- `plan`：节点为有限实数、权重为有限复数且不全为零、两者等长；`cauchy` 工厂生成节点 $k=-J\cdot h,\dots,J\cdot h$ 与权重 $w=h/[\pi(1+k^2)]$，核标记 `"finite_cauchy"`。
- `hamiltonian_function`：`(K, t) -> BlockEncoding` 的可替换协议，默认 `taylor_hamiltonian`（截断 Taylor 的非酉 BE）。
- `time`：非负实数。

返回 `StateOracle`（`target`/`signal`），模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"lchs_qode"` |
| `evolution_alpha` | 有限和 LCU 的归一化常数 $\sum_j|w_j|\,\alpha_{E_j}$ |
| `branch_count` | 求积节点数 |
| `quadrature_kernel` / `quadrature_nodes` | 计划核名与 JSON 节点表 |
| `input_assumption` | `"L>=0; autonomous; homogeneous"` |
| `remainder` | `"finite quadrature pending"` |

## 实现要点

逐节点构造 $K_j = H + k_j L$（`lcu([(1, parts.h), (k_j, parts.hermitian)])`，零系数项由 `lcu` 自动剔除，$\alpha_{K_j}=\alpha_H+|k_j|\,\alpha_L$），交给 `hamiltonian_function` 得到近似演化 $E_j$ 的 BE（必须保留 alpha），再把 $\{(w_j, E_j)\}$ 用 `lcu` 组合成 $\widetilde V$ 并经 `apply_be_to_state` 作用到初态。寄存器布局继承输入 oracle：target 宽度即 $L/H$ 宽度，signal 为 LCU 选择位与各分支信号位的拼接。成功子空间为 signal 全零，读出的是归一化方向 $\widetilde V|u_0\rangle/\alpha_V$；物理幅值还含初值范数，协议不自动恢复。

适用边界：$L\succeq 0$、自治、齐次均为输入声明；`remainder` 明示有限求和没有尾积分保证，不要把权重强行归一化后仍称为同一个近似算子。原文的时间依赖与非齐次结果不在当前接口内。求解前由 `operator_state_contract` 检查 $L$、$H$ 两个 BE 与初态的能力契约。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）：作用算子需在容差内逼近目标连续函数。`tests/core/test_differential.py` 定位为结构与开放输入差异测试（见其 docstring），数值对拍目前由输入模型侧承担。三层证据：

- 结构：`tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` 的 lchs 子测试——抽象 BE/SP 输入下开放槽（`input_A`/`input_b`）保留、`algorithm` 属性正确、RIR 序列化往返。
- 数值：`tests/core/test_sde.py:SolverContractTests.test_qode_problem_accepted_by_lchs`——4 点 OU 网格的 `QODEProblem`（`QuadraturePlan.cauchy(cutoff=0)`、Taylor 1 阶）经 `linear_qode("lchs")` 的 `check`/`solve` 全链路，输出宽度与耗散声明透传正确、序列化往返。
- 绑定：同类注册与契约检查；`FokkerPlanckProblem.qode_problem` 即以本协议为消费方（见 [Fokker–Planck 输入模型](fokker-planck.md)）。

## 已知缺口与计划阶段

统一的"解析可解 ODE 族"收敛基准缺失：同一解析可解问题过全部 ODE 求解器的批量对拍未建立，归入阶段 V2 的收敛性扫描框架；与 `validation-coverage.md` 的 lchs.py 行一致。

## 相关链接

- 源码：`src/pyqecclang/algorithms/lchs.py`
- API 参考：[LCHS](../../api/algorithms/lchs.rst)
- 相关页：[QODE 问题对象与协议](qode-problem.md) · [Carleman 线性化](carleman.md)（提升后接 LCHS）· [Fokker–Planck 输入模型](fokker-planck.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
