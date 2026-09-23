# CBMD 轮廓分解（CBMD）

> 类别 C2 · 模块 `pyqecclang.algorithms.qode.cbmd` · 阶段 V2

## 概述

基于轮廓分解的矩阵函数与演化组装：把 $u'=-Au$（$A=L+iH$，$L\succeq 0$）的 $e^{-At}$ 写成 Hermitian 分支 $H+qL$（$q$ 取实节点）的有限加权和，权重来自轮廓恒等式的留数离散。主级数与截断依据源码 `ContourPlan.metadata` 引用的 QST 11 035027 (2026), Eq.12–13（[arXiv:2511.10267v3](https://arxiv.org/abs/2511.10267)）。

实现立场是"显式记录遗漏"：辅助极点贡献的非 Hermitian 演化分支与无穷级数尾项当前**不生成**，两者在 `metadata()` 的 `omitted` 列表中明确保留；成立前提 `L>=0 and norm(integral L dt) <= 2*pi*a` 也一并写入元数据。

## 接口与输入模型

```python
cbmd_qode(model, time, *, plan=None, hamiltonian_function=taylor_hamiltonian)
ContourPlan(a=1.0, cutoff=2, poles=(2j, -1 + 1j, 1j, 1 + 1j))
cbmd_function(a, nodes, residue_weights, hermitian_function)
```

- `model`：`LinearODE(HermitianParts(L, H), initial)`，input model 为 ODE（与 LCHS 相同的输入面）。
- `ContourPlan`：`a > 0` 的轮廓参数、`cutoff >= 0` 的截断、有限互异非实且避开 $-i$ 的辅助极点。派生量：主级数节点 $q_k=k/a$（$k=-J..J$）、节点权重与辅助极点系数。
- `cbmd_qode` 的 `hamiltonian_function` 与 `time` 约定同 LCHS。
- `cbmd_function(a, nodes, residue_weights, hermitian_function)`：通用 $f(A)$ 组装点，`a` 为算子 BE，节点/留数权重由调用方按轮廓恒等式提供。

`cbmd_qode` 返回 `StateOracle`（`target`/`signal`），模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"cbmd_qode"` |
| `evolution_alpha` / `branch_count` | 有限和 LCU 归一化与节点数 |
| `contour_plan` | `ContourPlan.metadata()` 的 JSON（含 `omitted`、`assumption`、`source`） |
| `remainder` | `"auxiliary pole contribution and truncation pending"` |
| `input_assumption` | `"L>=0; autonomous; homogeneous"` |

`cbmd_function` 返回 `BlockEncoding`，属性 `algorithm="cbmd_matrix_function"`、`correctness="pending"`、`residue_sign_convention="caller supplies target-side weights from contour identity"`。

## 实现要点

主级数节点权重按闭式

$$
w_k=\frac{e^{-2\pi a}-1}{a\cdot 2\pi i\,(q_k+i)\prod_{p}\frac{q_k-p}{-i-p}}
$$

给出，辅助极点系数同理（分母换为 $(e^{-2\pi i p a}-1)\prod_{p'\neq p}\frac{p-p'}{-i-p'}$）。`cbmd_qode` 的生成与 LCHS 共用 `_lcu_dynamics`：逐节点 $K_k=H+q_kL$、可替换 `hamiltonian_function`、LCU 组合、`apply_be_to_state` 作用到初态；差异仅在计划类型与记录的元数据。

`cbmd_function` 的分支组合顺序为 `lcu([(node, parts.h), (1, parts.hermitian)])`，即 $qH+L$（与 `cbmd_qode` 的 $H+qL$ 约定不同）；它保持 Hermitian-function protocol 开放，不偷换为矩阵求逆，也不宣称完成了留数符号方向的验证。

适用边界：除 LCHS 的耗散/自治/齐次声明外，还需轮廓前提 $\lVert\int L\,dt\rVert\le 2\pi a$；省略的辅助极点分支与无穷级数尾意味着当前输出只是主级数的有限部分，不是完整的 $e^{-At}$ 近似。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。三层证据：

- 结构：`tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` 的 cbmd 子测试——抽象 BE/SP 输入下开放槽保留、`algorithm` 属性正确、序列化往返。
- 数值：`tests/core/test_differential.py:DifferentialStructureTests.test_cbmd_plan_records_omitted_terms`——`ContourPlan(cutoff=2)` 的权重数为 5、辅助系数与极点数一致，且 `metadata()` 的 `omitted` 列表包含 `auxiliary_nonhermitian_evolutions`。
- 绑定：同类注册与契约检查（共享 `_dynamics` 契约路径，见 [LCHS](lchs.md) 的契约链见证）。

## 已知缺口与计划阶段

统一的"解析可解 ODE 族"收敛基准缺失（阶段 V2）；辅助极点分支与无穷级数尾的实现属于算法层后续工作，当前以 `omitted` 元数据显式声明，不计入缺口列。与 `validation-coverage.md` 的 cbmd.py 行一致。

## 数值验证

论文级数值实验见 `tests/verification/verify_ode.py`（ode 组，本文件覆盖 `cbmd.py` 的部分）。经典参考全部独立：`ContourPlan` 派生量按 Eq.12–13 闭式直接重算、$t=0$ 轮廓恒等式残差、numpy 逐分支截断 Taylor 求和、`scipy.linalg.expm` 精确解；量子程序在真实后端 reference、rir-pysparq、OriginIR-ext 上执行。

**实验设计**：(a) 计划子结构——`ContourPlan(a=1, cutoff=2..8)` 的节点/权重/辅助极点系数对独立闭式重算，$t=0$ 恒等式 main+aux→1 的截断残差趋势；(b) 端到端——与 LCHS 同一非对易 $2\times2$ 问题（$L=[[1,0.3],[0.3,0.5]]$、$H=[[0.2,0.1],[0.1,-0.1]]$），`ContourPlan(a=1.0, cutoff=2)`、Taylor degree 3、$t=0.2$，后选择块对独立仿真与精确解；(c) 同题对拍——复用两方法的量子物理块比较解方向。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| cbmd-contour-plan-formula-and-identity | cutoff 2/4/8 | 计划对象 | 权重/辅助系数对闭式最大偏差 | 0.0 |
| 同上 | | | $t=0$ 恒等式残差（信息性） | 3.73e-2 / 3.20e-3 / 1.59e-4 |
| cbmd-parts-noncommuting | 15 量子位 | reference、originir | 实现误差（对独立仿真） | 1.1e-16 |
| 同上 | | | 方法误差（omitted 项，信息性） | 9.1e-3 |
| cbmd-vs-lchs-direction | 同题 | reference | 解方向误差：CBMD / LCHS | 3.0e-3 / 4.5e-2 |

恒等式残差随截断稳定下降，与元数据 `omitted` 中 `infinite_series_tail` 一致；端到端方法误差即文档声明省略的辅助极点分支贡献（本实例约 1%）。同题下 CBMD 主级数的解方向误差（3.0e-3）显著小于 LCHS 的 Cauchy 求积余项（4.5e-2）。

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

产物：`out/verification/ode.json`（`cbmd-*` 共 3 个案例）。

## 相关链接

- 源码：`src/pyqecclang/algorithms/qode/cbmd.py`
- API 参考：[CBMD](../../api/algorithms/qode/cbmd.rst)
- 相关页：[QODE 问题对象与协议](qode-problem.md) · [LCHS](lchs.md)（同一输入面与组装骨架）· [Carleman 线性化](carleman.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
