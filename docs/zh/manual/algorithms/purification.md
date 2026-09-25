# 纯化访问（Purification Access）

<a href="../../../en/index.html">English</a> · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.input_model.density`](../../api/algorithms/input_model/density.rst) · 阶段 V1

## 概述

密度矩阵（DM）input model 的核心抽象是**纯化访问**：对密度矩阵 $\rho$ 的访问定义为制备其纯化态的量子操作 $U$，$U|0\rangle = |\psi_\rho\rangle$ 作用在 system 与 environment 两个寄存器上，且对环境取偏迹后

$$
\operatorname{Tr}_{\mathrm{env}} |\psi_\rho\rangle\langle\psi_\rho| = \rho .
$$

任何密度矩阵都存在 environment 与 system 等宽的纯化，因此该视图是下游算法（B2 量子 SDP、[Gibbs 态制备](gibbs-state.md)）依赖的稳定接口。模块按仓库的三层范式给出三种构造：abstract 开放声明、gate 显式小矩阵见证、纯态适配。

## 接口与输入模型

```python
abstract_purification(name, width, environment_width=None, *, reversible=True)
gate_purification(rho, *, name=None)
maximally_mixed_purification(width, *, name=None)
PurificationAccess.from_state_preparation(preparation)
```

API 入口：{obj}`abstract_purification <oracq.algorithms.input_model.density.abstract_purification>`、{obj}`gate_purification <oracq.algorithms.input_model.density.gate_purification>`、{obj}`maximally_mixed_purification <oracq.algorithms.input_model.density.maximally_mixed_purification>`

- {obj}`abstract_purification <oracq.algorithms.input_model.density.abstract_purification>`：开放声明一个纯化访问槽（input model 为 DM 的 abstract 层），`environment_width` 缺省取 `width`；经 `linking.bind` 绑定见证实现后程序闭合。
- {obj}`gate_purification <oracq.algorithms.input_model.density.gate_purification>`：显式小密度矩阵的 gate 见证。`rho` 必须是维度为 2..16 内二的幂、迹为 1 的 Hermitian 半正定矩阵。
- {obj}`maximally_mixed_purification <oracq.algorithms.input_model.density.maximally_mixed_purification>`：最大混合态 $I/2^n$ 的纯化生成器（$n$ 对 Bell 对，宽度 ≤ 32）。
- `from_state_preparation`：纯态即平凡纯化，把 work 复净的 {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` 适配为纯化访问（work 扮演 environment）。

四个入口均返回 {obj}`PurificationAccess <oracq.algorithms.input_model.density.PurificationAccess>`（`oracle_kind = "purification_access"`，操作签名为 `("system", "environment")`）。模块属性与属性：

| 属性 | 含义 |
|---|---|
| `width` / `environment_width` | system / environment 寄存器宽度 |
| `operation` | 底层的 {obj}`Operation <oracq.infrastructure.builder.Operation>`，模块属性含 `density_model="purification_access"`、`zero_input=True` |
| {obj}`as_state_preparation() <oracq.algorithms.input_model.interfaces.as_state_preparation>` | 整体视为 system⊕environment 上的 `StatePreparation`（供 B2 组合） |
| `describe()` | oracle 描述（`main_qubit=width`、`anc_qubit=environment_width`） |

配套经典小矩阵工具：{obj}`partial_trace <oracq.algorithms.input_model.density.partial_trace>`（偏迹）、{obj}`trace_distance <oracq.algorithms.input_model.density.trace_distance>`（迹距离 $T(\rho,\sigma) = \lVert\rho-\sigma\rVert_1/2$）、{obj}`gibbs_state <oracq.algorithms.input_model.density.gibbs_state>`（经典参考 Gibbs 态）。

## 实现要点

gate 见证先对 $\rho$ 做循环 Jacobi 特征分解 $\rho = \sum_j p_j |v_j\rangle\langle v_j|$（复 Hermitian 矩阵，先对角相位旋转再实 Jacobi 消元，特征值降序），再取纯化态 $|\psi_\rho\rangle = \sum_j \sqrt{p_j}\,|v_j\rangle_{\mathrm{s}} |j\rangle_{\mathrm{e}}$，其幅度向量经 {obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>` 的多重旋转树在拼接寄存器上制备；基态下标约定为 system | (environment << system_width)。半正定性在生成期校验，负特征值低于 $-10^{-7}$ 时抛 `INPUT_PROMISE` 失败。

`maximally_mixed_purification` 逐位施加 H 与 CNOT 生成 $n$ 对 $|\Phi^+\rangle$；`from_state_preparation` 依赖制备契约的 `clean_work=True` 承诺，未承诺时拒绝适配。适用边界：gate 见证仅支持不超过 16 维的小密度矩阵（验证用小实例），大规模 DM 输入需自行提供 `abstract_purification` 声明的实现并绑定。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）：偏迹后的约化密度矩阵须与目标 $\rho$ 在容差内一致。三层证据：

- 结构：`tests/core/test_density.py:PurificationTests` 各用例的宽度与属性断言（如 `test_gate_purification_recovers_rho` 校验 `access.width == 1`）；`GibbsTests.test_invalid_inputs_fail_at_generation` 覆盖 `gate_purification` 的迹不为 1、非 Hermitian、维度非法，以及 `maximally_mixed_purification(0)`、`partial_trace` 长度不符、`trace_distance` 维度不一致等违例。
- 数值：`PurificationTests.test_gate_purification_recovers_rho` 模拟后对 environment 取偏迹与 2×2 混合态 `RHO` 逐元对拍（places = 7）；`test_maximally_mixed_purification` 对拍 $I/4$ 并校验 Schmidt 结构（非零幅度仅出现在 system == environment 基态）；`test_pure_state_adapter` 对拍纯态 $\sqrt{0.3}|0\rangle + \sqrt{0.7}|1\rangle$ 的约化矩阵；`test_classical_tools` 校验 `gibbs_state` 与 `trace_distance` 的解析值（places = 12）。
- 绑定：`PurificationTests.test_abstract_purification_binds`——abstract 声明经 {obj}`bind <oracq.infrastructure.linking.bind>` 绑定 `gate_purification` 见证后偏迹仍回到 `RHO`，覆盖 abstract / gate 两层一致性。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（偏迹对拍 + Schmidt 结构 + 绑定一致性）。

## 相关链接

- 源码：`src/oracq/algorithms/input_model/density.py`
- 同模块算法：[Gibbs 态制备](gibbs-state.md)
- API 参考：[密度矩阵输入模型与 Gibbs 态](../../api/algorithms/input_model/density.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

本节结果由 `tests/verification/verify_stateprep.py` 在真实后端上跑出（reference、rir-pysparq、adapter-pysparq、originir-ext + UniQC 四条独立路径）。

**实验设计**（5 个案例）：对每条路径制备出的纯化态，用 numpy 独立实现的偏迹（不调用本模块的 `partial_trace` / `trace_distance`）计算 system 约化密度矩阵，再以差矩阵特征值给出迹距离 $T(\rho_{\text{sim}}, \rho)$ 与逐元误差。实例：2×2 复 Hermitian 混合态（非对角元含虚部）、4×4 与 16×16 高斯随机密度矩阵（$G^\dagger G$ 归一化，确定性种子）、最大混合态 $I/8$ 的 3 对 Bell 纯化（附带 Schmidt 结构检查：非零幅度仅出现在 system == environment 基态）、以及 `from_state_preparation` 纯态适配（约化矩阵应等于 $|\psi\rangle\langle\psi|$）。

**关键指标**：

| 案例 | 规模（system+environment） | 路径 | trace_distance | max_element_error |
|---|---|---|---|---|
| purification-gate-rho2 | 1+1 qubit | 四路径 | 1.1e-16 | 1.1e-16 |
| purification-gate-rho4 | 2+2 qubit | 四路径 | 3.2e-11 | 1.9e-11 |
| purification-gate-rho16 | 4+4 qubit | 四路径 | 8.9e-11 | 1.4e-11 |
| purification-bell-w3 | 3+3 qubit | 四路径 | 2.2e-16 | 5.6e-17 |
| purification-pure-adapter | 1+0 qubit | 四路径 | 8.3e-17 | 1.1e-16 |

4×4 与 16×16 案例的残差（~1e-10）来自生成期循环 Jacobi 特征分解（容差 1e-13）的方法误差而非执行误差——后端之间两两偏差处于机器精度量级；Bell 纯化的 `off_schmidt_probability` 为 0，确认 Schmidt 结构。

**复现命令**：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_stateprep.py
```

**产物路径**：`out/verification/stateprep.json`。
