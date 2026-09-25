# QODE 问题对象与协议（QODE Problem and Protocol）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.qode.ode`](../../api/algorithms/qode/ode.rst) · 阶段 V2

## 概述

本模块不实现具体求解算法，而是为自治齐次线性 ODE

$$
u'(t) = G\,u(t),\qquad u(0)=u_0
$$

定义三个可组合构件：问题对象 {obj}`QODEProblem <oracq.algorithms.qode.ode.QODEProblem>`（生成元、初态、耗散声明与证据说明）、可替换协议 {obj}`QODEProtocol <oracq.algorithms.qode.ode.QODEProtocol>`（契约检查与声明透传），以及通用入口 {obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>`，把 $G$ 路由到 LCHS、CBMD 或 Schrödingerization 等具体方法。另提供 {obj}`make_euler_history_qode <oracq.algorithms.qode.ode.make_euler_history_qode>`，把隐式 Euler 时间推进写成一次量子线性系统求解。

设计立场是"声明优先"：`dissipative` 是数学声明 $\mathrm{Hermitian}(G)\le 0$，不是语言证明；`evidence` 字段强制记录声明来源。用户向的完整流程（输入范式适配、分批绑定）见[微分方程章](../differential-equations.md)。

## 接口与输入模型

```python
QODEProblem(generator, initial, dissipative=None, initial_norm=None,
            evidence="caller_declared_unverified")
linear_qode(method, *, hamiltonian_function=taylor_hamiltonian, **options)
make_euler_history_qode(qlss, *, steps=2)
```

API 入口：{obj}`QODEProblem <oracq.algorithms.qode.ode.QODEProblem>`、{obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>`、{obj}`make_euler_history_qode <oracq.algorithms.qode.ode.make_euler_history_qode>`

- `generator`：$G$ 的 {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`（input model 为 BE），需具备 adjoint/controlled 能力。
- `initial`：{obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>`（SP），零输入、干净 work，宽度与 `generator` 相同。
- `dissipative` / `initial_norm`：`bool | None` 耗散声明与非负初值范数。
- `evidence`：非空字符串，声明来源说明。
- `method`：`"lchs"` / `"cbmd"` / `"schrodingerization"`；`options` 仅接受 `plan`，类型须匹配所选方法。
- `qlss`：可调用的线性系统求解器（如 `qlss.py` 的 Costa 路线），签名 `(matrix, rhs)`。

`linear_qode` 返回 `QODEProtocol`（字段 `name`、`kernel`、`requires_dissipative`，其中 lchs/cbmd 为 `True`）。协议方法：

| 方法 | 行为 |
|---|---|
| `check(generator, initial=None, time=None)` | 契约报告：BE/SP 能力、同宽、时间非负、耗散前提 |
| `__call__(generator, initial, time)` | 兼容入口：check 通过后生成，返回 {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>` |
| `solve(problem, time)` | 问题级入口：校验 `QODEProblem` 后生成（推荐） |

`solve` 返回的 `StateOracle`（`target`/`signal`，宽度等于 `generator.width`）在方法自带属性之外透传问题声明：

| 属性 | 含义 |
|---|---|
| `qode_input_evidence` | 问题对象的 `evidence` |
| `qode_dissipative_promise` | 声明的 `dissipative`（仅当非 None） |
| `qode_initial_norm` | 声明的初值范数（仅当非 None） |

## 实现要点

`QODEProtocol.contract` 在共享的 {obj}`operator_state_contract <oracq.algorithms.input_model.interfaces.operator_state_contract>` 之上按方法附加假设：lchs/cbmd 追加 "Hermitian(G)<=0；问题级 solve 必须声明 dissipative=True"，schrodingerization 追加辅助窗口、Fourier 约定与恢复区域的应用层验证责任。`_generate` 强制输出为 `StateOracle`、宽度与输入一致且保持 adjoint/controlled 能力，违例在生成期抛 {obj}`ContractError <oracq.algorithms.input_model.contracts.ContractError>`。lchs/cbmd 路线内部先构造 {obj}`LinearODE(HermitianParts.from_operator(scale(-1, generator)), initial) <oracq.algorithms.qode.ode_models.LinearODE>`，即 $A=-G$ 的 $L+iH$ 分解；schrodingerization 直接消费 $G$。Hermitian/半正定性质均由调用方承担。

`make_euler_history_qode` 的历史态布局：时间寄存器 `nt = steps.bit_length()` 位在高位，物理寄存器在低位。系统矩阵为 LCU 符号和

$$
C = I - \Delta t\,(q\otimes G) - (S\otimes I),\qquad \Delta t = T/\text{steps},
$$

其中 `q = projector(nt, range(1, steps + 1))` 把 $G$ 限制在时间槽 $1..\text{steps}$，`S = truncated_shift(nt, steps)` 是槽间提升移位；展开后每个槽的方程恰为隐式 Euler $(I-\Delta t\,G)\,u_k = u_{k-1}$，槽 0 由 {obj}`extend_initial <oracq.algorithms.common.state_preparation.extend_initial>` 填入初态。`qlss(c, rhs)` 求解后经 {obj}`select_subspace <oracq.algorithms.common.state_preparation.select_subspace>` 选出时间槽 `steps` 的快照，选择条件并入 signal。

未实现部分：非齐次 $u'=Gu+f$、时间依赖 $G(t)$ 与统一的物理范数恢复对象都没有现成接口；`solve` 只保存声明来源，不添加虚构的范数恢复能力。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。本页覆盖组装与契约层，具体方法的数值见证见各方法页。三层证据：

- 结构：`tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles`——三种 `linear_qode` 方法在抽象 BE/SP 输入下生成后开放槽（`input_A`/`input_b`）原样保留、`algorithm` 属性正确、RIR 序列化往返。
- 数值：`tests/core/test_sde.py:SolverContractTests.test_qode_problem_accepted_by_lchs`——4 点 OU 网格的 `QODEProblem` 经 `check().require()` 与 `solve`，输出宽度、`qode_dissipative_promise` 透传与程序序列化往返均正确。
- 绑定：同类注册与契约检查；`make_euler_history_qode` 的组装经 catalog 的 QHAM m=1 案例（`applications/catalog.py`，`steps=1` 接 Costa QLSS）进入 `tests/core/test_workloads.py:WorkloadTests.test_every_catalog_case_has_a_closed_description` 的目录级往返、闭合与导出检查。

## 已知缺口与计划阶段

统一的"解析可解 ODE 族"收敛基准缺失（同一解析可解问题过全部求解器的批量对拍未建立），归入阶段 V2 的收敛性扫描框架；与 `validation-coverage.md` 的 ode.py 行一致。

## 数值验证

论文级数值实验见 `tests/verification/verify_ode.py`（ode 组，本文件覆盖 `ode.py` 与 `ode_models.py` 的部分）。`QODEProtocol`/`linear_qode` 路由与 `QODEProblem.solve` 的数值行为通过各方法端到端案例验证（详见 [LCHS](lchs.md)、[CBMD](cbmd.md)、[Schrödingerization](schrodingerization.md)、[Carleman 线性化](carleman.md) 的数值验证节）。

**实验设计**：(a) `linear_qode("lchs")` 路由——五种输入模型的端到端对拍（`HermitianParts.from_operator(scale(-1, G))` 的 $A=-G$ 对齐由独立 numpy 仿真确认）；(b) `QODEProblem.solve`——4 点 OU 零通量 Fokker–Planck 问题（`dissipative=True`、evidence 透传）经 `check().require()` 与 `solve`，输出宽度、`qode_dissipative_promise` 透传与数值解同时验证；(c) 声明透传——`qode_input_evidence`/`qode_dissipative_promise` 属性按声明出现在输出模块属性中。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| lchs-given-be-scalar-decay | 7 量子位 | 四路径 | 实现误差（路由+组装对独立仿真） | 1.2e-16 |
| lchs-fokker-planck-ou | 16 量子位 | reference、originir | 实现误差 / 经典参考互证 / `qode_dissipative_promise` | 1.1e-16 / 1.1e-16 / True 透传 |
| cbmd-parts-noncommuting | 15 量子位 | reference、originir | 实现误差（CBMD 路由） | 1.1e-16 |
| schrodingerization-* | 21 量子位 | reference、originir（decay-grid）/ reference（其余） | 组装保真度（Schrödingerization 路由） | 8.7e-19 |

`make_euler_history_qode` 的隐式 Euler 历史组装（LCU 符号和 $C=I-\Delta t(q\otimes G)-(S\otimes I)$）由 catalog 的 QHAM m=1 案例结构测试覆盖，本组未对其做数值端到端验证（需要完整 QLSS 内核，超出本组范围）。

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

产物：`out/verification/ode.json`。

## 相关链接

- 源码：`src/oracq/algorithms/qode/ode.py`
- 教程：[为同一个线性问题替换 QODE 方法](../../tutorials/differential-equations.md)
- API 参考：[QODE 组装接口](../../api/algorithms/qode/ode.rst)
- 同组算法页：[LCHS](lchs.md) · [CBMD](cbmd.md) · [Schrödingerization](schrodingerization.md) · [Carleman 线性化](carleman.md) · [Fokker–Planck 输入模型](fokker-planck.md)
- 概念：[算法自己的约定：从一个 gate 开始](../contracts.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
