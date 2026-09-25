# 导入路径迁移

<a href="../../en/manual/compatibility.html">English</a> · **简体中文**

0.8 将实现移到按职责划分的目录。根包中的常用导出继续可用；旧模块路径通过集中兼容表转发，引用同一个实现。新代码应使用规范路径。

| 旧路径 | 规范路径 |
|---|---|
| `oracq.ir` | [`oracq.infrastructure.ir`](../api/infrastructure/ir.rst) |
| `oracq.builder` | [`oracq.infrastructure.builder`](../api/infrastructure/builder.rst) |
| `oracq.validation` | [`oracq.infrastructure.validation`](../api/infrastructure/validation.rst) |
| `oracq.serialization` | [`oracq.infrastructure.serialization`](../api/infrastructure/serialization.rst) |
| `oracq.linking` | [`oracq.infrastructure.linking`](../api/infrastructure/linking.rst) |
| `oracq.execution` | [`oracq.infrastructure.execution`](../api/infrastructure/execution.rst) |
| `oracq.native` | [`oracq.infrastructure.native`](../api/infrastructure/native.rst) |
| `oracq.layout` | [`oracq.infrastructure.layout`](../api/infrastructure/layout.rst) |
| `oracq.readout` | [`oracq.infrastructure.readout`](../api/infrastructure/readout.rst) |
| `oracq.backends.originir` | [`oracq.infrastructure.backends.originir`](../api/infrastructure/backends/originir.rst) |
| `oracq.backends` | `oracq.infrastructure.backends` |
| `oracq.backends.pysparq` | [`oracq.infrastructure.backends.pysparq`](../api/infrastructure/backends/pysparq.rst) |
| `oracq.backends.basis` | [`oracq.infrastructure.backends.basis`](../api/infrastructure/backends/basis.rst) |
| `oracq.mathfunc.graph` | [`oracq.infrastructure.mathfunc.graph`](../api/infrastructure/mathfunc/graph.rst) |
| `oracq.mathfunc.frontend` | [`oracq.infrastructure.mathfunc.frontend`](../api/infrastructure/mathfunc/frontend.rst) |
| `oracq.mathfunc.numeric` | [`oracq.infrastructure.mathfunc.numeric`](../api/infrastructure/mathfunc/numeric.rst) |
| `oracq.mathfunc.lowering` | [`oracq.infrastructure.mathfunc.lowering`](../api/infrastructure/mathfunc/lowering.rst) |
| `oracq.mathfunc` | [`oracq.infrastructure.mathfunc`](../api/infrastructure/mathfunc.rst) |
| `oracq.mathfunc.roe_formulas` | [`oracq.applications.roe_formulas`](../api/applications/roe_formulas.rst) |
| `oracq.library` | [`oracq.algorithms.input_model.operators`](../api/algorithms/input_model/operators.rst) |
| `oracq.combinators` | [`oracq.algorithms.input_model.block_encoding`](../api/algorithms/input_model/block_encoding.rst) |
| `oracq.contracts` | [`oracq.algorithms.input_model.contracts`](../api/algorithms/input_model/contracts.rst) |
| `oracq.oracles` | [`oracq.algorithms.input_model.oracles`](../api/algorithms/input_model/oracles.rst) |
| `oracq.arithmetic` | [`oracq.algorithms.common.arithmetic`](../api/algorithms/common/arithmetic.rst) |
| `oracq.sparse_models` | [`oracq.algorithms.input_model.sparse`](../api/algorithms/input_model/sparse.rst) |
| `oracq.access` | [`oracq.algorithms.input_model.sparse`](../api/algorithms/input_model/sparse.rst) |
| `oracq.qlss` | [`oracq.algorithms.qlss.qlss`](../api/algorithms/qlss/qlss.rst) |
| `oracq.algorithms.costa` | [`oracq.algorithms.qlss.qlss`](../api/algorithms/qlss/qlss.rst) |
| `oracq.algorithms.cks` | [`oracq.algorithms.qlss.qlss`](../api/algorithms/qlss/qlss.rst) |
| `oracq.qode` | [`oracq.algorithms.qode.ode`](../api/algorithms/qode/ode.rst) |
| `oracq.qfvm` | [`oracq.applications.qfvm`](../api/applications/qfvm.rst) |
| `oracq.qfvm_sparse` | [`oracq.applications.qfvm`](../api/applications/qfvm.rst) |
| `oracq.flow_data` | [`oracq.applications.flow_data`](../api/applications/flow_data.rst) |
| `oracq.roe` | [`oracq.applications.roe`](../api/applications/roe.rst) |
| `oracq.applications` | `oracq.applications.legacy` |
| `oracq.workloads` | [`oracq.applications.catalog`](../api/applications/catalog.rst) |
| `oracq.qham.pde` | [`oracq.applications.qham.pde`](../api/applications/qham/pde.rst) |
| `oracq.qham.linearization` | [`oracq.applications.qham.linearization`](../api/applications/qham/linearization.rst) |
| `oracq.qham.reference` | [`oracq.applications.qham.reference`](../api/applications/qham/reference.rst) |
| `oracq.qham.quantum` | [`oracq.algorithms.input_model.qham`](../api/algorithms/input_model/qham.rst) |
| `oracq.qham` | `oracq.applications.qham` |
| `oracq.qham.examples` | [`oracq.applications.qham.examples`](../api/applications/qham/examples.rst) |
| `oracq.qham.report` | [`oracq.applications.qham.report`](../api/applications/qham/report.rst) |
| `oracq.qham.stencils` | [`oracq.applications.qham.stencils`](../api/applications/qham/stencils.rst) |

没有独立 API 页的包路径（如 `oracq.infrastructure.backends`、`oracq.applications.qham`）保持原样，可进入其子模块的 API 页查看。

## 拆分的旧入口

`algorithms.elementary` 中的查询、搜索、QFT、QPE 和矩阵变换分别进入 `oracle_algorithms`、`search`、`fourier`、`estimation` 和 `transforms`。

`algorithms.differential` 中的具体方法进入 `lchs`、`schrodingerization`、`cbmd` 和 `carleman`；通用入口在 `ode`，PDE 适配在 `pde`。

`algorithms.solvers` 的状态组合工具进入 `state_preparation`，Euler history 进入 `ode`，Trotter 进入 `hamiltonian`。旧工厂仍通过兼容入口可用。

## 文件格式

RIR 版本为 0.1，指令集合与执行语义不变；链接器新增的 `binding_captures`
是原有标量属性格式中的来源信息，Schema 与语义检查已同步说明其约束。
Python 类的模块位置与生成器组合发生变化后，模块哈希名可能改变；不要将哈希名字用作应用协议。

## 算法接口命名与资源台账

{obj}`AlgorithmContract <oracq.algorithms.input_model.contracts.AlgorithmContract>`、{obj}`QLSSSolver <oracq.algorithms.qlss.qlss.QLSSSolver>`、{obj}`QODESolver <oracq.algorithms.qode.ode.QODESolver>` 是新的规范名称；
{obj}`ProtocolContract <oracq.algorithms.input_model.contracts.ProtocolContract>`、{obj}`QLSSProtocol <oracq.algorithms.qlss.qlss.QLSSProtocol>`、{obj}`QODEProtocol <oracq.algorithms.qode.ode.QODEProtocol>` 保留为同类型别名。
问题构造器接受提供方协议，构造后的字段仍保存具体角色视图。

资源估计的 `rotations` 从列表改为紧凑只读序列。原来的长度、索引与迭代
读法继续可用；直接修改列表的代码应改为读取 `counts` 或使用报告。
开放分析中 `qubits` 为 None，已知下界放在 `qubits_lower_bound`。
