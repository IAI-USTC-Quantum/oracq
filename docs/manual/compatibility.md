# 导入路径迁移

0.8 将实现移到按职责划分的目录。根包中的常用导出继续可用；旧模块路径通过集中兼容表转发，引用同一个实现。新代码应使用规范路径。

| 旧路径 | 规范路径 |
|---|---|
| `pyqecclang.ir` | `pyqecclang.infrastructure.ir` |
| `pyqecclang.builder` | `pyqecclang.infrastructure.builder` |
| `pyqecclang.validation` | `pyqecclang.infrastructure.validation` |
| `pyqecclang.serialization` | `pyqecclang.infrastructure.serialization` |
| `pyqecclang.linking` | `pyqecclang.infrastructure.linking` |
| `pyqecclang.execution` | `pyqecclang.infrastructure.execution` |
| `pyqecclang.native` | `pyqecclang.infrastructure.native` |
| `pyqecclang.layout` | `pyqecclang.infrastructure.layout` |
| `pyqecclang.readout` | `pyqecclang.infrastructure.readout` |
| `pyqecclang.backends.originir` | `pyqecclang.infrastructure.backends.originir` |
| `pyqecclang.backends` | `pyqecclang.infrastructure.backends` |
| `pyqecclang.backends.pysparq` | `pyqecclang.infrastructure.backends.pysparq` |
| `pyqecclang.backends.basis` | `pyqecclang.infrastructure.backends.basis` |
| `pyqecclang.mathfunc.graph` | `pyqecclang.infrastructure.mathfunc.graph` |
| `pyqecclang.mathfunc.frontend` | `pyqecclang.infrastructure.mathfunc.frontend` |
| `pyqecclang.mathfunc.numeric` | `pyqecclang.infrastructure.mathfunc.numeric` |
| `pyqecclang.mathfunc.lowering` | `pyqecclang.infrastructure.mathfunc.lowering` |
| `pyqecclang.mathfunc` | `pyqecclang.infrastructure.mathfunc` |
| `pyqecclang.mathfunc.roe_formulas` | `pyqecclang.applications.roe_formulas` |
| `pyqecclang.library` | `pyqecclang.algorithms.operators` |
| `pyqecclang.combinators` | `pyqecclang.algorithms.block_encoding` |
| `pyqecclang.contracts` | `pyqecclang.algorithms.contracts` |
| `pyqecclang.oracles` | `pyqecclang.algorithms.oracles` |
| `pyqecclang.arithmetic` | `pyqecclang.algorithms.arithmetic` |
| `pyqecclang.sparse_models` | `pyqecclang.algorithms.sparse` |
| `pyqecclang.access` | `pyqecclang.algorithms.sparse` |
| `pyqecclang.qlss` | `pyqecclang.algorithms.qlss` |
| `pyqecclang.algorithms.costa` | `pyqecclang.algorithms.qlss` |
| `pyqecclang.algorithms.cks` | `pyqecclang.algorithms.qlss` |
| `pyqecclang.qode` | `pyqecclang.algorithms.ode` |
| `pyqecclang.qfvm` | `pyqecclang.applications.qfvm` |
| `pyqecclang.qfvm_sparse` | `pyqecclang.applications.qfvm` |
| `pyqecclang.flow_data` | `pyqecclang.applications.flow_data` |
| `pyqecclang.roe` | `pyqecclang.applications.roe` |
| `pyqecclang.applications` | `pyqecclang.applications.legacy` |
| `pyqecclang.workloads` | `pyqecclang.applications.catalog` |
| `pyqecclang.qham.pde` | `pyqecclang.applications.qham.pde` |
| `pyqecclang.qham.linearization` | `pyqecclang.applications.qham.linearization` |
| `pyqecclang.qham.reference` | `pyqecclang.applications.qham.reference` |
| `pyqecclang.qham.quantum` | `pyqecclang.algorithms.qham` |
| `pyqecclang.qham` | `pyqecclang.applications.qham` |
| `pyqecclang.qham.examples` | `pyqecclang.applications.qham.examples` |
| `pyqecclang.qham.report` | `pyqecclang.applications.qham.report` |
| `pyqecclang.qham.stencils` | `pyqecclang.applications.qham.stencils` |

## 拆分的旧入口

`algorithms.elementary` 中的查询、搜索、QFT、QPE 和矩阵变换分别进入 `oracle_algorithms`、`search`、`fourier`、`estimation` 和 `transforms`。

`algorithms.differential` 中的具体方法进入 `lchs`、`schrodingerization`、`cbmd` 和 `carleman`；通用入口在 `ode`，PDE 适配在 `pde`。

`algorithms.solvers` 的状态组合工具进入 `state_preparation`，Euler history 进入 `ode`，Trotter 进入 `hamiltonian`。旧工厂仍通过兼容入口可用。

## 文件格式

RIR 版本仍为 0.3，Schema 与指令语义没有改变。Python 类的模块位置与生成器组合发生变化后，模块哈希名可能改变；不要将哈希名字用作应用协议。
