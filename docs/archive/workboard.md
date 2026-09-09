# 范式实现工作面板

一般 QHAM 的本轮进度见 [QHAM 工作面板](qham-general-plan.md)。

QFVM / QLSS 输入模型修正见 [审查计划与结果](qfvm-qlss-plan.md)。

新增纯函数自动编译进度见 [函数编译面板](function-compiler-plan.md)。

后续第二阶段进度见 [寄存器级模拟与算法工作面板](stage2-board.md)。

更新时间：2026-09-09。本轮验收表达、开放 IR、绑定和后端描述可达性；数学正确性留待下一阶段。

| 编号 | 工作项 | 状态 | 证据 |
|---|---|---|---|
| P0 | 盘点 QECC.Lang 文档与用例，建立覆盖矩阵 | done | [coverage.json](coverage.json)、[coverage.md](coverage.md) |
| P1 | RIR 0.2：开放 oracle 声明、部分绑定、缺口报告和闭合导出 | done | [linking.py](../api/infrastructure/linking.rst)、[open-ir.md](../reference/open-ir.md)、[test_open_ir.py](../../tests/core/test_open_ir.py) |
| P2 | BE、XOR database、state-prep isometry、CKS sparse 等范式 | done | [oracles.py](../api/algorithms/oracles.rst)、[oracle-paradigms.md](oracle-paradigms.md) |
| P3 | Gate/QRAM 基础实现与组合适配器 | done | [access.py](../api/algorithms/sparse.rst)、[combinators.py](../api/algorithms/block_encoding.rst) |
| P4 | Deutsch–Jozsa、Grover 抽象与具体组装 | done | [elementary.py](../api/index.rst) |
| P5 | Costa QLSS walk、调度、filtering 及实例绑定 | done | [costa.py](../api/algorithms/qlss.rst) |
| P6 | QFVM 输入访问与 QLSS 组装 | done | [applications.py](../api/applications/legacy.rst)、[program.originir](../../out/catalog/qfvm_qram/program.originir) |
| P7 | QODE/QPDE 与 QHAM 提升系统组装 | done | [solvers.py](../api/index.rst)、[program.originir](../../out/catalog/qham_qode/program.originir) |
| P8 | 覆盖面板、开放/闭合 IR 和 OriginIR-ext 示例产物、结构验收 | done | [build_catalog.py](../../tools/build_catalog.py)、[test_workloads.py](../../tests/core/test_workloads.py)、[test_catalog.py](../../tests/integration/test_catalog.py)、[phase-validation.json](phase-validation.json) |

## 案例产物

每个目录提供 open.rir.json、partial.rir.json、closed.rir.json、bindings.json、memory.json 和 program.originir。它们是可重新生成的描述产物，保存在被 Git 忽略的 out/catalog。

| 案例 | 开放槽 | 绑定后模块 | OriginIR DEF | 原生解析 | 产物 |
|---|---:|---:|---:|---|---|
| bell | 0 | 1 | 1 | 通过 | [开放](../../out/catalog/bell/open.rir.json) / [部分绑定](../../out/catalog/bell/partial.rir.json) / [闭合](../../out/catalog/bell/closed.rir.json) / [OriginIR](../../out/catalog/bell/program.originir) / [报告](../../out/catalog/bell/report.json) |
| ghz | 0 | 1 | 1 | 通过 | [开放](../../out/catalog/ghz/open.rir.json) / [部分绑定](../../out/catalog/ghz/partial.rir.json) / [闭合](../../out/catalog/ghz/closed.rir.json) / [OriginIR](../../out/catalog/ghz/program.originir) / [报告](../../out/catalog/ghz/report.json) |
| python_generators | 0 | 3 | 7 | 通过 | [开放](../../out/catalog/python_generators/open.rir.json) / [部分绑定](../../out/catalog/python_generators/partial.rir.json) / [闭合](../../out/catalog/python_generators/closed.rir.json) / [OriginIR](../../out/catalog/python_generators/program.originir) / [报告](../../out/catalog/python_generators/report.json) |
| trotter_hamsim | 0 | 1 | 3 | 通过 | [开放](../../out/catalog/trotter_hamsim/open.rir.json) / [部分绑定](../../out/catalog/trotter_hamsim/partial.rir.json) / [闭合](../../out/catalog/trotter_hamsim/closed.rir.json) / [OriginIR](../../out/catalog/trotter_hamsim/program.originir) / [报告](../../out/catalog/trotter_hamsim/report.json) |
| dj_gate | 1 | 3 | 3 | 通过 | [开放](../../out/catalog/dj_gate/open.rir.json) / [部分绑定](../../out/catalog/dj_gate/partial.rir.json) / [闭合](../../out/catalog/dj_gate/closed.rir.json) / [OriginIR](../../out/catalog/dj_gate/program.originir) / [报告](../../out/catalog/dj_gate/report.json) |
| dj_qram | 1 | 3 | 3 | 通过 | [开放](../../out/catalog/dj_qram/open.rir.json) / [部分绑定](../../out/catalog/dj_qram/partial.rir.json) / [闭合](../../out/catalog/dj_qram/closed.rir.json) / [OriginIR](../../out/catalog/dj_qram/program.originir) / [报告](../../out/catalog/dj_qram/report.json) |
| grover_gate | 1 | 4 | 5 | 通过 | [开放](../../out/catalog/grover_gate/open.rir.json) / [部分绑定](../../out/catalog/grover_gate/partial.rir.json) / [闭合](../../out/catalog/grover_gate/closed.rir.json) / [OriginIR](../../out/catalog/grover_gate/program.originir) / [报告](../../out/catalog/grover_gate/report.json) |
| grover_qram | 1 | 5 | 6 | 通过 | [开放](../../out/catalog/grover_qram/open.rir.json) / [部分绑定](../../out/catalog/grover_qram/partial.rir.json) / [闭合](../../out/catalog/grover_qram/closed.rir.json) / [OriginIR](../../out/catalog/grover_qram/program.originir) / [报告](../../out/catalog/grover_qram/report.json) |
| stateprep_gate | 1 | 2 | 2 | 通过 | [开放](../../out/catalog/stateprep_gate/open.rir.json) / [部分绑定](../../out/catalog/stateprep_gate/partial.rir.json) / [闭合](../../out/catalog/stateprep_gate/closed.rir.json) / [OriginIR](../../out/catalog/stateprep_gate/program.originir) / [报告](../../out/catalog/stateprep_gate/report.json) |
| stateprep_qram | 1 | 2 | 2 | 通过 | [开放](../../out/catalog/stateprep_qram/open.rir.json) / [部分绑定](../../out/catalog/stateprep_qram/partial.rir.json) / [闭合](../../out/catalog/stateprep_qram/closed.rir.json) / [OriginIR](../../out/catalog/stateprep_qram/program.originir) / [报告](../../out/catalog/stateprep_qram/report.json) |
| sparse_gate | 3 | 10 | 10 | 通过 | [开放](../../out/catalog/sparse_gate/open.rir.json) / [部分绑定](../../out/catalog/sparse_gate/partial.rir.json) / [闭合](../../out/catalog/sparse_gate/closed.rir.json) / [OriginIR](../../out/catalog/sparse_gate/program.originir) / [报告](../../out/catalog/sparse_gate/report.json) |
| sparse_qram | 3 | 10 | 10 | 通过 | [开放](../../out/catalog/sparse_qram/open.rir.json) / [部分绑定](../../out/catalog/sparse_qram/partial.rir.json) / [闭合](../../out/catalog/sparse_qram/closed.rir.json) / [OriginIR](../../out/catalog/sparse_qram/program.originir) / [报告](../../out/catalog/sparse_qram/report.json) |
| costa_gate | 2 | 10 | 20 | 通过 | [开放](../../out/catalog/costa_gate/open.rir.json) / [部分绑定](../../out/catalog/costa_gate/partial.rir.json) / [闭合](../../out/catalog/costa_gate/closed.rir.json) / [OriginIR](../../out/catalog/costa_gate/program.originir) / [报告](../../out/catalog/costa_gate/report.json) |
| costa_qram | 2 | 10 | 20 | 通过 | [开放](../../out/catalog/costa_qram/open.rir.json) / [部分绑定](../../out/catalog/costa_qram/partial.rir.json) / [闭合](../../out/catalog/costa_qram/closed.rir.json) / [OriginIR](../../out/catalog/costa_qram/program.originir) / [报告](../../out/catalog/costa_qram/report.json) |
| costa_sparse_qram | 4 | 16 | 33 | 通过 | [开放](../../out/catalog/costa_sparse_qram/open.rir.json) / [部分绑定](../../out/catalog/costa_sparse_qram/partial.rir.json) / [闭合](../../out/catalog/costa_sparse_qram/closed.rir.json) / [OriginIR](../../out/catalog/costa_sparse_qram/program.originir) / [报告](../../out/catalog/costa_sparse_qram/report.json) |
| qfvm_gate | 9 | 24 | 49 | 通过 | [开放](../../out/catalog/qfvm_gate/open.rir.json) / [部分绑定](../../out/catalog/qfvm_gate/partial.rir.json) / [闭合](../../out/catalog/qfvm_gate/closed.rir.json) / [OriginIR](../../out/catalog/qfvm_gate/program.originir) / [报告](../../out/catalog/qfvm_gate/report.json) |
| qfvm_qram | 9 | 22 | 53 | 通过 | [开放](../../out/catalog/qfvm_qram/open.rir.json) / [部分绑定](../../out/catalog/qfvm_qram/partial.rir.json) / [闭合](../../out/catalog/qfvm_qram/closed.rir.json) / [OriginIR](../../out/catalog/qfvm_qram/program.originir) / [报告](../../out/catalog/qfvm_qram/report.json) |
| qham_qode | 3 | 45 | 89 | 通过 | [开放](../../out/catalog/qham_qode/open.rir.json) / [部分绑定](../../out/catalog/qham_qode/partial.rir.json) / [闭合](../../out/catalog/qham_qode/closed.rir.json) / [OriginIR](../../out/catalog/qham_qode/program.originir) / [报告](../../out/catalog/qham_qode/report.json) |
| qham_qpde | 3 | 45 | 89 | 通过 | [开放](../../out/catalog/qham_qpde/open.rir.json) / [部分绑定](../../out/catalog/qham_qpde/partial.rir.json) / [闭合](../../out/catalog/qham_qpde/closed.rir.json) / [OriginIR](../../out/catalog/qham_qpde/program.originir) / [报告](../../out/catalog/qham_qpde/report.json) |
| qpe | 1 | 6 | 8 | 通过 | [开放](../../out/catalog/qpe/open.rir.json) / [部分绑定](../../out/catalog/qpe/partial.rir.json) / [闭合](../../out/catalog/qpe/closed.rir.json) / [OriginIR](../../out/catalog/qpe/program.originir) / [报告](../../out/catalog/qpe/report.json) |
| qsvt | 1 | 4 | 4 | 通过 | [开放](../../out/catalog/qsvt/open.rir.json) / [部分绑定](../../out/catalog/qsvt/partial.rir.json) / [闭合](../../out/catalog/qsvt/closed.rir.json) / [OriginIR](../../out/catalog/qsvt/program.originir) / [报告](../../out/catalog/qsvt/report.json) |
| oaa | 1 | 4 | 5 | 通过 | [开放](../../out/catalog/oaa/open.rir.json) / [部分绑定](../../out/catalog/oaa/partial.rir.json) / [闭合](../../out/catalog/oaa/closed.rir.json) / [OriginIR](../../out/catalog/oaa/program.originir) / [报告](../../out/catalog/oaa/report.json) |
| lchs | 2 | 10 | 10 | 通过 | [开放](../../out/catalog/lchs/open.rir.json) / [部分绑定](../../out/catalog/lchs/partial.rir.json) / [闭合](../../out/catalog/lchs/closed.rir.json) / [OriginIR](../../out/catalog/lchs/program.originir) / [报告](../../out/catalog/lchs/report.json) |
| heat_qode | 2 | 21 | 42 | 通过 | [开放](../../out/catalog/heat_qode/open.rir.json) / [部分绑定](../../out/catalog/heat_qode/partial.rir.json) / [闭合](../../out/catalog/heat_qode/closed.rir.json) / [OriginIR](../../out/catalog/heat_qode/program.originir) / [报告](../../out/catalog/heat_qode/report.json) |
| poisson_qlss | 2 | 12 | 24 | 通过 | [开放](../../out/catalog/poisson_qlss/open.rir.json) / [部分绑定](../../out/catalog/poisson_qlss/partial.rir.json) / [闭合](../../out/catalog/poisson_qlss/closed.rir.json) / [OriginIR](../../out/catalog/poisson_qlss/program.originir) / [报告](../../out/catalog/poisson_qlss/report.json) |
| carleman_step | 2 | 30 | 60 | 通过 | [开放](../../out/catalog/carleman_step/open.rir.json) / [部分绑定](../../out/catalog/carleman_step/partial.rir.json) / [闭合](../../out/catalog/carleman_step/closed.rir.json) / [OriginIR](../../out/catalog/carleman_step/program.originir) / [报告](../../out/catalog/carleman_step/report.json) |
| schrodingerisation | 2 | 9 | 9 | 通过 | [开放](../../out/catalog/schrodingerisation/open.rir.json) / [部分绑定](../../out/catalog/schrodingerisation/partial.rir.json) / [闭合](../../out/catalog/schrodingerisation/closed.rir.json) / [OriginIR](../../out/catalog/schrodingerisation/program.originir) / [报告](../../out/catalog/schrodingerisation/report.json) |
| be_algebra | 2 | 18 | 23 | 通过 | [开放](../../out/catalog/be_algebra/open.rir.json) / [部分绑定](../../out/catalog/be_algebra/partial.rir.json) / [闭合](../../out/catalog/be_algebra/closed.rir.json) / [OriginIR](../../out/catalog/be_algebra/program.originir) / [报告](../../out/catalog/be_algebra/report.json) |
| arithmetic | 1 | 4 | 4 | 通过 | [开放](../../out/catalog/arithmetic/open.rir.json) / [部分绑定](../../out/catalog/arithmetic/partial.rir.json) / [闭合](../../out/catalog/arithmetic/closed.rir.json) / [OriginIR](../../out/catalog/arithmetic/program.originir) / [报告](../../out/catalog/arithmetic/report.json) |
| batch_qram | 1 | 3 | 3 | 通过 | [开放](../../out/catalog/batch_qram/open.rir.json) / [部分绑定](../../out/catalog/batch_qram/partial.rir.json) / [闭合](../../out/catalog/batch_qram/closed.rir.json) / [OriginIR](../../out/catalog/batch_qram/program.originir) / [报告](../../out/catalog/batch_qram/report.json) |
| banked_qram | 1 | 2 | 2 | 通过 | [开放](../../out/catalog/banked_qram/open.rir.json) / [部分绑定](../../out/catalog/banked_qram/partial.rir.json) / [闭合](../../out/catalog/banked_qram/closed.rir.json) / [OriginIR](../../out/catalog/banked_qram/program.originir) / [报告](../../out/catalog/banked_qram/report.json) |
| register_views | 0 | 1 | 1 | 通过 | [开放](../../out/catalog/register_views/open.rir.json) / [部分绑定](../../out/catalog/register_views/partial.rir.json) / [闭合](../../out/catalog/register_views/closed.rir.json) / [OriginIR](../../out/catalog/register_views/program.originir) / [报告](../../out/catalog/register_views/report.json) |
| measure_reset | 0 | 1 | 1 | 通过 | [开放](../../out/catalog/measure_reset/open.rir.json) / [部分绑定](../../out/catalog/measure_reset/partial.rir.json) / [闭合](../../out/catalog/measure_reset/closed.rir.json) / [OriginIR](../../out/catalog/measure_reset/program.originir) / [报告](../../out/catalog/measure_reset/report.json) |

## 覆盖与限制

- [旧案例逐项映射](coverage.md)包含 61 个正例和 16 个负例；这里不是旧源码/CLIR 的等价证明。
- QFVM 具体样例使用玩具物理查表核；完整 Roe 可以继续保留为开放声明。
- 此阶段 QHAM 为 m=1 特例；一般阶数的后续实现见 qham-general-plan.md。
- 当前绑定针对固定宽度和 alpha 的接口；更换这些常量需要重新运行 Python 生成器。
- 所有应用均标记 correctness=not_assessed；没有因为尚未完成精度或收敛证明而阻断组装。

## 下一阶段验证顺序

先检查普通 oracle 的矩阵和位语义，再检查 BE 与 state prep，随后是 Costa 初态/反射/filtering，最后验证 QODE/PDE、Roe 物理核与 QHAM 的条件输出和外层行为。
