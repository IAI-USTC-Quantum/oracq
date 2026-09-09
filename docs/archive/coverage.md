# 旧案例与新范式覆盖矩阵

本表覆盖 61 个正例和 16 个负例。映射证明有对应的表达和组装路径，不表示旧源码逐字迁移、旧 golden 等价或算法正确性认证。

| 旧用例 | 状态 | 新例子/证据 | 说明 |
|---|---|---|---|
| 00-primitives/adjoint-block.qec | paradigm_mapped | oaa | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 00-primitives/bell-pair.qec | paradigm_mapped | bell | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 00-primitives/control-nested.qec | paradigm_mapped | costa_gate | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 00-primitives/dagger-self-inverse.qec | paradigm_mapped | oaa | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 00-primitives/ghz.qec | paradigm_mapped | ghz | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 00-primitives/measure-reset.qec | paradigm_mapped | measure_reset | 由显式宿主 ReadoutAction 生成末端测量和重置；不混入酉 oracle 主体。 |
| 00-primitives/rot-angle.qec | paradigm_mapped | qsvt | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 00-primitives/slice-ops.qec | paradigm_mapped | register_views | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 01-registers/ancilla-clean.qec | paradigm_mapped | costa_gate | 工作区显式进入接口；不继承旧版 clean_ancilla 静态证明。 |
| 01-registers/fixed-point-add.qec | paradigm_mapped | arithmetic | 可逆算术接口与小型查表实现；RNE/溢出算法语义待下一阶段。 |
| 01-registers/fuse-basic.qec | paradigm_mapped | register_views | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 01-registers/qram-paged.qec | paradigm_mapped | batch_qram, banked_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 01-registers/register-array.qec | paradigm_mapped | batch_qram, banked_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 01-registers/reinterpret-fixed.qec | paradigm_mapped | arithmetic | 定点解释作为存储位模式和库元数据，不内建数值精度证明。 |
| 01-registers/split-basic.qec | paradigm_mapped | register_views | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 01-registers/split-in-fuse.qec | paradigm_mapped | register_views | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 02-generics-const/const-fn-iters.qec | paradigm_mapped | python_generators, be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 02-generics-const/const-recursion.qec | paradigm_mapped | python_generators, be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 02-generics-const/external-fn-phases.qec | paradigm_mapped | python_generators, be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 02-generics-const/generic-explicit.qec | paradigm_mapped | python_generators, be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 02-generics-const/generic-float-alpha.qec | paradigm_mapped | python_generators, be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 02-generics-const/generic-infer-width.qec | paradigm_mapped | python_generators, be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 02-generics-const/op-recursion.qec | paradigm_mapped | python_generators, be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 03-oracle/grover2.qec | paradigm_mapped | grover_gate, grover_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 03-oracle/phase-oracle-custom.qec | paradigm_mapped | grover_gate, grover_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 03-oracle/qram-lookup.qec | paradigm_mapped | dj_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 03-oracle/requires-controlled-ok.qec | paradigm_mapped | grover_gate, grover_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 03-oracle/sparse-entry-pair.qec | paradigm_mapped | sparse_gate, sparse_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 03-oracle/state-prep-isometry.qec | paradigm_mapped | stateprep_gate, stateprep_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 04-block-encoding/cfd-fuse-be.qec | paradigm_mapped | be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 04-block-encoding/lcu-select.qec | paradigm_mapped | be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 04-block-encoding/mybe-wrapper.qec | paradigm_mapped | be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 04-block-encoding/product-be.qec | paradigm_mapped | be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 04-block-encoding/sparse-be.qec | paradigm_mapped | sparse_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 04-block-encoding/sum-be.qec | paradigm_mapped | be_algebra | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 04-block-encoding/walk-operator.qec | paradigm_mapped | qpe | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 05-qsvt-qpe/eigen-filter.qec | paradigm_mapped | costa_gate, costa_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 05-qsvt-qpe/oaa-std.qec | paradigm_mapped | oaa | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 05-qsvt-qpe/qpe-walk.qec | paradigm_mapped | qpe | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 05-qsvt-qpe/qsp-phases-roundtrip.qec | paradigm_mapped | qsvt, python_generators | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 05-qsvt-qpe/qsvt-monomial.qec | paradigm_mapped | qsvt, python_generators | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 05-qsvt-qpe/walk-pow2.qec | paradigm_mapped | qpe | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 06-protocols/dependency-upload.qec | paradigm_mapped | costa_qram, qham_qode | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 06-protocols/hamsim-qsvt-bind.qec | paradigm_mapped | costa_qram, qham_qode | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 06-protocols/hamsim-trotter.qec | paradigm_mapped | trotter_hamsim | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 06-protocols/protocol-generic-oracle.qec | paradigm_mapped | costa_qram, qham_qode | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 06-protocols/qlss-nested-require.qec | paradigm_mapped | costa_qram, qham_qode | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 06-protocols/solver-wrapper-step.qec | paradigm_mapped | costa_qram, qham_qode | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 06-protocols/two-implementations-order1.qec | paradigm_mapped | dj_gate, dj_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 06-protocols/two-implementations-order2.qec | paradigm_mapped | dj_gate, dj_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 07-scientific/carleman-step.qec | paradigm_mapped | carleman_step | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 07-scientific/cfd-implicit-step.qec | paradigm_mapped | qfvm_qram | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 07-scientific/heat-equation.qec | paradigm_mapped | heat_qode | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 07-scientific/lchs-ode.qec | paradigm_mapped | lchs | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 07-scientific/poisson-qlss.qec | paradigm_mapped | poisson_qlss | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 07-scientific/schrodingerization.qec | paradigm_mapped | schrodingerisation | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 07-scientific/stateprep-vs-oracle.qec | paradigm_mapped | stateprep_gate, costa_gate | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 08-modules/cross-module-protocol/main.qec | paradigm_mapped | dj_gate, dj_qram, qham_qpde | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 08-modules/import-alias/main.qec | paradigm_mapped | dj_gate, dj_qram, qham_qpde | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 08-modules/qualified-binding-key/main.qec | paradigm_mapped | dj_gate, dj_qram, qham_qpde | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 08-modules/two-module-grover/main.qec | paradigm_mapped | dj_gate, dj_qram, qham_qpde | 映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。 |
| 20-errors/alias-violation.qec | structural_rejection | tests/core/test_language.py, tests/core/test_open_ir.py | 对应 RIR 的签名、类型、别名、常量或调用图验证。 |
| 20-errors/ancilla-not-clean.qec | intentional_design_change | tests/core/test_language.py, tests/core/test_open_ir.py | 所有工作区显式持有；不继承未实现的自动复净判定。 |
| 20-errors/cyclic-dependency.qec | structural_rejection | tests/core/test_language.py, tests/core/test_open_ir.py | 对应 RIR 的签名、类型、别名、常量或调用图验证。 |
| 20-errors/fused-source-use.qec | intentional_design_change | tests/core/test_language.py, tests/core/test_open_ir.py | 逻辑视图不进行持久冻结；同次调用的重叠和控制修改仍拒绝。 |
| 20-errors/generic-infer-fail.qec | structural_rejection | tests/core/test_language.py, tests/core/test_open_ir.py | 对应 RIR 的签名、类型、别名、常量或调用图验证。 |
| 20-errors/implements-mismatch-effect.qec | structural_rejection | tests/core/test_language.py, tests/core/test_open_ir.py | 对应 RIR 的签名、类型、别名、常量或调用图验证。 |
| 20-errors/implements-mismatch-signature.qec | structural_rejection | tests/core/test_language.py, tests/core/test_open_ir.py | 对应 RIR 的签名、类型、别名、常量或调用图验证。 |
| 20-errors/isometry-adjoint-cap.qec | intentional_design_change | tests/core/test_language.py, tests/core/test_open_ir.py | 等距角色可以提供 unitary 扩张；无声明逆能力时拒绝逆调用。 |
| 20-errors/isometry-nonzero-input.qec | intentional_design_change | tests/core/test_language.py, tests/core/test_open_ir.py | 零输入是接口契约，本阶段不证明零态。 |
| 20-errors/measure-in-if.qec | intentional_design_change | tests/core/test_language.py, tests/core/test_open_ir.py | 量子结果读出在宿主层，RIR 没有动态经典分支。 |
| 20-errors/oracle-takes-operation.qec | intentional_design_change | tests/core/test_language.py, tests/core/test_open_ir.py | 由 Python 高阶生成函数承担，旧语法限制不继承。 |
| 20-errors/param-in-static-position.qec | structural_rejection | tests/core/test_language.py, tests/core/test_open_ir.py | 对应 RIR 的签名、类型、别名、常量或调用图验证。 |
| 20-errors/partial-application.qec | intentional_design_change | tests/core/test_language.py, tests/core/test_open_ir.py | 允许 Python 绑定和闭包；生成后 IR 不保留 Python callback。 |
| 20-errors/require-in-program.qec | intentional_design_change | tests/core/test_language.py, tests/core/test_open_ir.py | require 文本语法不迁移，依赖通过 Python 参数和显式绑定表达。 |
| 20-errors/require-mid-body.qec | intentional_design_change | tests/core/test_language.py, tests/core/test_open_ir.py | require 文本语法不迁移。 |
| 20-errors/unbound-require-program.qec | intentional_design_change | tests/core/test_language.py, tests/core/test_open_ir.py | 现在允许保存开放 IR，在后端导出时报告缺口。 |

六组参考负载分别由 qfvm_gate/qram、be_algebra、arithmetic、costa/sparse、oracle 目录与 qham_qode/qpde 覆盖。Roe 物理核与规则化一般 QHAM 后续已实现，见 qham-general-implementation.md；严格 QSVT 相位、HAM/PDE 收敛证明及数值认证仍待后续工作。
