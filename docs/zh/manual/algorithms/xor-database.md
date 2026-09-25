# XOR 数据库（XOR Database）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C5 · 模块 [`oracq.algorithms.input_model.oracles`](../../api/algorithms/input_model/oracles.rst) · 阶段 V4

## 概述

经典数据量子访问的标准接口：给定查询表 `memory`，实现范式 `database_xor`

$$
|a,\,d\rangle \;\mapsto\; |a,\,d \oplus \mathtt{memory}[a]\rangle ,
$$

其中 `address` 与 `data` 均为 bits 寄存器。该 XOR 语义在任意初值下自逆，是函数表加载、矩阵元素查询、LCU 系数表等下游构造的公共底座（input model 词汇中的 QRAM）。模块提供 abstract 开放声明、gate 真值表与 QRAM 资源三种绑定层，另有把超过 64 位的宽字拆成多个 bank 的变体。

## 接口与输入模型

```python
abstract_database(name, address_width, data_width)
gate_database(address_width, data_width, table, *, name=None)
qram_database(address_width, data_width, *, name=None)
banked_database(address_width, data_width, *, name="BankedData", abstract=False)
```

API 入口：{obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>`、{obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>`、{obj}`qram_database <oracq.algorithms.input_model.oracles.qram_database>`、{obj}`banked_database <oracq.algorithms.input_model.oracles.banked_database>`

- {obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>`：开放声明（input model 为 QRAM 的 abstract 层），体为空，经 `linking.bind` 绑定实现后程序闭合。
- {obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>`：真值表见证。`table` 为稀疏字典或字序列（按地址枚举），地址与字在生成期校验范围；逐条目按地址受控翻转数据位。
- {obj}`qram_database <oracq.algorithms.input_model.oracles.qram_database>`：声明 {obj}`QRAM(address_width, data_width) <oracq.infrastructure.ir.QRAM>` 资源（名为 `table`），执行时数据以 memory 字典提供，不进入 IR。
- {obj}`banked_database <oracq.algorithms.input_model.oracles.banked_database>`：把 `data_width` 拆成不超过 64 位的 `data0, data1, ...` 寄存器与对应 QRAM bank（小端序），`abstract=True` 时给出带 `logical_data_width` 的开放声明。

四个入口均返回 {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>`（{obj}`OracleView <oracq.algorithms.input_model.contracts.OracleView>`，`xor_database()` 返回自身）。属性：

| 属性 | 含义 |
|---|---|
| `operation` | 底层 {obj}`Operation <oracq.infrastructure.builder.Operation>`（寄存器恰为 `address`、`data`） |
| `address_width` / `data_width` | 接口位宽 |
| `spec` / `describe()` | oracle 描述（类型 `database_xor`、能力、开放/闭合状态） |

模块属性：`implementation`（`gate_truth_table` / `qram`）、`implementation_status`；banked 版另有 `logical_data_width` 与 `word_order="little_endian_banks"`。

## 实现要点

gate 版对每个非零表字在地址控制下逐位施加 X，门数随表规模线性增长，仅适合小实例见证。QRAM 版把数据外置为执行期资源：IR 中只保留一条查询原语，嵌套调用时资源名按模块前缀提升（如 `db__table`）。数据表与门级子数据库都在生成时展开、不进入 IR JSON。

下游适配：{obj}`sparse_entry <oracq.algorithms.input_model.oracles.sparse_entry>` 把 `row|column` 拼接成 2w 位地址（见[稀疏访问](sparse-access.md)）；{obj}`diagonal_block_encoding <oracq.algorithms.input_model.oracles.diagonal_block_encoding>` 按查询字控制信号 Ry，编码角表决定的对角矩阵；{obj}`alias_prepare <oracq.algorithms.common.prepare_select.alias_prepare>` 经 `database` 参数注入 (keep, alt) 表（见 [Alias 采样制备](alias-preparation.md)）；QROM 系列返回同一 `XorDatabase` 接口（见 [QROM 查找](qrom-lookup.md)）。

## 验证方案

类别 C5（数据访问层，判定准则见 `../development/validation-plan.md` §2）：查询语义逐点正确 + 三层绑定一致 + 复净。三层证据：

- 结构：`tests/core/test_contracts.py:OracleContractTests.test_aggregate_report_and_kernel_not_executed`——`abstract_database` 输入不满足协议契约时内核不被执行（错误在运行算法之前定位）；`tests/core/test_open_ir.py:OpenIRTests.test_open_roundtrip_and_export_boundary`（开放声明经序列化往返仍保持未解析、导出被拒）。
- 数值：gate 版逐地址语义作为下游对拍基线：`tests/core/test_data_loading.py:DataLoadingTests.test_select_swap_matches_gate_database_per_address` 与 `test_result_invariant_across_partitions_and_qrom_lookup_baseline`（16 地址表逐点精确相等）；`tests/core/test_qlss_input_models.py:QLSSInputTests.test_signed_sparse_encoding_and_chebyshev` 以 `gate_database` 承载矩阵元素并逐元对拍块编码（places = 11）。
- 绑定：`tests/core/test_open_ir.py:OpenIRTests.test_partial_binding_and_lifted_resources`（`qram_database` 经 {obj}`Binding <oracq.infrastructure.linking.Binding>` 绑定后程序闭合、资源改名提升）与 `test_shared_capture_has_one_entry_resource`（同一实现绑定到两个槽位时共享单一资源）；契约层能力限制的存活见 `OracleContractTests.test_capability_restrictions_survive_annotation_and_wrapping`；宽字拆分见 `tests/core/test_workloads.py:WorkloadTests.test_banked_words_do_not_violate_register_storage_width`（96 位字拆为 64+32，闭合后恰两个资源）。

## 已知缺口与计划阶段

三绑定一致性参数化（V4 铺开）：同一 abstract 声明的 gate / QRAM（及下游变体）绑定尚未有统一的参数化对拍，目前各绑定独立见证。与 `validation-coverage.md` 的 oracles.py 行一致。

## 相关链接

- 源码：`src/oracq/algorithms/input_model/oracles.py`
- 教程：[给算法替换 oracle](../../tutorials/oracle-binding.md)
- 同组页面：[QROM 查找](qrom-lookup.md)、[Select-Swap QROM](select-swap.md)、[稀疏访问](sparse-access.md)、[态制备 Oracle](state-preparation.md)、[PREPARE–SELECT 分解](prepare-select.md)
- API 参考：[Oracle 声明与实现](../../api/algorithms/input_model/oracles.rst)
- 概念：[Oracle 与算子表示](../operators.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值验证脚本：`tests/verification/verify_oracles.py`（A 组案例），产物 `out/verification/oracles.json`。

实验设计：gate 真值表在两种模式下穷举全部输入域——**基态逐点**（全部 $2^{a+d}$ 个 $(a,d)$ 初态逐跑，断言输出恰为单一基态 $|a,\,d \oplus table[a]\rangle$）与**全叠加**（address/data 同时 Hadamard，一次运行覆盖全部分支，与经典置换逐振幅对拍）。QRAM 绑定以执行期 memory 字典走相同叠加穷举，并加非零 data 初值的基态抽点。双绑定一致性案例把同一 `abstract_database` 开放声明分别绑定 gate 与 QRAM 实现，三路执行结果两两对拍——即上文"已知缺口"中三绑定一致性参数化（V4）的首轮数值证据。后端路径：reference（内置参考执行器）、rir-pysparq（PySparQ 原生 RIR 解释器）、adapter-pysparq（oracq 事件适配器）、originir-ext（UniQC 全振幅态向量）。

| 案例 | 规模 (address×data) | 路径 | 指标值 |
|---|---|---|---|
| xor-gate-basis-2x3 | 2×3，32 初态 | reference, rir-pysparq, originir-ext | failures = 0 |
| xor-gate-basis-3x2 | 3×2，32 初态 | reference, rir-pysparq | failures = 0 |
| xor-gate-superposition-2x3 | 2×3，32 分支 | 全部四路径 | max_error = 5.6e-17 |
| xor-gate-superposition-3x4 | 3×4，128 分支 | 全部四路径 | max_error = 4.2e-17 |
| xor-gate-superposition-4x3 | 4×3，128 分支 | 全部四路径 | max_error = 4.2e-17 |
| xor-qram-superposition-3x4 | 3×4，128 分支 | 全部四路径 | max_error = 4.2e-17，basis_failures = 0 |
| xor-abstract-binding-consistency-3x2 | 3×2，gate/qram 双绑定 | reference, rir-pysparq, originir-ext | max_error = 5.6e-17 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_oracles.py
```

产物路径：`out/verification/oracles.json`。
