# XOR 数据库（XOR Database）

> 类别 C5 · 模块 `pyqecclang.algorithms.oracles` · 阶段 V4

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

- `abstract_database`：开放声明（input model 为 QRAM 的 abstract 层），体为空，经 `linking.bind` 绑定实现后程序闭合。
- `gate_database`：真值表见证。`table` 为稀疏字典或字序列（按地址枚举），地址与字在生成期校验范围；逐条目按地址受控翻转数据位。
- `qram_database`：声明 `QRAM(address_width, data_width)` 资源（名为 `table`），执行时数据以 memory 字典提供，不进入 IR。
- `banked_database`：把 `data_width` 拆成不超过 64 位的 `data0, data1, ...` 寄存器与对应 QRAM bank（小端序），`abstract=True` 时给出带 `logical_data_width` 的开放声明。

四个入口均返回 `XorDatabase`（`OracleView`，`xor_database()` 返回自身）。属性：

| 属性 | 含义 |
|---|---|
| `operation` | 底层 `Operation`（寄存器恰为 `address`、`data`） |
| `address_width` / `data_width` | 接口位宽 |
| `spec` / `describe()` | oracle 描述（类型 `database_xor`、能力、开放/闭合状态） |

模块属性：`implementation`（`gate_truth_table` / `qram`）、`implementation_status`；banked 版另有 `logical_data_width` 与 `word_order="little_endian_banks"`。

## 实现要点

gate 版对每个非零表字在地址控制下逐位施加 X，门数随表规模线性增长，仅适合小实例见证。QRAM 版把数据外置为执行期资源：IR 中只保留一条查询原语，嵌套调用时资源名按模块前缀提升（如 `db__table`）。数据表与门级子数据库都在生成时展开、不进入 IR JSON。

下游适配：`sparse_entry` 把 `row|column` 拼接成 2w 位地址（见[稀疏访问](sparse-access.md)）；`diagonal_block_encoding` 按查询字控制信号 Ry，编码角表决定的对角矩阵；`alias_prepare` 经 `database` 参数注入 (keep, alt) 表（见 [Alias 采样制备](alias-preparation.md)）；QROM 系列返回同一 `XorDatabase` 接口（见 [QROM 查找](qrom-lookup.md)）。

## 验证方案

类别 C5（数据访问层，判定准则见 `../development/validation-plan.md` §2）：查询语义逐点正确 + 三层绑定一致 + 复净。三层证据：

- 结构：`tests/core/test_contracts.py:OracleContractTests.test_aggregate_report_and_kernel_not_executed`——`abstract_database` 输入不满足协议契约时内核不被执行（错误在运行算法之前定位）；`tests/core/test_open_ir.py:OpenIRTests.test_open_roundtrip_and_export_boundary`（开放声明经序列化往返仍保持未解析、导出被拒）。
- 数值：gate 版逐地址语义作为下游对拍基线：`tests/core/test_data_loading.py:DataLoadingTests.test_select_swap_matches_gate_database_per_address` 与 `test_result_invariant_across_partitions_and_qrom_lookup_baseline`（16 地址表逐点精确相等）；`tests/core/test_qlss_input_models.py:QLSSInputTests.test_signed_sparse_encoding_and_chebyshev` 以 `gate_database` 承载矩阵元素并逐元对拍块编码（places = 11）。
- 绑定：`tests/core/test_open_ir.py:OpenIRTests.test_partial_binding_and_lifted_resources`（`qram_database` 经 `Binding` 绑定后程序闭合、资源改名提升）与 `test_shared_capture_has_one_entry_resource`（同一实现绑定到两个槽位时共享单一资源）；契约层能力限制的存活见 `OracleContractTests.test_capability_restrictions_survive_annotation_and_wrapping`；宽字拆分见 `tests/core/test_workloads.py:WorkloadTests.test_banked_words_do_not_violate_register_storage_width`（96 位字拆为 64+32，闭合后恰两个资源）。

## 已知缺口与计划阶段

三绑定一致性参数化（V4 铺开）：同一 abstract 声明的 gate / QRAM（及下游变体）绑定尚未有统一的参数化对拍，目前各绑定独立见证。与 `validation-coverage.md` 的 oracles.py 行一致。

## 相关链接

- 源码：`src/pyqecclang/algorithms/oracles.py`
- 同组页面：[QROM 查找](qrom-lookup.md)、[Select-Swap QROM](select-swap.md)、[稀疏访问](sparse-access.md)、[态制备 Oracle](state-preparation.md)
- API 参考：[Oracle 声明与实现](../../api/algorithms/oracles.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
