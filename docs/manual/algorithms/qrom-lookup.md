# QROM 查找（QROM Lookup）

> 类别 C5 · 模块 `pyqecclang.algorithms.data_loading` · 阶段 V4

## 概述

QROM 查找把经典表 `T` 实现为 XOR 数据库 $|a, d\rangle \mapsto |a, d \oplus T[a]\rangle$（input model 词汇中的 QRAM）。实现依据 Low–Kliuchnikov–Schaeffer 2018 与 Babbush et al. 2018（模块 docstring 引用）。本页的 `qrom_lookup` 是**逐地址受控 XOR 的一元迭代**基线：结构上即 `partitions = 1` 的 Select-Swap，作为查询深度与辅助比特权衡曲线的 λ = 1 端点；分区版见 [Select-Swap QROM](select-swap.md)。

## 接口与输入模型

```python
qrom_lookup(table, *, data_bits=None, name=None)
qrom_cost(n_addresses, data_bits, partitions=1)
```

- `table`：字序列或稀疏字典，缺失地址按 0 处理；地址与字必须是非负整数。
- `data_bits`：数据位宽，缺省取最大表字的位宽。
- `qrom_cost`：纯经典资源估算，不生成电路；`n_addresses` 为表字数 N，`partitions` 为分区数 λ（二的幂且不超过 $2^{\lceil\log_2 N\rceil}$）。

`qrom_lookup` 返回 `XorDatabase`（`implementation="qrom_unary_iteration"`），模块属性附带 `qrom_cost` 估算：`qrom_partitions`、`qrom_address_bits`、`qrom_data_bits`、`select_toffoli`、`swap_toffoli`、`t_count`、`round_trip_t_count`、`work_qubits`、`dirty_fanout_qubits`。

`qrom_cost` 返回 `QromCost` 快照：

| 字段 / 属性 | 含义 |
|---|---|
| `n_addresses` / `data_bits` / `partitions` / `address_bits` | 输入回显与地址位宽 |
| `select_toffoli` / `swap_toffoli` | select 段与交换段 Toffoli 计数（λ = 1 时后者为 0） |
| `t_count` | compute 单程 T 计数 = 4 ×（两段 Toffoli 之和） |
| `t_depth` | select 段串行深度加交换段逐级归并深度 |
| `round_trip_t_count` | 相干复净（compute + 伴随 uncompute）= 2 × `t_count` |
| `work_qubits` / `fanout_qubits` / `ancilla_qubits` | 数据窗口与低位扇出副本（后者可用 dirty qubit） |
| `to_dict()` | 供目录与后端报告使用的结构化字典 |

## 实现要点

基线结构复用 `oracles.gate_database`：逐地址受控翻转数据位，一元迭代树给出 `select_toffoli = 2^{k} - 1`（k 为地址位数），λ = 1 时 T 计数 $4(N-1)$（N 为二的幂时与公式一致）。与 `gate_database` 的差别只在标注与附带成本属性：`implementation="qrom_unary_iteration"` 使下游能区分"基线见证"与"资源优化实现"。

数据表在生成时展开为门级子数据库、不进入 IR JSON；资源估算由 `qrom_cost` 纯经典给出，生成程序的结构属性与其一致（双向检查，见验证方案）。适用边界：门数随 N 线性增长，大表应改用 `select_swap_qrom` 以辅助比特换深度。

## 验证方案

类别 C5（数据访问层，判定准则见 `../development/validation-plan.md` §2）。三层证据位于 `tests/core/test_data_loading.py:DataLoadingTests`：

- 结构：`test_generated_operations_carry_cost_attributes`——`implementation` / `qrom_partitions` / `t_count` 等属性与 `qrom_cost` 公式一致（基线 λ = 1）。
- 数值：`test_result_invariant_across_partitions_and_qrom_lookup_baseline`——16 地址表逐地址读出与表值精确相等，并作为 Select-Swap 各分区的公共基线；`test_cost_model_matches_formulas_and_tradeoff_curve` 校验 N = 16、b = 3 的公式值（λ = 1 时 t_count = 4·15）与曲线闭式 $4(N/\lambda - 1 + b(\lambda-1))$；`test_cost_report_is_structured` 校验 `to_dict`（N = 1024、b = 32、λ = 8：t_count = 4·(127 + 224)，辅助比特 277）。
- 绑定：gate / QRAM 两绑定由 `test_result_invariant_across_partitions_and_qrom_lookup_baseline` 间接覆盖（同一查询语义在不同实现 / 分区下结果不变）；三绑定一致性参数化待 V4。
- 负例：`test_invalid_inputs_fail_at_generation`——空表、负值、字宽越界、partitions 非二的幂等在生成期抛 `ValidationError`。

## 已知缺口与计划阶段

三绑定参数化（V4）：同一声明的 gate / QRAM 绑定统一参数化对拍未铺开。与 `validation-coverage.md` 的 data_loading.py 行一致。

## 相关链接

- 源码：`src/pyqecclang/algorithms/data_loading.py`
- 同组页面：[Select-Swap QROM](select-swap.md)、[XOR 数据库](xor-database.md)
- API 参考：[Select-Swap QROM 数据加载](../../api/algorithms/data_loading.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
