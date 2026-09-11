# Select-Swap QROM（Select-Swap QROM）

> 类别 C5 · 模块 `pyqecclang.algorithms.data_loading` · 阶段 V4

## 概述

Select-Swap QROM 是 QROM 查找的资源优化实现（Low–Kliuchnikov–Schaeffer 2018；Babbush et al. 2018，模块 docstring 引用）：把地址拆为高位 `h`（k 位）与低位 `y`（l 位），$\lambda = 2^l$ 个窗口分区共享高位地址并行加载子表 $T[h\cdot\lambda + i]$，再按 $y = i$ 受控交换归并，以 T 计数约 $4(2^k + \lambda b)$（b 为数据位宽）在查询深度与辅助比特之间权衡（input model 词汇中的 QRAM）。

## 接口与输入模型

```python
select_swap_qrom(table, *, partitions, data_bits=None, name=None)
qrom_cost(n_addresses, data_bits, partitions=1)
```

- `table`：字序列或稀疏字典，缺失地址按 0 处理。
- `partitions`：分区数 λ，必须是二的幂且不超过地址数 $2^{\lceil\log_2 N\rceil}$。
- `data_bits`：数据位宽，缺省取最大表字的位宽。

返回 `XorDatabase`（`implementation="select_swap"`），对外接口仍是 `address` / `data` 两个 bits 寄存器，XOR 语义与基线完全一致。模块属性附带 `qrom_cost` 估算（字段表见 [QROM 查找](qrom-lookup.md)），其中 `work_qubits = λ·b` 为窗口寄存器、`fanout_qubits = (λ-1)·l` 为低位扇出副本（可用 dirty qubit）。

## 实现要点

三个阶段：**select**——λ 个 clean 局部窗口寄存器 `window_i` 各自执行子表查询 `window_i ^= T[h·λ+i]`，全部窗口共享高位地址，因此一元迭代只支付一次 $2^k - 1$（k = 0 时子表为常数，退化为 X 门）；**merge**——按 `y == i` 对 `(window_i, data)` 施加受控 swap–xor–swap 复合，净效果是受控的 `data ^= window_i`，在 data 任意初值下保持 XOR 语义且窗口内容不被破坏；**uncompute**——伴随重放子查询把全部窗口复净为零。

以测试基线（N = 16、b = 3）为例，`qrom_cost` 给出的权衡曲线：

| λ | select | swap | t_count | t_depth | work | fanout | 辅助比特 |
|---|---|---|---|---|---|---|---|
| 1 | 15 | 0 | 60 | 15 | 3 | 0 | 3 |
| 2 | 7 | 3 | 40 | 8 | 6 | 1 | 7 |
| 4 | 3 | 9 | 48 | 6 | 12 | 6 | 18 |
| 8 | 1 | 21 | 88 | 8 | 24 | 21 | 45 |
| 16 | 0 | 45 | 180 | 15 | 48 | 60 | 108 |

t_count 在 λ = 4 处取谷（测试 `test_cost_model_matches_formulas_and_tradeoff_curve` 断言的两个不等式），t_depth 同样在中间分区最小；λ = 16 时高位耗尽、select 段退化为纯交换。适用边界：λ 增大降低 T 深度但线性增加辅助比特，最优 λ 依赖硬件的比特 / 深度折算。

## 验证方案

类别 C5（数据访问层，判定准则见 `../development/validation-plan.md` §2）。三层证据位于 `tests/core/test_data_loading.py:DataLoadingTests`（16 地址表 `TABLE16`，b = 3）：

- 结构：`test_generated_operations_carry_cost_attributes`——`select_toffoli` / `swap_toffoli` / `t_count` / `work_qubits` / `dirty_fanout_qubits` 等属性与 `qrom_cost(16, 3, 4)` 逐项一致。
- 数值：`test_select_swap_matches_gate_database_per_address`——λ ∈ {1, 2, 4, 8, 16} 下全部 16 个地址的读出与 `gate_database` 基线逐点相等；`test_xor_semantics_with_nonzero_data_register`（data 初值 5 时输出为 `T[a] XOR 5`）；`test_superposition_query_restores_clean_work`（地址均匀叠加上查询后边缘分布均匀，places = 12，窗口由 LocalExit 强制复净）；`test_cost_model_matches_formulas_and_tradeoff_curve`——λ = 4 处于曲线谷底（curve[4] < curve[1] = 4·15 且 < curve[16] = 4·45），并断言曲线与闭式 $4(N/\lambda - 1 + b(\lambda - 1))$ 逐点一致。
- 绑定：gate / QRAM 两绑定由 `test_result_invariant_across_partitions_and_qrom_lookup_baseline` 间接覆盖（结果不随分区与实现变化）；三绑定一致性参数化待 V4。
- 负例：`test_invalid_inputs_fail_at_generation`——partitions 非二的幂 / 越界、空表、负值、字宽越界等在生成期抛 `ValidationError`。

## 已知缺口与计划阶段

三绑定参数化（V4）：同一声明的 gate / QRAM 绑定统一参数化对拍未铺开；validation-plan §5 另建议把 `qrom_cost` 纳入目录报告属性。与 `validation-coverage.md` 的 data_loading.py 行一致。

## 相关链接

- 源码：`src/pyqecclang/algorithms/data_loading.py`
- 同组页面：[QROM 查找](qrom-lookup.md)、[XOR 数据库](xor-database.md)
- API 参考：[Select-Swap QROM 数据加载](../../api/algorithms/data_loading.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
