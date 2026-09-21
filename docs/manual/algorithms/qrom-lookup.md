# QROM 查找（QROM Lookup）

> 类别 C5 · 模块 `pyqecclang.algorithms.input_model.data_loading` · 阶段 V4

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
- API 参考：[Select-Swap QROM 数据加载](../../api/algorithms/input_model/data_loading.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

本节结果由 `tests/verification/verify_stateprep.py` 在真实后端上跑出（reference、rir-pysparq、adapter-pysparq、originir-ext + UniQC 四条独立路径）。

**实验设计**（14 个案例）：真值表穷举而非逐基态循环——对地址寄存器施加 Hadamard 后一次查询即并行读出整张表（期望态 $\frac{1}{\sqrt{N}}\sum_a |a, T[a]\rangle$），地址+数据全叠加进一步穷举 XOR 语义 $|a, d\rangle \mapsto |a, d \oplus T[a]\rangle$。实例：`qrom_lookup` 的 16 地址稠密表、缺失地址按 0 处理的稀疏字典表、64 地址宽表；`select_swap_qrom` 的 16 地址表在 $\lambda \in \{1, 2, 4, 8, 16\}$ 全分区曲线下逐振幅穷举，及 64 地址表 $\lambda \in \{2, 8\}$；`qram_database`（QRAM 资源绑定，memory 提供表数据）同样做两种叠加穷举。$\lambda = 8/16$ 与 64 地址 $\lambda = 8$ 的工作区分别为 31/42/55 qubit，超 OriginIR-ext 的 24 qubit 预算（由 `workspace_table` 预判），这些案例只走 pysparq 路径并在案例参数中注明。`qrom_cost` 与论文闭式公式 $4(N/\lambda - 1 + b(\lambda-1))$、T 深度 $N/\lambda - 1 + (\lambda - 1)$、辅助比特 $\lambda b + (\lambda-1)\log_2\lambda$ 在 $(N, b, \lambda)$ 网格上双向对拍，并校验生成操作的模块属性与估算一致。

**关键指标**：

| 案例 | 规模 | 路径 | max_error |
|---|---|---|---|
| qrom-lookup-truth-table-n16 | N=16, b=3 | 四路径 | 8.3e-17 |
| qrom-lookup-sparse-dict | N=8, b=3（稀疏表） | 四路径 | 5.6e-17 |
| qrom-lookup-wide-n64 | N=64, b=4 | 四路径 | 5.6e-17 |
| select-swap-truth-table-l1..l16 | N=16, b=3, λ∈{1,2,4,8,16} | ref+rir+adapter（λ≤4 加 originir） | ≤ 8.4e-17 |
| select-swap-xor-superposition-l4 | N=16, b=3，128 分支 | 四路径 | 4.2e-17 |
| select-swap-wide-n64-l2 / l8 | N=64, b=4 | ref+rir+adapter（l2 加 originir） | ≤ 5.6e-17 |
| qrom-cost-formula | 9 组 (N,b,λ) 网格 + 生成属性 | 纯经典 | failures = 0 |
| qram-database-truth-table / xor-superposition | N=4, b=3（QRAM 绑定） | 四路径 | ≤ 1.2e-16 |

所有分区的查询结果与真值表在机器精度量级逐振幅一致，窗口局部寄存器复净（模拟器在 LocalExit 强制检查），成本模型与闭式公式完全一致。

**复现命令**：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_stateprep.py
```

**产物路径**：`out/verification/stateprep.json`。
