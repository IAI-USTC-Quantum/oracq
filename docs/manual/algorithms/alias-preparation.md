# Alias 采样制备（Alias Sampling Preparation）

> 类别 C5/C1 · 模块 `pyqecclang.algorithms.common.prepare_select` · 阶段 V4

## 概述

Babbush et al. 2018 的 alias 采样 PREPARE（模块 docstring 引用；input model 词汇中的 CP→SP / QRAM）：把 LCU 概率 $p_i = |c_i|/\alpha$ 经 Vose alias 预处理打包成 `(keep, alt)` 表，量子侧只需均匀态、一次表加载、一次比较与一次受控交换，避免多路复用旋转树。`keep` 按 `precision` 位向下取整量化，采样分布与精确分布的总变差不超过 $2^{\text{selector}} \cdot 2^{-\text{precision}}$。

## 接口与输入模型

```python
alias_table(coefficients, *, precision=8)
alias_prepare(coefficients, *, precision=8, database=None)
```

- `alias_table`：纯经典预处理。系数先补零到 $2^{\text{selector}}$ 再按均值分裂 keep / alt（零概率槽也可工作），数据字打包为 `keep | (alt << precision)`；`precision` 2..32 且 selector 位宽 + precision ≤ 64（QRAM 字宽上限）。
- `alias_prepare`：量子侧组装。`database` 缺省为 `qram_database(width, precision + width)` 并返回默认 memory（键 `db__table`），也可注入 `gate_database` 等任何 `XorDatabase`（须 address = selector 位宽、data = precision + selector 位宽，否则生成期拒绝）。

返回对象：

| 对象 / 属性 | 含义 |
|---|---|
| `alias_table` → `AliasTable` | `probabilities` / `keep` / `alt` / `quantized` / `precision` / `table`；`distribution()` 给出量化 keep 与均匀抽取下的经典采样分布 |
| `alias_prepare` → `AliasPreparation` | `.preparation`（`StatePreparation`）、`.table`（`AliasTable`）、`.memory`（默认 QRAM 绑定数据）；`.state_preparation()` 返回句柄 |

`AliasPreparation.preparation` 的模块属性：`implementation="alias_sampling"`、`alias_precision`、`prepare_terms` / `prepare_alpha` / `selector_width`、`clean_work=False`。

## 实现要点

寄存器布局：`target = selector(width)`，`work = data(precision + width) | coin(precision)`；`data` 低 precision 位是 keep、高位是 alt。流程：对 target 与 coin 施加 H → 表查询 `data ^= (keep|alt)[target]` → 无符号比较 `coin < keep` 得 flag → flag = 0 时受控 `swap(target, alt)`（选中 alias 槽）→ 再次调用比较复净 flag。

data / coin 与 selector 纠缠，构成 $\sum_i \sqrt{p_i}\,|i\rangle|\mathrm{junk}_i\rangle$ 的**纯化 junk**：`clean_work=False` 诚实标注 work 不复净。下游只消费边缘分布或块编码 (0,0) 块时 junk 内积不产生影响，junk 由伴随 PREPARE 复净。比较器复用 `fixed_arithmetic("lt", FixedFormat(precision, 0, signed=False))`。

与 `lcu_prepare_select` 组合时，BE 的 signal 位宽为 selector 位宽加 work 宽度（`prepare_work` 属性）；以 4 系数、precision = 8 为例 target 宽 2、work 宽 18（data 10 + coin 8），默认 QRAM 路线的 memory 键为 `db__table`。gate 路线把整张表展开为门级子数据库，适合小实例见证。适用边界：数据表不进 IR；默认 QRAM 路线执行时须与电路一同提供 memory。

## 验证方案

类别 C5/C1：制备的 selector 边缘分布与经典采样分布闭式对拍，块编码组装的 (0,0) 块按 C1 与稠密矩阵逐点对拍（判定准则见 `../development/validation-plan.md` §2）。三层证据位于 `tests/core/test_prepare_select.py:PrepareSelectTests`：

- 结构：`test_alias_preparation_samples_target_distribution` 断言 `implementation == "alias_sampling"`、`clean_work` 为 False、返回的 `.table` 与独立构造一致；`test_bad_inputs_fail_at_generation` 覆盖 database 宽度不匹配的生成期拒绝。
- 数值：`test_alias_table_reproduces_normalized_distribution`——经典侧由 keep / alt 重构的概率与 $p_i$ 对拍（places = 12）、量化分布总变差 delta = 0.02、数据字不超过 `precision + selector` 位；`test_alias_preparation_samples_target_distribution`——gate_database 与默认 QRAM 两条路线的 selector 边缘分布均与 `table.distribution()` 对拍（places = 10）；`test_alias_block_encoding_approaches_hamiltonian`——alias PREPARE 组装的块编码 (0,0) 块乘 α 后与稠密哈密顿量对拍（delta = 0.02）。
- 绑定：数据表经 `XorDatabase` 接口注入，gate 表与 QRAM 资源两条路线在同一测试内对拍；三绑定一致性参数化待 V4。

## 已知缺口与计划阶段

alias 三绑定一致性（gate / QRAM 当前缺独立场景）：统一的参数化对拍未铺开，两条路线的分布对拍目前内嵌在数值见证中。阶段 V4，与 `validation-coverage.md` 的 prepare_select.py 行一致。

## 相关链接

- 源码：`src/pyqecclang/algorithms/prepare_select.py`
- 同组页面：[PREPARE–SELECT 分解](prepare-select.md)、[态制备 Oracle](state-preparation.md)、[XOR 数据库](xor-database.md)
- API 参考：[PREPARE-SELECT 分解](../../api/algorithms/common/prepare_select.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

本节结果由 `tests/verification/verify_stateprep.py` 在真实后端上跑出（reference 内置参考执行器与 rir-pysparq 原生 RIR 解释器；alias 程序工作区达 51 qubit，超出 OriginIR-ext 的 24 qubit 预算，故走稀疏模拟路径）。

**实验设计**（4 个案例）：对 selector（target）边缘分布做对拍，oracle 为两个独立构造的经典分布——精确分布 $p_i = |c_i|/\alpha$（由系数直接计算）与量化 alias 分布 $q_i \propto \text{quantized}_i + \sum_{j:\,\text{alt}_j=i}(2^p - \text{quantized}_j)$（由 keep/alt 表按定义独立展开，非调用 `AliasTable.distribution()`；后者与闭式的一致性作为信息性指标一并报告）。实例：4 系数（含复数）precision 8 的 gate_database 与默认 QRAM 双绑定、6 系数（含一个零概率槽，补零到 8）precision 10 的 gate 绑定、以及 keep 恰好全为 1 的均匀分布特例。量子侧一次运行给出全部 $2^{\text{selector}+\text{precision}}$ 个分支的边缘分布。

**关键指标**：

| 案例 | 规模 | 绑定 / 路径 | TVD（vs 量化分布） | TVD（vs 精确分布） | 理论界 $2^w \cdot 2^{-p}$ |
|---|---|---|---|---|---|
| alias-preparation-gate-w2-p8 | w=2, p=8 | gate / ref+rir | 8.3e-17 | 1.33e-3 | 1.56e-2 |
| alias-preparation-qram-w2-p8 | w=2, p=8 | QRAM / ref+rir | 8.3e-17 | 1.33e-3 | 1.56e-2 |
| alias-preparation-gate-w3-p10 | w=3, p=10（含零槽） | gate / ref+rir | 7.6e-17 | 1.63e-4 | 7.81e-3 |
| alias-preparation-uniform-p8 | w=2, p=8 均匀 | gate / ref+rir | 5.6e-17 | 5.6e-17 | 1.56e-2 |

量子边缘分布与量化 alias 分布在机器精度量级一致（gate 与 QRAM 两绑定的跨路径 TVD 为 0），与精确分布的 TVD 由 keep 量化主导且远小于理论界；均匀特例下量化无损，TVD 降到机器精度。`distribution_self_consistency` 全为 0，表明模块自带 `distribution()` 与独立闭式一致。

**复现命令**：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_stateprep.py
```

**产物路径**：`out/verification/stateprep.json`。
