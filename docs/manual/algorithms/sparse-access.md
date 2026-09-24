# 稀疏矩阵访问（Sparse Access）

> 类别 C5 · 模块 [`oracq.algorithms.input_model.oracles`](../../api/algorithms/input_model/oracles.rst) · 阶段 V4

## 概述

稀疏矩阵 oracle（input model 词汇中的 SO）按 CKS 的位置 / 元素分离约定组织成**两个查询操作的束**，而不是一个虚构的总酉：

- 位置 $P_A^{(1)}$（范式 `sparse_location_inplace`，原地）：$|\text{column}, i\rangle \mapsto |\text{column}, \nu(\text{column}, i)\rangle$，$\nu$ 枚举该列的非零行；
- 元素 $P_A^{(2)}$（范式 `sparse_entry_xor`）：$|\text{row}, \text{column}, d\rangle \mapsto |\text{row}, \text{column}, d \oplus A[\text{row},\text{column}]\rangle$，任意行列可查。

gate 实现要求每列给出完整置换扩张（把非零行的部分置换补全为 $[0, 2^w)$ 上的置换），QRAM 实现用正反两张查询表。

## 接口与输入模型

```python
abstract_sparse_access(name, width, value_width, sparsity, work_width=None)
sparse_location_gate(width, permutations, *, work_width=None, name=None)
sparse_location_qram(width)
sparse_entry(database, width)
```

API 入口：{obj}`abstract_sparse_access <oracq.algorithms.input_model.oracles.abstract_sparse_access>`、{obj}`sparse_location_gate <oracq.algorithms.input_model.oracles.sparse_location_gate>`、{obj}`sparse_location_qram <oracq.algorithms.input_model.oracles.sparse_location_qram>`、{obj}`sparse_entry <oracq.algorithms.input_model.oracles.sparse_entry>`

- {obj}`abstract_sparse_access <oracq.algorithms.input_model.oracles.abstract_sparse_access>`：开放声明，产出 `name_position` / `name_entry` 两个操作；位置操作带 `sparsity` 与 `full_permutation_extension=True` 属性；`work_width` 缺省取 `width`。
- {obj}`sparse_location_gate <oracq.algorithms.input_model.oracles.sparse_location_gate>`：`permutations` 按列给出 $2^w \times 2^w$ 的完整置换表。
- {obj}`sparse_location_qram <oracq.algorithms.input_model.oracles.sparse_location_qram>`：声明 `forward` / `inverse` 两张 {obj}`QRAM(2*width, width) <oracq.infrastructure.ir.QRAM>` 表。
- {obj}`sparse_entry <oracq.algorithms.input_model.oracles.sparse_entry>`：由 {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` 适配，要求 `address_width == 2*width`（`row|column` 拼接寻址）。

{obj}`SparseAccess(location, entry, width, value_width, sparsity) <oracq.algorithms.input_model.oracles.SparseAccess>` 为 frozen dataclass，约束 $1 \le \text{width} \le 64$、$1 \le \text{value\_width} \le 64$、$1 \le \text{sparsity} \le 2^{\text{width}}$，位置签名 `("column", "index", "work")`、元素签名 `("row", "column", "data")`。属性：

| 属性 | 含义 |
|---|---|
| `location` / `entry` | 两个底层 {obj}`Operation <oracq.infrastructure.builder.Operation>` |
| `width` / `value_width` / `sparsity` | 维度位数、值字宽、稀疏度 |
| `describe()` | {obj}`OracleSpec <oracq.algorithms.input_model.contracts.OracleSpec>`（类型 `cks_sparse`、`anc_qubit=None`、position / entry 组件描述、sparsity 参数） |

## 实现要点

gate 位置实现逐列受控执行置换：每个循环沿单比特差异路径分解为对换，对换再化为受控 X 的前向路径加回退路径，中间基态复原。门数随 $2^w$ 爆炸，只适合小实例见证。

QRAM 位置实现三步完成原地置换：`work ^= forward[column, index]`、`swap(index, work)`、`work ^= inverse[column, index]`——最后一步把 work 复净为零，index 持新值。正反两张表在 catalog 案例中经资源改名（`forward` / `inverse`）绑定。

`SparseAccess` 是束不是酉：整体 `anc_qubit=None`，避免把两个查询操作伪装成同一个 unitary。下游由 `qlss.SparseSystem` 消费（CKS Chebyshev / Costa walk 的稀疏输入模型），也可经 `sparse.py` 的稀疏→BE 适配转块编码；`applications/catalog.py` 与 `examples/ode_input_models.py` 给出组装样例。

## 验证方案

类别 C5（数据访问层，判定准则见 `../development/validation-plan.md` §2）。三层证据：

- 结构：`tests/core/test_contracts.py:OracleContractTests.test_sparse_is_a_bundle_not_a_fictitious_unitary`——`type == "cks_sparse"`、`anc_qubit` 为 `None`、组件描述与 sparsity 参数正确。
- 数值：`tests/core/test_qlss_input_models.py:QLSSInputTests.test_signed_sparse_encoding_and_chebyshev`——`sparse_location_gate` + `sparse_entry(gate_database(...))` 组装 2×2 对称稀疏系统后，稀疏→BE 适配的 (0,0) 块与稠密矩阵逐元对拍（places = 11，α = 2）；`test_matrix_probe_recovers_scalar_system_norm` 复用同一元素查询做范数探针（places = 10）。真实后端侧由 `tests/integration/test_qfvm_input_models.py:QfvmInputNativeTests` 消费（位置置换的完整性与 QRAM 原始字段）。
- 绑定：范式能力限制经 `OracleContractTests.test_capability_restrictions_survive_annotation_and_wrapping` 见证；QRAM 位置实现进入 catalog 案例（`tests/core/test_workloads.py:WorkloadTests.test_every_catalog_case_has_a_closed_description` 的闭合与序列化检查）。

## 已知缺口与计划阶段

三绑定一致性参数化（V4 铺开）：gate / QRAM 位置实现与 abstract 声明的统一参数化对拍未铺开，现有证据分散在各下游见证。与 `validation-coverage.md` 的 oracles.py 行一致。

## 相关链接

- 源码：`src/oracq/algorithms/input_model/oracles.py`
- 同组页面：[XOR 数据库](xor-database.md)
- API 参考：[Oracle 声明与实现](../../api/algorithms/input_model/oracles.rst)
- 概念：[Oracle 与算子表示](../operators.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值验证脚本：`tests/verification/verify_blockencoding.py`（真实后端执行，无 mock、无 skip；2026-09-16 共 73 个案例全部通过），产物 `out/verification/blockencoding.json`。本页对应 `sparse-access-layer` 与 `sparse-lookup-helpers`、`sparse-boolean-helpers`、`sparse-rotation-helpers` 共 4 个案例（访问层在块编码中的端到端表现见[稀疏矩阵块编码](sparse-block-encoding.md)数值验证节）。

实验设计（dim = 4、三对角结构位置、定点格式 3.1 有符号）：

1. 位置 oracle `sparse_location_gate`：全部 16 个 (column, index) 基态输入逐一读出，对照经典置换表，work 须复净为 0；
2. 元素 oracle `sparse_entry`：全部 16 个 (row, column) × 两种 data 初值（0 与非零 5）验证 XOR 语义 $d \oplus A_{rc}$；
3. QRAM 位置实现 `sparse_location_qram`：正/反两张表作为 memory 数据绑定，逐基态验证复净语义，reference 与 rir-pysparq 两后端逐振幅一致；
4. 访问层辅助：{obj}`reversible_lookup <oracq.algorithms.input_model.sparse.reversible_lookup>` 融合地址/数据视图 128 组输入、{obj}`batch_lookup <oracq.algorithms.input_model.sparse.batch_lookup>` 并发三路、{obj}`compare_words <oracq.algorithms.input_model.sparse.compare_words>`（eq/lt，位宽 1–4 全输入穷举）、{obj}`value_transposition <oracq.algorithms.input_model.sparse.value_transposition>`（3 位全 512 组输入）、{obj}`prefix_state <oracq.algorithms.input_model.sparse.prefix_state>` 前缀均匀叠加、{obj}`word_rotation <oracq.algorithms.input_model.sparse.word_rotation>` / {obj}`magnitude_rotation <oracq.algorithms.input_model.sparse.magnitude_rotation>` 的解析概率语义。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `sparse-access-layer`（gate 位置） | dim 4，16 基态 | reference | 置换失配数 | 0 |
| 同上（元素 XOR） | 16 × 2 组 | reference | 失配数 | 0 |
| 同上（QRAM 位置） | 16 基态 | reference + rir | 失配数 / 后端偏差 | 0 / 0 |
| `sparse-lookup-helpers` | 128 组融合视图 + 批量 | reference | 失配数 | 0 |
| `sparse-boolean-helpers` | eq/lt 穷举 + 转置 512 组 | reference | 失配数 / 前缀均匀性误差 | 0 / 1.1e-16 |
| `sparse-rotation-helpers` | 位宽 2–4 / 格式 4.1 | reference + rir | 概率语义误差 | 5.6e-16 / 1.1e-16 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_blockencoding.py
```

产物：`out/verification/blockencoding.json`。
