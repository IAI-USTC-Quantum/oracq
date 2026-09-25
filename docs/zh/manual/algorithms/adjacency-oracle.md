# 图邻接 oracle（Adjacency Oracle）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C1 · 模块 [`oracq.algorithms.input_model.graph_walks`](../../api/algorithms/input_model/graph_walks.rst) · 阶段 V3

## 概述

图行走框架的输入模型：以邻居表 oracle 访问图，查询语义为

$$
|v\rangle|j\rangle|c\rangle \mapsto |v\rangle|j\rangle|c \oplus N(v,j)\rangle,
$$

即给定顶点 $v$ 与出边下标 $j$，把第 $j$ 个邻居编号异或进 `neighbor` 寄存器。设 $g$ 为出边下标位宽、$D = 2^g$，每个顶点的出边表统一补齐到恰好 $D$ 列（通常用自环补齐），使邻居叠加态 $\frac{1}{\sqrt D}\sum_j |N(v,j)\rangle|j\rangle$ 无需知道具体度数即可归一化。该输入模型（FO，XOR 查询范式 `database_xor`）同时服务 [Szegedy 量子行走](szegedy-walk.md)与 [MNRS 量子行走搜索](mnrs-search.md)。

## 接口与输入模型

核心类型 {obj}`AdjacencyOracle <oracq.algorithms.input_model.graph_walks.AdjacencyOracle>`（冻结 dataclass，{obj}`OracleView <oracq.algorithms.input_model.contracts.OracleView>` 子类）封装一个 `(vertex, index, neighbor)` 签名的操作并复用 `database_xor` 范式，构造时校验签名及 `vertex`/`neighbor` 等宽。三个构造入口：

```python
abstract_adjacency(vertex_bits, degree_bits, *, name=None)   # 开放声明
gate_adjacency(neighbors, *, name=None)                      # 门级受控 X 网络
qram_adjacency(vertex_bits, degree_bits, *, name=None)       # QRAM 表
```

API 入口：{obj}`abstract_adjacency <oracq.algorithms.input_model.graph_walks.abstract_adjacency>`、{obj}`gate_adjacency <oracq.algorithms.input_model.graph_walks.gate_adjacency>`、{obj}`qram_adjacency <oracq.algorithms.input_model.graph_walks.qram_adjacency>`

- {obj}`abstract_adjacency <oracq.algorithms.input_model.graph_walks.abstract_adjacency>`：体为空的开放声明，`vertex_bits` 至多 32、`degree_bits` 取 0..32，可经 {obj}`bind <oracq.infrastructure.linking.bind>` 绑定 gate/QRAM 实现。
- {obj}`gate_adjacency <oracq.algorithms.input_model.graph_walks.gate_adjacency>`：小规模邻居表的门级实现；表必须是非空矩形且每项为合法顶点编号，位宽由表规模自动推出。
- {obj}`qram_adjacency <oracq.algorithms.input_model.graph_walks.qram_adjacency>`：QRAM 实现，表项按地址 `vertex | (index << vertex_bits)` 寻址。
- {obj}`as_adjacency(value) <oracq.algorithms.input_model.graph_walks.as_adjacency>` 把裸 {obj}`Operation <oracq.infrastructure.builder.Operation>` 适配为 `AdjacencyOracle`。

`AdjacencyOracle` 的属性：

| 属性 | 含义 |
|---|---|
| `vertex_bits` / `degree_bits` | 顶点与出边下标位宽 $v$ / $g$ |
| `vertices` / `degree` | $2^v$ 与 $2^g$ |
| `xor_database()` | 适配为 {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>`：`address = vertex‖index`（宽 $v+g$）、`data = neighbor` |

## 实现要点

`gate_adjacency` 枚举全部 $(v, j)$ 分支做受控 X 网络，仅对非零的目标位发射门，因此任意表上的查询语义与真值逐点一致。`qram_adjacency` 是单条 QRAM Load 查表。补齐是调用方的责任：框架只要求每行恰好 $D$ 列，行内重复（自环）不改变 XOR 语义，但决定 Szegedy 行走的有效转移。`xor_database()` 把图输入接入 `XorDatabase` 生态，可复用 QRAM 绑定与数据加载工具。

适用边界：`gate_adjacency` 的门数随表规模增长，面向小图与测试；更大的图用 `qram_adjacency` 或保持 abstract 声明延后绑定。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../../development/validation-plan.md` §2）：查询语义须与邻居表逐点相等。证据：

- 结构：`tests/core/test_graph_walks.py:AdjacencyOracleTests`；`test_bad_tables_fail_at_generation` 覆盖非矩形表、顶点越界、marked 集合不可达导致首达时间发散等生成期违例。
- 数值：`AdjacencyOracleTests.test_gate_adjacency_implements_neighbor_table`——对 Q3 超立方体邻居表（8 顶点 × 4 列，自环补齐到 $D=4$）穷举全部 32 个 $(v, j)$ 查询，读出幅度精确为 1（places = 10）；`test_xor_database_adapter` 校验适配后的地址/数据宽度与查询值。
- 绑定：`AdjacencyOracleTests.test_qram_adjacency_matches_gate_implementation`（gate 与 QRAM 两实现逐点对拍）与 `QuantumWalkSearchTests.test_abstract_handles_bind_to_gate_implementations`（abstract 声明绑定 gate 实现后端到端语义不变，见 [MNRS 量子行走搜索](mnrs-search.md)）共同构成邻接 oracle 的 gate/qram 双绑定证据。

## 已知缺口与计划阶段

与验证矩阵一致：Johnson 图（element distinctness 依赖）上尚无见证，需补新见证后接入，阶段 V3。

## 相关链接

- 同模块：[Szegedy 量子行走](szegedy-walk.md)、[MNRS 量子行走搜索](mnrs-search.md)、[周期格点硬币行走](coined-cycle-walk.md)
- 源码：`src/oracq/algorithms/input_model/graph_walks.py`
- API 参考：[图行走搜索](../../api/algorithms/input_model/graph_walks.rst)
- 概念：[Oracle 与算子表示](../operators.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

实验设计（端到端行走数值，邻接查询本身已由单测覆盖）：以超立方体 Q3 邻居表（8 顶点、度 4）的 `gate_adjacency` 为例，对 vertex/index 寄存器加 H 做一次叠加调用，穷举全部 32 个 $(v,j)$ 查询分支，四条后端路径（reference、rir-pysparq、adapter-pysparq、originir-ext）逐振幅核对 $\lvert v\rangle\lvert j\rangle\lvert N(v,j)\rangle$ 结构。另以 QRAM 绑定（`qram_adjacency` + 内存表）跑通 MNRS 搜索端到端（K4，steps=2，OriginIR-ext 不含 QRAM 资源，该实例仅 reference 与 rir-pysparq 两路径）。

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| adjacency-superposition-hypercube | 8 顶点 × 4 列 = 32 查询 | 4 路径 | 振幅最大误差 / 表外分支权重 | 5.6e-17 / 0.0 |
| mnrs-search-qram-k4 | N=4, D=4, steps=2 | reference, rir-pysparq | marked 概率 / 末态最大误差 | 0.578125 / 6.3e-16 |

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

产物：`out/verification/search_walks.json`（案例 `adjacency-superposition-hypercube`、`mnrs-search-qram-k4`）。
