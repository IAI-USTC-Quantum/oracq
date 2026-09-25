# Szegedy 量子行走（Szegedy Walk）

<a href="../../../en/index.html">English</a> · **简体中文**

> 类别 C1 · 模块 [`oracq.algorithms.input_model.graph_walks`](../../api/algorithms/input_model/graph_walks.rst) · 阶段 V3

## 概述

在二部行走空间 `(current, peer, index)` 上构造 Szegedy 量化行走步

$$
W = R_B R_A, \qquad R_A = A\,(2|0\rangle\langle 0| - I)\,A^\dagger,
$$

其中零反射作用在 `peer‖index` 联合寄存器上，邻居态制备 $A = \mathrm{Adj} \cdot H^{\otimes g}$ 把 $|v\rangle|0\rangle|0\rangle$ 映到 $\frac{1}{\sqrt D}\sum_j |v\rangle|N(v,j)\rangle|j\rangle$；$R_A$ 是关于 $\mathrm{span}\{A|v\rangle|0\rangle|0\rangle\}$ 的反射，$R_B$ 交换两端角色（等价于 SWAP 共轭的同一反射）。实现依据 Szegedy 2004 对行走算符谱的量化构造（模块 docstring 引用）；marked-vertex 搜索骨架见 [MNRS 量子行走搜索](mnrs-search.md)。

## 接口与输入模型

```python
szegedy_walk(adjacency, *, name=None)
szegedy_setup(adjacency, *, name=None)
```

API 入口：{obj}`szegedy_walk <oracq.algorithms.input_model.graph_walks.szegedy_walk>`、{obj}`szegedy_setup <oracq.algorithms.input_model.graph_walks.szegedy_setup>`

- `adjacency`：[图邻接 oracle](adjacency-oracle.md)（{obj}`AdjacencyOracle <oracq.algorithms.input_model.graph_walks.AdjacencyOracle>` 或经 {obj}`as_adjacency <oracq.algorithms.input_model.graph_walks.as_adjacency>` 适配的裸操作），input model 为 FO。
- {obj}`szegedy_walk <oracq.algorithms.input_model.graph_walks.szegedy_walk>` 返回行走步 {obj}`Operation <oracq.infrastructure.builder.Operation>`，寄存器为 `current: Bits(v)`、`peer: Bits(v)`、`index: Bits(g)`。属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"szegedy_walk"` |
| `vertex_bits` / `degree_bits` | 顶点与下标位宽 |
| `composition` | `"R_B.R_A"` |

- {obj}`szegedy_setup <oracq.algorithms.input_model.graph_walks.szegedy_setup>` 返回 {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>`，制备初态 $\frac{1}{\sqrt N}\sum_v A|v\rangle$：`target: Bits(2v+g)`（低位到高位依次为 current、peer、index）、`work: Bits(0)`；标注 `state_prep_isometry`、`zero_input=True`、`clean_work=True`、`implementation="uniform_vertex_plus_neighbor"`。

## 实现要点

每个反射按 Adj → H → {obj}`reflect_zero <oracq.algorithms.input_model.block_encoding.reflect_zero>`（$2|0\rangle\langle 0| - I$）→ H → Adj 的顺序组装：`Adj` 自逆且 $A$ 整体酉，该序列恰实现 $R = A(2|0\rangle\langle 0| - I)A^\dagger$；`index` 寄存器留在行走空间内，不作为工作区复净。$R_B$ 复用同一子程序、把 current/peer 角色对调。图不必正则：非正则时的有效转移由邻居表补齐方式决定（见邻接 oracle 页）。

`szegedy_setup` 的 `target` 布局与 `szegedy_walk` 的寄存器声明顺序一致，{obj}`quantum_walk_search <oracq.algorithms.input_model.graph_walks.quantum_walk_search>` 据此把行走寄存器映射到 `target` 的连续切片。适用边界：本模块只提供单步算符与初态制备，检测 marked 顶点需要 MNRS 骨架。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../../development/validation-plan.md` §2）：标准见证技术含酉性恒等式与不动点结构。证据：

- 结构：`tests/core/test_graph_walks.py:SzegedyWalkTests.test_walk_step_is_unitary`——W 后接 W† 在三个初态基矢上复原原基态（幅度 1，places = 10），覆盖 C1 的酉性恒等式 $W^\dagger W = I$。
- 数值：`SzegedyWalkTests.test_stationary_state_is_fixed_point`——偶环的交替边染色邻居表（满足对合性 $N(N(v,j),j) = v$）上，`szegedy_setup` 制备的均匀平稳态经一步行走逐幅度不变（places = 10）。
- 绑定：无独立绑定见证；abstract → gate 的端到端绑定证据由 `QuantumWalkSearchTests.test_abstract_handles_bind_to_gate_implementations` 覆盖（见 [MNRS 量子行走搜索](mnrs-search.md)）。

## 已知缺口与计划阶段

与验证矩阵一致：Johnson 图后需新见证（element distinctness），阶段 V3。偶环上搜索不放大属理论结果，讨论见 [MNRS 量子行走搜索](mnrs-search.md)。

## 相关链接

- 同模块：[图邻接 oracle](adjacency-oracle.md)、[MNRS 量子行走搜索](mnrs-search.md)、[周期格点硬币行走](coined-cycle-walk.md)
- 源码：`src/oracq/algorithms/input_model/graph_walks.py`
- API 参考：[图行走搜索](../../api/algorithms/input_model/graph_walks.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

实验设计：8 顶点的偶环（交替边染色邻居表，满足对合性 $N(N(v,j),j)=v$，$v=3, g=1$，行走空间 7 量子位）。两层对照，经典参考为独立 numpy 组装的 $W=R_B R_A$（$R_A, R_B$ 为邻居叠加态反射，由邻居表直接构造）：

1. 幺正层：OriginIR-ext 导出经 UniQC `Circuit.to_matrix` 得全幺正，与参考矩阵逐元素对比；
2. 态层：`szegedy_setup` 初态上作用 1 步与 3 步行走，四条后端路径（reference、rir-pysparq、adapter-pysparq、originir-ext）逐振幅对比。

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| szegedy-walk-cycle8 | 7 qubits, steps=1,3 | 4 路径 + to_matrix | 幺正矩阵最大误差 | 4.1e-16 |
| 〃 | 〃 | 〃 | 演化态最大误差 | 5.3e-16 |

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

产物：`out/verification/search_walks.json`（案例 `szegedy-walk-cycle8`；MNRS 端到端搜索见 [MNRS 量子行走搜索](mnrs-search.md) 的数值验证一节）。
