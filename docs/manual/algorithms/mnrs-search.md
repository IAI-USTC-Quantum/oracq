# MNRS 量子行走搜索（MNRS Quantum Walk Search）

> 类别 C1 · 模块 `pyqecclang.algorithms.graph_walks` · 阶段 V3

## 概述

marked-vertex 搜索的量子行走骨架：`setup` 制备初态 $\frac{1}{\sqrt N}\sum_v A|v\rangle$，随后迭代 `Repeat{ walk; marked 相位翻转 }`。实现依据 Magniez–Nayak–Roland–Santha 2011（Search via quantum walk，hitting time 与步数标度，模块 docstring 引用）：步数取经典随机游走平均首达时间 $H_{\text{avg}}$ 的标度 $\lceil \frac{\pi}{4}\sqrt{H_{\text{avg}}}\rceil$；完全图上 $H_{\text{avg}} = N/|M|$，步数恢复 Grover 标度 $\frac{\pi}{4}\sqrt{N/|M|}$。

## 接口与输入模型

```python
quantum_walk_search(setup, walk, marked, steps, *, name=None)
transition_matrix(neighbors)
hitting_times(transition, marked)
suggest_steps(transition, marked)
```

- `setup`：初态制备，须满足 `StatePreparation` 协议（如 [Szegedy 量子行走](szegedy-walk.md)的 `szegedy_setup` 产物）。
- `walk`：行走步句柄，其寄存器按声明顺序映射到 `setup.target` 的连续切片；可为开放声明。
- `marked`：作用于行走空间前若干位的相位 oracle（`target[, work]` 接口，由 `phase_marks(width, marked)` 构造）；对 `szegedy_walk` 而言作用的就是 `current` 寄存器。
- `steps`：非负迭代次数，可参照 `suggest_steps` 选择。
- 三个句柄都允许是开放声明，生成的程序经 `bind` 分批绑定实现后执行。输入模型整体为 FO：图经[图邻接 oracle](adjacency-oracle.md) 给出，marked 集合经相位 oracle 给出。

`quantum_walk_search` 返回 `Operation`，寄存器为 `target: Bits(width)`（width 为行走空间总宽）与 `work: Bits(prep.work_width + marked_work)`。属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"quantum_walk_search"` |
| `framework` | `"MNRS"` |
| `steps` / `walk_width` | 迭代步数与行走空间宽度 |

经典侧工具：`transition_matrix` 由邻居表得 $P[v][u] = |\{j: N(v,j) = u\}| / D$；`hitting_times` 解线性方程 $(I - P_{\text{free}})h = \mathbf{1}$ 给出各顶点首次到达 marked 集合的期望步数（marked 顶点为 0，不可达时抛 `ValidationError`）；`suggest_steps` 返回 $\max(1, \lceil \frac{\pi}{4}\sqrt{H_{\text{avg}}}\rceil)$，全图 marked 时为 0。

## 实现要点

生成结构为 setup 调用后一个 `Repeat(steps)` 循环：循环体内先把行走步的各寄存器按声明顺序映射到 `target` 的连续切片并调用，再对 `target` 的前 `marked_width` 位调用 marked 相位翻转。宽度一致性在生成期校验（walk 总宽等于 setup 宽度、marked 不超出行走空间、steps 非负），违例抛 `ValidationError`。

适用边界：步数由调用方给出，`suggest_steps` 只是参照——它是经典 hitting time 的函数，不感知量子谱；偶环（交替边染色表）上该骨架不放大，见已知缺口一节。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../../development/validation-plan.md` §2）：输出分布须与精确构造对拍。证据：

- 结构：`tests/core/test_graph_walks.py:QuantumWalkSearchTests` / `HittingTimeTests`；`test_search_rejects_mismatched_handles` 覆盖 steps 为负、setup 宽度不匹配等生成期违例。
- 数值：`QuantumWalkSearchTests.test_search_amplifies_marked_vertex_on_hypercube`——Q3 超立方体、marked {0}、`suggest_steps` 给 3 步：基线（0 步）marked 概率 0.125，3 步后 0.78125（places = 9）；`test_search_recovers_grover_limit_on_complete_graph`——K4 完全图 2 步后概率精确 1.0（places = 10）；`HittingTimeTests.test_cycle_hitting_times_match_closed_form`——8 环各顶点首达时间逐点对拍闭式 $h_k = d(8-d)$、$d = \min(k, 8-k)$（places = 9）；`test_suggest_steps_matches_grover_scaling` 与 `test_all_marked_needs_no_steps` 校验步数公式。
- 绑定：`QuantumWalkSearchTests.test_abstract_handles_bind_to_gate_implementations`——邻接 oracle 与 marked oracle 均为开放声明（`unresolved` 为 {G, M}），绑定 gate 实现后 marked 概率与直接生成一致（0.78125，places = 9）；`test_walk_handle_may_stay_abstract` 见证 setup/walk 句柄同样可以延后绑定。

## 已知缺口与计划阶段

与验证矩阵一致：Johnson 图上的 element distinctness 见证待补，阶段 V3。另一已确认事实：偶环（交替边染色邻居表）上 marked 概率恒停留在均匀基线 $1/N$、不随步数放大（参考模拟复核 8 环多组步数均为 0.125），已证实为理论结果而非实现缺陷（`validation-plan.md` §3 记录该口径）；可放大行为的见证由超立方体与完全图实例承担。

## 相关链接

- 同模块：[图邻接 oracle](adjacency-oracle.md)、[Szegedy 量子行走](szegedy-walk.md)
- 源码：`src/pyqecclang/algorithms/graph_walks.py`
- API 参考：[图行走搜索](../../api/algorithms/graph_walks.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
