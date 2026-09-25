# 资源估计（Toffoli+Clifford+T+QRAM）

[English](../../manual/resource-estimation.html) · **简体中文**

> 模块 [`oracq.infrastructure.estimate`](../api/infrastructure/estimate.rst) 与 [`oracq.infrastructure.backends.strict`](../api/infrastructure/backends/strict.rst)

## 概述

资源估计回答「这个程序在容错级别要花多少」：把 RIR 程序沿
OriginIR-ext → basis 的真实降低链路计数到 **Toffoli + Clifford + T + QRAM**
级别，输出量子位数、Toffoli 数、Clifford 原子数、精确 T 数、待合成旋转
原子数（及其 Ross–Selinger 型 T 开销估计）和 **QRAM 查询数**。估计器是
组合式的：按 `(module, 控制数)` 记忆化聚合，{obj}`Repeat <oracq.infrastructure.ir.Repeat>` 以符号计数相乘——
`2^40` 次重复的程序不需要展开即可估计，这正是 RIR「结构保持」设计直接
买来的能力。

同一降低链路上的 {obj}`export_strict(program) <oracq.infrastructure.backends.strict.export_strict>` 给出**可检查的编译产物**：模块
保持的严格网表（TOFFOLI/CZ/H/S/SDG/T/TDG/X/Y/Z 原子 + `RZ`/`RY` 待合成
旋转 + QRAM 行），`DEF` 调用、`Repeat` 符号与 QRAM 声明都不展开。

## 接口

开放程序使用 {obj}`estimate_resources(program, require_closed=False) <oracq.infrastructure.estimate.estimate_resources>`。报告把已知
线路成本和 `oracle_calls` 分开列出；每项记录槽名、控制数、伴随标记、资源实参
及调用次数。绑定后重新估计即可得到具体实现成本。`complete=False` 表示仍有
未实现依赖，已知门数不等于全程序总成本。

开放程序的 `qubits=None`；`qubits_lower_bound` 是已知寄存器和私有工作区
下界，`unknown_workspace` 列出未实现模块。MCX 辅助位也只覆盖已知实现。
完整程序保持原有整数位数结果。QRAM 读写按入口逻辑资源分账，沿调用边替换
形参，不把不同 bank 合并到内部同名参数。同一分析也有命令行入口，
见[开放与闭合资源分析](cli.md#开放与闭合资源分析)。

`rotations` 现在是紧凑的 {obj}`RotationCounts <oracq.infrastructure.estimate.RotationCounts>`：`rotations.total` 给出任意大小的
总数，`rotations.counts` 按 `(轴, 角度)` 保存重数。保留只读的 `len`、索引和
惰性迭代，小型程序可继续按序列读取；它不是可修改列表，也不是执行顺序。
巨大 Repeat 的计数和 JSON 报告不会展开旋转列表。

`epsilon` 表示每个旋转的合成配置，不是整个科学计算算法的总误差保证。
资源报告不为 QRAM 硬件、数据装载或物理纠错成本补上默认零值。

```python
from oracq import estimate_resources, export_strict

estimate = estimate_resources(program)          # 或 operation.estimate()
estimate.qubits            # 入口寄存器 + 模块工作区（不含 MCX 池位）
estimate.mcx_ancilla       # basis 池位（max(0, 最大控制数-1)）
estimate.toffoli           # Toffoli 总数
estimate.clifford          # H/S/SDG/X/Y/Z/CZ 原子总数
estimate.t_exact           # 精确 T/TDG 数（π/4 网格上的旋转）
estimate.rotations         # 待合成旋转原子 [(axis, angle), ...]
estimate.t_total(1e-10)    # 精确 T + 旋转 × ceil(3·log2(1/eps))
estimate.qram_queries      # 逐资源查询数；estimate.qram_total 为总和
estimate.to_dict()         # JSON 友好台账
```

计数模型与导出器逐行对齐（这不是理想化成本模型，而是**实际发射电路**
的账本）：X 型多控走 `basis.mcx` 配方（c≥3 时 2c−3 Toffoli）；受控单比特
门走 `basis.controlled_u3` 配方（c≥2 时 2(c−1) Toffoli 梯子加中间网络）；
零值控制位在每个发射行组两侧各翻转一次；每次 {obj}`Load <oracq.infrastructure.ir.Load>` 记 1 次 QRAM 查询。
角度恰为 π/4 整数倍的旋转按精确 Clifford+T 计入（S/T/Z 系列），其余记为
待合成旋转。全局相位不计（与 strict 导出一致）。

## 适用边界与已知差距

- 估计器按**当前编译器的发射配方**计数：常量加法走 O(n³) 阶梯（文献有
  O(n) 可逆加法器）、{obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>`/{obj}`qrom_lookup <oracq.algorithms.input_model.data_loading.qrom_lookup>` 走逐分支真值表网络
  （理想化 unary-iteration QROM 为 N−1 Toffoli，见
  {obj}`qrom_cost <oracq.algorithms.input_model.data_loading.qrom_cost>` 模型与 babbush2018encoding）。真实发射与理想化模型之间的
  差距是编译器 lowering 的优化空间，估计器把它变成可测量的量。
- 旋转的 T 开销是前导项估计（≈3·log2(1/ε)，Ross–Selinger 型），不含
  具体合成器实现；物理层成本（T 工厂、蒸馏、布线、表面码周期）不在
  本层范围内。
- QRAM 查询数把 QRAM 当作第一类资源：它回答「算法向内存系统发起多少次
  查询」，与「若用门合成同一表需要多少 Toffoli」互为对照（见下表）。

## 数值验证

实验设计（`tests/core/test_estimate.py` 手算计数与 strict 网表逐行对拍；
`tools/build_resource_estimates.py` 十组缩放实验，产物
`out/resource-estimates/*.json`，已注册进 `tools/check_project.py`）：

- 手算计数：阶梯加法、受控门配方、零值控制翻转、符号 Repeat（2^40）
  相乘、QRAM 逐资源查询，15 例全过；`estimate_resources` 与
  `export_strict` 的逐原子计数在平铺程序上逐项一致。
- 缩放实验（log-log 拟合斜率 vs 理论预期）：

| 组 | 规模 | 关键读数 | 拟合/预期 |
|---|---|---|---|
| `add_const` | n=4..64 | n=64 时 Toffoli 81 375 | 阶梯 ≈3.5（立方区）；文献最优 O(n) |
| `fixed_mul` | n=4..32 | Toffoli 112→6 356 | 1.94 vs O(n²) |
| {obj}`qft <oracq.algorithms.common.fourier.qft>` | n=3..12 | 旋转 3→165，精确 T 6→33 | ≈2.85（π/4 网格精确旋转在小 n 处掉出） |
| `qpe_add_const` | p=3..12 | p=12 时门数 ~5.5万 | 符号 Repeat：~2^p × cost(add) |
| {obj}`grover <oracq.algorithms.common.search.grover>` | n=4..14 | 迭代 ⌊π/4·2^{n/2}⌋ | 总门数 ~2^{n/2}·O(n) |
| `state_prep_dense` | n=2..10 | 旋转 17→5 114 | O(2^n) |
| `qram_lookup` | n=4..16 | **查询 = 1**，门 = 0 | 与规模无关 |
| `qrom_lookup_gate` | n=4..12 | n=12 时 Toffoli 344 064 | 真实发射 vs `qrom_cost` 理想模型 4 095（≈84×） |
| `diagonal_be_qram` | n=2..12 | **查询 = 2**，旋转 = 2 | 常数 |
| `diagonal_be_gate` | n=2..8 | n=8 时 Toffoli 11 388 | 同一输入模型的门级对照 |
| {obj}`roe_face <oracq.applications.roe.roe_face>` / `math_polynomial` | w=6..16 | Roe 面 w=16 时 Toffoli 446 076 | ≈1.8，字长多项式 |

- 暴露的文档级差距：`qrom_lookup` 的 `qrom_unary_iteration` 成本属性
  描述的是理想化模型，而当前发射为真值表网络——已记录为 lowering
  优化项（与 OAA、Schrödingerization 两处记录并列，属验证驱动的发现）。

复现：`PYTHONPATH=src python tools/build_resource_estimates.py`；
`python -m unittest tests.core.test_estimate -v`。
