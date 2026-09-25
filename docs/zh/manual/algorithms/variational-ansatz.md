# 硬件高效拟设（Hardware-Efficient Ansatz）

<a href="../../../en/index.html">English</a> · **简体中文**

> 类别 C4 · 模块 [`oracq.algorithms.optimization.variational`](../../api/algorithms/optimization/variational.rst) · 阶段 V3

## 概述

硬件高效拟设是一类参数化电路：每层先对每个量子比特施加独立的 Ry/Rz 单比特旋转，再用相邻比特间的 CNOT 链产生纠缠。角度向量由经典优化器选择，使目标泛函（能量、损失等）在拟设表达的态上最优。本模块生成给定角度的电路，参数优化与重复执行由调用方组织。

## 接口与输入模型

```python
hardware_efficient_ansatz(width, layers)
```

API 入口：{obj}`hardware_efficient_ansatz <oracq.algorithms.optimization.variational.hardware_efficient_ansatz>`

- `width`：target 位宽，范围 1..64。
- `layers`：角度张量，形状为 `[layer][qubit][Ry, Rz]`——非空、每层恰有 width 个角度对、每对恰有两个有限实数角度（弧度）。

变分角度是经典数据直接参数化，input model 为 CP。返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`，寄存器只有 `target`（width 位）。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"hardware_efficient_ansatz"` |
| `layers` | 层数 |

## 实现要点

每层按"先旋转、后纠缠"两段生成：对每个比特依次施加 `Ry(θ)` 与 `Rz(φ)`，然后对相邻比特对 (i, i+1) 逐个施加 CNOT（线路原语为 XOR）。纠缠结构固定为线性链，没有层间重排或纠缠策略开关；非链式连通的硬件由后端降低阶段插入交换。参数形状违例（层数为零、每层长度不等于 width、角度对长度不是 2、非有限角度）在生成期抛 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`。

零角度时整条电路为恒等，这是当前数值见证的锚点。适用边界：拟设只负责生成电路；表达力与可训练性（如 barren plateau）不在验证范围，也没有内置参数初始化策略。

## 验证方案

类别 C4（启发式/优化语义，判定准则见 `../../development/validation-plan.md` §2）。三层证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_bad_inputs_fail_at_generation` 覆盖参数布局违例（如每层角度对数量不足）。
- 数值：`AlgorithmExpansionTests.test_ansatz_zero_angles_and_walk_one_step`——全零角度的单层拟设作用后态保持 $\lvert 0\ldots0\rangle$（幅度精确为 1），钉死旋转与 CNOT 链在零参数下的恒等性（同测试的后半段见证 `walks.py` 的 cycle walk，见该模块）。
- 绑定：本模块无独立绑定见证（输入为经典参数）。

## 已知缺口与计划阶段

与 `validation-coverage.md` 的 `variational.py` 行一致：该行登记的唯一缺口是"达到小实例最优割"的强见证（落在 QAOA 侧，V3 接入）；拟设电路自身无登记缺口。

## 相关链接

- 源码：`src/oracq/algorithms/optimization/variational.py`
- 同模块页面：[MaxCut QAOA](qaoa-maxcut.md)、[VQE 测量电路](vqe.md)
- API 参考：[变分算法电路](../../api/algorithms/optimization/variational.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_misc_algorithms.py`（misc_algorithms 组），全部在真实后端上执行。

**实验设计**：width = 3、两层共 12 个非平凡角度（每层每比特 Ry/Rz 对）的硬件高效拟设；经典预言机为逐门构造的 numpy 态向量（独立实现 Ry/Rz/CNOT 的 exactly-same 层序与低位比特约定）。后端路径：`reference`、`rir-pysparq`、`adapter-pysparq`、`originir-ext` 四路径全振幅对拍。

**关键指标**：

| 案例 | 规模 | 路径 | 最大幅度误差 | 保真度 | 跨后端偏差 |
|---|---|---|---|---|---|
| ansatz-parameter-intent | 3 比特、2 层 | 四路径全振幅对拍 | 0 | 1 − 4e-16 | 0 |

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

产物：`out/verification/misc_algorithms.json`（24 个案例全过，本页对应 `ansatz-parameter-intent` 案例）。
