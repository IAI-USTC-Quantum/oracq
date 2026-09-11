# VQE 测量电路（VQE Pauli Measurements）

> 类别 C4 · 模块 `pyqecclang.algorithms.variational` · 阶段 V3

## 概述

变分量子本征求解器（VQE）用参数化态 $\lvert\psi\rangle$ 逼近哈密顿量 $H=\sum_k c_k P_k$ 的基态能量，其中 $P_k$ 是 Pauli 字（I/X/Y/Z 的逐位张量积）、$c_k$ 为实系数。量子侧只需对每个 $P_k$ 反复测量 $\langle P_k\rangle$，经典侧加权求和得能量 $E=\sum_k c_k\langle P_k\rangle$ 并驱动参数优化。

本模块实现测量电路族：`pauli_measurement` 把制备态旋转到指定 Pauli 测量基，`vqe_measurements` 为整组 Hamiltonian 项批量生成电路。能量汇总与参数优化在经典侧进行，模块本身不含优化器。

## 接口与输入模型

```python
vqe_measurements(preparation, terms)
pauli_measurement(preparation, word)
```

- `preparation`：参数已实例化的 ansatz 或其他态制备（`StatePreparation`），input model 为 SP；必须满足零输入、干净工作区契约（`checked_state_preparation` 检查，违例抛 `ContractError`）。
- `terms`：`(实系数, Pauli 字)` 列表；系数须为有限实数，Pauli 字与态同宽（I/X/Y/Z），空列表抛 `ValidationError`。
- `word`：同宽的 I/X/Y/Z 字符串，第一个字符对应最低位。

`vqe_measurements` 返回元组，元素为 `(系数, 测量 Operation)`。每个测量电路的寄存器为 `target`（态宽度）与 `work`（制备的工作区）。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"pauli_measurement"` |
| `pauli_word` | 被测的 Pauli 字 |
| `readout_register` | `"target"` |

## 实现要点

生成链：先调用制备 oracle 填充 `target`（work 语义继承制备），再做逐位基变换——Y 位先施加角度 $-\pi/2$ 的相位门（$S^\dagger$）再 Hadamard，X 位只做 Hadamard，I 位不动。变换后对 `target` 做标准 Z 测量，读取非 I 位的 Z 奇偶性即可估计该 Pauli 字的期望。

每个 Pauli 字生成一条独立电路，重复测量与统计由调用方组织；同名项可复用同一电路的采样，模块本身不做 qubit-wise commuting 分组。适用边界：只提供测量电路，期望值的置信区间与参数优化策略不在本模块内。

## 验证方案

类别 C4（启发式/优化语义，判定准则见 `../../development/validation-plan.md` §2）。三层证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言；同类 `test_bad_inputs_fail_at_generation` 覆盖本模块入口的生成期违例（Pauli 字宽/字符校验在源码中生成期执行，无独立负例测试）。
- 数值：`AlgorithmExpansionTests.test_vqe_measurement_of_y_eigenstate`——制备 Y 本征态 $\lvert y_+\rangle$（Hadamard 后接 $\pi/2$ 相位门），对项 `(2, "Y")` 生成测量电路，系数回显为 2，`target` 读出 0 的概率为 1，即 $\langle Y\rangle=+1$ 恰被恢复。
- 绑定：本模块无独立绑定见证（输入已是具体态制备，无抽象槽位）。

## 已知缺口与计划阶段

与 `validation-coverage.md` 的 `variational.py` 行一致：该行登记的唯一缺口是"达到小实例最优割"的强见证（落在 QAOA 侧，V3 接入）；VQE 测量电路自身无登记缺口。

## 相关链接

- 源码：`src/pyqecclang/algorithms/variational.py`
- 同模块页面：[MaxCut QAOA](qaoa-maxcut.md)、[硬件高效拟设](variational-ansatz.md)
- API 参考：[变分算法电路](../../api/algorithms/variational.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
