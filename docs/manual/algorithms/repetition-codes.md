# 三位重复码（Repetition Codes）

> 类别 C1 · 模块 `pyqecclang.algorithms.qec.error_correction` · 阶段 V1

## 概述

把一位逻辑态 $\alpha|0\rangle + \beta|1\rangle$ 编码为三个物理位的重复码 $\alpha|000\rangle + \beta|111\rangle$，并从任意单个指定类型的错误中相干恢复。`error="bit"` 纠正单个 X 错误；`error="phase"` 是同一编码的 H 共轭版本，纠正单个 Z 错误。编码与恢复都是纯酉电路，不含测量与重置，因此逻辑幅度（包括相对相位）在"编码—错误—恢复"循环前后逐点保持。

## 接口与输入模型

```python
repetition_encode(*, error="bit")
repetition_recover(*, error="bit")
```

- `error`：`"bit"` 或 `"phase"`，编码器与恢复器必须取同一值；其他取值在生成期抛出 `ValidationError`。
- 输入模型为 CP：没有 oracle 输入，错误类型与码结构都以经典参数给出。
- 两个入口的寄存器相同：`target: Bits(1)`（逻辑位）与 `syndrome: Bits(2)`；编码要求 `syndrome` 输入为零。三个物理位按 `target`、`syndrome[0]`、`syndrome[1]` 排列。

模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"repetition_encode"` / `"repetition_recover"` |
| `error_kind` | `"bit"` 或 `"phase"` |
| `syndrome_policy` | 仅恢复器：`"retained; host reset required before reuse"` |

## 实现要点

编码器做两条 CNOT（`target` → `syndrome[0]`、`target` → `syndrome[1]`）得到三位重复；`"phase"` 变体随后对三个物理位整体做 H，把码字转到相位翻转基。恢复器先把 `target` 异或进两个 syndrome 位得到校验子：错误落在 `target` 位时两位同时为 1，受控 X（以 `syndrome == 3` 为控制条件）把 `target` 翻回逻辑值；错误落在某个 syndrome 位时校验子不全 1，`target` 本就未受影响。`"phase"` 变体在恢复前先整体 H 回到位翻转基，复用同一套逻辑。

设计约束：恢复后 syndrome 确定地保留错误位置（无错误时为 0），但**不被复净**——复用寄存器前宿主必须复位，这正是 `syndrome_policy` 属性的口径。适用边界：只保证纠正至多一个与 `error` 一致的错误；两个及以上同类型错误或混合错误超出三位码的能力。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../../development/validation-plan.md` §2）：编码—错误—恢复的复合作用须与恒等算子在任意逻辑幅度上逐点相等。证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言；`test_bad_inputs_fail_at_generation` 覆盖 `repetition_encode(error="unknown")` 的生成期违例。
- 数值：`AlgorithmExpansionTests.test_repetition_codes_preserve_arbitrary_logical_amplitudes`——对两种 `error` × 三个物理位置，用 Ry(0.73)/Rz(0.29) 制备任意逻辑态，插入单错误并恢复后，`target` 幅度精确复原为 $e^{-0.145j}\cos 0.365$ 与 $e^{0.145j}\sin 0.365$（places = 10），且 syndrome 收敛到单一确定值。
- 绑定：本算法无独立绑定见证（无开放声明入口）。

## 已知缺口与计划阶段

无已知缺口（验证矩阵缺口列为 —），阶段 V1。

## 相关链接

- 源码：`src/pyqecclang/algorithms/qec/error_correction.py`
- API 参考：[重复码与错误恢复](../../api/algorithms/qec/error_correction.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_misc_algorithms.py`（misc_algorithms 组），全部在真实后端上执行。

**实验设计**：任意逻辑态 $\mathrm{Ry}(0.73)\mathrm{Rz}(0.29)\lvert 0\rangle$（numpy 独立给出精确幅度），对两种 `error` 类型 × 4 个注入位置（无错误 / target / syndrome[0] / syndrome[1]）执行"编码—注入—恢复"：检查 target 幅度逐点复原与 syndrome 的确定性取值（位置映射 0/3/1/2）；另对无错误的编码—恢复复合作用经 `originir-ext + UniQC Circuit.to_matrix` 提取 syndrome = 0 输入块与泄漏。后端路径：`reference`、`rir-pysparq`、`adapter-pysparq`、`originir-ext`（3 量子位全振幅对拍）。

**关键指标**：

| 案例 | 规模 | 路径 | 幅度复原误差 | syndrome | 幺正块误差 / 泄漏 |
|---|---|---|---|---|---|
| repetition-bit-flip-injection | 3 比特 × 4 位置 | 四路径 | 0 | 全部确定且等于错误位置 ✓ | — |
| repetition-phase-flip-injection | 3 比特 × 4 位置 | 四路径 | ≤ 4.6e-16 | 全部确定且等于错误位置 ✓ | — |
| repetition-bit-encode-recover-unitary | 3 比特 | to_matrix | — | — | 0 / 0 |
| repetition-phase-encode-recover-unitary | 3 比特 | to_matrix | — | — | 4.4e-16 / 2.4e-17 |

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

产物：`out/verification/misc_algorithms.json`（24 个案例全过，本页对应 `repetition-*` 四个案例）。
