# Oblivious 振幅放大（Oblivious Amplitude Amplification）

> 类别 C2 · 模块 `pyqecclang.algorithms.transforms` · 阶段 V1

## 概述

对块编码 $U$ 做 oblivious 振幅放大（OAA）：不依赖初态知识的迭代式提升零信号块幅度。记 $\Pi = |0\rangle\langle 0|_{\text{signal}} \otimes I_{\text{target}}$，单次迭代在时间顺序上依次作用 $U$、反射 $R = I - 2\Pi$、$U^\dagger$、$R$，算子为

$$
R\,U^\dagger R\,U \;=\; (2\Pi - I)\,U^\dagger\,(2\Pi - I)\,U
$$

（$R$ 相对 $2\Pi - I$ 的两个全局负号相消）。典型用途是 LCU 组合出的"两支各占一半"块编码——其零信号块约为目标算子的 $1/2$ 倍——单次迭代即可把幅度恢复到 $O(1)$。该构造是 qubitization 框架的标准配套组件（框架文献 Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)；行走算子背景见 [量子化行走](qubitization-walk.md)）。

## 接口与输入模型

```python
oblivious_amplification(a, iterations=1)
```

- `a`：`BlockEncoding`，被放大的块编码（input model 为 BE）。
- `iterations`：迭代次数，默认 1。

返回 `Operation`（裸操作而非 `BlockEncoding`），寄存器为 `target`（宽度 `a.width`）与 `signal`（宽度 `a.signal_qubits`）。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"oaa"` |

BE 归一化不经属性传播，放大后的幅度语义由调用方解释。

## 实现要点

每次迭代为四次结构调用：`invoke` 挂载输入 BE → `reflect_zero(signal)`（默认 `positive=False`，信号零分支取 $-1$、其余分支取 $+1$，即 $R = I - 2\Pi$）→ 伴随 `invoke` → 再次 `reflect_zero(signal)`。迭代次数经 IR 的 `Repeat` 结构表达，不在生成或 JSON 序列化阶段无条件展开。

适用边界：输入必须是块编码。放大保证依赖"零信号块接近某个 $1/2$ 缩放的部分等距算子"这一结构假设，本函数不检查该假设；归一化远离 $1/2$ 的块编码迭代不会按 OAA 语义收敛，需先经缩放或 LCU 组合调整。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../../development/validation-plan.md` §2）。当前证据（与验证覆盖矩阵 `transforms.py` 行一致）：

- 结构：`tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` 与 `test_trotter_only_input_does_not_need_block_encoding` 覆盖同模块组装链的协议与寄存器契约（含 Repeat 结构保持的同类断言）。
- 数值：`tests/core/test_qsvt.py:PhaseSynthesisTests.test_convention_matches_qsvt_sequence` 钉死本模块共享的信号/反射约定（delta = 1e-10）；OAA 与 [量子化行走](qubitization-walk.md) 共用 `invoke` + `reflect_zero` 组件，由此间接覆盖。
- 绑定：无独立绑定见证（输入已要求具体 BE）。

## 已知缺口与计划阶段

无已知缺口（与验证覆盖矩阵 `transforms.py` 行一致），阶段 V1。OAA 迭代自身的放大语义（$1/2$ 缩放块编码经一次迭代恢复）暂无独立数值见证，由共享组件的约定与协议测试间接覆盖。

## 相关链接

- 源码：`src/pyqecclang/algorithms/transforms.py`
- API 参考：[矩阵变换序列](../../api/algorithms/transforms.rst)
- 同族页面：[量子化行走](qubitization-walk.md)、[QSVT 相位序列](qsvt-sequence.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_fourier.py`（真实后端执行，无模拟替身）。输入 BE 为 gate 级绑定小实例。组装正确性：迭代 $it = 1, 2, 3$ 的全幺正经 OriginIR-ext + UniQC `to_matrix` 取出，与 $[R\,U^\dagger R\,U]^{it}$（$R = I - 2\Pi$）逐元素对比，与文档给出的算子公式一致。行为量测：构造零信号块恰为 $X/2$ 的门级最小实例（$\sin\theta = 1/2$ 的"$V/2$"夹具），一次迭代后实测零信号块为 $-I/2$，与实现算子的代数结果 $\Pi R U^\dagger R U \Pi = 2B^\dagger B - I$ 精确一致；但与标准三查询 OAA 语义（$\sin\theta \to \sin 3\theta$，$\theta = \pi/6$ 时幅度应恢复为 1、零块回到 $X$）的偏差为 1.0——即本模块实现的偶数次迭代是 qubitization 迭代的伴随对称形式（零块按 $2B^\dagger B - I$ 演化），不具备概述所述"1/2 缩放块编码经单次迭代恢复到 $O(1)$"的放大效应，该叙事与实现算子的差距在此记录（信息性指标，详见验证脚本注释与最终报告）。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `oaa-unitary-it1` | 8 维，1 次迭代 | originir-ext + UniQC to_matrix | max_error | 2.7e-16 |
| `oaa-unitary-it2` | 8 维，2 次迭代（Repeat 结构） | originir-ext + UniQC to_matrix | max_error | 3.8e-16 |
| `oaa-unitary-it3` | 8 维，3 次迭代（Repeat 结构） | originir-ext + UniQC to_matrix | max_error | 3.9e-16 |
| `oaa-half-block-behavior` | 零块 = $X/2$ 夹具 | originir-ext + UniQC to_matrix | implemented_operator_error | 2.2e-16 |
| 同上 | — | — | standard_oaa_deviation（信息性指标） | 1.0 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_fourier.py
```

产物：`out/verification/fourier.json`。

### 补充验证：一般块的迭代恒等式与标准序列正面恢复（verify_hamiltonian.py）

`tests/verification/verify_hamiltonian.py`（27 案例全 PASS）给出两组互补证据，与上文 fourier 组的量测相互印证：其一，对**一般非等距块**（$2\times2$ 矩阵 BE，$\Pi U\Pi = A/\alpha$）在 reference / rir-pysparq / originir-ext 三条路径上验证库迭代体 $S = R\,U^\dagger R\,U$ 的代数恒等式 $\Pi S\Pi = 2(\Pi U\Pi)^\dagger(\Pi U\Pi) - \Pi$（max_error 1.8e-16）——对 $V/2$ 型输入该恒等式退化为标量 $-I/2$，确认库迭代体本身不执行文献 OAA 放大；其二，用库公开组件（`invoke` + `reflect_zero` + `adjoint`）按文献标准序列 $U\,R\,U^\dagger R\,U$（比库迭代体多闭合一次 $U$ 调用）组装，对 $V/2$ 夹具（$V = R_z(0.4)R_y(0.9)$）的零信号块精确恢复 $-V$（max_error 1.1e-16），目标幅度由 0.4502 放大到 0.9004（恰好 ×2.0）。结论：放大语义需要闭合的奇数次 $U$ 调用，库 `oblivious_amplification` 的迭代体缺少末次 $U$，疑似实现缺陷（库未修改，在此与验证产物中记录）。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `oaa-iterate-block-identity` | 一般块，1 次迭代 | reference、rir-pysparq、originir-ext | 恒等式 max_error | 1.8e-16 |
| `oaa-standard-sequence-amplification` | $V/2$ 夹具，$URU^\dagger RU$ | reference、rir-pysparq、originir-ext | 与 $-V$ 的 max_error | 1.1e-16 |
| 同上 | — | — | 幅度放大（0.5 → 1.0 的模长比） | 0.4502 → 0.9004（×2.0） |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_hamiltonian.py
```

产物：`out/verification/hamiltonian.json`。
