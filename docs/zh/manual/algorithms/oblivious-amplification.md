# Oblivious 振幅放大（Oblivious Amplitude Amplification）

[English](../../../index.html) · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.common.transforms`](../../api/algorithms/common/transforms.rst) · 阶段 V1

## 概述

对块编码 $U$ 做 oblivious 振幅放大（OAA）：不依赖初态知识的迭代式提升零信号块幅度。记 $\Pi = |0\rangle\langle 0|_{\text{signal}} \otimes I_{\text{target}}$、$R = I - 2\Pi$，$m$ 次迭代的算子为

$$
W_m \;=\; U\,\big(R\,U^\dagger R\,U\big)^{m}
$$

（$m=1$ 即文献标准三查询形式 $U\,R\,U^\dagger R\,U$；$R$ 相对 $2\Pi - I$ 的两个全局负号相消）。对零信号块为 $V/2$（$V$ 部分等距）的输入，$W_m$ 的零信号块为 $(-1)^m \sin\!\big((2m+1)\theta\big)\, V$（$\sin\theta = 1/2$）——幅值按标准 OAA 叙事 $\sin\theta \to \sin 3\theta \to \cdots$ 演化，一次迭代即可把幅度恢复到 $O(1)$。一般块编码满足切比雪夫恒等式 $\Pi W_1 \Pi = B\,(4B^\dagger B - 3I)$。该构造是 qubitization 框架的标准配套组件（框架文献 Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)；行走算子背景见 [量子化行走](qubitization-walk.md)）。

> 历史缺陷记录：2026-09 之前的版本装配 $[R\,U^\dagger R\,U]^m$（缺收尾 $U$），零信号块退化为 $2B^\dagger B - I$、不执行放大；验证轮发现后修复为上式，并以下文"库算子 = 标准序列"的回归案例钉住，详见 [数值验证](#数值验证)。

## 接口与输入模型

```python
oblivious_amplification(a, iterations=1)
```

API 入口：{obj}`oblivious_amplification <oracq.algorithms.common.transforms.oblivious_amplification>`

- `a`：{obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`，被放大的块编码（input model 为 BE）。
- `iterations`：迭代次数，默认 1。

返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`（裸操作而非 `BlockEncoding`），寄存器为 `target`（宽度 `a.width`）与 `signal`（宽度 `a.signal_qubits`）。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"oaa"` |

BE 归一化不经属性传播，放大后的幅度语义由调用方解释。

## 实现要点

程序时间顺序为：先 {obj}`invoke <oracq.algorithms.input_model.oracles.invoke>` 挂载输入 BE（初始 $U$），随后每次迭代四次结构调用：{obj}`reflect_zero(signal) <oracq.algorithms.input_model.block_encoding.reflect_zero>`（默认 `positive=False`，信号零分支取 $-1$、其余分支取 $+1$，即 $R = I - 2\Pi$）→ 伴随 `invoke` → 再次 `reflect_zero(signal)` → 收尾 `invoke`。迭代次数经 IR 的 {obj}`Repeat <oracq.infrastructure.ir.Repeat>` 结构表达，不在生成或文本序列化阶段无条件展开。

适用边界：输入必须是块编码。放大保证依赖"零信号块接近某个 $1/2$ 缩放的部分等距算子"这一结构假设，本函数不检查该假设；归一化远离 $1/2$ 的块编码迭代不会按 OAA 语义收敛，需先经缩放或 LCU 组合调整。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../../development/validation-plan.md` §2）。当前证据（与验证覆盖矩阵 `transforms.py` 行一致）：

- 结构：`tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` 与 `test_trotter_only_input_does_not_need_block_encoding` 覆盖同模块组装链的协议与寄存器契约（含 Repeat 结构保持的同类断言）。
- 数值：`tests/core/test_qsvt.py:PhaseSynthesisTests.test_convention_matches_qsvt_sequence` 钉死本模块共享的信号/反射约定（delta = 1e-10）；OAA 与 [量子化行走](qubitization-walk.md) 共用 `invoke` + `reflect_zero` 组件，由此间接覆盖。
- 绑定：无独立绑定见证（输入已要求具体 BE）。

## 已知缺口与计划阶段

无已知缺口（与验证覆盖矩阵 `transforms.py` 行一致），阶段 V1。放大语义有独立数值见证：$V/2$ 夹具经一次迭代幅值 $0.5 \to 1.0$（$\sin 3\theta$），一般块的切比雪夫恒等式与"库算子 = 标准序列"回归钉见 [数值验证](#数值验证)。

## 相关链接

- 源码：`src/oracq/algorithms/common/transforms.py`
- API 参考：[矩阵变换序列](../../api/algorithms/common/transforms.rst)
- 同族页面：[量子化行走](qubitization-walk.md)、[QSVT 相位序列](qsvt-sequence.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_fourier.py`（真实后端执行，无模拟替身）。输入 BE 为 gate 级绑定小实例。组装正确性：迭代 $it = 1, 2, 3$ 的全幺正经 OriginIR-ext + UniQC `to_matrix` 取出，与 $U\,[R\,U^\dagger R\,U]^{it}$（$R = I - 2\Pi$）逐元素对比。行为量测：构造零信号块恰为 $X/2$ 的门级最小实例（$\sin\theta = 1/2$ 的"$V/2$"夹具），一次迭代后零信号块幅值恢复 $1.0$（$-\!X$，即 $(-1)^1\sin 3\theta \cdot X$），$it = 2, 3$ 分别为 $0.5$、$0.5$（$\sin 5\theta$、$\sin 7\theta$）——幅值与标准三查询 OAA 叙事逐次一致，全局符号 $(-1)^{it}$ 作为信息性指标记录（幺正层面可观测、测量层面不可区分）。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `oaa-unitary-it1` | 8 维，1 次迭代 | originir-ext + UniQC to_matrix | max_error | 2.3e-16 |
| `oaa-unitary-it2` | 8 维，2 次迭代（Repeat 结构） | originir-ext + UniQC to_matrix | max_error | 4.4e-16 |
| `oaa-unitary-it3` | 8 维，3 次迭代（Repeat 结构） | originir-ext + UniQC to_matrix | max_error | 4.6e-16 |
| `oaa-half-block-it1` | 零块 = $X/2$ 夹具 | originir-ext + UniQC to_matrix | 零块 vs $-\sin 3\theta\,X$ | 6.1e-17 |
| `oaa-half-block-it2/3` | 同上（Repeat 结构） | originir-ext + UniQC to_matrix | 零块 vs $\pm\sin 5\theta\,X$、$\sin 7\theta\,X$ | 8.9e-16 / 1.1e-15 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_fourier.py
```

产物：`out/verification/fourier.json`。

### 补充验证：一般块的切比雪夫恒等式与标准序列回归钉（verify_hamiltonian.py）

`tests/verification/verify_hamiltonian.py`（27 案例全 PASS）给出两组互补证据，与上文 fourier 组的量测相互印证：其一，对**一般非等距块**（$2\times2$ 矩阵 BE，$B = \Pi U \Pi = A/\alpha$）在 reference / rir-pysparq / originir-ext 三条路径上验证库算子 $W_1 = U\,R\,U^\dagger R\,U$ 的切比雪夫放大恒等式 $\Pi W_1 \Pi = B\,(4B^\dagger B - 3I)$（max_error 3.6e-16）；其二，用库公开组件（`invoke` + `reflect_zero` + `adjoint`）按文献标准序列独立组装，对 $V/2$ 夹具（$V = R_z(0.4)R_y(0.9)$）的零信号块精确恢复 $-V$（max_error 1.1e-16），目标幅度由 0.4502 放大到 0.9004（恰好 ×2.0），且**库算子与脚本组装的标准序列在全部三条路径上逐振幅一致**（library_vs_script_error = 0.0）——这是"库 = 标准三查询序列"的回归钉，防止缺收尾 $U$ 的历史缺陷回退（该缺陷由验证轮发现并修复，见上文概述的历史记录）。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `oaa-iterate-block-identity` | 一般块，1 次迭代 | reference、rir-pysparq、originir-ext | 切比雪夫恒等式 max_error | 3.6e-16 |
| `oaa-standard-sequence-amplification` | $V/2$ 夹具，$URU^\dagger RU$ | reference、rir-pysparq、originir-ext | 与 $-V$ 的 max_error | 1.1e-16 |
| 同上 | — | — | 库 vs 脚本组装逐振幅偏差 | 0.0 |
| 同上 | — | — | 幅度放大（×2.0） | 0.4502 → 0.9004 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_hamiltonian.py
```

产物：`out/verification/hamiltonian.json`。
