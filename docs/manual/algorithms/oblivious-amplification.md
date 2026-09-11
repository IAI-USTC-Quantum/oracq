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
