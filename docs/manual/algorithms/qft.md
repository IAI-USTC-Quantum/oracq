# 量子 Fourier 变换（Quantum Fourier Transform）

> 类别 C1 · 模块 `pyqecclang.algorithms.fourier` · 阶段 V1

## 概述

$n$ 位正号离散 Fourier 变换，矩阵元为

$$
F_{y x} = \frac{1}{\sqrt{2^n}} \exp\!\bigl(2\pi i \cdot x y / 2^n\bigr),
$$

即 $|x\rangle \mapsto \frac{1}{\sqrt{2^n}} \sum_y \exp(2\pi i \cdot x y / 2^n) |y\rangle$。实现采用标准的 Hadamard + 受控相位结构并附带末尾位序交换，位序约定为 little endian。QFT 是 [Fourier 加法](fourier-addition.md)与相位估计等模块的基础组件。

## 接口与输入模型

```python
qft(width)
qft_with_work(width)
inverse_qft(width)
```

- `qft(width)`：`width` 范围 1..64，非正整数或超界在生成期抛出 `ValidationError`。返回 `Operation`，只有 `target: Bits(width)` 寄存器。
- `qft_with_work(width)`：保留早期零宽 `work` 接口的适配，寄存器为 `target: Bits(width)` 与 `work: Bits(0)`，内部直接调用 `qft`。
- `inverse_qft(width)`：QFT 的伴随操作（逆变换），以 adjoint 上下文包裹模块调用实现。

三个入口都不消费 oracle 输入（无 input model），位宽即全部参数。

## 实现要点

`qft` 的生成顺序：从最高位向最低位遍历 `high`，先对 `target[high]` 加 H，再对每个 `low < high` 施加以 `target[low]` 控制的相位门，角度 $\pi / 2^{\text{high} - \text{low}}$；随后用 $\lfloor n/2 \rfloor$ 个 swap 交换镜像位，把变换末端的位序翻转为 little endian。`inverse_qft` 用 `with b.adjoint()` 包裹对 `qft` 的调用，变换本身保持为模块调用、不在生成或序列化阶段展开，逆与正变换共享同一子模块定义。

适用边界：`width` 上限 64 由寄存器位宽契约给出；输出位序是 little endian，与外部位序对接时需要注意末尾 swap 已经完成。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../development/validation-plan.md` §2）：作用酉须与 Fourier 矩阵逐点相等。三层证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects` 钉死历史导入名 `pyqecclang.algorithms.elementary.qft` 指回本模块的规范 `qft` 对象（矩阵口径）。
- 数值：`AlgorithmExpansionTests.test_qft_matches_positive_fourier_matrix`——`qft(3)` 对全部 8 个基态列逐一模拟，每个输出幅度与正号 Fourier 矩阵元 $\exp(2\pi i \cdot x y / 8)/\sqrt{8}$ 对拍，精度 places = 10。
- 绑定：本算法无独立绑定见证（矩阵口径为 —；无开放声明入口）。

## 已知缺口与计划阶段

无已知缺口（验证矩阵缺口列为 —），阶段 V1。

## 相关链接

- 同模块：[Fourier 加法](fourier-addition.md)
- 源码：`src/pyqecclang/algorithms/fourier.py`
- API 参考：[Fourier 变换与算术](../../api/algorithms/fourier.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
