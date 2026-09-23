# 量子 Fourier 变换（Quantum Fourier Transform）

> 类别 C1 · 模块 `pyqecclang.algorithms.common.fourier` · 阶段 V1

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
- 源码：`src/pyqecclang/algorithms/common/fourier.py`
- API 参考：[Fourier 变换与算术](../../api/algorithms/common/fourier.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_fourier.py`（真实后端执行，无模拟替身）。实验设计：幺正级用 OriginIR-ext 导出经 UniQC `Circuit.to_matrix` 取全幺正（n ≤ 4，无工作区故 `effective_block` 即全矩阵、泄漏为 0），与独立构造的正号 DFT 矩阵 $F_{yx} = \exp(2\pi i x y/2^n)/\sqrt{2^n}$ 逐元素对比；寄存器级用 rir-pysparq 在 n = 6/8/12 上做两类逐点对拍——基态 $|x\rangle$ 经 QFT 的全部 $2^n$ 个振幅对照 DFT 行，以及仅用 H 与单比特相位门独立制备的 Fourier 模式 $\sum_x \omega^{-kx}|x\rangle/\sqrt{2^n}$ 经 QFT 后对基态索引 $k$ 的聚焦；n = 4 富相位态另做 reference / rir-pysparq / adapter-pysparq / originir-ext 四路径交叉对拍（含 `qft_with_work` 适配等价性）；逆 QFT 往返恒等在幺正级（n = 3/4 复合幺正 = $I$）与寄存器级（n = 8/12 全叠加恢复均匀态、采样基态往返）两级验证。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `qft-unitary-dft` | n = 1/2/3/4 | originir-ext + UniQC to_matrix | max_error | 8.7e-17 / 2.7e-16 / 1.4e-15 / 3.8e-15 |
| `inverse-qft-unitary` | n = 1/2/3/4 | originir-ext + UniQC to_matrix | max_error | 8.7e-17 / 2.7e-16 / 1.4e-15 / 3.8e-15 |
| `inverse-qft-roundtrip-unitary` | n = 3/4 | originir-ext + UniQC to_matrix | max_error | 4.4e-16 / 7.8e-16 |
| `qft-basis-row-pointwise` | n = 6/8/12 | rir-pysparq | max_error | 4.2e-15 / 9.5e-15 / 5.4e-14 |
| `qft-fourier-mode-focus` | n = 6/8/12 | rir-pysparq | min_success_probability | ≥ 1 − 4.4e-15（泄漏振幅 0） |
| `qft-cross-path` | n = 4 | 四路径对拍 | max_pairwise_deviation | 0.0 |
| `inverse-qft-roundtrip-wide` | n = 8/12 | rir-pysparq | uniform_max_error | 1.2e-16 / 4.5e-17（基态往返失败 0） |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_fourier.py
```

产物：`out/verification/fourier.json`（39 个案例全过，总运行约 80 秒）。
