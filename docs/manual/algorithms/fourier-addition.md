# Fourier 加法（Fourier Addition）

> 类别 C1 · 模块 `pyqecclang.algorithms.fourier` · 阶段 V1

## 概述

在 [QFT](qft.md) 域中实现的模 $2^n$ 加法器：$|a, b\rangle \mapsto |a,\ (a+b) \bmod 2^n\rangle$。实现依据 Draper 的 QFT 加法构造（[arXiv:quant-ph/0008033](https://arxiv.org/abs/quant-ph/0008033)，见算法目录"实现依据"节）。整数加法在 Fourier 域中退化为单比特相位旋转：对 $b$ 做 QFT 后，每个 Fourier 基分量携带的相位线性于被加数，因此只需以 $a$ 的各位为控制、对 $b$ 的 Fourier 位施加固定角度相位门，再逆变换回计算基，全程不需要进位寄存器。

## 接口与输入模型

```python
fourier_add(width)
```

- `width`：`a` 与 `b` 的共同位宽，范围 1..64，非正整数或超界在生成期抛出 `ValidationError`。

返回 `Operation`，寄存器为 `a: Bits(width)` 与 `b: Bits(width)`：`a` 作为控制保留原值，`b` 更新为 $(a+b) \bmod 2^{\text{width}}$，无额外公开工作区。本入口不消费 oracle 输入（无 input model）。

## 实现要点

生成序列为：对 `b` 调用 `qft(width)` → 双重循环加受控相位：以 `a[i]` 控制、对 `b[j]` 施加角度 $2\pi / 2^{\text{width} - i - j}$ 的相位门（$0 \le i, j$ 且 $i + j < \text{width}$）→ 对 `b` 调用 `inverse_qft(width)`。QFT 与逆 QFT 均保留为模块调用，不在生成或序列化阶段展开。角度下标 `width - i - j` 的来源是 little endian 位权：`a` 的第 $i$ 位贡献 $a_i 2^i$，对 `b` 的第 $j$ 个 Fourier 位需要旋转 $2\pi a_i 2^i / 2^{\text{width} - j}$ 的整数倍相位。

适用边界：结果是模 $2^{\text{width}}$ 回绕的加法，不提供溢出标志；需要带符号或检测溢出的算术时应使用[可逆定点算术](fixed-point-arithmetic.md)的 `add`（含 `status` 标志）。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../development/validation-plan.md` §2）：作用酉须与真值表逐点相等。三层证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests`（与 QFT 同组的构造与属性断言，矩阵口径）。
- 数值：`AlgorithmExpansionTests.test_fourier_add_all_basis_inputs`——`fourier_add(3)` 对全部 $8 \times 8 = 64$ 对基态输入 $(a, b)$ 逐一模拟，输出幅度集中在 $(a,\ (a+b) \bmod 8)$ 上且等于 1，精度 places = 10，即穷举小实例全输入对拍。
- 绑定：本算法无独立绑定见证（矩阵口径为 —；无开放声明入口）。

## 已知缺口与计划阶段

无已知缺口（验证矩阵缺口列为 —），阶段 V1。

## 相关链接

- 同模块：[量子 Fourier 变换](qft.md)
- 源码：`src/pyqecclang/algorithms/fourier.py`
- API 参考：[Fourier 变换与算术](../../api/algorithms/fourier.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
