# 量子相位估计（Quantum Phase Estimation）

> 类别 C3 · 模块 `pyqecclang.algorithms.estimation` · 阶段 V2

## 概述

给定酉算子 $U$ 的本征态 $|\psi\rangle$，$U|\psi\rangle = e^{2\pi i\varphi}|\psi\rangle$，量子相位估计（QPE）用 $p$ 位相位寄存器估计 $\varphi\in[0,1)$：读出值近似 $2^p\varphi$，即把相位编到 $2^p$ 等分的整数栅格上。电路为标准教科书构造——Hadamard 层、受控幂 $U^{2^k}$、逆 QFT。

QPE 是估计模块的读出核：数论模块的 `order_finding` 在其上组装求阶，本模块的 `amplitude_estimation` 则对 Grover 迭代做 QPE 完成振幅估计。

## 接口与输入模型

```python
phase_estimation(operation, *, precision=2)
```

- `operation`：支持受控调用的完整 `Operation`（input model 为 UO，酉算子 oracle）；输入态由调用方准备，QPE 自身不做制备。
- `precision`：相位寄存器位数，范围 1..63。

返回 `Operation`：保留输入的全部公开寄存器并追加 `phase`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"qpe"` |

生成期抛 `ValidationError` 的情形：精度越界、被调接口占用了 `phase` 参数名、输入不支持所需受控调用。

读出示例（`applications/gallery.py` 的 `phase_estimation` 条目）：对 $|1\rangle$ 上的相位门 $e^{i\pi/2}$（$\varphi=1/4$）取 `precision=3`，`phase` 读出为 2。

## 实现要点

寄存器布局：输入公开寄存器的名称、位宽和视图原样保留到后端降低阶段，`phase`（`Bits(precision)`）追加在外。各次幂以 `repeat(1 << bit)` 保存——Repeat 节点按幂次计数，生成与 JSON 序列化阶段不展开门序列。逆变换用 `qft(precision)` 的伴随调用实现（正号 Fourier 约定，见 `fourier.py`）。

读出与统计在宿主侧完成（模块 docstring 口径）：$\varphi$ 恰落在栅格点上时相位分布峰唯一；否则质量分布于最近栅格点的邻域，有限精度读出可能对应多个近似值。

## 验证方案

类别 C3（判定准则见 `docs/development/validation-plan.md` §2：输出分布等于闭式期望）。三层证据与 `docs/development/validation-coverage.md` 的 `estimation.py` 行一致：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言。
- 数值：QPE 无独立数值见证条目，经由两个下游消费方间接覆盖——`AlgorithmExpansionTests.test_amplitude_estimation_half_probability`（QAE 内部对 Grover 迭代调用 `phase_estimation`，相位分布 places=10 对拍）与 `number_theory` 行的 `AlgorithmExpansionTests.test_modular_multiplication_total_permutation_and_order`（`order_finding` 组装 QPE，相位分布 $p[0]=p[2]=0.5$）。
- 绑定：模块行登记为"—"（输入已要求具体 `Operation`，无开放槽位）。

见证技术为精确态矢量模拟取相位边际分布后与闭式对拍，不做采样断言。

## 已知缺口与计划阶段

本算法无单独缺口。`estimation.py` 在 validation-coverage 中登记的唯一缺口属于 QAE 的置信区间声明（见[振幅估计](qae.md)），模块整体因此处于阶段 V2。

## 相关链接

- 源码：`src/pyqecclang/algorithms/estimation.py`
- API 参考：[相位、振幅与重叠估计](../../api/algorithms/estimation.rst)
- 同组页面：[振幅估计](qae.md)、[量子计数](quantum-counting.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

`tests/verification/verify_estimation.py` 在真实后端上对本接口做分布级数值验证（QPE 部分共 13 个案例，全部通过）。实验设计：

- 已知相位幺正：单比特对角相位门 $U|1\rangle=e^{2\pi i\varphi}|1\rangle$，栅格上相位（$\varphi=0.5$ 与 $0.625$，$p=2..6$）与栅格外相位（$\varphi=0.3$，$p=3..6$）；相位寄存器直方图对照 Dirichlet 核闭式 $D(\delta)=\sin^2(\pi 2^p\delta)/(4^p\sin^2(\pi\delta))$。
- 非对角幺正：2 比特 `add_const(1)` 以其 Fourier 本征态 $|\tilde\varphi_j\rangle$（$j=1,3$）为输入；期望本征相位 $(-j/4)\bmod 1$ 由 numpy 独立构造置换矩阵与 Fourier 态求本征值得到，不复用库内实现。
- 幺正层面：$p=3$ 线路经 OriginIR-ext → UniQC `Circuit.to_matrix` 取全幺正，与 numpy 独立组装（H 层、受控 $U^z$、逆 DFT）逐元素对比。
- 后端路径：reference（内置参考执行器）、rir-pysparq（PySparQ 原生 RIR 解释器）、adapter-pysparq（PySparQ 适配器）、originir-ext（UniQC 全振幅态向量）。

| 案例 | 规模 | 路径 | 指标值 |
|---|---|---|---|
| 栅格相位读出 | $p=2..6$，$\varphi=0.5/0.625$ | 全部四条 | 峰值 $=2^p\varphi$，峰概率 $=1$，max TVD $\le 4.4\times10^{-16}$ |
| 栅格外相位读出 | $p=3..6$，$\varphi=0.3$ | 全部四条 | 峰值 $=\mathrm{round}(2^p\varphi)$，峰概率 $0.573/0.876$（$\ge 4/\pi^2$），max TVD $\le 7.3\times10^{-16}$ |
| Fourier 本征态 | $j=1,3$，$p=3,4$ | 全部四条 | 峰值 $=2^p\varphi$（$\varphi=0.75/0.25$），峰概率 $=1$，max TVD $\le 7.8\times10^{-16}$ |
| 全幺正对比 | $p=3$，$16\times16$ | originir-ext+to_matrix | max_error $=1.0\times10^{-15}$ |

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_estimation.py
```

产物路径：`out/verification/estimation.json`（案例名前缀 `qpe-`）。
