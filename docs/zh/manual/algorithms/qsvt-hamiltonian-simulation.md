# QSVT 哈密顿模拟（QSVT Hamiltonian Simulation）

[English](../../../index.html) · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.common.qsvt`](../../api/algorithms/common/qsvt.rst) · 阶段 V2

## 概述

对厄米矩阵 $A$ 的块编码构造 $e^{itA/\alpha}$ 的块编码（$\alpha$ 为输入 BE 归一化）。实现依据 QSVT 框架（Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)，模块 docstring 引用）与 $e^{itx}$ 的 Jacobi–Anger 展开 $e^{itx} = \cos(tx) + i\sin(tx)$：偶支近似 $\cos(tx)$、奇支近似 $\sin(tx)$，两支分别经 QSP 相位合成（见 [QSP 相位合成](qsp-phase-synthesis.md)），各用 $(U_\Phi + U_{-\Phi})/2$ 提取实部，最后按系数 1 与 $i$ 做 LCU 组合。返回 BE 的零信号块约为 $e^{itA/\alpha}/\text{sim\_scale}$。

## 接口与输入模型

```python
qsvt_hamiltonian_simulation(a, t, *, error=0.01)
```

API 入口：{obj}`qsvt_hamiltonian_simulation <oracq.algorithms.common.qsvt.qsvt_hamiltonian_simulation>`

- `a`：{obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`，被模拟哈密顿量的块编码（input model 为 BE）。
- `t`：演化时间，必须是非零有限实数。
- `error`：多项式近似误差，必须在 $(0, 1)$ 内。

返回 `BlockEncoding`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"qsvt_hamiltonian_simulation"` |
| `be_alpha` | 输出 BE 的归一化（1 与 $i$ 两支 LCU，为 2.0） |
| `time` / `error` | 调用参数回显 |
| `qsp_degree` | Bessel 截断度数 $K$（两支实际合成度数为不超过 $K$ 的偶 / 奇度数） |
| `sim_scale` | 输出缩放 $2s$（零信号块约为 $e^{itA/\alpha}/\text{sim\_scale}$） |

## 实现要点

截断度数自适应：$J_k(t)$ 用幂级数（纯 Python）计算，$K$ 取 Bessel 尾项绝对值之和不超过 `error`/4 的最小度数，上界封顶 40。两支统一缩放 $s = 1.5\max(\lVert f_c\rVert_\infty, \lVert f_s\rVert_\infty, 10^{-3})$，为虚部补全留出余量，$\text{sim\_scale} = 2s$。奇偶性校验（cos 支为偶函数、sin 支为奇函数）违例抛 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`；虚部补全候选族（cos 支为常数加单偶次幂、sin 支为单奇次幂）逐个尝试，取首个可合成的候选。

适用边界：$t = 0$ 被拒绝（平凡情形）；$|t|$ 很大时所需截断度数会先于误差条件触及 40 上限，此时 `error` 界可能不被满足而不报错，调用方需自行核对。输出以 `sim_scale` 缩放，读出后需除回。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../../development/validation-plan.md` §2）：作用算子需在容差内逼近目标连续函数。三层证据：

- 结构：`tests/core/test_qsvt.py:TransformWitnessTests.test_transform_input_validation`（t = 0、error 越界等负例）与 `PhaseSynthesisTests.test_input_validation`（底层合成的输入校验）。
- 数值：`tests/core/test_qsvt.py:TransformWitnessTests.test_hamiltonian_simulation_block`——2×2 对角矩阵（本征值 0.5、$-0.25$，$\alpha = 0.5$，谱变量 1.0 与 $-0.5$），t = 0.7、error = 0.01，逐列读出零信号块幅度与 $e^{itx}/\text{sim\_scale}$ 对拍，delta = 5e-4。
- 绑定：本算法无独立绑定见证（输入已要求具体 BE）。

## 已知缺口与计划阶段

与验证覆盖矩阵 `qsvt.py` 行一致：收敛性扫描缺失（误差随 t / 度数的批量曲线未自动化，待阶段 V2 的收敛性扫描框架接入）；病态负例（度数 ≳ 16 或精度过低）的参数化测试已通过 input_validation 局部覆盖，单独的"病态输入必须抛错"批量参数化待补。

## 相关链接

- 源码：`src/oracq/algorithms/common/qsvt.py`
- API 参考：[QSVT 标准变换](../../api/algorithms/common/qsvt.rst)
- 同族页面：[QSP 相位合成](qsp-phase-synthesis.md)、[QSVT 矩阵求逆](qsvt-matrix-inversion.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_hamiltonian.py`（真实后端执行，无模拟替身）。实验设计：输入 BE 为 {obj}`matrix_pauli_encoding <oracq.algorithms.input_model.block_encoding.matrix_pauli_encoding>` 编码的非对角厄米矩阵 $A = \begin{pmatrix} 0.5 & 0.2 \\ 0.2 & -0.3 \end{pmatrix}$（$\alpha = 0.7$），$t \in \{0.7, 2.0\}$、`error = 0.01`；在 reference / rir-pysparq / originir-ext 三条路径上逐列读出完整 $2\times2$ 零信号块，与 `scipy.linalg.expm` 计算的 $e^{itA/\alpha}/\text{sim\_scale}$ 对拍。实现误差（电路块 vs $e^{itx}$ 参考）与方法误差（Jacobi–Anger 截断尾部上界 $2\sum_{j>K}|J_j(t)|$，由 `scipy.special.jv` 独立求值）分列报告。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `qsvt-hamsim-t0.7` | $K = 4$，sim_scale ≈ 3.0 | reference、rir-pysparq、originir-ext | impl_error | 2.72e-5 |
| 同上 | — | — | 方法尾部上界 / 总误差界 | 9.10e-5 / 1.18e-4 |
| `qsvt-hamsim-t2.0` | $K = 6$ | reference、rir-pysparq、originir-ext | impl_error | 5.42e-5 |
| 同上 | — | — | 方法尾部上界 / 总误差界 | 4.00e-4 / 4.54e-4 |

两个时间点上的总误差界（1.2e-4 与 4.5e-4）均远小于请求的 `error = 0.01`，且实现误差显著低于方法误差，表明截断度数选择是精度的主导项、相位合成与电路组装处于噪声量级。

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_hamiltonian.py
```

产物：`out/verification/hamiltonian.json`。
