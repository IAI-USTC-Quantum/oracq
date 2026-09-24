# QSVT 矩阵求逆（QSVT Matrix Inversion）

> 类别 C2 · 模块 [`oracq.algorithms.common.qsvt`](../../api/algorithms/common/qsvt.rst) · 阶段 V2

## 概述

给定厄米矩阵 $A$ 的块编码（block encoding），构造近似 $A^{-1}$ 的块编码。实现依据 Gilyén 等人的量子奇异值变换框架（Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)）：对奇异值 $x$ 施加多项式变换，目标多项式取奇扩展求逆多项式

$$
J_b(x) = \frac{1-(1-x^2)^b}{x}, \qquad f(x) = c \cdot J_b(x),
$$

其中奇次幂系数 $(-1)^m \binom{b}{m+1}$ 解析给出，缩放常数 $c$ 使 $\lVert f \rVert_\infty \le 1/3$。在 $|x| \ge 1/\kappa$ 上，$f(x)$ 以不超过 `error` 的相对误差逼近 $c/x$；截断参数 $b$ 由条件数 $\kappa$ 与 `error` 决定，多项式度数 $d = 2b - 1$。

## 接口与输入模型

```python
qsvt_matrix_inversion(a, kappa, *, error=0.05)
```

API 入口：{obj}`qsvt_matrix_inversion <oracq.algorithms.common.qsvt.qsvt_matrix_inversion>`

- `a`：{obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`，被求逆矩阵的块编码（input model 为 BE）。
- `kappa`：条件数 $\kappa$，必须是不小于 1 的有限数。
- `error`：相对近似误差，必须在 $(0, 1)$ 内。

返回 `BlockEncoding`，其零信号块约为 `inverse_scale · A⁻¹`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"qsvt_matrix_inversion"` |
| `be_alpha` | 输出 BE 的归一化（本算法为 1.0） |
| `kappa` / `error` | 调用参数回显 |
| `qsp_degree` | 合成的多项式度数 $d = 2b - 1$ |
| `inverse_scale` | 输出缩放 $c \cdot \alpha$（$\alpha$ 为输入 BE 归一化） |

## 实现要点

度数由 $\kappa$ 与 `error` 闭式推出：$b = \max\!\bigl(1,\ \lceil \log(1/\varepsilon) / -\log(1 - 1/\kappa^2) \rceil\bigr)$，$\kappa = 1$ 时 $b = 1$；所需度数超过合成上限 40 时抛出 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`，调用方需放宽 `error`。

目标多项式 $f$ 端点未饱和（$|f(\pm 1)| < 1$），必须借助非零虚部补全 $h$ 才能合成相位（见模块 docstring 的可实现条件）。合成后用 $(U_\Phi + U_{-\Phi})/2$ 的 LCU 组合提取实部——$-Φ$ 恰好实现共轭多项式 $\bar{P}$——所得块编码即 $f(A/\alpha)$。相位合成采用补多项式求根加逐层剥离（layer stripping），每次合成后经 {obj}`qsp_response <oracq.algorithms.common.qsvt.qsp_response>` 往返自检，病态输入直接拒绝。

适用边界：输入必须是块编码而非稀疏 oracle 或 QRAM；矩阵未先块编码时需先经 `block_encoding.py` / `sparse.py` / `lowrank.py` 等适配。合成度数受 40 上限约束，$\kappa$ 大且 `error` 小的组合会在生成期报错而非静默降级。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）：作用算子需在容差内逼近目标连续函数。三层证据：

- 结构：`tests/core/test_qsvt.py:PhaseSynthesisTests.test_input_validation`（奇偶性违例、上界违例、端点未饱和、非实系数、虚部奇偶）与 `TransformWitnessTests.test_transform_input_validation`（κ < 1、error 越界、所需度数超限等）。
- 数值：`TransformWitnessTests.test_matrix_inversion_block`——对 2×2 矩阵（本征值 0.4、0.8，κ = 2，error = 0.15）在参考模拟器上逐列读出零信号块幅度，与 `inverse_scale` 乘精确逆 $A^{-1} = \begin{pmatrix} 1.875 & 0.625 \\ 0.625 & 1.875 \end{pmatrix}$ 对拍，相对容差为期望值的 35%。
- 绑定：本算法无独立绑定见证（输入已要求具体 BE）。

底层序列约定由 `PhaseSynthesisTests.test_convention_matches_qsvt_sequence` 钉死（随机相位下电路零信号块与 `qsp_response` 逐点一致，delta = 1e-10），实部提取所依赖的共轭性质由 `test_negated_phases_conjugate_polynomial` 见证（delta = 1e-12）。

## 已知缺口与计划阶段

收敛性扫描缺失：误差随度数 / $\kappa$ 的批量曲线未自动化（阶段 V2 的收敛性扫描框架落地后接入）。病态负例（度数 ≳ 16 或精度过低）已由 `test_input_validation` / `test_transform_input_validation` 局部覆盖，单独的"病态输入必须抛错"批量参数化测试待补。

## 相关链接

- 源码：`src/oracq/algorithms/common/qsvt.py`
- 同族页面：[特征态过滤](eigenstate-filtering.md)、[QSP 相位合成](qsp-phase-synthesis.md)、[QSVT 哈密顿模拟](qsvt-hamiltonian-simulation.md)、[块编码组合代数](block-encoding-algebra.md)（BE 组合入口）
- API 参考：[QSVT 标准变换](../../api/algorithms/common/qsvt.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_hamiltonian.py`（真实后端执行，无模拟替身）。实验设计：三组矩阵——2×2（本征值 0.4/0.8，$\kappa = 2$，`error` ∈ {0.15, 0.10}，对应度数 $d = 13, 17$）、2×2（0.2/0.6，$\kappa = 3$，`error = 0.4`，$d = 15$）、4×4 对角（0.3/0.4/0.6/0.9，$\kappa = 3$，`error = 0.4`，$d = 15$）；在 reference / rir-pysparq / originir-ext 路径上读出完整零信号块并除以 `inverse_scale`，与 `numpy.linalg.inv` 的精确逆对拍，报告相对谱误差、恢复逆的条件数与理论区间 $\kappa\cdot(1\pm e)/(1\mp e)$（每个奇异值带 ≤ `error` 的相对偏差），以及各列成功概率。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `qsvt-inversion-kappa2-error0.15` | 2×2，$d = 13$ | reference、rir-pysparq | 相对谱误差 | 0.1335（≤ 0.15） |
| 同上 | — | — | $\kappa$ 恢复 vs 精确（理论区间） | 1.733 vs 2.0（[1.478, 2.706]） |
| `qsvt-inversion-kappa2-error0.1` | 2×2，$d = 17$ | reference、rir-pysparq | 相对谱误差 | 0.0751（≤ 0.10） |
| 同上 | — | — | $\kappa$ 恢复 vs 精确 | 1.850 vs 2.0（[1.636, 2.444]） |
| `qsvt-inversion-kappa3-1q` | 2×2，$d = 15$ | reference、rir-pysparq、originir-ext | 相对谱误差 | 0.3897（≤ 0.40） |
| 同上 | — | — | $\kappa$ 恢复 vs 精确 | 1.831 vs 3.0（[1.286, 7.0]） |
| `qsvt-inversion-kappa3-2q` | 4×4，$d = 15$ | reference、rir-pysparq、originir-ext | 相对谱误差 | 0.3897（≤ 0.40） |
| 同上 | — | — | $\kappa$ 恢复 vs 精确；各列成功概率 | 1.863 vs 3.0；0.031–0.109 |

相对谱误差均落在请求 `error` 之内，恢复条件数均落在理论区间之内，且随 `error` 收紧而改善（0.1335 → 0.0751）。数值边界实测：求逆多项式在 $d \ge 19$（如 $\kappa = 3$、`error` ≤ 0.3）被相位合成自检拒绝（`ValidationError`，与"病态输入拒绝而非静默降级"的文档承诺一致）；$d \le 17$ 的可合成区域全部通过。

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_hamiltonian.py
```

产物：`out/verification/hamiltonian.json`。
