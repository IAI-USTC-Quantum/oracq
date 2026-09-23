# Gibbs 态制备（Gibbs State Preparation）

> 类别 C2 · 模块 `pyqecclang.algorithms.input_model.density` · 阶段 V1

## 概述

给定哈密顿量 $H$ 的块编码与逆温 $\beta \ge 0$，制备 Gibbs 态 $\rho = e^{-\beta H}/Z$（$Z = \operatorname{Tr} e^{-\beta H}$）的近似纯化。实现走 QSVT 纯化路线（Chowdhury–Somma 2017、van Apeldoorn–Gilyén 2019、Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)）：先在 system 与 environment 上制备最大混合态的纯化（$n$ 对 Bell 对，见[纯化访问](purification.md)），再对 system 作用零信号块正比于 $g(H/\alpha)$ 的 QSVT 块编码，其中谱变量 $x = \lambda/\alpha \in [-1,1]$ 上

$$
g(x) = \exp\!\bigl(-\tfrac{\beta\alpha}{2}(x+1)\bigr) \in (0, 1], \qquad c = \beta\alpha/2 .
$$

$g$ 按奇偶分解为 $e^{-c}\cosh(cx)$ 与 $-e^{-c}\sinh(cx)$ 两支，各以修正 Bessel 截断逼近；后置选择 signal == 0 后，system 的约化密度矩阵正比于 $g(H/\alpha)^2 = e^{-\beta H}$（至多相差归一化）。

## 接口与输入模型

```python
gibbs_purification(hamiltonian, beta, *, error=0.01)
```

- `hamiltonian`：`BlockEncoding`，约定 $H$ 的谱含于 $[-\alpha, \alpha]$（input model 为 BE + DM；$\alpha$ = `be_alpha`）。
- `beta`：逆温 $\beta \ge 0$；$\beta = 0$ 时退化为最大混合态纯化（`qsp_degree = 0`、`gibbs_scale = 1.0`）。
- `error`：多项式一致逼近误差，必须在 $(0, 1)$ 内；每支截断尾部按 $\le$ `error`/8 控制。

返回 `ApproximatePurification`（`oracle_kind = "approximate_purification"`，操作签名为 `("system", "environment", "signal")`）。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"gibbs_purification"` |
| `beta` / `error` | 调用参数回显 |
| `qsp_degree` | 偶/奇两支截断度数的较大者 |
| `gibbs_scale` | 零信号块归一化 $2s$（与配分函数无关，不影响约化态） |
| `success_condition` | `"signal == 0"` |
| `width` / `environment_width` / `signal_qubits` | 三个寄存器的宽度 |

经典参考 `gibbs_state(hamiltonian, beta)`（小矩阵 $e^{-\beta H}/\operatorname{Tr}$）与 `trace_distance` 供见证与下游复用。

## 实现要点

生成链为 Bell 对制备 → 偶/奇支截断 → 相位合成 → LCU 相加。偶/奇支均为凸函数，端点匹配的虚部补全恒满足单位圆盘约束（$f^2(x)$ 不超过端点连线），相位合成必然可行；每支经虚部补全合成相位后用 $(U_\Phi + U_{-\Phi})/2$ 提取实部，两支再经 `linear_combination` 相加。缩放 $s = 1.5\max(\lVert g_{\mathrm{even}}\rVert_\infty, \lVert g_{\mathrm{odd}}\rVert_\infty, 10^{-3})$，两支合计一致误差不超过 `error`/4。

寄存器布局：`system(n) | environment(n) | signal(gibbs_be.signal_qubits)`，其中 environment 承载 Bell 对的另一半。适用边界：$\beta\alpha$ 过大时多项式度数超过合成上限（40），生成期抛出 `ValidationError`（提示减小 $\beta$ 或先缩小 $H$ 的谱尺度），不静默降级；输入必须是块编码而非稀疏 oracle。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）：后置选择后的约化密度矩阵须在容差内逼近经典参考 Gibbs 态。三层证据：

- 结构：`tests/core/test_density.py:GibbsTests.test_beta_zero_is_maximally_mixed`（$\beta=0$ 退化路径的属性与 $I/2$ 对拍）与 `test_invalid_inputs_fail_at_generation`（$\beta < 0$、`error` 为 0 或 $\ge 1$、hamiltonian 非 BlockEncoding 等）。
- 数值：`GibbsTests.test_diagonal_hamiltonian_gibbs` 与 `test_non_diagonal_hamiltonian_gibbs`——对 2×2 对角 / 非对角哈密顿量（$\beta = 0.6 / 0.8$、`error = 0.05`）在参考模拟器上取 signal == 0 分支，偏迹重归一后与 `gibbs_state` 经典参考对拍，迹距离 < 0.1。
- 绑定：本算法无独立绑定见证（输入已要求具体 BE）。

收敛见证（V1 新增）：`test_error_convergence_decreases` 固定 $\beta = 0.8$ 扫 `error` ∈ {0.4, 0.2, 0.1}，断言各档迹距离单调不增且每档 ≤ `error`——实测口径 5.6e-4 / 8.7e-5 / 8.7e-5（`error` 0.2 与 0.1 落到同一截断度数，距离相同，故断言不增而非严格下降）；`test_error_bound_uniform_in_beta` 固定 `error = 0.1` 扫 $\beta$ ∈ {0.2, 0.5, 1.0}，各点迹距离 ≤ `error`，实测 2.7e-6 / 1.5e-5 / 1.9e-4，随 $c = \beta\alpha/2$ 增长但始终远小于 `error`。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（对角/非对角对拍 + error 三档单调不增扫描 + β 网格一致性）。

## 相关链接

- 源码：`src/pyqecclang/algorithms/input_model/density.py`
- 同模块算法：[纯化访问](purification.md)
- API 参考：[密度矩阵输入模型与 Gibbs 态](../../api/algorithms/input_model/density.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_hamiltonian.py`（真实后端执行，无模拟替身）。实验设计：非对角哈密顿量 $H = \begin{pmatrix} 0.5 & 0.2 \\ 0.2 & -0.3 \end{pmatrix}$（$\alpha = 0.7$，$\beta = 1.2$）与对角哈密顿量 $H = \operatorname{diag}(0.8, -0.4)$（$\beta = 0.6$），`error = 0.02`；另有 $\beta = 0$ 退化例。在 reference / rir-pysparq / originir-ext 路径上执行纯化程序，取 signal == 0 分支，用 numpy 独立实现偏迹并重归一，与 `scipy.linalg.expm` 计算的经典 Gibbs 态 $e^{-\beta H}/Z$ 比较迹距离（`numpy.linalg.eigvalsh`），并报告后选择成功概率。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `gibbs-trace-distance-nondiag` | 13 量子位，$\beta = 1.2$ | reference、rir-pysparq、originir-ext | 迹距离 | 5.10e-5 |
| 同上 | — | — | 成功概率 | 0.0953 |
| `gibbs-trace-distance-diag` | 对角 $H$，$\beta = 0.6$ | reference、rir-pysparq | 迹距离 | 1.28e-5 |
| 同上 | — | — | 成功概率 | 0.0992 |
| `gibbs-beta0-maximally-mixed` | $\beta = 0$ 退化路径 | reference、rir-pysparq、originir-ext | 迹距离 | 0.0（精确 $I/2$） |

三例迹距离均远小于 `3×error = 0.06` 的判据（实测最高 5.1e-5），与页面所载收敛见证（迹距离随 `error` 单调不增且远小于 `error`）一致。

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_hamiltonian.py
```

产物：`out/verification/hamiltonian.json`。
