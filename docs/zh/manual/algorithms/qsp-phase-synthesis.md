# QSP 相位合成（QSP Phase Synthesis）

[English](../../../index.html) · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.common.qsvt`](../../api/algorithms/common/qsvt.rst) · 阶段 V2

## 概述

量子信号处理（QSP）相位合成：给定实系数目标多项式，求出使 [QSVT 相位序列](qsvt-sequence.md) 实现该多项式的相位序列 $\Phi$。模块导出一对互逆例程：{obj}`qsp_response(x, phases) <oracq.algorithms.common.qsvt.qsp_response>` 正向求值，返回相位序列在反射约定下的顶层左块 $p(x)$（$x \in [-1, 1]$）；{obj}`qsp_phases(coeffs, imag=None) <oracq.algorithms.common.qsvt.qsp_phases>` 逆向合成，返回时间正序相位（长度 $d + 1$，$d$ 为目标度数）。

合成采用 Gilyén et al. 2019（[arXiv:1806.01838](https://arxiv.org/abs/1806.01838)，模块 docstring 引用）的补多项式求根加逐层剥离（layer stripping）。端点未饱和的实目标（如 $1/x$ 的截断近似）必须借助非零虚部补全 $h$：合成 $P = f + i h$ 后，再用 $(U_\Phi + U_{-\Phi})/2$ 的 LCU 组合提取实部——$-\Phi$ 恰好实现共轭多项式 $\bar{P}$，块编码即为 $f(A/\alpha)$。

## 接口与输入模型

```python
qsp_response(x, phases)
qsp_phases(coeffs, imag=None)
```

API 入口：{obj}`qsp_response <oracq.algorithms.common.qsvt.qsp_response>`、{obj}`qsp_phases <oracq.algorithms.common.qsvt.qsp_phases>`

- `x`：标量谱变量（$x \in [-1, 1]$）；`phases`：相位序列。返回复数 $p(x)$。
- `coeffs`：目标多项式 $f$ 的升幂实系数（常数项在前）；`imag`：可选虚部补全 $h$，同为升幂实系数。

两个入口都是纯数值例程，不生成量子程序，因此没有 input model；产物（相位元组）供 {obj}`qsvt_sequence <oracq.algorithms.common.transforms.qsvt_sequence>` 与标准变换族消费。可实现条件（违反时抛 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`）：$f$ 与 $h$ 的奇偶性均为 $d \bmod 2$；$[-1,1]$ 上 $f^2 + h^2 \le 1$；端点饱和 $f(\pm 1)^2 + h(\pm 1)^2 = 1$；补多项式 $R = (1 - f^2 - h^2)/(1 - x^2)$ 非负且满足谱分解的根重数条件。$d = 0$ 的目标必须是单位模常数，此时返回单相位。

## 实现要点

合成管线：输入校验（实系数、奇偶性、上界、端点饱和、$R$ 非负）→ $R$ 求根（纯 Python Durand–Kerner 迭代，按相对距离聚类重根）→ 根结构分类（实根成对、纯虚根成对、复根共轭-反号四元组）构造 $Q$ 使 $Q\bar{Q} = R$ → layer stripping 逆向恢复相位（双精度复数）→ 用 `qsp_response` 在 512 点网格上往返自检，误差超过 1e-5 抛出 `ValidationError`。

`qsp_response` 的求值是 2×2 矩阵递推：从恒等阵出发，每步先左乘 $W(x)$（时间上先于本步相位的作用，首相位除外）、再左乘 $S(\varphi)$，返回累积矩阵的 $(0,0)$ 元。除独立调用外，它还是合成器内部的自检工具：`qsp_phases` 与 {obj}`fixed_point_search_phases <oracq.algorithms.common.qsvt.fixed_point_search_phases>` 在剥离结束后都用它做网格往返校验。

数值边界：合成度数上限 40；度数几十以内、补多项式根分离良好时往返误差典型在 1e-9 量级。度数较高（≳ 16）或补多项式近重根时求根与剥离数值不稳定，自检会拒绝病态输入而非输出错误相位。适用边界：只接受实系数目标，复目标需拆分为实部加虚部补全；返回相位的时间顺序约定与 `qsvt_sequence` 电路约定逐步对齐（见 [QSVT 相位序列](qsvt-sequence.md)）。标准变换族（[矩阵求逆](qsvt-matrix-inversion.md)、[特征态过滤](eigenstate-filtering.md)、[哈密顿模拟](qsvt-hamiltonian-simulation.md)、[定点搜索](fixed-point-search.md)）的相位全部经本管线合成。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../../development/validation-plan.md` §2）。三层证据：

- 结构：`tests/core/test_qsvt.py:PhaseSynthesisTests.test_input_validation`——奇偶性违例、上界违例、端点未饱和、非实系数、虚部奇偶性五类负例全部在合成期抛 `ValidationError`。
- 数值：`test_roundtrip_chebyshev`（Chebyshev 多项式 $T_n$，$n = 1..6$，51 点网格往返一致 delta = 1e-8，相位长度 $n + 1$）；`test_roundtrip_with_explicit_imaginary_part`（$P = (0.5 + i\sqrt{0.75})x$，$|P|^2 = x^2$、$Q = 1$ 的精确可实现例，delta = 1e-9）；`test_negated_phases_conjugate_polynomial`（$-\Phi$ 实现 $\bar{P}$，实部提取所依赖，delta = 1e-12）；`test_convention_matches_qsvt_sequence`（`qsp_response` 与电路零信号块逐点一致 delta = 1e-10）。
- 绑定：纯数值例程，无绑定见证。

## 已知缺口与计划阶段

与验证覆盖矩阵 `qsvt.py` 行一致：收敛性扫描缺失（误差随度数的批量曲线未自动化，待阶段 V2 的收敛性扫描框架接入）；病态负例（度数 ≳ 16 或精度过低）的参数化测试已通过 `test_input_validation` 局部覆盖，单独的"病态输入必须抛错"批量参数化待补。

## 相关链接

- 源码：`src/oracq/algorithms/common/qsvt.py`
- API 参考：[QSVT 标准变换](../../api/algorithms/common/qsvt.rst)
- 同族页面：[QSVT 相位序列](qsvt-sequence.md)（相位消费端）、[QSVT 矩阵求逆](qsvt-matrix-inversion.md)（虚部补全的典型用例）
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_hamiltonian.py`（真实后端执行，无模拟替身）。实验设计：目标多项式取 Chebyshev $T_1,\dots,T_6$ 与精确可实现例 $P = (0.5 + i\sqrt{0.75})x$，由 `qsp_phases` 合成相位后在 401 点均匀网格上用 numpy 按本文档 $W(x)$ / $S(\varphi)$ 约定**独立实现**的 2×2 矩阵递推响应求值（不调用库内 `qsp_response`，两者互为独立实现），与目标多项式的独立 Horner 求值逐点对拍。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `qsp-phase-synthesis-roundtrip` | 7 个目标，401 点网格 | numpy 独立响应（合成器不生成量子程序） | max_error | 3.8e-15 |
| 同上（Chebyshev $T_1$–$T_6$） | — | — | per-target 误差 | 0 – 3.8e-15 |
| 同上（显式虚部目标） | — | — | 误差 | 2.2e-16 |

合成相位在电路层面的消费端对拍（四后端路径、随机相位、矩阵多项式块）见 [QSVT 相位序列](qsvt-sequence.md) 的数值验证节。

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_hamiltonian.py
```

产物：`out/verification/hamiltonian.json`。
