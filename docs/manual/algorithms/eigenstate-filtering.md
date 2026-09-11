# 特征态过滤（Eigenstate Filtering）

> 类别 C2 · 模块 `pyqecclang.algorithms.qsvt` · 阶段 V2

## 概述

特征态过滤：对矩阵 $A$ 的块编码构造"尖峰"多项式块编码——谱变量落在过滤中心附近通带内的分量被保留（$|f|$ 接近 1），通带外的分量被压制到 $1/T_d(r)$ 以下。目标为 Lin–Tong 型过滤多项式（源码 docstring 引用的构造依据）：

$$
f(x) = \frac{T_d\!\bigl(g(x^2)\bigr)}{T_d(r)}, \qquad
g(y) = \frac{2(y - \Delta^2)}{1 - \Delta^2} - 1, \qquad
r = \frac{1 + \Delta^2}{1 - \Delta^2},
$$

其中 $\Delta$ 为过滤宽度、$d$ 为 Chebyshev 度数。$f$ 在 $x = 0$ 处饱和（$|f(0)| = 1$），在 $|x| \ge \Delta$ 上 $|f(x)| \le 1/T_d(r)$，压制率随 $d \cdot \operatorname{arcosh}(r)$ 指数下降。相位合成框架同 [QSP 相位合成](qsp-phase-synthesis.md)（Gilyén et al. 2019, [arXiv:1806.01838](https://arxiv.org/abs/1806.01838)）。

## 接口与输入模型

```python
eigenstate_filter(a, gap, degree, *, center=0.0)
```

- `a`：`BlockEncoding`，被过滤矩阵的块编码（input model 为 BE；谱变量 $x = \lambda/\alpha$，厄米情形即归一化本征值）。
- `gap`：过滤宽度 $\Delta$，必须在 $(0, 1)$ 内。
- `degree`：Chebyshev 度数 $d$，必须是正整数（浮点被拒绝）。
- `center`：过滤中心，必须在 $(-1, 1)$ 内；非零时先平移谱（见下）。

返回 `BlockEncoding`，零信号块在通带外被压到 `suppression` 以下。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"eigenstate_filter"` |
| `be_alpha` | 输出 BE 的归一化（实部提取 LCU 各半，为 1.0） |
| `gap` / `filter_degree` / `center` | 调用参数回显 |
| `qsp_degree` | 合成的多项式度数 $2d$ |
| `suppression` | 通带外上界 $1/T_d(r)$ |

## 实现要点

合成度数为 $2d$（$f$ 是 $x^2$ 的 $d$ 次 Chebyshev 组合，单项式度数 $2d$），受合成上限 40 约束（即 $d \le 20$），超限抛 `ValidationError`。`center` 非零时先经 BE 线性组合 `linear_combination(1.0, a, -center, identity(a.width))` 平移谱：构造 $A - \text{center}\cdot I$ 的块编码（归一化变为 $\alpha + |\text{center}|$），过滤作用在平移后的谱变量 $x' = (\lambda - \text{center})/(\alpha + |\text{center}|)$ 上，`gap` 以该变量计量。

$f$ 在 $x = 0$ 饱和但 $|f(\pm 1)| = 1/T_d(r) < 1$ 端点未饱和，需借助虚部补全 $h = \text{sat}\cdot x^2$（$\text{sat} = \sqrt{1 - f(1)^2}$）才能合成相位；合成后经 $(U_\Phi + U_{-\Phi})/2$ 提取实部（机制见 [QSP 相位合成](qsp-phase-synthesis.md)）。适用边界：输入必须是块编码；平移后谱值应落在 $[-1, 1]$ 内，超出部分失去上界保证。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../../development/validation-plan.md` §2）：作用算子需在容差内逼近目标连续函数。三层证据：

- 结构：`tests/core/test_qsvt.py:TransformWitnessTests.test_transform_input_validation`（gap 越界、度数非正整数等负例）与 `PhaseSynthesisTests.test_input_validation`（底层合成的输入校验）。
- 数值：`tests/core/test_qsvt.py:TransformWitnessTests.test_eigenstate_filter_isolates_eigenvalue`——2×2 对角 BE（本征值 0.05 与 0.5，$\alpha = 0.5$），gap = 0.2、$d = 8$、center = 0.1：`suppression` < 0.08；通带本征值平移后 $\approx -0.083$（$|x'| < \Delta$），零信号幅度 > 0.7；阻带本征值平移后 $\approx 0.667$，幅度 < 0.1。
- 绑定：本算法无独立绑定见证（输入已要求具体 BE）。

## 已知缺口与计划阶段

与验证覆盖矩阵 `qsvt.py` 行一致：收敛性扫描缺失（压制率随度数的指数衰减曲线未自动化，待阶段 V2 的收敛性扫描框架接入）；病态负例（度数 ≳ 16 或精度过低）的参数化测试已通过 input_validation 局部覆盖，单独的"病态输入必须抛错"批量参数化待补。

## 相关链接

- 源码：`src/pyqecclang/algorithms/qsvt.py`
- API 参考：[QSVT 标准变换](../../api/algorithms/qsvt.rst)
- 同族页面：[QSP 相位合成](qsp-phase-synthesis.md)、[QSVT 矩阵求逆](qsvt-matrix-inversion.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
