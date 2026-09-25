# 定点搜索（Fixed-Point Search）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.common.qsvt`](../../api/algorithms/common/qsvt.rst) · 阶段 V2

## 概述

定点振幅放大（fixed-point amplitude amplification）：不需要预知初始振幅、也不会因迭代次数过多而过冲的搜索放大。构造依据 Yoder–Low–Chuang 的闭式补多项式（源码 docstring 引用）：记 $L$ 为 BE 调用次数（必须为奇数）、$c = T_{1/L}(1/\delta)$，补多项式 $Q(x) = \delta c\, R(c^2(1 - x^2))$（$R(u) = T_L(\sqrt{u})/\sqrt{u}$ 为解析多项式），$P$ 由 $1 - (1 - x^2)Q^2$ 的求根谱分解得到。相位序列实现的成功概率恰为

$$
P_S(x) = 1 - \delta^2\, T_L^2\!\bigl(c\sqrt{1 - x^2}\bigr),
$$

$|x| \ge \sqrt{1 - 1/c^2}$ 时 $P_S \ge 1 - \delta^2$，且阈值随 $L$ 单调下降趋于 0（不动点性质）。

## 接口与输入模型

```python
fixed_point_search_phases(delta, degree)
fixed_point_search(a, delta, degree)
```

API 入口：{obj}`fixed_point_search_phases <oracq.algorithms.common.qsvt.fixed_point_search_phases>`、{obj}`fixed_point_search <oracq.algorithms.common.qsvt.fixed_point_search>`

- {obj}`fixed_point_search_phases <oracq.algorithms.common.qsvt.fixed_point_search_phases>`：纯数值例程，无 input model。`delta` 为误差 $\delta \in (0, 1)$；`degree` 为 $L$，必须是正奇数且不超过 20。返回时间正序相位元组（长度 $L + 1$）。
- {obj}`fixed_point_search <oracq.algorithms.common.qsvt.fixed_point_search>`：`a` 为 {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`（input model 为 BE），`delta` / `degree` 同上。返回 `BlockEncoding`，零信号块为复多项式 $P(A/\alpha)$，成功概率 $|P(x)|^2$ 满足上述 YLC 保证。

模块属性（`fixed_point_search`）：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"fixed_point_search"` |
| `be_alpha` | 输出 BE 的归一化（本算法为 1.0） |
| `delta` / `qsp_degree` | 误差参数回显 / $L$ |
| `threshold` | 放大阈值 $\sqrt{1 - 1/c^2}$：$\|x\| \ge$ 该值时成功概率 $\ge 1 - \delta^2$ |

## 实现要点

$Q$ 由 Chebyshev 奇次系数加 $(1 - x^2)$ 幂解析展开；$P$ 的谱分解用 Durand–Kerner 求根、重根聚类与共轭-反号四元组因子，因子计数自检 $1 + 2\times(\text{因子数}) = L$；随后 layer stripping 恢复相位（管线同 [QSP 相位合成](qsp-phase-synthesis.md)），自检 $|P(x)|^2$ 与 $1 - (1 - x^2)Q(x)^2$ 在 512 点网格逐点一致（误差超过 1e-5 抛 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`）。

本算法不做实部提取：$P$ 本身即目标（$|P| \le 1$ 保证块编码有效）。度数上限 20（合成管线 40 上限的一半）。适用边界：与 `search.py` 中作用于 marked oracle 的 Grover / 振幅放大不同，本入口作用于块编码的奇异值结构，谱变量需落在 $[-1, 1]$ 内。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../../development/validation-plan.md` §2）。三层证据：

- 结构：`tests/core/test_qsvt.py:TransformWitnessTests.test_transform_input_validation`（$\delta \ge 1$、偶数度数等负例）。
- 数值：`tests/core/test_qsvt.py:TransformWitnessTests.test_fixed_point_search_phase_properties`——$\delta = 0.4$、$L \in \{3, 5, 7\}$：相位长度 $L + 1$；400 点扫描 $|P(x)|^2 \le 1 + 10^{-9}$ 且与 YLC 闭式 $1 - \delta^2 T_L^2(c\sqrt{1-x^2})$ 逐点一致（delta = 1e-5）；阈值以上 $P_S \ge 1 - \delta^2$；阈值随 $L$ 单调下降。`test_fixed_point_search_circuit_amplifies`——2×2 对角 BE（$\alpha = 0.6$，谱变量 1.0 与 $1/6$），$\delta = 0.4$、$L = 5$：`threshold` 属性与 $\sqrt{1 - 1/c^2}$ 对拍（places = 3）；电路在 $x = 1$ 列的零信号幅度平方 $\ge 1 - \delta^2 - 10^{-9}$。
- 绑定：本算法无独立绑定见证（输入已要求具体 BE）。

## 已知缺口与计划阶段

与验证覆盖矩阵 `qsvt.py` 行一致：收敛性扫描缺失（阈值—误差权衡随 $L$ 的批量曲线未自动化，待阶段 V2 的收敛性扫描框架接入）；病态负例（度数 ≳ 16 或精度过低）的参数化测试已通过 input_validation 局部覆盖，单独的"病态输入必须抛错"批量参数化待补。

## 相关链接

- 源码：`src/oracq/algorithms/common/qsvt.py`
- API 参考：[QSVT 标准变换](../../api/algorithms/common/qsvt.rst)
- 同族页面：[QSP 相位合成](qsp-phase-synthesis.md)、[QSVT 相位序列](qsvt-sequence.md)、[VTAA-CKS 变时求解器](vtaa-cks.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

实验设计：取 $\delta=0.3$、$L=5$（阈值 $\sqrt{1-1/c^2}\approx0.3582$），用对角块编码把标量 $x$ 编码进零信号块（{obj}`diagonal_block_encoding <oracq.algorithms.input_model.oracles.diagonal_block_encoding>`，两基态对角元均为 $x$），对 `fixed_point_search` 的输出 BE 在四条后端路径（reference、rir-pysparq、adapter-pysparq、originir-ext）上测零信号成功概率，对照 YLC 闭式 $P_S(x)=1-\delta^2 T_L^2(c\sqrt{1-x^2})$。两个谱点分别位于阈值两侧：$x=\cos(\pi/6)\approx0.8660$（阈值上，应满足 $P_S\ge1-\delta^2=0.91$）与 $x=0.3$（阈值下）。

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| fixed-point-search-above-threshold | x=0.8660, δ=0.3, L=5 | 4 路径 | 成功概率（闭式 0.991311）/ 误差 | 0.991311 / 8.9e-16 |
| fixed-point-search-below-threshold | x=0.3, δ=0.3, L=5 | 4 路径 | 成功概率（闭式 0.772035）/ 误差 | 0.772035 / 1.1e-16 |

阈值上实例满足 $1-\delta^2$ 保证（0.991311 ≥ 0.91），阈值下实例仍与闭式逐点一致，说明相位序列实现的成功概率曲线在全谱区间符合理论。

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_search_walks.py
```

产物：`out/verification/search_walks.json`（案例 `fixed-point-search-above-threshold`、`fixed-point-search-below-threshold`）。
