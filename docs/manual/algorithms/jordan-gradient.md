# Jordan 量子梯度估计（Jordan Quantum Gradient Estimation）

> 类别 C3 · 模块 `pyqecclang.algorithms.gradient` · 阶段 V1

## 概述

估计 $d$ 维实函数 $f$ 在给定点处的梯度 $\nabla f$。实现依据 Jordan 2005（PRL 95, 050501），缩放约定与 Gilyén–Arunachalam–Wiebe 2019 的推广一致：相位 oracle 实现

$$
O\,|x\rangle = e^{2\pi i \cdot N \cdot f(x)}\,|x\rangle, \qquad N = 2^{m},
$$

其中每个坐标取 $m$ 位定点网格（整数值 $k$ 对应网格点 $x = k/N$）。$f$ 在网格上近似线性时，单次 oracle 调用把 $N \cdot \partial f/\partial x_i$ 写入第 $i$ 个坐标寄存器的 Fourier 基，逐坐标逆 QFT 后一次读出全部 $d$ 个分量；经典确定性评估同一梯度需要 $O(d)$ 次函数查询。

## 接口与输入模型

```python
gradient_estimation(oracle, *, dimension, grid_bits)
```

- `oracle`：`PhaseOracle` 或具 `target` 签名的相位 oracle 操作，input model 为 FO（相位 oracle）；`phase_scale` 属性必须等于 $2^{\text{grid\_bits}}$，不符时按 INPUT_PROMISE 违例抛错。
- `dimension`：网格维数 $d$，范围 1..16。
- `grid_bits`：每个坐标的位数 $m$，范围 1..32；总宽度 $d \cdot m$ 不超过 64。

oracle 的三个构造入口：`abstract_phase_oracle(name, width, *, phase_scale)`（开放声明，实现留待 `bind`）、`gate_phase_oracle(width, angles, *, phase_scale, name=None)`（显式相位表，长度须为 $2^{\text{width}}$）、`function_phase_oracle(source, *, dimension, grid_bits, fmt=None, scale=None, ...)`（由 mathfunc 算术编译 $f$ 后做相位踢回）。

返回 `Operation`，寄存器为 `target`（宽 $d \cdot m$）。读出后用 `gradient_from_readout(value, *, dimension, grid_bits)` 解码：第 $i$ 个分量把位段 $[i \cdot m, (i+1) \cdot m)$ 按二进制补码解释后除以 $2^m$。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"jordan_gradient"` |
| `readout_register` / `decoder` | `"target"` / `"gradient_from_readout"` |
| `dimension` / `grid_bits` | 维数 $d$ 与每坐标位数 $m$ |
| `oracle_queries` / `classical_queries` | `1` / `O(dimension)` |

## 实现要点

生成链为：全部坐标寄存器制备均匀叠加 → 单次调用相位 oracle → 对每个坐标位段施加 `inverse_qft(grid_bits)`。寄存器布局即单一的 `target`，第 $i$ 个坐标占位段 $[i \cdot m, (i+1) \cdot m)$；无 work 寄存器，`function_phase_oracle` 内部的算术工作位在相位踢回后逆调用复原。

设计决策：相位缩放约定（`phase_scale == 2**grid_bits`）在生成期核对而非隐式假设，oracle 宽度也必须等于 $d \cdot m$，两者不符均在生成期抛 `ValidationError`。`function_phase_oracle` 默认定点格式 `FixedFormat(grid_bits+12, grid_bits+8)`，要求带符号且 fraction ≥ grid_bits 以精确表示网格坐标。

适用边界：读出分量只在 mod 1 意义下可分辨（相位 $e^{2\pi i N a j}$ 对 $a$ 与 $a \pm 1$ 相同），见证因此取 $|\partial f/\partial x_i| < 1/2$ 的分量；大梯度需调用方预缩放 $f$。$f$ 偏离线性时读出峰围绕真值展宽，失败概率随网格细化按近似二次率衰减（见验证方案）。

## 验证方案

类别 C3（概率分布语义，判定准则见 `../development/validation-plan.md` §2）：输出分布须等于闭式期望。三层证据：

- 结构：`tests/core/test_gradient.py:GradientTests.test_invalid_inputs_fail_at_generation` 覆盖维度 / 位数越界、oracle 宽度不符、相位表长度不符、phase_scale 非正、读出值越界等全部入口违例。
- 数值：`test_linear_function_exact_two_dimensions`——线性函数（分量 $3/16$、$-2/16$，$m = 4$）读出分布确定性落在编码值上（概率 1.0，places = 12），解码精确等于梯度；`test_function_phase_oracle_from_mathfunc`——mathfunc 路径 $f(x) = 0.25x$ 精确恢复分量 $0.25$；`test_perturbed_linear_concentrates_with_grid_bits`——$f(x) = ax + x^2/N^2$（$a = 3/8$）在 $m = 3, 4, 5$ 下峰位恒为真值、正确概率单调上升，并断言失败概率衰减率 $q_{m+1} \le 0.34 \cdot q_m$（相位扰动 $o(1/N)$ 下峰外泄漏近似二次收敛；实测比率 0.292 / 0.268）。
- 绑定：`test_abstract_oracle_binds_to_gate_implementation`——abstract 声明经 `bind` 绑定 gate 相位表后读出解码与直接构造一致，覆盖声明/实现两层一致性。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（线性精确读出 + mathfunc 路径 + 失败概率衰减率断言 $q_{m+1} \le 0.34 \cdot q_m$，实测 0.292 / 0.268 + 绑定一致性）。

## 相关链接

- 源码：`src/pyqecclang/algorithms/gradient.py`
- API 参考：[量子梯度估计](../../api/algorithms/gradient.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

`tests/verification/verify_estimation.py` 在真实后端上对本接口做读出级数值验证（共 11 个案例，全部通过）。实验设计：

- 线性函数精确读出：`gate_phase_oracle` 相位表，$d=1$（$m=3,4,5$，梯度分量 $3/8$、$-5/16$、$9/32$）、$d=2$（$m=4$）、$d=3$（$m=3$）；读出分布应确定性落在编码值上。
- mathfunc 相位 oracle：`function_phase_oracle` 编译 $f(x)=0.25x$（$d=1,m=4$）与 $f(x_0,x_1)=0.25x_0-0.125x_1$（$d=2,m=3$），同时对照中心有限差分（线性函数时与真值完全一致，差距 $0$）。该路径工作区按 `workspace_table` 预算为 $1768$ / $1638$ 量子位，超出 OriginIR 24 位预算，只走 reference / rir-pysparq / adapter-pysparq 三条寄存器级路径。
- 扰动线性收敛：$f(x)=ax+x^2/N^2$（$a=3/8$），$m=3,4,5$；峰位恒为真值，失败概率衰减率对照近似二次收敛判据；信息性指标给出 $x=1/2$ 处中心有限差分（$f$ 非线性时与线性系数差 $O(1/N^2)$）。
- 后端路径：reference、rir-pysparq、adapter-pysparq、originir-ext（gate 相位表案例，总宽 $\le 9$ 量子位）。

| 案例 | 规模 | 路径 | 指标值 |
|---|---|---|---|
| 线性读出（gate 表） | $d=1$，$m=3,4,5$ | 全部四条 | 解码误差 $=0$，峰概率 $=1$ |
| 线性读出（gate 表） | $d=2$，$m=4$；$d=3$，$m=3$ | 全部四条 | 解码误差 $=0$，峰概率 $=1$ |
| mathfunc oracle | $d=1,m=4$；$d=2,m=3$ | ref/rir/adp | 解码误差 $=0$，峰概率 $=1$，与有限差分一致 |
| 扰动收敛 | $d=1$，$m=3,4,5$ | 全部四条 | 成功概率 $0.9588/0.9880/0.9968$，衰减率 $0.2922/0.2678$（$<0.34$）；有限差分 $0.3906/0.3789/0.3760\to 0.375$ |

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_estimation.py
```

产物路径：`out/verification/estimation.json`（案例名前缀 `jordan-`）。
