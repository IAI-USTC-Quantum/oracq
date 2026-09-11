# QPCA 主成分分析（Quantum Principal Component Analysis）

> 类别 C3 · 模块 `pyqecclang.algorithms.qpca` · 阶段 V1

## 概述

读出密度矩阵 $\rho$ 的谱：对酉步 $e^{-i\rho\Delta t}$ 做相位估计，本征态 $|\lambda\rangle$ 对应的相位读出为

$$
\varphi = -\frac{\lambda\,\Delta t}{2\pi} \pmod 1,
$$

解码即得本征值 $\lambda$。酉步由 LMR 密度矩阵指数化提供（见[密度矩阵指数化](density-matrix-exponentiation.md)），一次调用共消耗 $2^{\text{precision}} - 1$ 份 $\rho$ 拷贝。实现依据 Lloyd–Mohseni–Rebentrost 2014（"Quantum principal component analysis", Nature Physics 10, 631）。

## 接口与输入模型

```python
qpca(preparation, *, precision, step_time, system=None, swap_width=None, name=None)
```

- `preparation`：$\rho$ 拷贝的态制备（`StatePreparation`），input model 为 SP + QRAM：纯态用 `gate_state_prep(amplitudes)`，混合态经 `density.gate_purification(rho)` 的 `PurificationAccess.as_state_preparation()` 适配；制备必须零 work。
- `precision`：相位寄存器位数，范围 1..6；拷贝数随精度指数增长（precision = 6 时已达 63 份，与 `density_matrix_exponentiation` 的 copies 上限 63 一致）。
- `step_time`：单步演化时间 $\Delta t$，必须为正；须满足 $\lambda\Delta t \ll 2\pi$ 以免读出混叠。
- `system`：可选的系统输入态制备（缺省 $|0\rangle$），宽度必须等于 swap_width；输入不同的本征态，读出对应不同的本征值峰。
- `swap_width`：参与交换的前缀位宽，语义同 `density_matrix_exponentiation`（纯化场景取 $\rho$ 的系统位宽）。

返回 `Operation`，寄存器为 `system`、`copies`、`phase`。读出 `phase` 后用 `eigenvalue_from_phase(value, precision, step_time)` 解码本征值：按 $\lambda \in [0, \pi/\Delta t)$ 分支取两补码相位，$\lambda\Delta t$ 超出该范围时混叠由调用方负责。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"qpca"` |
| `readout_register` / `decoder` | `"phase"` / `"eigenvalue_from_phase"` |
| `copies` | 总拷贝数 $2^{\text{precision}} - 1$ |
| `step_time` | 单步演化时间 $\Delta t$ |
| `reference` | `"Lloyd-Mohseni-Rebentrost 2014, Nature Physics 10, 631"` |

## 实现要点

生成链为：可选系统初态制备 → 全部拷贝从全零由制备 oracle 填充 → `phase` 制备均匀叠加 → 相位估计的受控幂展开：`phase` 第 $b$ 位控制 $2^b$ 次连续部分交换（拷贝按固定游标顺序逐份消费，恰好用完 $2^p - 1$ 份）→ 对 `phase` 施加 `inverse_qft(precision)`。寄存器布局 `system(swap_width) | copies((2^p-1)\cdot w) | phase(p)`，$w$ 为制备宽度。

部分交换的逐位精确分解（XX+YY+ZZ = 2·SWAP − I）与指数化原语共用，细节见[密度矩阵指数化](density-matrix-exponentiation.md)。

设计决策与适用边界：大 $\Delta t$ 的 LMR 一阶误差表现为峰展宽而不移峰位——峰位（众数）解码仍精确、分布的圆周均值仍在容差内（见验证方案），因此见证允许直接按峰位读出。解码分支限定 $\lambda\Delta t \in [0, \pi)$，更大的本征值需调用方预缩放 $\Delta t$。

## 验证方案

类别 C3（概率分布语义，判定准则见 `../development/validation-plan.md` §2）：读出分布的峰位须与 $\rho$ 的闭式谱对拍。三层证据：

- 结构：`tests/core/test_qpca.py:QpcaTests.test_invalid_inputs_fail_at_generation` 覆盖 precision 越界（0 / 7）、step_time = 0、swap_width 越界、system 宽度不符、解码值越界等全部入口违例。
- 数值：`test_pure_state_eigenvalues`——$\rho = |+\rangle\langle+|$（本征值 1、0）：系统输入 $|+\rangle$（拷贝态本身，部分交换作用在 SWAP 对称本征态上对任意 $\Delta t$ 精确）时取 precision = 3、$\Delta t = \pi/4$（$\varphi = 7/8$ 落在 QPE 栅格上），读出确定（places = 9）、解码 $\lambda = 1$（places = 12）；输入 $|-\rangle$ 时峰展宽，分布圆周均值仍以 $\lambda = 0$ 为中心（delta = 0.15）。`test_mixed_state_eigenvalue_via_purification`——$\rho = \mathrm{diag}(0.75, 0.25)$ 经纯化适配，系统输入 $|0\rangle$ 读出主本征值：峰位解码恰为 0.75（places = 12，峰高 > 0.4），圆周均值 delta = 0.15。
- 绑定：混合态的 DM 输入经 `gate_purification(...).as_state_preparation()` 适配为 SP 进入电路，由 `test_mixed_state_eigenvalue_via_purification` 覆盖。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（纯态/混合态特征值峰对拍 + 指数化侧 $\Delta t$ 一阶率与 copies 标度见证）。按 validation-plan §5 的 V3 计划，QPCA 待登记 catalog 案例并接入真实后端对拍。

## 相关链接

- 源码：`src/pyqecclang/algorithms/qpca.py`
- 同模块页面：[密度矩阵指数化](density-matrix-exponentiation.md)
- API 参考：[QPCA 量子主成分分析](../../api/algorithms/qpca.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
