# 密度矩阵指数化（Density Matrix Exponentiation）

> 类别 C3 · 模块 `pyqecclang.algorithms.qml.qpca` · 阶段 V1

## 概述

给定密度矩阵 $\rho$ 的拷贝供应，在系统寄存器上近似实现 $e^{-i\rho t}$。实现依据 Lloyd–Mohseni–Rebentrost 2014（"Quantum principal component analysis", Nature Physics 10, 631）的 LMR 协议：利用 SWAP 算符的本征值结构，对一份 $\rho$ 拷贝与系统做部分交换 $e^{-i\Delta t\cdot\mathrm{SWAP}}$，丢弃拷贝后系统上的有效通道即 $e^{-i\rho\Delta t}$ 的一阶近似，单步误差 $O(\Delta t^2)$；`copies` 份拷贝串联给出总时间 $t = \mathrm{copies}\cdot\Delta t$，一阶误差 $O(t\cdot\Delta t)$ 随拷贝数线性下降。

该原语是 [QPCA 主成分分析](qpca.md)相位估计的酉步来源。

## 接口与输入模型

```python
density_matrix_exponentiation(preparation, *, time, copies, swap_width=None, name=None)
```

- `preparation`：$\rho$ 拷贝的态制备（`StatePreparation`），input model 为 SP + QRAM：纯态用 `gate_state_prep(amplitudes)`，混合态经 `density.gate_purification(rho)` 的 `PurificationAccess.as_state_preparation()` 适配（target = system ⊕ environment，零 work）；制备必须零 work，非零 work 在生成期抛 `ValidationError`。
- `time`：总演化时间 $t$，必须为正的有限实数。
- `copies`：拷贝数，范围 1..63；每步 $\Delta t = t/\mathrm{copies}$。
- `swap_width`：参与交换的前缀位宽，缺省为制备的整个 target；纯化场景应取 $\rho$ 的系统位宽，环境位留在拷贝中不参与交换（效果等同于取偏迹）。

返回 `Operation`，寄存器为 `system`（swap_width 位）与 `copies`（copies × 制备宽度）。系统的输入态由调用方准备；拷贝寄存器从全零由制备 oracle 填充。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"density_matrix_exponentiation"` |
| `copies` | 拷贝数 |
| `step_time` | 单步时长 $\Delta t = t/\mathrm{copies}$ |
| `error_scaling` | `"O(time * step_time)"`（一阶误差标度） |
| `reference` | `"Lloyd-Mohseni-Rebentrost 2014, Nature Physics 10, 631"` |

## 实现要点

部分交换按恒等式 $\mathrm{XX}+\mathrm{YY}+\mathrm{ZZ} = 2\cdot\mathrm{SWAP} - I$ 精确分解：三个轴的 Pauli 对旋转 $e^{-i(\Delta t/2)\,P\otimes P}$ 相互可交换、逐轴逐位施加；单轴旋转基于奇偶校验（XOR–Rz–XOR），$y$ 轴经 $R_x(\pi/2)$ 共轭切换基，$I$ 分量以全局相位 $-\Delta t/2$ 修正。

生成链为：第 $k$ 份拷贝从全零由制备 oracle 填充 → `system` 与该拷贝的前 swap_width 位逐位施加角度 $\Delta t$ 的部分交换。寄存器布局 `system | copies`，第 $k$ 份拷贝占位段 $[k\cdot w, (k+1)\cdot w)$（$w$ 为制备宽度）。

设计决策：算法只要求拷贝可重复制备——多拷贝即对制备 oracle 的多次调用，$\rho$ 的矩阵元素无需显式给出，这是 QRAM 系 QML 的资源前提。丢弃拷贝即对拷贝寄存器取偏迹，纯化的环境位保持纠缠留在拷贝中。适用边界：一阶近似，需要更高精度时增加 copies（减小 $\Delta t$），误差标度见 `error_scaling` 属性。

## 验证方案

类别 C3（概率分布语义，判定准则见 `../development/validation-plan.md` §2）：丢弃拷贝后的系统态须在迹距离上逼近闭式期望。三层证据：

- 结构：`tests/core/test_qpca.py:QpcaTests.test_invalid_inputs_fail_at_generation` 覆盖指数化入口的 `time = 0`、`copies = 0`、非 `StatePreparation` 输入等生成期违例。
- 数值：`test_small_step_first_order_accurate`——$\rho = |+\rangle\langle+|$、单步 $\Delta t = 0.05$ 作用在 $|0\rangle$ 上，丢弃拷贝后与精确 $e^{-it\,|+\rangle\langle+|}$ 作用的迹距离 $< 0.6\cdot\Delta t^2$（$\Delta t$ 一阶收敛率，实测系数约 0.53）；`test_error_halves_with_copies`——固定总时间 $t = 0.4$、copies = 1, 2, 4，迹距离逐级下降且每级 $< 0.6\times$ 前级（copies 翻倍即 $\Delta t$ 减半，误差近似减半）。
- 绑定：混合态的 DM 输入经 `gate_purification(...).as_state_preparation()` 适配为 SP 进入电路，由 `QpcaTests.test_mixed_state_eigenvalue_via_purification` 见证（见 [QPCA 主成分分析](qpca.md)）。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（$\Delta t$ 一阶率 + copies 翻倍误差减半标度 + 纯化适配路径）。按 validation-plan §5 的 V3 计划，QPCA 模块待登记 catalog 案例并接入真实后端对拍。

## 相关链接

- 源码：`src/pyqecclang/algorithms/qpca.py`
- 同模块页面：[QPCA 主成分分析](qpca.md)
- API 参考：[QPCA 量子主成分分析](../../api/algorithms/qml/qpca.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_misc_algorithms.py`（misc_algorithms 组），全部在真实后端上执行。

**实验设计**：两个密度矩阵指数化案例——(a) 纯态 $\rho = \lvert+\rangle\langle+\rvert$（`gate_state_prep`，总时间 $t = 0.4$）；(b) 混合态 $\rho = \mathrm{diag}(0.75, 0.25)$（`gate_purification` 纯化适配、`swap_width = 1`，$t = 0.5$）。系统初态均为 $\lvert 0\rangle\langle 0\rvert$，拷贝数 copies = 1/2/4（最多 9 量子位）。经典预言机为 numpy 特征分解的精确演化 $e^{-i\rho t}\sigma e^{i\rho t}$（与被测实现完全独立）。后端路径：`reference`（内置参考执行器）、`rir-pysparq`（PySparQ 原生 RIR）、`adapter-pysparq`（PySparQ 适配器）、`originir-ext`（UniQC 全振幅态向量），逐振幅与偏迹后密度矩阵两级对拍。

**关键指标**：

| 案例 | 规模 | 路径 | 迹距离（copies = 1/2/4） | 误差比率 | 跨后端偏差 |
|---|---|---|---|---|---|
| dm-exponentiation-pure-convergence | t = 0.4，≤5 量子位 | 四路径全振幅对拍 | 8.55e-2 / 4.34e-2 / 2.20e-2 | 0.508 / 0.507 | ≤ 2.1e-16 |
| dm-exponentiation-mixed-convergence | t = 0.5，≤9 量子位 | 四路径全振幅对拍 | 5.75e-2 / 2.97e-2 / 1.52e-2 | 0.516 / 0.512 | ≤ 3.3e-16 |

copies 翻倍（即 $\Delta t$ 减半）时迹距离比率稳定在约 0.51，直接观测到 LMR 协议的一阶收敛标度 $O(t\cdot\Delta t)$；四条后端路径在机器精度内一致。

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

产物：`out/verification/misc_algorithms.json`（24 个案例全过，本页对应 `dm-exponentiation-*` 两个案例）。
