# Costa 行走线性系统求解器（Costa Walk QLSS）

<a href="../../../manual/algorithms/costa-walk.html">English</a> · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.qlss.qlss`](../../api/algorithms/qlss/qlss.rst) · 阶段 V3

## 概述

基于离散参数化量子行走的 QLSS（Costa 等人, arXiv:2111.08152，输入模型讨论见 [QFVM 输入模型审查](../../reference/qfvm-input-models.md)）：把 $A^{-1}b$ 的求解改写为插值哈密顿量 $H(s)$ 的行走算子序列（RHS 零态反射、$R(s)$ 旋转、受控 $U_A$ / $U_A^\dagger$ 与信号零反射），再用 Laurent 多项式 filter 在行走幂上读取解。行走算子（{obj}`costa_walk <oracq.algorithms.qlss.qlss.costa_walk>`）、调度（{obj}`schedule <oracq.algorithms.qlss.qlss.schedule>`）、Dolph–Chebyshev filter（{obj}`dolph_chebyshev_plan <oracq.algorithms.qlss.qlss.dolph_chebyshev_plan>` + {obj}`lcu_filter <oracq.algorithms.qlss.qlss.lcu_filter>`）与问题层入口（{obj}`make_costa_qlss <oracq.algorithms.qlss.qlss.make_costa_qlss>`）分开组装。输入模型为 BE 加 SP；稀疏输入经显式适配（[稀疏矩阵块编码](sparse-block-encoding.md)）转换，不违反论文的输入假设。

## 接口与输入模型

```python
costa_walk(a, bprep, fs)                          # 单步行走算子，返回 Operation
schedule(s, kappa, power=1.5)                      # 调度点 s∈[0,1] → f(s)
CostaConfig(steps=2, kappa=4.0, schedule_power=1.5,
            filter_degree=2, filter_attenuation=0.2)
dolph_chebyshev_plan(degree=2, attenuation=0.2)    # → FilterPlan
lcu_filter(walk, plan)                             # 相干 LCU filter，返回 Operation
unary_weight_preparation(weights)                  # unary 前缀叠加制备
costa_qlss(a, bprep, config=None, *, filtering=None)  # 内核，返回 StateOracle
make_costa_qlss(config=None)                       # → QLSSProtocol（input_model="block_encoding"）
```

API 入口：{obj}`costa_walk <oracq.algorithms.qlss.qlss.costa_walk>`、{obj}`schedule <oracq.algorithms.qlss.qlss.schedule>`、{obj}`CostaConfig <oracq.algorithms.qlss.qlss.CostaConfig>`、{obj}`dolph_chebyshev_plan <oracq.algorithms.qlss.qlss.dolph_chebyshev_plan>`

- `a`：{obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`（BE）；`bprep`：{obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>`（SP）；两者同宽，`fs ∈ [0,1]`。
- `CostaConfig.kappa` 是**编码矩阵** $A/\alpha$ 的逆谱界（要求 $\sigma_{\min}(A/\alpha) \ge 1/\kappa$），不是任意尺度下的 cond(A)。问题层入口按谱声明推导：`replace(config, kappa=system.inverse_norm_bound)`，避免 alpha 与 σ_min 声明脱节。
- 问题输入：{obj}`LinearSystem(block=BlockSystem(encoding, rhs, spectrum)) <oracq.algorithms.qlss.qlss.LinearSystem>`，或稀疏输入自动经 CKS $T^\dagger S T$ 适配；返回 {obj}`SolveResult <oracq.algorithms.qlss.qlss.SolveResult>`（属性与范数探针契约见 [CKS 求解器](cks.md)）。

{obj}`costa_qlss <oracq.algorithms.qlss.qlss.costa_qlss>` 输出属性：`algorithm = "costa_qlss"`、`input_alpha`、`encoded_inverse_bound`、`normalization_assumption = "sigma_min(A / alpha) >= 1 / kappa"`、`steps`、`filtering`、`kernel_status = "prototype; finite-step accuracy and success probability require validation"`、`success_condition = "signal == 0; probability and solution accuracy unverified"`。

## 实现要点

`costa_walk` 信号布局为 `enc(a.signal_qubits) | bw(RHS work) | a1 | a2 | a3 | a4`。电路按论文结构组装：$U_b^\dagger$ → RHS 零态反射（同时要求 target 与 RHS work 位为零——$U_b$ 是 target+work 上的酉扩张，投影对象必须完整）→ $U_b$；$R(s)$ 旋转写成 $R_y(2\arctan\frac{f}{1-f}) \cdot Z$ 作用在 a2；受控 $U_A$ / $U_A^\dagger$ 与 a2 零反射构成行走主体；末尾对 `enc | a2 | a3` 做正反射并附 $\pi/2$ 全局相位。输入契约经 {obj}`operator_state_contract <oracq.algorithms.input_model.interfaces.operator_state_contract>` 检查。

调度 $f(s) = \frac{\kappa}{\kappa-1}\bigl(1 - (1 + s(\kappa^{p-1} - 1))^{1/(1-p)}\bigr)$（$\kappa = 1$ 时退化为 $s$）。`costa_qlss` 在 $s_i = (i+1)/\text{steps}$ 处生成各步行走，在最后一步 walk 上施加 Dolph–Chebyshev filter：{obj}`unary_weight_preparation <oracq.algorithms.qlss.qlss.unary_weight_preparation>` 在 clock 寄存器制备权重，逐 clock 位受控重复行走（负 offset 幂经 adjoint），相干叠加出 Laurent 多项式。问题层 `solve()` 再统一做物理通道选择与独立范数探针。

适用边界：`kernel_status` 标注有限步精度与成功概率仍需验证，默认两步只演示程序构造，不保证求解精度；alpha 与 σ_min 声明冲突在 check 阶段报 `INPUT_SPECTRUM`，不静默生成错误调度的电路。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。与验证覆盖矩阵 `qlss.py` 行一致：

- 结构：`tests/core/test_qlss_input_models.py:QLSSInputTests`（构造与契约断言）。
- 数值：`test_protocols_consume_different_models_and_scale_kappa`——声明 `kappa=999` 的配置在问题层被替换为 `encoded_inverse_bound = 8`（$\alpha = 2$、$\sigma_{\min}$ 下界 0.25），底层 `costa_qlss` 模块属性同步，且与 CKS 路线同宽、同 alpha、程序无未绑定槽位；`test_matrix_probe_recovers_scalar_system_norm`（CKS 范数探针恢复 $\lVert x\rVert = 4$，places = 10）。同类的 `test_qfvm_preserves_sparse_input_for_both_solvers` 见证 QFVM 稀疏问题两条路线输出同宽解态、`encoded_inverse_bound = 72`，仅 Roe oracle 槽位保持开放。
- 绑定：`test_costa_rhs_reflection_is_independent_of_unitary_extension`——同一 $|b\rangle$ 的两个酉扩张（直接基态制备 vs swap 进 work 位）生成的 `costa_walk(identity(1), ·, 0.3)` 在全部 target/signal 基态输入下幅度逐点一致（places = 11），见证 RHS 零态反射对 work 位的处理不依赖具体扩张。

## 已知缺口与计划阶段

一般谱下的有限步误差界及按目标精度自动选步数/filter 参数尚未完成。下述二维案例只验证指定系统，不认证任意谱、应用成功成本、层析或硬件性能。增加步数本身不构成单调收敛证书。

## 相关链接

- 源码：`src/oracq/algorithms/qlss/qlss.py`
- 同族页面：[CKS Chebyshev 求解器](cks.md)、[稀疏矩阵块编码](sparse-block-encoding.md)、[VTAA-CKS 变时求解器](vtaa-cks.md)
- API 参考：[量子线性系统](../../api/algorithms/qlss/qlss.rst)
- 输入模型审查：[QFVM 输入模型审查](../../reference/qfvm-input-models.md)
- 概念：[算法自己的约定：从一个 gate 开始](../contracts.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

实验设计：最小实例（1 量子位目标）：`a` 为对角块编码（{obj}`diagonal_block_encoding <oracq.algorithms.input_model.oracles.diagonal_block_encoding>`，angle_scale $=\pi/3$），`bprep` 为基态制备 {obj}`basis_state(1, 0) <oracq.algorithms.input_model.oracles.basis_state>`，调度点 $f_s=0.5$，信号共 6 位（enc 2 + a1..a4），总计 7 量子位。验证两个可判定性质：

1. 幺正性：OriginIR-ext 导出经 UniQC `Circuit.to_matrix` 得 128 维矩阵 $W$，计算 $\|W^\dagger W - I\|_{\max}$；
2. 后端一致性：零输入态上的执行结果在 reference、rir-pysparq、adapter-pysparq、originir-ext 四路径间逐振幅对拍。

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| costa-walk-unitarity | 7 qubits, fs=0.5 | to_matrix + 4 路径 | $\|W^\dagger W-I\|_{\max}$ | 5.6e-16 |
| 〃 | 〃 | 〃 | 跨路径态最大误差 | 3.1e-17 |

行走算子幺正且各后端一致，说明电路组装（反射、受控 $U_A$/$U_A^\dagger$、旋转调度与末尾正反射）在算符层面自洽。该历史单步行走检查不覆盖求解成功通道；下面的独立回归覆盖初始化、读出和逆解数值。

复现命令：

```bash
PYTHONPATH=src python tests/verification/verify_search_walks.py
```

产物：`out/verification/search_walks.json`（案例 `costa-walk-unitarity`）。


## 初始本征态与物理解读出（0.3.1）

本实现采用 $W(0)=i(2\Pi-I)U_H(0)$，其中 $\Pi$ 将 `enc | a2 | a3`
投影到零。设 $v$ 为信号全零的 RHS 态，令 $w=U_H(0)v$；此时
`a4=a3=1`、`a2=+`、`a1=0`。有 $W(0)v=-iw$、$W(0)w=iv$，
故 $(v-iw)/\sqrt{2}$ 才是 $+1$ 本征态。内核用 `H(a4)`、精确
`S†(a4)`（`Rz(-pi/2)` 加 `-pi/4` 全局相位），然后以 `a4` 控制
`H(a2)` 与 `X(a3)` 制备。保留全局相位以保持受控调用的相干语义。

最后插值点的逆解在本仓 $A_f$ 约定下位于 `a1=1`。filter 后的
`X(a1)` 将该物理分支映射到公开的 signal 全零成功条件，所有辅助位
仍参与后选择。这修正了 0.3.0 内核只从 $v$ 出发且选择相反末态分支的问题。

`tests/core/test_costa_channels.py` 覆盖 Hermitian/复非对称矩阵、奇偶步数、
带干净 work 的酉扩张、受控 adjoint，以及 Costa 问题层独立范数探针。
相同成功通道案例在原 0.3.0 确实失败；并非只修改标签，也未使用 mock。

`tests/integration/test_costa_channels_native.py` 在真实 PySparQ 逐门执行
$A=L\operatorname{diag}(0.4,0.9)R$，$L=R_y(0.61)$，分别取
$R=R_y(-0.61)$ 或 $R=R_z(0.83)$，以及
$b=(\sqrt{0.3},e^{0.4i}\sqrt{0.7})$。512 步、$\kappa=2.5$、默认
filter 下，两例成功概率约 0.499974，对独立 NumPy 逆解的相位对齐条件
$L_2$ 误差约 $4.24\times10^{-4}$（要求 $10^{-3}$）。复数例还排除了
误算 $A^\dagger$ 逆解；不使用 native 算术捷径。这是正确性见证，不是加速证明。

```bash
PYTHONPATH=src python -m unittest discover -s tests/core -p test_costa_channels.py -v
# 需真实可选 PySparQ 后端及 NumPy：
PYTHONPATH=src python -m unittest discover -s tests/integration -p test_costa_channels_native.py -v
```

推导参考：[Costa 等人 arXiv:2111.08152](https://arxiv.org/abs/2111.08152)
第 IV 节与附录 E；上述方程采用本仓辅助位约定。
