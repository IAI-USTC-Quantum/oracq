# VTAA-CKS 变时线性系统求解器（VTAA-CKS QLSS）

> 类别 C2 · 模块 `pyqecclang.algorithms.qlss.vtaa_cks` · 阶段 V3

## 概述

Childs–Kothari–Somma（arXiv:1511.02306, SIAM J. Comput. 2017）第 5 节的变时幅度放大（VTAA）线性系统求解器：在 [CKS 第 4 节基础求解器](cks.md) 的 $\kappa^2$ 路线之上叠加变时层，把查询复杂度降到 $O(d\,\kappa\,\mathrm{polylog}(d\kappa/\epsilon))$。截至实现时，全网没有 VTAA 或 CKS §5 的公开实现（Qiskit/PennyLane/Qrisp 等仅有固定时间幅度放大或经典成本模型），本模块是该算法路线的首个电路级实现。

变时层按论文结构组装：

- **时钟位 $C_j$**（Lemma 22 GPE）：每一步对 block encoding 的 qubitization walk 施加判决多项式，把"本征值是否大到可以在本频带求逆"相干写入时钟位；判 1 振幅即 $|P(\lambda/\alpha)|$。
- **分频带逆 LCU**（Lemma 23 $W(\lambda,\delta)$）：受控于 $C_j=1$ 时施加本频带的截断 Chebyshev 逆多项式，式 (98) 的旋转把各频带成功振幅均匀压到 $1/\alpha_{\max}$。
- **时钟感知嵌套放大**：Ambainis（arXiv:1010.4458）的 VTAA 级联，操作符形式取 Low–Su（arXiv:2410.18178）式 (47)–(53)；正确性不依赖放大日程，日程只影响成功率。
- **$A'$ 反计算**（式 (99)–(112)）：把 $W_j$ 替换为纯旗标翻转的 $A'$ 取逆，抹除 GPE 时钟与垃圾，解态留在 target，成功条件为 `signal == 0`。

GPE 判决采用 Low–Su Prop 23 的确定性路线（与 CKS Lemma 22 的 PEA + 多数投票同阶 $O((\alpha/\theta_j)\log(1/\epsilon))$ 次查询）：复用 `qsvt.fixed_point_search_phases` 已验证的 Yoder–Low–Chuang 定点多项式合成，对 walk 施加 `qsvt_sequence`。fire 侧（$|x| \geq \theta_j$）有硬界 $\epsilon$；谱低端的过渡带（CKS 未承诺带）行为确定、可精确计算，但两侧都不落在验证区时只有混合语义保证。

## 接口与输入模型

```python
make_vtaa_cks_qlss(config=None)   # 返回 QLSSProtocol（input_model="sparse"）
vtaa_cks(system, config=None)     # 底层内核，接收 SparseSystem
VTAAConfig(order=2, terms=None, clock_steps=None, marker_epsilon=0.02,
           degree_cap=40, rounds=None)
gapped_phase_estimation(a, threshold, x_edge, *, epsilon=0.02, degree_cap=40)
band_inverse_step(a, coefficients, alpha_max)
tunable_rounds(stage_amplitudes, thresholds=None)  # Low–Su 式 (52)–(53) 日程
```

问题输入与 `make_cks_qlss` 相同的 `LinearSystem(sparse=...)`（Hermitian、非负对角声明）。`clock_steps` 缺省由声明的物理条件数 $\kappa_{\mathrm{phys}} = $ `norm_upper / sigma_min_lower` 推导（$\lceil\log_2\kappa\rceil+1$），要求 $2^{m-1} \geq \kappa_{\mathrm{phys}}$ 覆盖最细频带。`rounds` 是各阶段放大轮数（缺省全零，即纯变时层加后选）；生产部署应按 Ambainis 算法 2 的振幅估计或 `tunable_rounds` 的确定性日程选择。

编码归一化需相对谱上界留出松弛（`norm_upper / alpha < 1`）：walk 相位 $\arccos(\lambda/\alpha)$ 的判决几何由此确定，过紧的编码（如对角谱配 `entry_bound = max|A|`）会被拒绝，请放宽 `entry_bound`。

内核输出属性：`algorithm = "vtaa_cks"`、`clock_steps`、`fire_thresholds`、`band_orders`、`band_lcu_normalizations`、`alpha_max`、`marker_degrees`、`marker_epsilon`、`rounds`、`implementation_scope`（引用三条论文线索）、`kernel_status = "prototype; band polynomial accuracy and VTAA schedule pending"`、`success_condition = "signal == 0"`。协议层沿用 `QLSSProtocol.solve()` 的物理子空间选择与独立范数探针。

## 实现要点

生成链：`real_symmetric_sparse_encoding`（$\alpha = s\cdot a_{\max}$）→ 每频带 `gapped_phase_estimation`（`qsvt_sequence` 判决 + `signal==0` 受控翻转时钟位，P_j 垃圾按论文保留）与 `band_inverse_step`（`chebyshev_block` 奇次幂 LCU + 式 (98) 均匀化旋转）→ `vtaa_variable_step`（前缀全 0 受控 GPE、$C_j=1$ 受控 $W_j$）逐级嵌入 `vtaa_prefix` / `vtaa_amplified_stage`（反射 $R_f$ 翻转 stopped∧失败 分支相位、$R_s$ 经前缀逆 + 全零反射构造）→ 顶层调用放大链后以 `vtaa_uncompute_step`（GPE 重放 + 纯旗标翻转）的逆抹除。模块调用与 Repeat 全部符号化保留，判决与求逆的每步信号寄存器独立分配。

与论文的两点差异（都记录在 `implementation_scope`）：GPE 用确定性 QSP 判决替代 PEA + 多数投票（Low–Su Prop 23，同一引理的现代实现，避免概率性判决分布并使参考模拟开销与多项式度数线性相关）；频带逆多项式的阶数日程为原型的几何放大（$2^{j-1}$），对 $1/x$ 的逼近精度未做理论标定。

## 验证方案

类别 C2（近似连续语义）。与验证覆盖矩阵 `vtaa_cks.py` 行一致：

- 结构：`tests/core/test_vtaa_cks.py:VTAAStructureTests`（协议输出属性、序列化往返、放大轮数的 Repeat 保留、符号化资源估计、配置拒绝路径）。
- 数值：`GappedPhaseEstimationTests.test_marker_circuit_matches_qsp_response_exactly`（判决电路与 QSP 响应在定点网格特征值上逐点一致，places = 11）；`BandInverseTests.test_band_inverse_step_matches_chebyshev_polynomial`（频带逆 LCU 的成功分支与 Chebyshev 多项式对拍，places = 11）；`VTAAEndToEndTests.test_uniform_spectrum_single_band_is_exact`（单频带端到端精确携带 $|b\rangle$ 方向）；`test_variable_time_clock_separates_bands`（双频带变时结构：逐本征值路径幅值由判决响应精确给出，端到端幅值落在归一化耦合区间内，且 band1 直接求逆分量占优）。
- 绑定：`test_amplification_schedule_preserves_conditional_solution`（同一条件解态在放大轮数改变后不变——日程无关性见证）。

## 已知缺口与计划阶段

频带逆多项式对 $1/x$ 的精度未标定（阶段 V3）：当前系数复用 `CKSConfig` 闭式并以几何阶数放大，与论文要求的逐频带 $\widetilde O(2^j)$ 度数日程的对应关系待核验；过渡带混合语义的端到端对拍只到幅值区间。VTAA 放大日程未接振幅估计通道（`tunable_rounds` 已给出 Low–Su 公式，缺阶段范数估计电路）。度数上限 40（继承 `qsvt.py` 合成上限）限制了可覆盖的 $\kappa$ 声明范围。与 HHL 论文参考值的端到端对拍同 `qlss.py` 行缺口。

## 相关链接

- 源码：`src/pyqecclang/algorithms/vtaa_cks.py`
- 论文：[CKS arXiv:1511.02306](https://arxiv.org/abs/1511.02306) §5、[Ambainis arXiv:1010.4458](https://arxiv.org/abs/1010.4458)、[Low–Su arXiv:2410.18178](https://arxiv.org/abs/2410.18178)
- 同族页面：[CKS Chebyshev 求解器](cks.md)、[Costa 行走求解器](costa-walk.md)、[稀疏矩阵块编码](sparse-block-encoding.md)、[定点搜索](fixed-point-search.md)
- API 参考：[VTAA-CKS 变时线性系统求解器](../../api/algorithms/qlss/vtaa_cks.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
