# Heinrich 量子求和（Heinrich Quantum Summation）

[English](../../../index.html) · **简体中文**

> 类别 C3 · 模块 [`oracq.algorithms.common.integration`](../../api/algorithms/common/integration.rst) · 阶段 V1

## 概述

估计非负整数函数 $f$ 在均匀网格上的均值 $E[f] = \frac{1}{N}\sum_i f(i)$，并以此组装一维数值积分。实现依据 Heinrich 2002（"Quantum Summation with an Application to Integration", J. Complexity 18(1)，另见 Novak 2001 的函数类量子求积率）。

读出采用比较器构造：阈值寄存器取均匀叠加后与函数值比较，好状态（flag = 1）概率**恰为** $E[v]/2^w$——对 $v$ 严格线性，无需小角度近似；对该标记做标准振幅估计，查询复杂度 $O(1/\varepsilon)$，相对经典 Monte Carlo 的 $O(1/\varepsilon^2)$ 呈二次改进。

## 接口与输入模型

```python
quantum_sum(database, *, precision=4, name=None)
```

API 入口：{obj}`quantum_sum <oracq.algorithms.common.integration.quantum_sum>`

- `database`：函数值加载器（{obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>`，address = index、data = value），input model 为 FO + QRAM，与仓库的三层绑定（abstract / gate / qram）直接兼容；通常由 {obj}`table_loader(values, data_width=None, backend="gate"|"qram") <oracq.algorithms.common.integration.table_loader>` 构造。
- `precision`：相位寄存器位数，范围 1..63；估计误差量级 $O(1/2^{\text{precision}})$。

返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`，寄存器为 `target`、`work`、`phase`。读出 `phase` 后用 {obj}`mean_from_phase(value, precision, data_width) <oracq.algorithms.common.integration.mean_from_phase>` 解码均值估计。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"quantum_sum"` |
| `readout_register` / `decoder` | `"phase"` / `"mean_from_phase"` |
| `value_bits` / `index_bits` | 值字宽 $w$ 与下标位数 $n$ |
| `query_complexity` / `classical_query_complexity` | `O(1/epsilon)` / `O(1/epsilon**2)` |

相关入口：{obj}`quantum_integral <oracq.algorithms.common.integration.quantum_integral>`（一维积分，均值乘区间长度，解码用 {obj}`integral_from_phase <oracq.algorithms.common.integration.integral_from_phase>`）与 {obj}`heinrich_rate(smoothness, dimension) <oracq.algorithms.common.integration.heinrich_rate>`（函数类最优收敛率：确定性 $s/d$、随机化 $s/d + 1/2$、量子 $s/d + 1$）。

## 实现要点

生成链为 {obj}`sum_preparation <oracq.algorithms.common.integration.sum_preparation>`（均匀 index + 函数值加载 + 阈值比较）→ {obj}`sum_iterate <oracq.algorithms.common.integration.sum_iterate>`（标记由制备 target 的 flag 位驱动的 Grover 迭代）→ {obj}`phase_estimation <oracq.algorithms.common.estimation.phase_estimation>`。寄存器布局：`target = index(n) | threshold(w) | flag(1)`，`work = value(w)`；value 字与 index 纠缠留在 work（比较器读出不需要复净它，制备标注 `clean_work=False`），调用方按 flag 标记后应逆调用制备复原。

适用边界：函数值必须量化为 $w$ 位非负整数（`table_loader` 逐值校验字宽）；均值估计的精度由 QAE 栅格决定，`precision` 每加 1 位栅格密度翻倍。积分路线（`quantum_integral`）的总误差 = 离散化误差（由网格与光滑性决定，见 `heinrich_rate`）+ QAE 估计误差。

## 验证方案

类别 C3（概率分布语义，判定准则见 `../development/validation-plan.md` §2）：输出分布须等于闭式期望。见证已齐，三层证据：

- 结构：`tests/core/test_integration.py:SumPreparationTests` / `QuantumSumTests` / `RateTests` 的构造与属性断言；`RateTests.test_invalid_inputs_fail_at_generation` 覆盖全部入口的参数违例（空表、负值、字宽越界、非法 backend、precision / interval / 光滑性参数越界等）。
- 数值：`SumPreparationTests.test_flag_probability_matches_mean` 与 `test_constant_table_exact` 对拍线性恒等式 $P(\text{flag}=1) = E[v]/2^w$（places = 12）；`QuantumSumTests.test_mean_on_qae_grid_is_exact` 在均值恰落 QAE 栅格时要求全部非零概率读出精确等于真值（places = 9）；`test_ramp_mean_within_qae_resolution`（delta = 0.5）与 `test_quantum_integral_trapezoid_scale`（对拍量化均值 delta = 0.06、积分真值 0.5 处 delta = 0.09）覆盖栅格外的分辨率界；`RateTests.test_heinrich_rate_known_values` 校验收敛率闭式值。
- 绑定：`SumPreparationTests.test_qram_binding_matches_gate`（gate / qram 两绑定的 flag 概率逐点对拍，places = 12）与 `test_abstract_loader_binds`（abstract 声明经 {obj}`bind <oracq.infrastructure.linking.bind>` 绑定后概率不变）共同覆盖三层一致性。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（闭式均值对拍 + `heinrich_rate` 校验 + 三层一致性）。按 validation-plan §5 的 V3 计划，Heinrich 积分待登记 catalog 案例并接入真实后端对拍。

## 相关链接

- 源码：`src/oracq/algorithms/common/integration.py`
- 同族页面：[Heinrich 量子积分](heinrich-integration.md)
- API 参考：[量子求和与积分](../../api/algorithms/common/integration.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_nt_qlss_sde.py`（nt_qlss_sde 组），全部在真实后端执行；经典预言机（经典求和 $E[v]$、QAE 双峰 Dirichlet 理论分布）独立构造，不经过被测实现的辅助函数。

**实验设计**：(a) 比较器恒等式——常数表、斜坡表、伪随机 16 值表三个实例，在 reference、rir-pysparq、adapter-pysparq 与 OriginIR-ext 上验证好状态概率**恰为** $E[v]/2^w$；(b) 均值恰落 QAE 栅格的常数表（$E[v]/2^w = 1/2$），要求所有非零概率读出解码后等于真值；(c) QAE 读出分布对照独立理论 $\frac{1}{2}\left[D^2(y - y_\theta) + D^2(y + y_\theta)\right]$（$y_\theta = 2^p\theta/\pi$，$\sin^2\theta = E[v]/2^w$），并做 precision 3→5 的期望误差收敛扫描；(d) gate 与 qram 两种加载器绑定在制备全振幅与求和 phase 分布上逐点对拍（qram 程序的资源带嵌套前缀 `prep__db__table`、`qpe__u__prep__db__table`，按入口资源名绑定同一数据表）；(e) `heinrich_rate` 收敛率闭式值。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| sum-preparation-flag（3 张表） | $n \le 4$、$w \le 4$ | 四路径 | $\lvert P(\text{flag}) - E[v]/2^w\rvert$ | 1.1e-16 |
| quantum-sum-on-grid | $w=2$、precision 4 | 三路径 | 读出解码最大偏差 / 跨后端 TVD | 4.4e-16 / 0.0 |
| quantum-sum-qae-p3/p4/p5 | $w=2$ | 双路径 | TVD vs QAE 理论 | 3.5e-15 / 6.2e-15 / 1.3e-14 |
| quantum-sum-qae-convergence | precision 3→5 | reference | 期望绝对误差 $E\lvert\hat\mu - 1.5\rvert$ | 0.7205 → 0.4752 → 0.2267（严格递减） |
| table-loader-qram-vs-gate | $w=3$、precision 3 | 双路径 | 制备振幅 / flag 概率 / 求和分布偏差 | < 1e-9 一致 |
| heinrich-rate-closed-form | 4 组 $(s,d)$ | 经典 oracle | 与 $s/d$、$+1/2$、$+1$ 的偏差 | < 1e-15 |

QAE 期望误差按 $\Theta(\log M / M)$ 收缩（Dirichlet 核重尾的一阶矩），单调递减验证了 $O(1/\varepsilon)$ 查询复杂度的实际表现；众数估计在各 precision 下均落在一阶分辨率界 $2^w\pi/2^p$ 内。

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

产物：`out/verification/nt_qlss_sde.json`（`sum-preparation-*`、`quantum-sum-*`、`table-loader-*`、`heinrich-rate-*` 共 10 个案例）。
