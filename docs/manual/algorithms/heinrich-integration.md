# Heinrich 量子积分（Heinrich Quantum Integration）

> 类别 C3 · 模块 `pyqecclang.algorithms.common.integration` · 阶段 V1

## 概述

估计一维定积分 $\int_0^{L} f(x)\,dx$：把区间均匀划分为 $N$ 个网格点，函数值量化为 $w$ 位非负整数，先以 Heinrich 量子求和估计网格均值 $E[v]$，再乘区间长度得到积分（复合矩形法则）。实现依据 Heinrich 2002（"Quantum Summation with an Application to Integration", J. Complexity 18(1)，另见 Novak 2001 的函数类量子求积率）。

积分估计是量子求和的直接组装：读出侧仍以比较器构造保证好状态概率对函数值严格线性，对标记做标准振幅估计，查询复杂度 $O(1/\varepsilon)$，相对经典 Monte Carlo 的 $O(1/\varepsilon^2)$ 呈二次改进。求和原语本身的细节见 [Heinrich 量子求和](heinrich-summation.md)。

## 接口与输入模型

```python
quantum_integral(database, *, precision=4, interval=1.0, name=None)
```

- `database`：函数值加载器（`XorDatabase`，address = index、data = value），input model 为 FO + QRAM；通常由 `table_loader(values, data_width=None, backend="gate"|"qram")` 构造，函数值按 `v/full_scale` 量化（`full_scale` 缺省取 $2^w - 1$）。
- `precision`：相位寄存器位数，范围 1..63；QAE 估计误差量级 $O(1/2^{\text{precision}})$。
- `interval`：区间长度 $L$，必须为正的有限实数。

返回 `Operation`，寄存器为 `target`、`work`、`phase`（与 `quantum_sum` 相同）。读出 `phase` 后用 `integral_from_phase(value, precision, data_width, interval=1.0, full_scale=None)` 解码：积分估计 $= E[v]/\text{full\_scale} \times L$。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"quantum_integral"` |
| `readout_register` / `decoder` | `"phase"` / `"integral_from_phase"` |
| `interval` | 区间长度 $L$ 回显 |
| `value_bits` / `index_bits` | 值字宽 $w$ 与下标位数 $n$ |
| `query_complexity` / `classical_query_complexity` | `O(1/epsilon)` / `O(1/epsilon**2)` |

相关入口：`heinrich_rate(smoothness, dimension)` 给出函数类最优收敛率（确定性 $s/d$、随机化 $s/d + 1/2$、量子 $s/d + 1$），用于按网格规模估计离散化误差。

## 实现要点

`quantum_integral` 不生成新电路：它调用 `quantum_sum` 得到同一模块，仅把 `algorithm`、`decoder` 改写为积分版本并追加 `interval` 属性（经 `dataclasses.replace` 重建模块属性表）。因此生成链、寄存器布局（`target = index(n) | threshold(w) | flag(1)`，`work = value(w)`）与 Grover 迭代结构与求和完全一致。

设计决策：区间长度不进电路而只进解码器，避免在 IR 中引入浮点缩放；量化尺度 `full_scale` 同理留给 `integral_from_phase`，使同一电路可服务不同区间与量程的重解释。总误差 = 离散化误差（由网格密度与光滑性决定，收敛率见 `heinrich_rate`）+ QAE 估计误差（由 `precision` 决定）。

适用边界：仅一维、均匀网格、非负量化函数值；多维积分与自适应网格未实现。`interval` 必须为正，非正值在生成期抛 `ValidationError`。

## 验证方案

类别 C3（概率分布语义，判定准则见 `../development/validation-plan.md` §2）：输出分布须等于闭式期望。三层证据：

- 结构：`tests/core/test_integration.py:QuantumSumTests.test_quantum_integral_trapezoid_scale` 断言模块属性 `algorithm == "quantum_integral"`；`RateTests.test_invalid_inputs_fail_at_generation` 覆盖 `interval = 0.0` 等参数违例在生成期抛错。
- 数值：`test_quantum_integral_trapezoid_scale`——$f(x) = x$ 在 $[0,1]$ 上取 8 个中点网格（$w = 4$、precision = 4），读出解码后与量化均值闭式值对拍（delta = 0.06），并与积分真值 $0.5$ 对拍（delta = 0.09）；`RateTests.test_heinrich_rate_known_values` 校验收敛率闭式值（$s/d$ 系列与量子率的二次间隔）。
- 绑定：底层加载器的三层一致性由求和侧见证覆盖（`SumPreparationTests.test_qram_binding_matches_gate` 与 `test_abstract_loader_binds`，gate / qram 绑定 flag 概率逐点对拍 places = 12），积分复用同一电路无需独立绑定见证。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（积分闭式对拍 + `heinrich_rate` 校验 + 三层一致性复用求和侧）。按 validation-plan §5 的 V3 计划，Heinrich 积分待登记 catalog 案例并接入真实后端对拍。

## 相关链接

- 源码：`src/pyqecclang/algorithms/common/integration.py`
- 同模块页面：[Heinrich 量子求和](heinrich-summation.md)
- API 参考：[量子求和与积分](../../api/algorithms/common/integration.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_nt_qlss_sde.py`（nt_qlss_sde 组），全部在真实后端执行；经典预言机为独立经典求和与 QAE 双峰 Dirichlet 理论分布。

**实验设计**：(a) 4 点中点网格实例——$f(x)=x$ 在 $[0,1]$ 上取 4 个中点、$w=3$（量化均值恰为 0.5）、precision 4，在 reference、rir-pysparq、adapter-pysparq 上验证完整 phase 分布对照 QAE 理论，并以 `interval=2.0` 重解码同一分布验证区间缩放恒等式（区间长度不进电路，只进解码器）；(b) 8 点中点规范实例（核心测试同款，$w=4$、precision 4、14 位寄存器）在 rir-pysparq 上运行并对照 QAE 理论（reference 交叉由同结构 4 点例的三后端对拍覆盖）。两个实例的判据均含与积分真值 $0.5$ 的偏差。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| quantum-integral-midpoint4 | 4 点、$w=3$、p 4 | 三路径 | TVD vs QAE 理论 / 跨后端 TVD | 8.0e-15 / 7.1e-17 |
| 同上 | | | 众数估计 vs 量化真值 0.5（分辨率界 0.2244） | 0.5714（偏差 0.0714） |
| 同上 | | | interval=2 缩放偏差 | 0.0 |
| quantum-integral-midpoint8 | 8 点、$w=4$、p 4 | rir-pysparq | TVD vs QAE 理论 | 9.9e-15 |
| 同上 | | | 众数估计 vs 量化真值 / 积分真值 0.5 | 0.5333（偏差 0.0333 / 0.0333） |

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

产物：`out/verification/nt_qlss_sde.json`（`quantum-integral-*` 共 2 个案例）。
