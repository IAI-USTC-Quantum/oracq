# 解码量子干涉（Decoded Quantum Interferometry）

<a href="../../../en/index.html">English</a> · **简体中文**

> 类别 C4 · 模块 [`oracq.algorithms.optimization.dqi`](../../api/algorithms/optimization/dqi.rst) · 阶段 V1

## 概述

DQI 求解 GF(2) max-XORSAT：给定约束组 $Bx=v$（$B$ 为 $m\times n$ 的 0/1 矩阵的稀疏行表示），寻找满足约束数最多的赋值。实现依据 Jordan et al. 2024（[arXiv:2408.08292](https://arxiv.org/abs/2408.08292)，图 4 的单权重版本）的线路骨架：在 $m$ 比特 error 寄存器制备权重 $l$ 的 Dicke 态，施加右端项相位 $(-1)^{v\cdot y}$，把 $B^T y$ 可逆计算进 $n$ 比特 syndrome 寄存器，再用可逆经典译码器把 error 寄存器卸载回 $\lvert 0\rangle$，最后对 syndrome 做 Hadamard 变换并测量。后选 error 为零的分支后，测得赋值 $x$ 的概率正比于 $K_l(u(x))^2$——$u(x)$ 是未满足约束数、$K_l$ 是 Krawtchouk 多项式，采样因此偏向满足更多约束的赋值。

## 接口与输入模型

```python
dqi(instance, decoder, *, weight)
XorSatInstance(rows, rhs, num_variables)
abstract_decoder(name, syndrome_width, error_width)
table_decoder(syndrome_width, error_width, table, *, name=None)
bruteforce_decoder(instance, *, max_weight=None, name=None)
dicke_state(m, weight)
```

API 入口：{obj}`dqi <oracq.algorithms.optimization.dqi.dqi>`、{obj}`XorSatInstance <oracq.algorithms.optimization.dqi.XorSatInstance>`、{obj}`abstract_decoder <oracq.algorithms.optimization.dqi.abstract_decoder>`、{obj}`table_decoder <oracq.algorithms.optimization.dqi.table_decoder>`

- `instance`：{obj}`XorSatInstance <oracq.algorithms.optimization.dqi.XorSatInstance>`，约束的稀疏行表示（每行是参与该约束的变量下标，行内不重复；`rhs` 取 0 或 1）。`num_variables` 即 syndrome 位宽 $n$，约束数 `num_constraints` 即 error 位宽 $m$。约束实例是经典数据直接参数化，input model 为 CP（编码约束）。
- `decoder`：{obj}`DecoderOracle <oracq.algorithms.optimization.dqi.DecoderOracle>`，语义为 `|syndrome, error> → |syndrome, error XOR D(syndrome)>`，是 FO 类开放输入（可逆经典函数，`reversible_function` 范式）。{obj}`abstract_decoder <oracq.algorithms.optimization.dqi.abstract_decoder>` 声明槽位后经 {obj}`bind <oracq.infrastructure.linking.bind>` 分批绑定；{obj}`table_decoder <oracq.algorithms.optimization.dqi.table_decoder>` 用显式查询表做 gate 见证（未列出的综合征映射到零错误）；{obj}`bruteforce_decoder <oracq.algorithms.optimization.dqi.bruteforce_decoder>` 枚举全部 $2^m$ 个错误模式给出不超过 `max_weight` 的最轻错误，只接受 $m \le 16$。
- `weight`：Dicke 态权重 $l$，范围 0..m；译码半径需覆盖该权重，译码器位宽必须与实例的 n/m 一致。

返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`，寄存器为 `error`(m) 与 `syndrome`(n)。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` / `field` | `"dqi"` / `"GF(2)"` |
| `num_constraints` / `num_variables` | m / n |
| `dicke_weight` | l |
| `readout_register` | `"syndrome"` |
| `postselection` | `"error_zero"` |

{obj}`dicke_state(m, weight) <oracq.algorithms.optimization.dqi.dicke_state>` 单独导出，制备 $\lvert D_l^m\rangle$（显式幅度构造，$m \le 16$）；`XorSatInstance.satisfied_count(assignment)` 提供经典侧对拍用的满足数统计。

## 实现要点

生成链：Dicke 制备（{obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>` 显式幅度，以 syndrome 的零宽视图作 work）→ 右端项相位（对 `rhs[i] = 1` 的 error 位施加 Z）→ 综合征计算（对每条约束的每个变量下标施加 `xor(error_i → syndrome_j)`，合计即 $B^T y$）→ 译码器调用 → syndrome 上的 Hadamard。译码成功的分支 error 回到 $\lvert 0\rangle$；后选（丢弃 error 非零的分支）由调用方按 `postselection` 属性完成。

设计决策：译码器是输入模型的一部分而非算法内部细节——无法译码的综合征可以映射到任意错误模式（对应分支在后选中被淘汰），高效经典译码（如 belief propagation 的可逆实现）经 `abstract_decoder` 加 `bind` 接入，与仓库的开放声明三层范式一致。适用边界：当前只支持 GF(2)，GF(q) 情形需要 q 元离散 Fourier 变换与广义 Dicke 态、寄存器按 $\log_2 q$ 分子组织，留作扩展；显式幅度 Dicke 与穷举译码只服务 $m \le 16$ 的小实例见证，大实例应换成专用 Dicke 线路（如 Bartschi–Eidenbenz 的 O(l·m) 构造，经 `state_prep_isometry` 开放声明接入）与高效译码器。

## 验证方案

类别 C4（判定准则见 `../../development/validation-plan.md` §2：严格优于随机基线，且小实例达到已知最优/理论分数）。三层证据位于 `tests/core/test_dqi.py:DqiTests`：

- 结构：`test_invalid_inputs_fail_at_generation` 覆盖 12 类生成期违例（越界/重复下标、rhs 非 0/1、空约束、Dicke 权重越界、译码表越界、位宽不匹配、decoder 类型错误、负权重等）。
- 数值：`test_planted_instance_beats_random_guessing`——7 约束 3 变量的植入实例（由赋值 x\* = 0b101 植入右端项），权重 1 时 syndrome 分布与 Krawtchouk 闭式 $K_l(u(x))^2$ 逐点对拍（places = 10），期望满足数 6.5 严格优于随机基线 3.5；`test_identity_instance_with_weight_two`——B 取单位阵的 m = n = 5 实例权重 2，分布对拍后期望满足数回到随机基线 m/2（断言用 `assertGreaterEqual` 加 1e-9 容差）；`test_dicke_state_weight_and_uniformity`——$\lvert D_2^5\rangle$ 恰有 C(5,2) = 10 个等幅非零分量（1/√10，places = 12）且 work 复净。
- 绑定：`test_abstract_decoder_binds_to_witness`——抽象译码器是程序唯一的未解析槽，经 `bind` 绑定穷举实现后无未解析声明，syndrome 分布与 gate 见证逐点一致（places = 10）。

共享的 `check_distribution` 还钉死译码成功时 error 寄存器确定性复净（P(error = 0) = 1，places = 12）。

## 已知缺口与计划阶段

与 `validation-coverage.md` 一致：无未决缺口，阶段 V1。曾出现的 identity 边界问题——译码器不提供优势时期望恰等于随机基线，严格不等断言在边界失败——已按 validation-plan §4 的确定性策略修复为 ≥/≤ 加数值容差（`test_identity_instance_with_weight_two` 即该修复的见证）。

## 相关链接

- 源码：`src/oracq/algorithms/optimization/dqi.py`
- API 参考：[DQI 解码量子干涉优化](../../api/algorithms/optimization/dqi.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_misc_algorithms.py`（misc_algorithms 组），全部在真实后端上执行。

**实验设计**：(a) 植入实例（7 约束 3 变量、右端项由 $x^* = \mathtt{0b101}$ 植入、Dicke 权重 1、穷举译码器），syndrome 分布对照论文的 Krawtchouk 闭式 $K_l(u(x))^2$（`math.comb` 独立计算），优化质量对照经典蛮力枚举的全部 8 个赋值；(b) 抽象译码器经 `bind` 绑定穷举见证后与直接见证逐振幅对拍；(c) Dicke 态 $\lvert D_2^5\rangle$ 的幅度均匀性。后端路径：`reference`、`rir-pysparq`、`adapter-pysparq`、`originir-ext`（实例 10 量子位，四路径分布对拍）。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| dqi-planted-krawtchouk | m = 7、n = 3，10 量子位 | 四路径 | 分布 TVD / 逐点误差 | 4.0e-16 / 6.7e-16 |
| 同上 | — | — | P(error = 0)（译码复净） | 1.0000 |
| 同上 | — | — | 期望满足数（随机基线 3.5） | 6.5000 |
| 同上 | — | — | 概率峰值赋值 = 蛮力最优（满足 7/7） | $x^*$ = 5 ✓ |
| dqi-abstract-decoder-bind | 同上 | rir-pysparq | bind 与直接见证最大振幅偏差 | 0 |
| dicke-state-uniformity | m = 5、l = 2 | reference + originir-ext | 支撑 / 幅度误差 | 10 基态 / 5.6e-17 |

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

产物：`out/verification/misc_algorithms.json`（24 个案例全过，本页对应 `dqi-*` 与 `dicke-state-*` 三个案例）。
