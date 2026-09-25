# Deutsch–Jozsa 查询（Deutsch–Jozsa）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C1 · 模块 [`oracq.algorithms.basics.oracle_algorithms`](../../api/algorithms/basics/oracle_algorithms.rst) · 阶段 V1

## 概述

给定承诺为常量或平衡的布尔函数 $f: \{0,1\}^n \to \{0,1\}$（常量指 $f$ 恒为 0 或恒为 1，平衡指一半输入取 0、一半取 1），用一次 oracle 查询判定属于哪种情形。电路利用相位反冲：把单比特 answer 制备到 $|-\rangle$，使 XOR oracle 的作用等价于相位 $(-1)^{f(x)}$；对均匀叠加的 input 做 Hadamard 干涉后，常量函数的 input 全部塌缩回 $|0^n\rangle$，平衡函数的读出必然非零。经典确定性判定在最坏情形需要 $2^{n-1}+1$ 次查询。

## 接口与输入模型

```python
deutsch_jozsa(function: XorDatabase)
```

API 入口：{obj}`deutsch_jozsa <oracq.algorithms.basics.oracle_algorithms.deutsch_jozsa>`

- `function`：单结果位 XOR database（input model 为 FO：布尔函数真值表），`data_width` 必须为 1，否则在生成期抛出 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`。{obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` 支持三层构造：{obj}`gate_database <oracq.algorithms.input_model.oracles.gate_database>`（门级真值表）、{obj}`qram_database <oracq.algorithms.input_model.oracles.qram_database>`（QRAM 查找）与 {obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>`（开放声明，事后绑定）。

返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`，寄存器为 `input: Bits(address_width)` 与 `answer: Bits(1)`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"deutsch_jozsa"` |
| `readout_register` | `"input"`（判定只读 input，answer 不读） |
| `validation_stage` | `"paradigm"`（见证级别标注） |

## 实现要点

生成序列为：answer 加 X 再加 H（制备 $|-\rangle$）→ input 整体 H（均匀叠加）→ 经 {obj}`invoke <oracq.algorithms.input_model.oracles.invoke>` 以 `function` 前缀调用 oracle 资源 → input 整体 H。oracle 保留为模块调用（带 `function__` 前缀的资源绑定），不在生成或序列化阶段展开。常量/平衡承诺由调用者承担：读出 input 全零判常量，非零判平衡；对不满足承诺的函数输出不作任何保证。

适用边界：只接受 `data_width == 1` 的 XOR database；相位反冲路线要求 answer 初态可制备 $|-\rangle$，电路内部已固定处理，调用方无需（也不应）预处理 answer。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../development/validation-plan.md` §2）：输出分布须与承诺情形的精确构造逐点相等。三层证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects`（模块规范对象与历史导入名的一致性，矩阵口径）。
- 数值：本算法无专属 DJ 数值测试（矩阵口径），但 [Bernstein–Vazirani](bernstein-vazirani.md) 的见证 `AlgorithmExpansionTests.test_bernstein_vazirani_recovers_secret_with_affine_bias` 运行的正是同一电路——{obj}`bernstein_vazirani <oracq.algorithms.basics.oracle_algorithms.bernstein_vazirani>` 直接包装 {obj}`deutsch_jozsa <oracq.algorithms.basics.oracle_algorithms.deutsch_jozsa>` 的产物，在 3 位全部秘密串上以 places = 10 逐分布对拍，DJ 的干涉结构由此间接覆盖。
- 绑定：tests/core 无独立绑定见证（矩阵口径为 —）。`applications/catalog.py` 的 `dj_gate` / `dj_qram` 案例在 L3 层把开放声明 `BooleanFunction` 分别绑到门级真值表与 QRAM，作为绑定路径的目录演示。

## 已知缺口与计划阶段

无已知缺口（验证矩阵缺口列为 —），阶段 V1。

## 相关链接

- 同模块：[Bernstein–Vazirani 秘密串读出](bernstein-vazirani.md)、[Simon 采样](simon.md)
- 源码：`src/oracq/algorithms/basics/oracle_algorithms.py`
- API 参考：[Oracle 查询算法](../../api/algorithms/basics/oracle_algorithms.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值验证脚本：`tests/verification/verify_oracles.py`（C 组案例），产物 `out/verification/oracles.json`。

实验设计：宽度 2–5 bit，每个宽度五个承诺函数——常量 0、常量 1、平衡 parity（$s\cdot x$，$s$ 全 1）、平衡 MSB、平衡阈值函数（非线性），oracle 均为 `gate_database` 真值表。判定指标为 $P(\mathrm{input}=0)$：常量情形须精确为 1、平衡情形须精确为 0；常量与线性平衡两类有闭式末态 $(-1)^c|0^n\rangle|-\rangle$ / $|s\rangle|-\rangle$，另做相位敏感的逐振幅对拍。另设 {obj}`BooleanNetwork <oracq.algorithms.common.arithmetic.BooleanNetwork>` 编译 oracle 的端到端案例（parity3 / const3）：网络经 `operation()` 编译后包装为 XOR database 接口，先在叠加态下对 `net.evaluate` 的 16 个分支逐振幅穷举，再做 DJ 判定（见[布尔网络](boolean-networks.md)数值验证节）。后端路径：reference、rir-pysparq、originir-ext。

| 案例 | 规模 | 路径 | 指标值 |
|---|---|---|---|
| dj-decision-w2 | 2 bit，5 函数 | reference, rir-pysparq, originir-ext | p_zero_max_error = 7.8e-16，decision_errors = 0，闭式 max_error = 2.2e-16 |
| dj-decision-w3 | 3 bit，5 函数 | 同上 | p_zero_max_error = 1.1e-15，decision_errors = 0，闭式 max_error = 3.3e-16 |
| dj-decision-w4 | 4 bit，5 函数 | 同上 | p_zero_max_error = 1.4e-15，decision_errors = 0，闭式 max_error = 4.4e-16 |
| dj-decision-w5 | 5 bit，5 函数 | 同上 | p_zero_max_error = 1.8e-15，decision_errors = 0，闭式 max_error = 5.6e-16 |
| dj-boolean-network-parity3 | 3 bit 平衡（网络编译） | rir-pysparq, originir-ext | max_error = 8.3e-17，p_zero_error = 0.0 |
| dj-boolean-network-const3 | 3 bit 常量（网络编译） | rir-pysparq, originir-ext | max_error = 8.3e-17，p_zero_error = 1.1e-15 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_oracles.py
```

产物路径：`out/verification/oracles.json`。
