# Bernstein–Vazirani 秘密串读出（Bernstein–Vazirani）

> 类别 C1 · 模块 `pyqecclang.algorithms.oracle_algorithms` · 阶段 V1

## 概述

给定承诺为仿射布尔函数 $f(x) = s \cdot x \oplus c$（$s \in \{0,1\}^n$ 为秘密串，$c \in \{0,1\}$ 为偏置，点积与异或在 GF(2) 上进行），用一次 oracle 查询完整读出 $s$。与 [Deutsch–Jozsa](deutsch-jozsa.md) 相同的相位反冲电路上，相位 $(-1)^{f(x)} = (-1)^c (-1)^{s \cdot x}$ 的线性部分在 Hadamard 干涉后把 input 塌缩到 $|s\rangle$，常数因子 $(-1)^c$ 不影响读出。

## 接口与输入模型

```python
bernstein_vazirani(function)
affine_boolean_oracle(width, secret, *, bias=0)
```

- `bernstein_vazirani(function)`：`function` 为单结果位 XOR database（input model 为 FO：仿射布尔函数真值表），内部先经 `deutsch_jozsa` 校验 `data_width == 1`；仿射承诺由调用者承担。oracle 可以保留为开放声明（`abstract_database`），生成后再绑定 gate 或 QRAM 实现。
- `affine_boolean_oracle(width, secret, *, bias=0)`：构造 $f(x) = s \cdot x \oplus c$ 的门级 oracle，返回 `XorDatabase`。`width` 范围 1..64，`secret` 范围 $0 .. 2^{\text{width}}-1$，`bias` 取 0 或 1；`secret` 的第 $i$ 位对应地址第 $i$ 位（little endian）。

`bernstein_vazirani` 返回 `Operation`，寄存器与 DJ 相同：`input: Bits(address_width)` 与 `answer: Bits(1)`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"bernstein_vazirani"` |
| `input_promise` | `"f(x)=dot(s,x) XOR c over GF(2)"`（承诺声明） |
| `readout_register` | `"input"`（读取 input 直接得到 $s$） |
| `validation_stage` | `"paradigm"`（自 DJ 电路继承） |

## 实现要点

实现复用 DJ 电路：调用 `deutsch_jozsa(function)` 后，仅以 `dataclasses.replace` 重标注模块名与属性（`algorithm` / `input_promise`），依赖与电路结构完全不变。读出满足承诺时 `input` 精确等于 $s$；对非仿射函数输出不作保证。oracle 保留为模块调用；开放声明路线（`abstract_database` → `bind`）在运行展示目录中有端到端示例：先以 `bernstein_vazirani(abstract_database(...))` 生成含未绑定槽位的程序，再用具体 oracle 的 `operation` 绑定。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../development/validation-plan.md` §2）：读出分布须逐点等于 $|s\rangle$。三层证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects`（矩阵口径）。
- 数值：`AlgorithmExpansionTests.test_bernstein_vazirani_recovers_secret_with_affine_bias`——`width = 3` 下全部 8 个秘密串 × 偏置 $\{0, 1\}$ 共 16 个组合，input 分布集中于 $s$ 的断言精度 places = 10。
- 绑定：tests/core 无独立绑定见证（矩阵口径为 —）；开放声明 → 门级绑定的用法由运行展示目录与 `applications/catalog.py` 的 `dj_gate` / `dj_qram` 案例（同一 XOR database 绑定机制）演示。

## 已知缺口与计划阶段

无已知缺口（验证矩阵缺口列为 —），阶段 V1。

## 相关链接

- 同模块：[Deutsch–Jozsa 查询](deutsch-jozsa.md)、[Simon 采样](simon.md)
- 源码：`src/pyqecclang/algorithms/oracle_algorithms.py`
- API 参考：[Oracle 查询算法](../../api/algorithms/oracle_algorithms.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
