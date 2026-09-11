# Deutsch–Jozsa 查询（Deutsch–Jozsa）

> 类别 C1 · 模块 `pyqecclang.algorithms.oracle_algorithms` · 阶段 V1

## 概述

给定承诺为常量或平衡的布尔函数 $f: \{0,1\}^n \to \{0,1\}$（常量指 $f$ 恒为 0 或恒为 1，平衡指一半输入取 0、一半取 1），用一次 oracle 查询判定属于哪种情形。电路利用相位反冲：把单比特 answer 制备到 $|-\rangle$，使 XOR oracle 的作用等价于相位 $(-1)^{f(x)}$；对均匀叠加的 input 做 Hadamard 干涉后，常量函数的 input 全部塌缩回 $|0^n\rangle$，平衡函数的读出必然非零。经典确定性判定在最坏情形需要 $2^{n-1}+1$ 次查询。

## 接口与输入模型

```python
deutsch_jozsa(function: XorDatabase)
```

- `function`：单结果位 XOR database（input model 为 FO：布尔函数真值表），`data_width` 必须为 1，否则在生成期抛出 `ValidationError`。`XorDatabase` 支持三层构造：`gate_database`（门级真值表）、`qram_database`（QRAM 查找）与 `abstract_database`（开放声明，事后绑定）。

返回 `Operation`，寄存器为 `input: Bits(address_width)` 与 `answer: Bits(1)`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"deutsch_jozsa"` |
| `readout_register` | `"input"`（判定只读 input，answer 不读） |
| `validation_stage` | `"paradigm"`（见证级别标注） |

## 实现要点

生成序列为：answer 加 X 再加 H（制备 $|-\rangle$）→ input 整体 H（均匀叠加）→ 经 `invoke` 以 `function` 前缀调用 oracle 资源 → input 整体 H。oracle 保留为模块调用（带 `function__` 前缀的资源绑定），不在生成或序列化阶段展开。常量/平衡承诺由调用者承担：读出 input 全零判常量，非零判平衡；对不满足承诺的函数输出不作任何保证。

适用边界：只接受 `data_width == 1` 的 XOR database；相位反冲路线要求 answer 初态可制备 $|-\rangle$，电路内部已固定处理，调用方无需（也不应）预处理 answer。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../development/validation-plan.md` §2）：输出分布须与承诺情形的精确构造逐点相等。三层证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects`（模块规范对象与历史导入名的一致性，矩阵口径）。
- 数值：本算法无专属 DJ 数值测试（矩阵口径），但 [Bernstein–Vazirani](bernstein-vazirani.md) 的见证 `AlgorithmExpansionTests.test_bernstein_vazirani_recovers_secret_with_affine_bias` 运行的正是同一电路——`bernstein_vazirani` 直接包装 `deutsch_jozsa` 的产物，在 3 位全部秘密串上以 places = 10 逐分布对拍，DJ 的干涉结构由此间接覆盖。
- 绑定：tests/core 无独立绑定见证（矩阵口径为 —）。`applications/catalog.py` 的 `dj_gate` / `dj_qram` 案例在 L3 层把开放声明 `BooleanFunction` 分别绑到门级真值表与 QRAM，作为绑定路径的目录演示。

## 已知缺口与计划阶段

无已知缺口（验证矩阵缺口列为 —），阶段 V1。

## 相关链接

- 同模块：[Bernstein–Vazirani 秘密串读出](bernstein-vazirani.md)、[Simon 采样](simon.md)
- 源码：`src/pyqecclang/algorithms/oracle_algorithms.py`
- API 参考：[Oracle 查询算法](../../api/algorithms/oracle_algorithms.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
