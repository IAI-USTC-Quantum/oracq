# Simon 采样（Simon Sampling）

> 类别 C1 · 模块 `pyqecclang.algorithms.oracle_algorithms` · 阶段 V1

## 概述

给定承诺为带非零 XOR 周期 $s$ 的二对一函数 $f: \{0,1\}^n \to \{0,1\}^m$（即 $f(x) = f(y)$ 当且仅当 $y = x \oplus s$），求周期 $s$。量子侧一次采样给出满足 $y \cdot s = 0$（GF(2) 点积）的均匀随机向量 $y$：电路对 input 做 Hadamard 后查询 oracle，output 上纠缠掉函数值，再对 input 做 Hadamard，input 即塌缩到与 $s$ 正交的约束上。收集约 $n - 1$ 个线性无关样本后在经典侧做 GF(2) 消元——本仓库与算法目录的口径一致，消元放在经典侧，量子线路只负责采样。

## 接口与输入模型

```python
simon_sample(function)
simon_nullspace(samples, width)
```

- `simon_sample(function)`：`function` 必须是 `XorDatabase` 实例（input model 为 FO：二对一函数真值表），其他类型在生成期抛出 `ValidationError`；二对一承诺由属性 `input_promise` 声明、由调用者承担。
- `simon_nullspace(samples, width)`：纯经典后处理。`samples` 为按 little endian 编码的整数样本列表，`width`（1..64）为秘密串位宽，返回样本行空间在 GF(2) 上零空间的基向量元组。样本恰好撑满 $s^\perp$ 时零空间一维，其非零元即周期 $s$；样本不足时保留多个基向量，调用方继续采样。

`simon_sample` 返回 `Operation`，寄存器为 `input: Bits(address_width)` 与 `output: Bits(data_width)`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"simon_sample"` |
| `input_promise` | `"two-to-one XOR period"`（二对一周期承诺） |
| `readout_register` | `"input"`（只读 input 即得一个正交约束） |

## 实现要点

生成序列为：input 整体 H → 经 `invoke` 以 `function` 前缀调用 oracle → input 整体 H。output 寄存器无需在电路中提前测量——读出 input 已经是完整的一次采样结果，output 上的函数值纠缠不进入判定。oracle 保留为模块调用。消元实现为逐位选主元的高斯消去：对每一位找当前行中该位为 1 的行作主元并交换，再对其他行做异或清位；自由变量逐个回代得到零空间基。

适用边界：对不满足二对一承诺的函数，读出分布与消元结果均不作保证；`simon_nullspace` 输入样本越界（不在 $0 .. 2^{\text{width}}-1$）时在生成期抛出 `ValidationError`。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../development/validation-plan.md` §2）：采样分布与消元结果须与闭式构造逐点相等。三层证据：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects`（矩阵口径）。
- 数值：`AlgorithmExpansionTests.test_simon_sampling_and_rank_aware_postprocess`——对周期 $s = 3$ 的二对一函数 `gate_database(2, 1, [0, 1, 1, 0])`，input 分布逐点对拍 $p(0) = p(3) = 0.5$；消元三组闭式对拍：`simon_nullspace([3], 2) == (3,)`（撑满零空间）、`simon_nullspace([], 3)` 给出 3 个基向量（无样本）与 `simon_nullspace([1, 2, 4], 3) == ()`（样本已满秩，零空间只有零向量）。
- 绑定：本算法无独立绑定见证（矩阵口径为 —）。

## 已知缺口与计划阶段

无已知缺口（验证矩阵缺口列为 —），阶段 V1。

## 相关链接

- 同模块：[Deutsch–Jozsa 查询](deutsch-jozsa.md)、[Bernstein–Vazirani 秘密串读出](bernstein-vazirani.md)
- 源码：`src/pyqecclang/algorithms/oracle_algorithms.py`
- API 参考：[Oracle 查询算法](../../api/algorithms/oracle_algorithms.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
