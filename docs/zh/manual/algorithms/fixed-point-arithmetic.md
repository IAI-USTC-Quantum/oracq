# 可逆定点算术（Fixed-Point Arithmetic）

[English](../../../index.html) · **简体中文**

> 类别 C1 · 模块 [`oracq.algorithms.common.arithmetic`](../../api/algorithms/common/arithmetic.rst) · 阶段 V1

## 概述

定点格式下的可逆算术族：加、减、乘（以及取负、绝对值、除法、倒数、开方、比较与按位逻辑等）。定点数由 {obj}`FixedFormat(width, fraction, signed) <oracq.algorithms.common.arithmetic.FixedFormat>` 描述——`width` 位字长、低 `fraction` 位为小数、`signed` 时最高位为补码符号位。输出采用 XOR 语义并附带两位 `status` 标志：`status[0]` 标记定义域失效（除零、负数开方），`status[1]` 标记结果越出字长（而非精度界内的舍入）；取整口径为向零舍入加模回绕（`rounding = "toward_zero; modular_wrap"`）。线路由[布尔网络](boolean-networks.md)的 compute/XOR/uncompute 模式生成，全程无脏工作区。

## 接口与输入模型

```python
FixedFormat(width=8, fraction=3, signed=True)
fixed_arithmetic(kind, fmt=DEFAULT_FIXED_FORMAT)
```

API 入口：{obj}`FixedFormat <oracq.algorithms.common.arithmetic.FixedFormat>`、{obj}`fixed_arithmetic <oracq.algorithms.common.arithmetic.fixed_arithmetic>`

- `FixedFormat`：要求 $2 \le \text{width} \le 64$ 且 $0 \le \text{fraction} < \text{width} - \text{signed}$（小数位不占用符号位），违例抛 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`；{obj}`encode(value) <oracq.infrastructure.serialization.encode>` / {obj}`decode(value) <oracq.infrastructure.serialization.decode>` 在实数与补码整数间换算。`DEFAULT_FIXED_FORMAT` 即默认的 (8, 3, 有符号)。
- {obj}`fixed_arithmetic(kind, fmt) <oracq.algorithms.common.arithmetic.fixed_arithmetic>`：`kind` 取 `add` / `sub` / `neg` / `abs` / `mul` / `div` / `reciprocal` / `sqrt` / `lt` / `eq` / `select` / `and` / `or` / `xor` 之一，未知值抛 `ValidationError`；结果按 `(kind, fmt)` 以 `lru_cache` 缓存复用。

返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`，寄存器为输入 `a`（一元 `kind` 只有 `a`；`select` 另有 1 位 `select`）加输出 `out`（`lt` / `eq` 为 1 位，其余 `width` 位）与 `status: Bits(2)`。本入口不消费 oracle 输入（CP：定点格式与操作种类直接参数化）。模块属性：

| 属性 | 含义 |
|---|---|
| `arithmetic_kind` / `fixed_width` / `fixed_fraction` / `fixed_signed` | 调用参数回显 |
| `rounding` | `"toward_zero; modular_wrap"` |
| `arithmetic_network` | 布尔网络 payload（原生执行入口使用） |

## 实现要点

生成策略是先在布尔 SSA 上搭建经典组合逻辑（见[布尔网络](boolean-networks.md)），再整体编译为可逆线路。有符号乘除走幅度域路线：两侧取绝对值，在加宽的幅度域完成乘（$2n$ 位移位累加后取 `[fraction:]` 段）或带恢复余数的除法，再按异或出来的符号位回贴符号并取负；`status[0]` 由分母是否全零（或开方输入符号位）驱动，失效时输出 mux 到全零。溢出判据 `status[1]`：无符号加减用进位/借位，有符号加减用符号一致性判据（同号相加变号 / 异号相减变号），乘除与开方检查幅度是否侵入符号位。除法与倒数共用 `1 / b` 的除法核（倒数的分子固定为 $2^{2f}$）。

适用边界：结果在字长内截断并模回绕，`status` 位是调用方判断溢出的唯一途径；精度由 `fraction` 事先决定，线路内无增长字长的机制。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../development/validation-plan.md` §2）：作用算子须与经典求值逐点相等。三层证据：

- 结构：`tests/core/test_stage2.py:Stage2StructureTests`——`test_arithmetic_construction_and_basis` 对全部 14 种 `kind` × `FixedFormat(4, 1)` 逐一检查：模块含私有工作区、程序 JSON 往返相等、Toffoli/U3/CZ 基导出只含三种基门且无 `controlled_by`。
- 数值：`Stage2StructureTests.test_gate_network_is_executable_small_example`——无符号 2 位 `add` 在 `a = 1, b = 2` 上模拟，幅度精确集中在 `(a, b, out, status) = (1, 2, 3, 0)`；同一模块的 `arithmetic_network` payload 经 `from_payload` 复原后经典求值 `evaluate(a=1, b=2) == {"out": 3, "status": 0}`，电路与经典语义对拍。其余 `kind` 由结构见证与 payload 经典求值一致性覆盖，电路级数值对拍集中在 `add` 小实例。
- 绑定：本算法无独立绑定见证（矩阵口径为 —）。真实后端上的原生执行（{obj}`arithmetic_native_registry <oracq.algorithms.common.arithmetic.arithmetic_native_registry>`）另由 `tests/integration/test_stage2_native.py:Stage2NativeTests` 在 PySparQ 上与参考模拟器对拍（L4 冒烟，非矩阵口径）。

## 数值验证

实验设计（`tests/verification/verify_arithmetic.py`，产物 `out/verification/arithmetic.json`，20 案例）：原语算术走「OriginIR-ext 全振幅（1–8 bit）+ UniQC `Circuit.to_matrix` 幺正（3–4 bit）+ PySparQ RIR 宽寄存器（16/32/64 bit）」三条路径；定点编译函数（{obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>` 的 Boolean SSA 降低）用叠加态一次穷举全部输入编码，对照独立经典语义逐分支核对；初等函数把「实现误差（对照模块属性 `math_approximation` 中的 Chebyshev 系数）」与「方法误差（系数多项式 vs 真函数）」分开报告。复现：`PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_arithmetic.py`。

| 案例 | 规模 | 路径 | 指标 |
|---|---|---|---|
| `add_const` 叠加穷举 | w=1..8 × 4 常量 | originir-ext + reference/rir/adapter 对拍 | max_error ≤ 1.2e-16 |
| `add_const` 幺正 | w=3/4 | originir-ext + to_matrix | 与经典置换矩阵逐元素差 = 0.0 |
| `add_const` 宽寄存器 | w=16/32/64 抽样基态 | rir-pysparq | failures = 0；w=12 全叠加 max_error = 1.6e-17 |
| xor/swap/受控算术 | 4+4+2 bit 结构化程序 | 四路径两两对拍 | max_pairwise_deviation = 0.0 |
| 编译多项式 `x*x+0.5` | FixedFormat(4,1)/(6,2)/(8,3) 全域 | rir-pysparq | 无旗标输出误差 ≤ 1 量子（判据 2 量子）；status 越界集合与预测一致 |
| 编译 `sin(x)`（degree=3） | FixedFormat(8,4)，区间 [−1,1] | rir-pysparq | 实现误差 0.055 ≤ 2 量子；方法误差 0.0453（信息性）；misflagged = 0 |
| 编译函数叠加对拍 | FixedFormat(4,1) 全输入 | 四路径 | max_pairwise_deviation = 0.0 |

## 已知缺口与计划阶段

无已知缺口（验证矩阵缺口列为 —），阶段 V1。

## 相关链接

- 同模块：[布尔网络](boolean-networks.md)、[Fourier 加法](fourier-addition.md)
- 源码：`src/oracq/algorithms/common/arithmetic.py`
- API 参考：[可逆算术](../../api/algorithms/common/arithmetic.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
