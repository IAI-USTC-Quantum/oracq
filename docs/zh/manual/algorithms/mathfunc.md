# 数学函数可逆编译（Math Function Compilation）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C5 · 模块 [`oracq.infrastructure.mathfunc`](../../api/infrastructure/mathfunc.rst) · 阶段 V2

## 概述

把纯 Python 数学函数（`math` / `cmath` 子集）编译为定点可逆量子线路（input model 词汇中的 FO，函数 oracle）。公开契约是 `reversible_function` 范式：输入保持不变，输出与 status 按位 XOR 更新，所有临时寄存器自动复净——线路实现的是有限字长下的函数 $\tilde F$，不声称精确表示连续数学函数。编译链为 Python AST → MIR 0.1（有类型 SSA 图）→ 模块化 RIR；初等函数用 Chebyshev–Clenshaw 多项式在配置区间上近似。使用教程见[普通数学函数自动生成](../math-functions.md)，MIR 格式见 [math-ir](../../reference/math-ir.md)。

## 接口与输入模型

```python
compile_function(function, *, fmt=None, inputs=None, constants=None, helpers=None,
                 output_names=None, config=None, max_unroll=128, entry=None)
lower_math_ir(program, *, fmt=None, config=None, output_names=None)
```

API 入口：{obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>`、{obj}`lower_math_ir <oracq.infrastructure.mathfunc.lower_math_ir>`

- `function`：普通函数对象或 `def` 源码字符串；源码只允许 `math` / `cmath` 导入与纯函数定义（模块 docstring 与 `__future__` 导入被忽略），`entry` 指定多函数源码的入口（缺省最后一个）。
- `inputs`：`{参数名: "real" | "complex" | "bool" | Index(width)}`；缺省按注解推断（`float` / `int` → real、`complex` → complex、`bool` → bool）。
- `constants`：生成期常量（含默认参数）；`helpers`：辅助函数命名空间（不支持递归）；`max_unroll`：静态 `range` 循环展开上限。
- `fmt`：{obj}`FixedFormat <oracq.algorithms.common.arithmetic.FixedFormat>`（缺省 `(12, 6)`，必须有符号）；`config`：{obj}`MathConfig(degree=6, intervals=()) <oracq.infrastructure.mathfunc.numeric.MathConfig>`，degree 为数学核多项式阶数（1..32），intervals 逐函数替换近似区间。
- {obj}`lower_math_ir <oracq.infrastructure.mathfunc.lower_math_ir>`：由 MIR 单独降低（可换 fmt / config），供 JSON 往返后重建。

返回 {obj}`CompiledFunction <oracq.infrastructure.mathfunc.lowering.CompiledFunction>`（`oracq.compile_function` 顶层再导出）。属性：

| 属性 | 含义 |
|---|---|
| `operation` / `program()` | 生成的 {obj}`Operation <oracq.infrastructure.builder.Operation>` / 完整 RIR 程序 |
| `math_ir` | 入口 {obj}`MathProgram <oracq.infrastructure.mathfunc.graph.MathProgram>`（含全部 helper 函数） |
| `fmt` | 定点格式回显 |
| `input_layout` / `output_layout` | 参数到物理寄存器的展开（复数拆 `_real` / `_imag`，bool 1 位，{obj}`Index <oracq.infrastructure.mathfunc.graph.Index>` `width` 位） |

寄存器为各输入、输出（缺省 `out` / `out_i`，可由 `output_names` 命名）加 `status`（2 位）。模块属性：`oracle_paradigm="reversible_function"`、`math_function`、`math_ir_version="0.1"`、`fixed_width` / `fixed_fraction`、`math_config`、`update_semantics`、`correctness="pending"`；入口模块附完整 `math_ir` JSON。

## 实现要点

前端**不执行**源码：副作用语句、`while`、递归、lambda 在编译期拒绝；静态 `range` 循环展开为直线上；`if` / 三元表达式 / 布尔短路编译为 `select` 节点——两分支电路都被计算，未走分支的算术错误经 status 掩码（如 `0 if x == 0 else 1/x` 在 x = 0 处不置位）。捕获的数值型全局名字与 `constants` 都成为生成期常量。

MIR 是顺序 SSA 图（节点引用的下标严格小于自身），`MathProgram.validate()` 校验操作元数与类型约束，{obj}`dumps <oracq.infrastructure.serialization.dumps>` / {obj}`loads <oracq.infrastructure.serialization.loads>` JSON 往返。降低器 {obj}`NumericEmitter <oracq.infrastructure.mathfunc.numeric.NumericEmitter>` 全部使用局部寄存器并在出口 XOR 输出后伴随重放整帧，实现自动 compute / XOR / uncompute；重叠的只读参数按无别名调用 ABI 复制。复数四则与初等函数按恒等式分解为实核（如 $e^{x+iy} = e^x(\cos y + i\sin y)$，$\arctan$ 系经 log / sqrt 恒等式）；整数指数幂用平方乘（上限 128），一般幂走 $\exp(b\log a)$；`atan2` 由幅值比加象限修正合成。

初等核 {obj}`elementary_kernel <oracq.infrastructure.mathfunc.numeric.elementary_kernel>`：在区间上取 degree+1 个 Chebyshev 节点经典采样、量子侧 Clenshaw 递推求值，配方（区间、阶数、系数）写入模块属性 `math_approximation`；越出近似区间置 status 位 1，`log` / `asin` / `atanh` 等的定义域违例置位 0（位 0 = 定义域失效、位 1 = 值域 / 字长越界）。核按 `(函数名, fmt, config)` 缓存复用。{obj}`Index <oracq.infrastructure.mathfunc.graph.Index>` 输入按无符号整数嵌入定点字的 fraction 偏移处（要求 `width + fraction < fmt.width`）。适用边界：无自动区间推导、误差证明与最优算术电路选择，近似区间必须覆盖各个中间数学核的输入。

## 验证方案

类别 C5（数据访问层，判定准则见 `../development/validation-plan.md` §2）：查询语义逐点正确 + 输入保留 / XOR 更新 / 复净。三层证据位于 `tests/core/test_mathfunc.py:MathFunctionTests`：

- 结构：`test_helper_module_reuse_and_roundtrip`——helper 去重（同一调用点只留一个 call 节点）、MIR JSON 往返、按重建 MIR 降低得到的 RIR 与直接编译逐字节一致；`test_source_is_not_executed_and_effects_rejected`（副作用 / 递归 / `while` / 文件操作被拒且源码不被执行）。
- 数值：`test_all_elementary_families_construct`——`math` / `cmath` 的全部 16 个实初等族在实 / 复输入下构造出闭合程序，`polar` / `rect` / `conjugate` 多返回值路径亦覆盖；`test_roe_is_a_compiled_pure_function`——Euler ROE 面通量整链编译，helper 族标签齐全。另有 `test_nonzero_output_xor_superposition_and_inverse`（叠加输入与非零输出初值下的 XOR 语义、双调用自逆）与 `test_unused_branch_error_is_masked`（未走分支的除零被 status 掩码）。
- 绑定：无独立绑定见证——编译产物是闭合模块、不含抽象槽位，与 `validation-coverage.md` 的 mathfunc 行一致。

## 已知缺口与计划阶段

定点量化误差界的显式见证未补：字长、阶数、区间与系数已进模块属性，但"给定 fmt / degree 下输出偏差不超过多少"没有自动化断言。登记为阶段 V2，与 `validation-coverage.md` 的 mathfunc 行一致。

## 数值验证

论文级数值实验见 `tests/verification/verify_mathfunc.py`（真实后端执行，无模拟替身），直接补齐上一节的误差界缺口。实验设计：覆盖 `examples/math_functions.py` 全部 5 个函数与 `applications/roe_formulas.frozen_roe_face`（6 实输入 + 2 个 Index(2) + 2 个生成期常量 + 双输出），FixedFormat(6,2)/(8,3) 两种格式；rir-pysparq 叠加穷举（单输入全域 2^6/2^8 分支，多输入逐轴纤维 + 联合立方体，Index 全域 16 组联合）。oracle 为按 {obj}`fixed_arithmetic <oracq.algorithms.common.arithmetic.fixed_arithmetic>` 语义独立实现的定点仿真机 `_Fx`（mul/div/sqrt 幅度向零截断、add/sub 模 wrap、status 位 0 = 定义域失效、位 1 = 值域/字长越界、`select` 掩码未选分支而 helper 调用合并全部参数旗标、helper 形参不参与折叠），初等核另按 `math_approximation` 系数逐比特复现 Clenshaw 递推——因此"实现误差"对照的是**量化后的系数配方**，"方法误差"（配方 vs 真函数）单独报告。结果：除 guarded_reciprocal 外全部案例与 _Fx 仿真**逐比特一致**（max_error = 0，覆盖数千个无旗标分支），guarded_reciprocal 全域误差 ≤ 1 量子（除法截断界）且全域无旗标；status 旗标（rho=0 除零、rho<0 平方根、c2≤0、z=±i、核区间越界、中间量字长越界含 wrap 后除零）逐分支符合预测（misflagged = 0）；三后端（rir-pysparq / reference / adapter-pysparq）振幅两两偏差 0.0；基态确定性（输入保持、输出 XOR、locals 复净）6 函数全过。工作区实测 358–2178 比特，OriginIR-ext 24 比特预算不适用。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| polynomial 全域穷举 | 6.2 / 8.3，64 / 256 分支 | rir-pysparq | max_error | 0 / 0（exact_fraction 1.0） |
| guarded_reciprocal 全域穷举 | 6.2 / 8.3，64 / 256 分支 | rir-pysparq | max_error | 0.94 / 0.98 量子（全域 status = 0） |
| pressure 纤维 + 立方体 | 两格式，192+64 / 768+64 分支 | rir-pysparq | max_error / 常量量化偏差 | 0 / 1.4（6.2）、0 / 0.578（8.3） |
| roe_speed 纤维 + 立方体 | 两格式，256+16 / 1024+16 分支 | rir-pysparq | max_error | 0 / 0 |
| phase_response 纤维 + 立方体 | 两格式，192+16 / 512+16 分支 | rir-pysparq | max_error / 方法误差 | 0 / 1.770（6.2）、0 / 1.349（8.3） |
| frozen_roe_face 索引联合 + 纤维 | 16 + 96 分支 × 两格式 | rir-pysparq | max_error | 0 / 0（left/right 双输出） |
| 核方法误差（信息性） | 4001 点网格，degree=3 | 经典对照 | exp / sin / cos | 0.148 / 0.195 / 0.364 |
| 跨后端对拍（6 案例） | 16–64 分支 | rir-pysparq / reference / adapter-pysparq | max_pairwise_deviation | 0.0 |
| 基态确定性 | 6 函数 | rir-pysparq | failures | 0 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_mathfunc.py
```

产物：`out/verification/mathfunc.json`（28 个案例全过，总运行约 185 秒）。

## 相关链接

- 源码：`src/oracq/infrastructure/mathfunc/`（`frontend.py` / `graph.py` / `numeric.py` / `lowering.py`）
- 使用教程：[普通数学函数自动生成可逆量子模块](../math-functions.md)
- API 参考：[数学函数编译入口](../../api/infrastructure/mathfunc.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
