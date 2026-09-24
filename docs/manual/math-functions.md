# 普通数学函数自动生成可逆量子模块

先编写纯 Python 数学函数，再调用 {obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>` 生成可逆量子模块。函数仍可用于经典计算；量子侧由编译器处理临时寄存器、别名复制、结果 XOR 和反算。[QFVM](qfvm.md) 的 Roe face 已采用这条路径。

```python
from oracq import FixedFormat, compile_function, export_toffoli_u3_cz

def pressure(rho, momentum, energy, gamma=1.4):
    velocity = momentum / rho
    return (gamma - 1) * (energy - 0.5 * momentum * velocity)

compiled = compile_function(
    pressure,
    fmt=FixedFormat(width=12, fraction=6),
    constants={"gamma": 1.4},
)
operation = compiled.operation
program = compiled.program()
originir = export_toffoli_u3_cz(program).text
```

该模块的接口为 rho、momentum、energy、out、status；gamma 是生成期参数。调用方式与其他 {obj}`Operation <oracq.infrastructure.builder.Operation>` 相同：

```python
from oracq import Builder, Bits

b = Builder("flow_pressure", {
    "rho": Bits(12), "momentum": Bits(12), "energy": Bits(12),
    "value": Bits(12), "flags": Bits(2),
})
# 输入通常先来自 QRAM。下面的 Call 不会在 RIR 中自动展开。
b.call(operation, rho=b["rho"], momentum=b["momentum"],
       energy=b["energy"], out=b["value"], status=b["flags"])
application = b.finish().program()
```

输入保持不变，结果 XOR 写入目标寄存器，是本功能的公开契约。它计算有限字长下的函数实现 F_tilde；不会声称有限电路精确表示任意连续数学函数。

## 已实现的数学范围

- 实数：四则运算、负号、绝对值、比较、条件选择、整数常量幂及一般幂；现有恢复除法和逐双位平方根电路复用。
- 实函数：sqrt、exp、log、log10、sin、cos、tan、asin、acos、atan、sinh、cosh、tanh、asinh、acosh、atanh。
- 复数：加减乘除、绝对值、实部/虚部、conjugate，以及上述 cmath 函数的复数分解；另有 phase、polar、rect。
- 补充分解：atan2 的幅值比与象限选择、hypot、多项式 Clenshaw 递推、复数代数组合、按分支选择的状态电路。
- 源码组合：纯 helper、静态有界 range、多个返回值和静态 tuple。helper 保留独立 MIR/RIR 模块。

```python
import cmath
from oracq import MathConfig

def response(z: complex):
    return cmath.exp(1j * z) / (1 + z * z)

compiled = compile_function(
    response,
    fmt=FixedFormat(16, 8),
    config=MathConfig(
        degree=8,
        intervals=(("exp", -2.0, 2.0), ("sin", -3.14159, 3.14159)),
    ),
)
print(compiled.input_layout)   # z_real, z_imag
print(compiled.output_layout)  # out_real, out_imag
```

普通无注解参数默认实数；可用 inputs 显式指定 real/complex/bool。{obj}`Index <oracq.infrastructure.mathfunc.graph.Index>`(width) 为 QFVM 行列索引等无符号整数提供较短的公开寄存器，进入计算时自动转换到定点表示。常量默认参数视为生成期参数；显式放入 inputs 后也可成为量子输入。多结果可用 output_names 指定名称，例如 Roe 的 left/right。

Python 的 cmath 在分支切线上区分有符号零；当前定点编码没有该信息。因此这里提供函数族和可替换的近似实现，不承诺完整浮点/分支兼容。[对应的 cmath 行为见官方文档](https://docs.python.org/3/library/cmath.html)。

## 近似与状态

{obj}`MathConfig <oracq.infrastructure.mathfunc.numeric.MathConfig>`.degree 控制实函数的 Chebyshev 多项式阶数；intervals 可替换每个数学核的区间。默认区间见 [numeric.py](../api/infrastructure/mathfunc/numeric.rst) 的 BOUNDS。配置必须覆盖各个**中间数学核输入**，不只是最外层函数的输入。状态会标记超出区间或发生定义域/字长问题的路径。当前没有自动区间推导、误差证明或最优算术电路选择。

函数系数通过 degree+1 个经典采样点生成，然后执行量子乘加递推，没有预先枚举每个输入对应的函数值。字长、近似阶数、区间和系数进入生成配置/模块属性，eps 仍不进入语言核心。

```python
def safe_inverse(x):
    return 0.0 if x == 0 else 1.0 / x
```

这个函数在 x=0 时不会因未选中分支的除零而设置最终 domain 标志。两分支电路可以被计算，但其结果和状态按量子条件选择，所有中间位最后反算。

## 中间表示与后端

完整链路为 Python 纯函数 → [MIR 0.1](../reference/math-ir.md) → RIR 0.3 → 模块化 OriginIR-ext / PySparQ。MIR 可独立 JSON 往返，再由 {obj}`lower_math_ir <oracq.infrastructure.mathfunc.lower_math_ir>` 按其他配置降低。最终 RIR 不含 Python callback；带数学核和 helper 的调用继续保留为 {obj}`Module <oracq.infrastructure.ir.Module>`/{obj}`Call <oracq.infrastructure.ir.Call>`。

现有 {obj}`arithmetic_native_registry <oracq.algorithms.common.arithmetic.arithmetic_native_registry>` 会识别生成模块内部的算术/布尔实现。PySparQ 在这些模块边界执行真实自定义 C++ 算子，跳过内部 Boolean 工作区。原生路径与门级路径使用同一组算术网络，并已通过实际执行检查。这里没有用 Python 原函数的直接求值冒充量子模拟。

```python
from oracq import arithmetic_native_registry, run_pysparq
state = run_pysparq(
    application,
    native_registry=arithmetic_native_registry(application),
    memory=memory,  # 若入口声明了 QRAM
)
```

QFVM 仍查询原始守恒量；[roe_formulas.py](../api/applications/roe_formulas.rst) 是普通经典公式，{obj}`roe_face <oracq.applications.roe.roe_face>` 通过 compile_function 生成原 ABI 的模块。经典 Riemann 更新和 QRAM 数据结构继续与量子矩阵元计算分离。

## 命令行与案例

```bash
oracq compile-function examples/math_functions.py \
  --function pressure --width 12 --fraction 6 \
  --constants '{"gamma": 1.4}' \
  --mir-output out/pressure.mir.json -o out/pressure.rir.yaml

oracq emit out/pressure.rir.yaml --basis toffoli-u3-cz -o out/pressure.originir

PYTHONPATH=src .venv/bin/python tools/build_math_functions.py
```

--inputs 接收 JSON 类型映射，例如 '{"z":"complex"}' 或 '{"row":{"index":2},"x":"real"}'。不支持的源语句会给出 {obj}`FunctionCompileError <oracq.infrastructure.mathfunc.frontend.FunctionCompileError>`；这不是可以编译任意 Python 程序的工具。

out/math-functions/ 包含 pressure、roe_speed、phase_response、guarded_reciprocal、polynomial、自动 Roe face 与接入后的 QFVM，共七组产物。见 [实施面板](../development/contributing.md) 与 [验证记录](../archive/function-compiler-validation.json)。

数学精度、复杂分支切线、Roe 数值结果和资源优化仍待核验。现有测试覆盖编译、类型、可逆更新、模块复用和真实后端消费。

## 数值验证

论文级数值实验见 `tests/verification/verify_mathfunc.py`（真实后端执行，无模拟替身）。实验设计：本页 5 个函数（pressure、roe_speed、phase_response、guarded_reciprocal、polynomial）与 `roe_formulas.frozen_roe_face` 经 `compile_function` 编译为 {obj}`FixedFormat <oracq.algorithms.common.arithmetic.FixedFormat>`(6,2)/(8,3) 两种格式（gamma、order、entropy_delta 等默认参数作生成期常量，roe_face 的 row/col 为 Index(2)、输出 left/right 双寄存器），在 rir-pysparq 上以叠加态一次穷举输入域——单输入函数全域 2^6/2^8 分支；多输入函数逐轴纤维穷举加联合立方体；roe_face 另做 row×col 全 16 组联合（含越界索引 3）与六实轴全幅值 16 点网格。期望值由独立的定点语义逐比特仿真（_Fx：mul/div/sqrt 幅度向零截断、add/sub 模 wrap、status 位 0 = 定义域失效、位 1 = 值域/字长越界、helper 调用合并全部参数旗标）与 float64 原式双层给出；phase_response 的初等核按模块属性 `math_approximation` 的 Chebyshev 系数逐比特复现（degree=3），另报系数配方与真函数的方法误差。status 旗标逐分支核对；reference / adapter-pysparq 在代表性程序上振幅级三方对拍。编译函数工作区实测 358–2178 量子比特，远超 OriginIR-ext 的 24 比特预算，故不走态向量路径。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `polynomial-exhaustive-6.2/8.3` | 64 / 256 分支全域 | rir-pysparq | max_error | 0 / 0（逐比特一致；旗标 44 / 197 分支全部符合值域越界预测） |
| `guarded-reciprocal-exhaustive-6.2/8.3` | 64 / 256 分支全域 | rir-pysparq | max_error | 0.235 / 0.123（0.94 / 0.98 量子，即除法截断界 1 量子内；全域无旗标，x=0 守护生效） |
| `pressure-fibers-6.2/8.3` | 3 轴 × 64 / 256 分支 | rir-pysparq | max_error / method_error | 0 / 0（常量量化偏差 1.4 / 0.578；rho=0 除零旗标与越界旗标逐分支符合） |
| `roe-speed-fibers-6.2/8.3` | 4 轴 × 64 / 256 分支 | rir-pysparq | max_error / method_error | 0 / 0（含 sqrt 截断仿真；rho≤0 定义域旗标逐分支符合） |
| `phase-response-fibers-6.2/8.3` | 2–3 轴 × 64 / 256 分支 | rir-pysparq | max_error / method_error | 0 / 0（对照系数配方；方法误差 1.770 / 1.349 为 degree=3 核固有近似误差；z=±i 除零与核区间越界旗标逐分支符合） |
| `phase-response-kernel-method` | 4001 点稠密网格 | 经典对照 | method_error（exp/sin/cos） | 0.148 / 0.195 / 0.364 |
| `roe-face-index-joint-6.2/8.3` | row×col 16 组 | rir-pysparq | max_error | 0 / 0（含越界索引 3 → 0.0；双输出 left/right） |
| `roe-face-fibers-6.2/8.3` | 6 轴 × 16 点全幅值网格 | rir-pysparq | max_error | 0 / 0（rho≤0、c2≤0 定义域旗标逐分支符合） |
| `cross-backend-*`（6 案例） | 16–64 分支 | rir-pysparq / reference / adapter-pysparq | max_pairwise_deviation | 0.0（三方振幅完全一致） |
| `basis-determinism-6.2` | 6 函数各 1 物理点 | rir-pysparq | error / failures | 0 / 0（单基态入 → 单基态出，输入保持） |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_mathfunc.py
```

产物：`out/verification/mathfunc.json`（28 个案例全过，总运行约 185 秒，VERIFY_WORKERS 控制并行度）。
