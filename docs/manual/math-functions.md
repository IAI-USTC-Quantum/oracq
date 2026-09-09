# 普通数学函数自动生成可逆量子模块

先编写纯 Python 数学函数，再调用 `compile_function` 生成可逆量子模块。函数仍可用于经典计算；量子侧由编译器处理临时寄存器、别名复制、结果 XOR 和反算。QFVM 的 Roe face 已采用这条路径。

```python
from pyqecclang import FixedFormat, compile_function, export_toffoli_u3_cz

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

该模块的接口为 rho、momentum、energy、out、status；gamma 是生成期参数。调用方式与其他 Operation 相同：

```python
from pyqecclang import Builder, Bits

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
from pyqecclang import MathConfig

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

普通无注解参数默认实数；可用 inputs 显式指定 real/complex/bool。Index(width) 为 QFVM 行列索引等无符号整数提供较短的公开寄存器，进入计算时自动转换到定点表示。常量默认参数视为生成期参数；显式放入 inputs 后也可成为量子输入。多结果可用 output_names 指定名称，例如 Roe 的 left/right。

Python 的 cmath 在分支切线上区分有符号零；当前定点编码没有该信息。因此这里提供函数族和可替换的近似实现，不承诺完整浮点/分支兼容。[对应的 cmath 行为见官方文档](https://docs.python.org/3/library/cmath.html)。

## 近似与状态

MathConfig.degree 控制实函数的 Chebyshev 多项式阶数；intervals 可替换每个数学核的区间。默认区间见 numeric.py 的 BOUNDS。配置必须覆盖各个**中间数学核输入**，不只是最外层函数的输入。状态会标记超出区间或发生定义域/字长问题的路径。当前没有自动区间推导、误差证明或最优算术电路选择。

函数系数通过 degree+1 个经典采样点生成，然后执行量子乘加递推，没有预先枚举每个输入对应的函数值。字长、近似阶数、区间和系数进入生成配置/模块属性，eps 仍不进入语言核心。

```python
def safe_inverse(x):
    return 0.0 if x == 0 else 1.0 / x
```

这个函数在 x=0 时不会因未选中分支的除零而设置最终 domain 标志。两分支电路可以被计算，但其结果和状态按量子条件选择，所有中间位最后反算。

## 中间表示与后端

完整链路为 Python 纯函数 → [MIR 0.1](../reference/math-ir.md) → RIR 0.3 → 模块化 OriginIR-ext / PySparQ。MIR 可独立 JSON 往返，再由 lower_math_ir 按其他配置降低。最终 RIR 不含 Python callback；带数学核和 helper 的调用继续保留为 Module/Call。

现有 arithmetic_native_registry 会识别生成模块内部的算术/布尔实现。PySparQ 在这些模块边界执行真实自定义 C++ 算子，跳过内部 Boolean 工作区。原生路径与门级路径使用同一组算术网络，并已通过实际执行检查。这里没有用 Python 原函数的直接求值冒充量子模拟。

```python
from pyqecclang import arithmetic_native_registry, run_pysparq
state = run_pysparq(
    application,
    native_registry=arithmetic_native_registry(application),
    memory=memory,  # 若入口声明了 QRAM
)
```

QFVM 仍查询原始守恒量；[roe_formulas.py](../api/applications/roe_formulas.rst) 是普通经典公式，roe_face 通过 compile_function 生成原 ABI 的模块。经典 Riemann 更新和 QRAM 数据结构继续与量子矩阵元计算分离。

## 命令行与案例

```bash
pyqecclang compile-function examples/math_functions.py \
  --function pressure --width 12 --fraction 6 \
  --constants '{"gamma": 1.4}' \
  --mir-output out/pressure.mir.json -o out/pressure.rir.json

pyqecclang emit out/pressure.rir.json --basis toffoli-u3-cz -o out/pressure.originir

PYTHONPATH=src .venv/bin/python tools/build_math_functions.py
```

--inputs 接收 JSON 类型映射，例如 '{"z":"complex"}' 或 '{"row":{"index":2},"x":"real"}'。不支持的源语句会给出 FunctionCompileError；这不是可以编译任意 Python 程序的工具。

out/math-functions/ 包含 pressure、roe_speed、phase_response、guarded_reciprocal、polynomial、自动 Roe face 与接入后的 QFVM，共七组产物。见 [实施面板](../development/contributing.md) 与 [验证记录](../archive/function-compiler-validation.json)。

数学精度、复杂分支切线、Roe 数值结果和资源优化仍待核验。现有测试覆盖编译、类型、可逆更新、模块复用和真实后端消费。
