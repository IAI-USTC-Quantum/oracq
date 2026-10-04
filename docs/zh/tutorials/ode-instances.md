# 组合线性 ODE 输入与求解实例

<a href="../../tutorials/ode-instances.html">English</a> · **简体中文**

`ODEProblem` 描述自治齐次系统 `u' = G u`、物理初值和终止时间。
`qode_solve` 选择算法，返回尚未准备的实例。与 [PDE 求解实例](pde-instances.md)
相同，`prepare()` 检查输入契约并缓存模块化 RIR；导出、资源估算和后端执行
都使用同一条线路。

## 准备、检查与执行

矩阵直接表示生成元 **G**，包含它的符号。例如对角元 `-1` 和 `-2` 表示衰减：

```{testcode}
import oracq

problem = oracq.ODEProblem(
    generator=oracq.MatrixInput([[-1, 2], [0, -2]]),
    initial=oracq.ArrayInput([2, -1]),
    components=("u", "v"),
    final_time=0.1,
)
instance = oracq.qode_solve(
    problem, method="taylor", config=oracq.ODEConfig(taylor_degree=2),
).prepare()
circuit = instance.circuit()
assert instance.prepared.report.ok
assert oracq.loads(circuit.dumps()) == circuit.program
assert circuit.resource_estimate().complete
assert "DEF " in circuit.originir_ext().text
assert instance.prepare().circuit() is circuit
```

准备与导出不导入 PySparQ 或 UnifiedQuantum。安装真实可选后端后，可以执行：

```python
result = instance.run_pysparq()
result = instance.run_originir_ext()
print(result.success_probability)
print(result.physical_amplitudes())
print(result.physical_values())  # approximately (1.64, -0.82) for this polynomial
```

`physical_amplitudes()` 默认返回物理成功通道的归一化振幅；
`normalize=False` 返回未归一化成功振幅。`physical_values()` 利用已知 LCU
归一化和初态范数恢复有限 Taylor 多项式的数值。此例计算
`(I + t G + t² G² / 2) u(0)`，仍有独立的 Taylor 截断误差；结果来自量子线路振幅。

输出向量按声明的分量顺序排列。标量使用一位 target，其他维度使用
`ceil(log2(size))` 位，最少一位。填充地址不参与读出，线性方法在 G 两侧
施加物理地址投影，防止通过填充分量产生额外路径。零初值产生零物理值和
零成功概率，此时没有归一化条件态。

## 独立替换算子与初态输入

| 算子输入 | 访问模型 | 范围 |
| --- | --- | --- |
| `MatrixInput(values)` | Pauli 门展开 | 有限实数或复数矩阵，最多五位 target |
| `MatrixInput(values, encoding="qram")` | 稠密实数矩阵角度表 | 实数条目，绝对值不超过 `value_scale` |
| `QRAMMatrixInput(bank, config)` | 已有矩阵角度表 | 显式维度、偏移、缩放及不可变内存 |
| `OracleGeneratorInput(encoding, size, qrams=...)` | 已有 block-encoding provider | 调用方提供算子语义及底层内存快照 |

这些选择与初态的 `UniformInput`、`ArrayInput`、`QRAMInput`、`OracleInput`
相互独立；它们也是 PDE 使用的[共用初态输入](../api/algorithms/input_model/initial.rst)。
`ArrayInput` 接受物理分量的一维数组；具名映射为每个分量提供一个样本，
如 `{"u": [2], "v": [-1]}`。门编码支持有符号及复数初态；内置 QRAM
旋转树支持非负实数振幅。自定义制备可以通过 `OracleInput` 提供物理范数和内存。

已有 block encoding 可替换显式矩阵：

```{testcode}
from dataclasses import replace

encoded = problem.generator.prepare_generator(oracq.ODELayout(2))
oracle_problem = replace(
    problem,
    generator=oracq.OracleGeneratorInput(encoded.encoding, size=2),
)
assert oracq.qode_solve(oracle_problem).prepare().prepared.report.ok
```

注册 QRAM 本身只提供字查询。下例通过显式适配器实现 G 的 block encoding：

```{testcode}
matrix_config = oracq.QRAMMatrixConfig(size=2, value_scale=2, offset=4)
words = matrix_config.encode_matrix([[-1, 1], [0.5, -2]], word_length=4)
bank = oracq.register_qram(words, address_length=4, word_length=4,
                           name="ode_generator")
qram_problem = replace(
    problem,
    generator=oracq.QRAMMatrixInput(bank, matrix_config),
    initial=oracq.ArrayInput([1, 1], encoding="qram", angle_width=3),
)
qram_instance = oracq.qode_solve(
    qram_problem, config=oracq.ODEConfig(taylor_degree=1),
).prepare()
assert qram_instance.circuit().memory
assert qram_instance.circuit().resource_estimate().qram_total > 0
```

矩阵条目地址是 `offset + row*2**width + column`，字代表
`value_scale*cos(pi*word/2**word_length)`。零字表示正 `value_scale`，
零矩阵条目需要半圈字。`encode_matrix` 填充所有补零行列，并饱和处理负端点。
每个条目的量化误差界是 `value_scale*pi/2**word_length`；每次编码调用查询两次，
归一化为 `2**width*value_scale`。稠密 QRAM 存储本身不带来稀疏访问加速。

内存快照留在 RIR 之外。准备阶段沿模块调用解析资源别名，不展开 `Repeat`。
缺失快照、同名冲突、未绑定 oracle 或位宽不匹配都会在执行前失败。

## 选择现有 QODE 算法

| 算法 | 配置 | 前提与物理值读出 |
| --- | --- | --- |
| `taylor`（默认） | `taylor_degree` | 有限多项式，已知归一化支持 `physical_values()` |
| `lchs` | `plan` 中的 `QuadraturePlan` | 需要 `dissipative=True`；求积和 Hamiltonian 模拟误差单独评估 |
| `cbmd` | `plan` 中的 `ContourPlan` | 需要 `dissipative=True`；辅助极点贡献及无穷级数尾项仍未生成 |
| `schrodingerization` | `plan` 中的 `SchrodingerPlan` | 调用方验证辅助网格和物理恢复区域 |

例如复用 LCHS：

```{testcode}
from oracq.algorithms.qode.lchs import QuadraturePlan

decay = oracq.ODEProblem(
    generator=oracq.MatrixInput([[-1, 0], [0, -2]]),
    initial=oracq.UniformInput(), final_time=0.01,
    dissipative=True, evidence="diagonal real nonpositive generator",
)
lchs_instance = oracq.qode_solve(
    decay, method="lchs",
    config=oracq.ODEConfig(plan=QuadraturePlan.cauchy(cutoff=1)),
).prepare()
assert lchs_instance.prepared.report.ok
```

`dissipative=True` 是数学声明；准备阶段检查声明，不证明 `Hermitian(G) <= 0`，
也不会自动平移生成元。`ODEConfig.hamiltonian_function` 可以替换 LCHS、CBMD
和 Schrödingerization 内部的 Hamiltonian 模拟协议。

这三种方法当前支持条件振幅读出，尚未提供完整的物理幅值恢复契约，因此
`physical_values()` 会报错。任意外部 `QODESolver` 作为 `method` 传入时也遵守
此边界。自定义 `ODEMethod.prepare(problem) -> PreparedODE` 可以提供验证过的
恢复尺度，并复用同一生命周期；自定义算法的配置直接传给算法对象。

此入口覆盖自治齐次**线性**系统。时变生成元、强迫项及非线性 ODE 需要显式适配。
现有低层 `QODEProblem`、`QODESolver` 和 `linear_qode` 仍可用于线路组合。

完整示例：

```bash
python examples/ode_instances.py
python examples/ode_instances.py --generator qram --initial qram --degree 1 --backend pysparq
python examples/ode_instances.py --backend originir_ext
python examples/ode_instances.py --method lchs
```

参见 [ODE 输入 API](../api/algorithms/input_model/ode.rst)、
[求解实例 API](../api/algorithms/qode/solver.rst) 和
[已有 QODE 组装教程](differential-equations.md)。
