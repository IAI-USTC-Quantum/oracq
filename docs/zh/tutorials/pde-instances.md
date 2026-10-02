# 组合 PDE 输入与求解实例

<a href="../../en/tutorials/pde-instances.html">English</a> · **简体中文**

`PDEProblem` 描述方程、网格、物理初值和终止时间。`qpde_solve` 选择算法，
`prepare()` 检查算法的输入要求并生成可复用的模块化 RIR 线路。
导出、资源估计和真实后端执行都使用同一份已准备线路。
创建问题和准备线路均不导入量子后端。

## 准备、检查与执行

```{testcode}
import oracq

problem = oracq.PDEProblem(
    type="heat",
    grid=oracq.UniformGrid1D(0, 2, 2),
    initial=oracq.ArrayInput([0.2, 0.1]),
    viscosity=0.1,
    final_time=0.01,
)
instance = oracq.qpde_solve(
    problem, method="QHAM",
    config=oracq.QHAMConfig(order=1, eta=-1, taylor_degree=1),
)
instance.prepare()
circuit = instance.circuit()
assert instance.prepared.report.ok
assert oracq.loads(circuit.dumps()) == circuit.program
assert circuit.resource_estimate(basic_gates="clifford+t+qram").complete
assert "DEF " in circuit.originir_ext().text
```

`prepare()` 缓存生成结果；准备前检查或执行线路会报错。
更换配置时新建实例。初态位宽错误、缺少导数表、oracle 未绑定或缺少内存快照，
均在准备阶段报错，而不是等到后端执行。

安装真实可选后端后，可以执行任一路径：

```python
results = instance.run_pysparq()
results = instance.run_originir_ext()
print(results.success_probability)
print(results.physical_amplitudes())
print(results.physical_values())
```

两者均返回 `PDEResult`，包括公开寄存器态及按场分组的物理读出。
`physical_amplitudes()` 返回成功通道的条件归一化振幅，
`physical_values()` 恢复已生成近似的物理幅值。读出排除补齐地址。
零物理解的成功概率为零，可以恢复零值，但没有归一化的条件态。
`run_originir_ext()` 通过 UnifiedQuantum 执行导出产物；
`originir_ext()` 只导出文本。
执行有状态数、指令数或比特数预算；线路可以生成和估计，并不代表它适合某个模拟器的预算。

默认 QHAM 配置用提升线性生成元的有限 Taylor 多项式演化。
`physical_values()` 恢复该多项式已知的 LCU 归一化因子，不承诺已收敛到连续 PDE。
空间离散误差、HAM 截断与收敛、角度量化误差、线性求解器误差需要分别判断。

## 独立替换初态输入

| 输入 | 提供的访问模型 | 物理范数 |
| --- | --- | --- |
| `UniformInput(value)` | 物理节点和场上的常数振幅 | 自动计算 |
| `ArrayInput(values, encoding="gates")` | 实数或复数样本的门制备 | 自动计算 |
| `ArrayInput(values, encoding="qram")` | 非负实数样本的量化旋转树制备 | 角度量化前计算 |
| `QRAMInput(bank, config)` | 已存储的旋转树角度字，支持基地址偏移 | 显式提供 |
| `OracleInput(preparation, norm, qrams=...)` | 已有的可逆振幅制备 oracle | 显式提供 |

同一种初态可以组合均匀网格或 QRAM 网格。
数组只描述物理节点，每个场对应连续样本；耦合 PDE 可使用
`ArrayInput({"u": u_samples, "v": v_samples})`，寄存器补齐自动插入。
自定义的带符号或复数 QRAM 加载器，可以连同内存快照一起包装成 `OracleInput`。

数值访问 `|i,0> -> |i,data[i]>` 本身不等于振幅态
`sum_i data[i]|i>` 的制备。因此 `register_qram` 返回主机侧内存对象，
而不是振幅 oracle；其 `database()` 方法提供 XOR 数值查询。
初态适配器显式选择制备线路及物理范数。

下面用外部保存的根节点角度制备均匀物理振幅，地址银行为四位，目标态为一位：

```{testcode}
import math

bank = oracq.register_qram({3: 16}, address_length=4, word_length=6,
                           name="initial_angles")
initial = oracq.QRAMInput(
    bank, oracq.QRAMInputConfig(state_width=1, norm=math.sqrt(0.08), offset=3),
)
prepared = initial.prepare_initial(oracq.InitialLayout(problem.grid, ("u",)))
state = oracq.simulate(prepared.preparation.operation.program(),
                      {bank.name: bank.snapshot()})
assert abs(state.amplitudes[(0, 0)] - 1 / math.sqrt(2)) < 1e-12
assert abs(state.amplitudes[(1, 0)] - 1 / math.sqrt(2)) < 1e-12
```

旋转树节点地址为 `offset + 2**depth - 1 + prefix`，
数据字以 `2*pi/2**word_length` 为单位表示 Ry 角度。
旋转树需要表示预期物理态，包括补齐地址上的零振幅。
物理范数是单独声明的输入承诺，与归一化振幅不同。

## 选择网格表示

`UniformGrid1D(0, 10, 10)` 有十个物理采样点和四位地址寄存器。
周期网格坐标覆盖 `[0, 10)`，间距为 `1`；
零 Dirichlet 网格使用内部节点，间距为 `(stop-start)/(points+1)`。

二次幂大小的周期网格使用移位线路。
其他大小和零 Dirichlet 边界当前使用最多五位空间寄存器的 Pauli 展开。
这支持十个点，但不是高效的大规模网格实现；仅展开空间导数算子，
不展开整个非线性端口或提升矩阵。应用可以实现 `PDEGrid`，直接提供更大网格的编码。

无结构网格需要先选择几何离散方法和边界处理，再存储**离散导数算子**。
只有坐标或连接关系还不能确定算子。内置布局将实矩阵存为角度银行，支持非对称算子：

```{testcode}
d2 = (("x", 2),)
layout = oracq.UnstructuredGridConfig(
    size=2, derivative_offsets={d2: 0}, value_scale=3,
)
data = layout.encode_matrices({d2: [[-2, 2], [2, -2]]}, word_length=8)
mesh_bank = oracq.register_qram(data, address_length=2, word_length=8, name="mesh")
mesh = oracq.UnstructuredGrid(mesh_bank, layout)
mesh_problem = oracq.PDEProblem(type="heat", grid=mesh,
                                initial=oracq.UniformInput(0.2), final_time=0.01)
mesh_instance = oracq.qpde_solve(mesh_problem).prepare()
assert mesh_instance.circuit().memory
assert mesh_instance.circuit().resource_estimate().qram_total > 0
```

每张导数表使用 `base + row*2**spatial_width + column` 寻址。
数据字表示 `value_scale*cos(pi*word/2**word_length)`：
零矩阵元素使用半圈角度字，零数据字则表示正的 `value_scale`。
`encode_matrices` 插入补齐，并对负端点采用饱和编码，避免整圈回绕改变符号。
每个元素的误差上界为 `value_scale*pi/2**word_length`。
相干编码使用两次查询，归一化因子为 `2**spatial_width*value_scale`；
QRAM 存储本身不意味着稀疏访问加速。
邻接表等其他访问模型可以通过自己的已验证块编码实现 `PDEGrid`。

内存快照是 RIR 外部的不可变主机数据。
准备时沿模块图的资源实参追踪入口绑定，包括共享银行的多个别名，且不展开 `Repeat`。
银行名称冲突或位宽不兼容会明确报错。

## 算法声明自己的要求

QHAM 要求各多线性端口提供可逆、可控的块编码，
初态制备具有匹配位宽、伴随和受控调用能力、零输入语义以及干净工作位。
合同报告可在 `instance.prepared.report` 查看，提升模型位于 `instance.prepared.model`。

已有 PDE 表达式可作为 `PDEProblem.type` 传入 `PolynomialPDE`；
字符串预设目前支持 `burgers` 和 `heat`。
使用 `QHAMConfig(linear_solver=...)` 可独立选择已有的三参数 QODE 函数。
该函数的 Hermitian 或耗散性要求仍然适用；实例不会隐式平移生成元。
此类函数不会自动向实例提供幅值恢复信息，所以 `physical_values()` 报错，
条件振幅读出仍然可用。

其他 PDE 算法实现 `PDEMethod.prepare(problem) -> PreparedPDE`，
自行定义网格限制、oracle 适配、输入合同报告和数学承诺。
把算法对象传给 `qpde_solve` 即可共用缓存生命周期、RIR 导出、资源估计、后端执行和物理读出。
QHAM 当前是唯一内置的具名 PDE 方法。

`clifford+t+qram` 估计包含七 T 的 Toffoli 分解；
`toffoli+clifford+t+qram` 则将 Toffoli 保留为单独基本门。
任意角旋转仍然需要合成，`t_total(epsilon)` 估计其 T 成本。
QRAM 查询次数与实际 QRAM 硬件开销分开计量。

完整示例默认只生成和估计，真实执行可选：

```bash
python examples/pde_instances.py
python examples/pde_instances.py --grid qram --initial qram --backend pysparq
python examples/pde_instances.py --backend originir_ext
```

另见 [QHAM 表达式组装](qham.md)、[QHAM 实现](../manual/qham.md)、
[PDE 输入 API](../api/algorithms/input_model/pde.rst) 和
[求解实例 API](../api/algorithms/qpde/solver.rst)。
