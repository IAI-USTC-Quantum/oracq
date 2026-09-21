# pyqecclang 指南与案例集

0.7 的新入口是 [算法自己的约定](../manual/contracts.md)：从一个 gate 的多个角色讲起，使用 Python 结构协议定义输入与输出。以下保留 0.6 九个案例的背景与运行说明。

本指南面向第一次接触 pyqecclang 的开发者，用九个可以完整运行的案例串起整门语言：从最小的 QRAM 查询模块，一直到从偏微分方程自动推导量子线性系统。每个案例的代码与输出都在当前版本（0.6.0）上实际运行过；涉及算法数学正确性的地方，文中会如实标注哪些已经验证、哪些仍在待核验清单上。规范性的完整定义见 [pyqecclang-spec.md](../reference/language.md)（下称《规范》），本指南按其章节编号引用。

专门实现或替换微分方程求解器时，配合阅读 [QPDE/QODE 实作指南](../manual/differential-equations.md)。该指南以 LCHS、Schrödingerization 和 Carleman 说明多种 input paradigm，并提供 11 组可运行的 oracle 绑定与 protocol 替换示例。

## 1. 这门语言解决什么问题

写量子算法的人通常面临一个别扭的选择：要么直接操作量子位门，把寄存器、算术和数据访问全部手工展平；要么用高层描述，却在导出时丢失全部结构，得到一张巨大的门列表。pyqecclang 的立场是两条路都不必走。你用普通 Python 函数描述算法结构，生成器产出一种**寄存器级的模块化中间表示（RIR）**：一个作用在整个 8 位寄存器上的 Hadamard 门就是一条指令；调用别的模块就是一条调用指令，不会把对方复制进来；重复一千次就是重复计数为一千的单个节点。模块结构一直保留到导出与执行阶段，只有具体后端在需要时才按自己的能力展开。

第二个核心立场是**允许未完成的实现**。一个量子算法往往依赖若干 oracle（数据库查询、态制备、块编码），它们各有物理实现方案，成熟度各不相同。pyqecclang 允许先把算法骨架生成为带有空槽的合法程序，保存、检查、报告缺口，之后再逐个绑定实现。空槽不会被空线路悄悄冒充，导出闭合程序时未绑定的槽会被点名拒绝。

最后是诚实的边界。语言核心验证的是结构正确、序列化确定、后端可消费；算法层面的数值精度、成功概率与复杂度优势不属于当前承诺，所有产物都带着待核验标记。本指南会在每个案例末尾说明该案例验证了什么、没有验证什么。

## 2. 安装与五分钟上手

语言核心零第三方运行依赖，要求 Python 3.11 或更新版本。推荐用 uv：

```bash
uv sync --all-extras
uv run python -m unittest discover -s tests/core -v
uv run ruff check src tests examples tools
```

没有 uv 也可以直接用标准库跑：

```bash
PYTHONPATH=src python3 -m unittest discover -s tests/core -v
PYTHONPATH=src python3 examples/qram_modules.py
```

真实后端（uniqc 与 pysparq）只在集成测试和原生执行时需要，且必须真实安装，测试不以缺少后端为理由跳过：

```bash
PYTHONPATH=src /path/to/backend/python -m unittest discover -s tests/integration -v
```

## 3. 概念地图

在进入案例之前，先建立三张思维地图。

**生成阶段与运行阶段。** 你写的 Python 函数在生成阶段执行完毕：求常量、选实现、组装模块。生成结束后的 RIR 里没有任何 Python 回调、没有运行期分支，只有寄存器声明、指令、模块调用与属性。所以"换个 gamma 重新算"永远是回到 Python 重新生成，而不是改 JSON 里的一个数字。

**五层表示。** 从书写到执行依次是：Python 生成层、可选的生成性中间表示（数学函数图 MIR、PDE 与 QCL plan）、RIR 0.3、描述后端（OriginIR-ext 文本或 Toffoli/U3/CZ 严格门集文本）、执行后端（纯标准库参考执行器、uniqc、PySparQ）。只有中间三层可序列化，它们之间的每次转换都保持模块边界。

**开放与闭合。** 模块主体为空的模块是开放声明（带 `oracle_paradigm` 属性标记接口约定）；主体写好的模块是闭合实现。`unresolved` 列出缺口，`bind` 逐个补上，`export_originir` 与全部执行入口要求入口可达的部分已经闭合。绑定核对签名形状、范式与能力，不核对数学内容。

## 4. 案例一：模块化 QRAM 查询程序

第一个例子来自 `examples/qram_modules.py`，它展示三件事：带 QRAM 资源形参的模块、跨模块的资源重命名、以及查询的异或语义。

```python
from pyqecclang import QRAM, Builder, UInt, dumps, export_originir, simulate

def make_lookup(address_width=2, data_width=3):
    b = Builder(
        f"lookup_{address_width}_{data_width}",
        {"address": UInt(address_width), "data": UInt(data_width)},
        {"table": QRAM(address_width, data_width)},
    )
    b.qram("table", b["address"], b["data"])
    return b.finish()

def make_program():
    lookup = make_lookup()
    wrapper = Builder("wrapped_lookup", {"a": UInt(2), "d": UInt(3)}, {"memory": QRAM(2, 3)})
    wrapper.call(lookup, address=wrapper["a"], data=wrapper["d"], resources={"table": "memory"})
    wrapped = wrapper.finish()
    b = Builder("demo", {"address": UInt(2), "data": UInt(3)}, {"values": QRAM(2, 3)})
    b.h(b["address"])
    b.x(b["data"][0])
    b.call(wrapped, a=b["address"], d=b["data"], resources={"memory": "values"})
    return b.finish().program()

p = make_program()
print(export_originir(p).text)
print(simulate(p, {"values": [1, 2, 4, 7]}).amplitudes)
```

`lookup_2_3` 模块声明了一个名为 `table` 的 QRAM 形参（2 位地址、3 位数据），主体只有一条查询指令。`wrapped_lookup` 包装它并把资源改名为 `memory`；入口模块 `demo` 再把它改名为 `values`。QRAM 内存内容不属于程序，执行时按入口资源名提供：这里 `values = [1, 2, 4, 7]`。

导出的 OriginIR-ext 文本如下（节选），可以看到模块边界完整保留：三个 DEF 定义，入口调用在最外层，QRAM 以 `QRAMDECL` 声明。

```
QRAMDECL ram_values 2,3
QINIT 5
CREG 0
DEF m_lookup_2_3_54524a0d1220cf4bd723f2aa(v_address[2], v_data[3])
ram_values v_address[0], v_address[1], v_data[0], v_data[1], v_data[2]
ENDDEF
DEF m_wrapped_lookup_16dc443ca6e4d78bd9ef03a6(v_a[2], v_d[3])
m_lookup_2_3_54524a0d1220cf4bd723f2aa(v_a[0], v_a[1], v_d[0], v_d[1], v_d[2])
ENDDEF
DEF m_demo_b89906731f4fe79582ed7931(v_address[2], v_data[3])
H v_address[0]
H v_address[1]
X v_data[0]
m_wrapped_lookup_16dc443ca6e4d78bd9ef03a6(v_address[0], v_address[1], v_data[2]...)
ENDDEF
m_demo_b89906731f4fe79582ed7931(q[0], q[1], q[2], q[3], q[4])
```

注意 `h(address)` 在 RIR 中是一条指令，只有导出到 OriginIR 时才降低为逐比特的 `H v_address[0]` 与 `H v_address[1]`；模块调用在这一步仍然是调用。

参考执行器的输出是四个基矢分量，每个幅度为 1/2：

```python
{(0, 0): (0.4999999999999999+0j), (2, 5): (0.4999999999999999+0j),
 (1, 3): (0.4999999999999999+0j), (3, 6): (0.4999999999999999+0j)}
```

读这个结果需要记住两件事。地址寄存器经过 Hadamard 处于均匀叠加，所以四个地址各占 1/2 的概率。数据寄存器初始被置为 1（`x(data[0])`），而 QRAM 查询的语义是**异或**而非覆盖：data 的新值是 1 XOR M[address]。对内存 [1, 2, 4, 7]，四个地址分别得到 1^1=0、1^2=3、1^4=5、1^7=6。异或语义使查询自逆，这也是后面一切"载入后反算"模式的来源。

配套的固化产物在 `examples/qram.rir.json`（规范 JSON）、`examples/qram.originir` 与 `examples/memory.json`，CLI 一章会直接消费它们。

## 5. 案例二：寄存器视图、控制与重复

第二个案例演示视图操作与三个结构块。目标很朴素：在标志位为零的条件下，把一个 8 位字加 7 加一千次，最后对低两位做一次逆旋转。

```python
from pyqecclang import Builder, Bits, UInt, simulate

b = Builder("arithmetic", {"word": UInt(8), "flag": Bits(1)})
with b.control(b["flag"], 0):
    with b.repeat(1000):
        b.add_const(b["word"], 7)
with b.adjoint():
    b.rz(b["word"][:2], 0.3)
operation = b.finish()
program = operation.program()
state = simulate(program)
print(state.registers)
print(state.amplitudes)
print(len(operation.module.body))
```

运行结果：

```python
(word=uint/8, flag=bits/1)
{(88, 0): (0.9553364891256059+0.29552020666133955j)}
2
```

三个数字都值得停下来看。一千次加 7 等于加 7000，而 `add_const` 在 8 位字上按模 256 回绕，7000 mod 256 = 88——这就是第一个寄存器的终值。相位 0.9553 + 0.2955i 恰好是 e^{i·0.3}：伴随块里的 `rz(0.3)` 被逆转为 −0.3，作用在低位为零的基矢上给出两倍的 +0.15 相移。最后一个 2 是要点所在：一千次重复和全部控制逻辑在模块主体里只占**两条顶层指令**（一个 Control 包着一个 Repeat）。IR 不随重复次数膨胀，导出时才会用对数规模的辅助 DEF 表达（2^40 次重复产出不到 45 个 DEF）。

这个案例还顺带说明了受保护控制位的含义：`flag` 在整个控制块内不可被任何操作数触碰，即便模块主体看似不修改它。切片 `b["word"][:2]` 是零成本的逻辑视图，不产生任何交换网络；跨寄存器拼接用 `fuse`，改存储解释用 `reinterpret`，同样都是生成期的元数据操作（《规范》3.3）。

## 6. 案例三：开放 oracle 与绑定工作流

现在看这门语言最有特色的工作流：先生成带空槽的算法骨架，再绑定实现。例子是 Deutsch–Jozsa 算法，它的 XOR 数据库 oracle 先保持抽象。

```python
from pyqecclang import bind, export_originir, unresolved, ValidationError
from pyqecclang.algorithms.basics.oracle_algorithms import deutsch_jozsa
from pyqecclang.algorithms.input_model.oracles import abstract_database, gate_database

abstract_fn = abstract_database("Function", 2, 1)
open_program = deutsch_jozsa(abstract_fn).program()

for req in unresolved(open_program):
    print(req.name, req.paradigm, req.path)
```

输出是：

```python
Function database_xor ('deutsch_jozsa_95da411eed915249aba2', 'Function')
```

`abstract_database` 声明了一个 `database_xor` 范式的开放槽；`unresolved` 报告槽名、范式与从入口出发的最短调用路径。这份带空槽的程序是完全合法的 RIR：可以 `dumps` 保存、可以 `validate`。但导出会被拒绝，并点名缺口：

```python
export_originir(open_program)
# ValidationError: 未绑定 oracle：Function (database_xor; deutsch_jozsa_95da... -> Function)
```

绑定一个门级真值表实现（表 [0, 1, 1, 0] 是平衡函数）：

```python
concrete = gate_database(2, 1, [0, 1, 1, 0])
closed = bind(open_program, {"Function": concrete.operation})
print(len(closed.modules))
text = export_originir(closed).text
```

绑定后程序闭合为 3 个模块，OriginIR 中出现三份 DEF：真值表实现、被改写的 `Function` 槽（主体变为一条对实现的调用，属性记 `oracle_bound_to`），以及算法入口。绑定核对的是签名形状（寄存器类型逐位相等）、范式一致性与能力需求；换成 `qram_database(2, 1)` 的实现也完全合法，只是运行时需要提供内存（这是目录中 `dj_gate` 与 `dj_qram` 两个案例的差别）。

同一流程可以完全走命令行。以目录里的 Costa 案例为例：

```bash
uv run python tools/build_catalog.py                 # 生成 out/catalog/ 全部案例产物
uv run pyqecclang requirements out/catalog/costa_qram/open.rir.json
uv run pyqecclang bind out/catalog/costa_qram/open.rir.json \
    --bindings out/catalog/costa_qram/bindings.json -o out/costa-bound.rir.json
uv run pyqecclang emit out/costa-bound.rir.json -o out/costa.originir
```

`--bindings` 清单把每个槽指到一个实现程序文件与可选的资源映射；分批绑定是正常用法，绑定一部分后其余槽位保持开放。QRAM 捕获提升（《规范》8.3）处理实现私需资源的情况：链接器把资源贯通调用链提升到入口，共享同一资源的槽只占一个入口资源。

## 7. 案例四：块编码代数

块编码库维护信号位视角的算子代数，组合时归一化常数 alpha 随之传播并保存在模块属性 `be_alpha` 里。README 中的经典片段如下，断言全部成立：

```python
from pyqecclang import identity, pauli_x, product, linear_combination, scale

a = scale(2, identity(1))          # alpha = 2
b = scale(3, pauli_x(1))           # alpha = 3
c = linear_combination(-0.5, a, 2, b)
assert c.alpha == 7                # |-0.5|*2 + |2|*3
d = product(a, b)
assert d.alpha == 6                # 2 * 3
```

`linear_combination(-0.5, a, 2, b)` 的内部结构是教科书式的 LCU：用一个选择位做 Ry 制备（角度由两分支的权重比决定），按分支调用各自的块编码并施加系数相位，最后逆制备。测试直接断言矩阵语义 ⟨0|U†(αH)|0⟩ = c·H，含负系数的实例也已覆盖。alpha 经过 JSON 往返仍然存在，绑定开放块编码槽时它参与逐值核对。

需要强调的边界是：alpha 是构造时写下的归一化事实，不是可随意修改的成本标签；目标矩阵的近似误差不进入语言核心。未导出但常被算法层使用的组合子还有 `lcu`（N 项加权和）、`tensor`、`direct_sum`、`kronecker_sum`、`projector`、`reflect_zero` 等（《规范》11），目录案例 `be_algebra` 把它们组合成一个演示程序。

## 8. 案例五：算法原型——Grover 与它的同族

`algorithms/elementary.py` 提供一组可直接组装的范式算法：Deutsch–Jozsa、Grover、QFT、qubitization walk、相位估计、QSVT 序列与 oblivious amplitude amplification。以 Grover 为例，用一个具体的相位 oracle 跑一步振幅放大：

```python
from pyqecclang import simulate
from pyqecclang.algorithms.common.search import grover
from pyqecclang.algorithms.input_model.oracles import phase_marks

marked = phase_marks(2, [3])            # 给基矢 |3> 加相位
program = grover(marked, 2).operation.program()
for values, amp in sorted(simulate(program).amplitudes.items()):
    print(values, abs(amp) ** 2)
```

输出只有一行：`(3, 0) 0.9999999999999991`。两比特、单标记态在一步迭代后概率回到约 1，教科书行为。

这个例子的意义不在算法本身，而在组装方式：`grover` 内部把相位 oracle 调用、逆制备与绕零反射放进一个 `repeat(iterations)` 块，全部是显式 IR 节点；oracle 是普通参数，可以是这里的门级实现，也可以是 `declare` 出来的开放槽（目录案例 `grover_gate` 与 `grover_qram` 正是同一骨架的两种绑定，后者由 QRAM 数据库构造相位 oracle）。同族的 `phase_estimation` 在 walk 算子上做受控幂次调用加逆 QFT；`qsvt_sequence` 按相位序列交替正调用与伴随调用；这些案例的开放槽都绑定了小型门级实现并全部通过真实 OriginIR 解析。验证范围是结构与后端可达性——严格 QSVT 相位设计等数学问题在待核验清单上。

## 9. 案例六：可逆定点算术与严格门集

`arithmetic.py` 提供 14 类定点算术的自动可逆合成：`add, sub, neg, abs, mul, div, reciprocal, sqrt, lt, eq, select, and, or, xor`。生成的模块长这样：

```python
from pyqecclang import FixedFormat, fixed_arithmetic

mul = fixed_arithmetic("mul", FixedFormat(8, 3))
print([(r.name, r.type.width) for r in mul.module.registers])
print([(r.name, r.type.width) for r in mul.module.locals])
```

输出：

```python
[('a', 8), ('b', 8), ('out', 8), ('status', 2)]
[('ssa_0', 64), ('ssa_1', 64), ('ssa_2', 64), ('ssa_3', 64), ('ssa_4', 64), ('ssa_5', 64), ('ssa_6', 64), ('ssa_7', 17)]
```

接口是两个操作数、XOR 输出 `out` 与两位状态字（`status[0]` 定义域失效，`status[1]` 越出字长）。内部是一张 Boolean 网络（加法进位链、乘法移位累加、恢复除法、逐双位开方），降低到 RIR 时 SSA 节点打包进 64 位私有 bank——也就是上面那些 `ssa_i` 局部寄存器——算完把结果异或进输出，再用一个 `Adjoint` 把整段前向计算反算掉。模块属性里的 `workspace_contract="zero_in_zero_out"` 描述的正是这个零进零出契约，参考执行器与 OriginIR 后端都会在模块边界实际检查。

同一张网络有三种消费方式。降到严格门集：

```python
from pyqecclang import export_toffoli_u3_cz
text = export_toffoli_u3_cz(mul.program()).text
```

得到 3368 行文本，门统计为 U3×1970、CZ×961、TOFFOLI×430，字母表确实只剩三种门（多控制位经 Toffoli 梯子借池位完成）。第三种方式是 PySparQ 原生执行：`arithmetic_native_registry` 从模块属性里的网络 JSON 识别算术模块，`BooleanCppFactory` 把同一张网络编译成真实的 C++ 自定义算子，在模块边界跳过内部门——两条路径用同一组算术网络，没有用 Python 求值冒充量子模拟。

## 10. 案例七：把普通 Python 数学函数编译成可逆模块

`compile_function` 接受受限的纯数学函数（`examples/math_functions.py` 里的五个函数都是合法输入），产出满足"输入保持、输出异或、状态标志、工作区复净"契约的模块。编译理想气体压强公式并实际执行：

```python
from pyqecclang import FixedFormat, MathConfig, compile_function, simulate

fmt = FixedFormat(12, 6)
source = '''
def pressure(rho, momentum, energy, gamma=1.4):
    velocity = momentum / rho
    return (gamma - 1) * (energy - 0.5 * momentum * velocity)
'''
cf = compile_function(source, fmt=fmt, config=MathConfig(degree=6),
                      constants={"gamma": 1.4})
print([(r.name, r.type.width) for r in cf.operation.module.registers])

program = cf.operation.program()
state = simulate(program, initial={"rho": fmt.encode(1.0),
                                   "momentum": fmt.encode(0.5),
                                   "energy": fmt.encode(2.0)})
(values,) = list(state.amplitudes)
print(values)
print("decoded out:", fmt.decode(values[3]))
print("classic:", 0.4 * (2.0 - 0.5 * 0.5 * 0.5 / 1.0))
```

运行输出：

```python
[('rho', 12), ('momentum', 12), ('energy', 12), ('out', 12), ('status', 2)]
(64, 32, 128, 46, 0)
decoded out: 0.71875
classic: 0.75
```

寄存器值读法：输入被定点编码成 64/32/128（即 1.0/0.5/2.0 × 64）且执行后保持不变——输入保持契约的直接体现；`out` 是 46，`status` 为 0。解码得 0.71875，而同一公式的浮点值是 0.75。差异来自定点常量 0.4 编码时的截断（25.6 落到 25）与乘法后的向零舍入，这正是属性里 `rounding="toward_zero; modular_wrap"` 的行为。这个差距值得直面：电路精确实现的是某个有限布尔函数 F̃，F̃ 与理想数学函数的接近程度由字长、区间与阶数配置决定，语言不替你证明。把 `FixedFormat(16, 8)` 之类的配置调宽，误差相应缩小。

前端接受的东西够写实用公式：四则、比较、条件表达式与结构化 if（编译成 select）、静态有界 for 循环、纯 helper、`math`/`cmath` 白名单（sqrt/exp/log/三角/双曲及复数分解）、`Index(width)` 型整数输入、多返回值。拒绝的东西同样明确：递归、lambda、动态循环、异常处理、任何副作用——遇到就抛 `FunctionCompileError`。非多项式核用配置区间的 Chebyshev 采样加 Clenshaw 递推实现，采样数只依赖阶数，不枚举输入真值表。

命令行同样可以完成编译：`pyqecclang compile-function examples/math_functions.py --function pressure --width 12 --fraction 6 --mir-output mir.json -o pressure.rir.json`。MIR 数学图独立于 RIR 序列化，换一份配置可以从同一张图重新降低。这条链路最有分量的应用在 QFVM 里：Roe 通量公式就是由 `frozen_roe_face` 纯函数自动编译的（案例八）。

## 11. 案例八：QFVM 中替换 QLSS——CKS 与 Costa 两条路线

这个案例展示库里最深的一条管线：从 QRAM 里的原始流场出发，组装 frozen-Roe 线性系统，然后在两种量子线性求解器协议之间替换。完整审查记录在 [qfvm-qlss-input-model-review.md](../manual/qfvm.md)，工具 `tools/build_qlss_comparison.py` 一键生成两份对照产物。

管线的前半段与求解器无关：`roe_qfvm_inputs` 声明密度、动量、能量三个场数据库与几何表（几何表只存邻居、槽位、分量索引等结构信息，不含矩阵元素）；`roe_entry` 模块在量子侧取三单元守恒量、调用两次自动编译的 Roe face、按谱带选择，算出矩阵元素；`qfvm_sparse_access` 把九个结构位置补全成完整置换，给出 CKS 形式的位置与元素 oracle。问题的谱性质由调用者以 `SpectralPromise` 显式声明，这是硬性参数。

替换发生在问题对象与协议之间：

```python
from pyqecclang import SpectralPromise
from pyqecclang.applications.qfvm import roe_qfvm_inputs, roe_qfvm_problem
from pyqecclang.algorithms.qlss.qlss import CKSConfig, make_cks_qlss
from pyqecclang.algorithms.qlss.qlss import CostaConfig, make_costa_qlss

inputs = roe_qfvm_inputs()
problem = roe_qfvm_problem(
    inputs,
    spectrum=SpectralPromise(norm_upper=10.0, sigma_min_lower=0.5),
    amax=4.0,
)
cks_state = make_cks_qlss(CKSConfig(order=2)).solve(problem)      # 直接消费 SparseSystem
costa_state = make_costa_qlss(CostaConfig(steps=1)).solve(problem)  # 触发稀疏→BE 显式适配
```

配置取值有讲究：`amax` 是条目幅值上界的声明（参与谱界与 PTheta 角表的换算），Costa 的 `steps` 与调度参数必须与问题尺寸相称——默认的两步调度在这个四单元系统上会被调度点校验拒绝，这正是"配置是问题侧声明"的体现。

同一个问题对象交出两条路线（以下数字来自刚生成的 `out/qlss-audit/` 产物）：

| | CKS Chebyshev | Costa general walk |
| --- | --- | --- |
| 输入模型 | sparse | block_encoding（经 CKS T†ST 适配） |
| alpha / 编码逆谱界 | 36 / 72 | 36 / 72 |
| 物理输出宽度 / 信号宽度 | 4 / 20 | 4 / 25 |
| 闭合模块数 | 50 | 49 |

两个求解器都返回带范数探针的 `SolveResult`，物理范数恢复公式统一为 `‖x‖ = ‖r‖ / (alpha·sqrt(p_joint/p_solver))`。值得注意的是"替换"的准确含义：输入模型与输出契约已经打通，适配过程显式可查（Costa 的产物属性里记录了完整的 adapter trace）；但两种求解器数值等价的认证、以及完整 CFD 经典-量子闭环，都明确不在已完成之列。替换协议后应重新生成上层寄存器布局，因为两条路线的信号位签名不同（20 对 25）。

## 12. 案例九：一般 QHAM——从 PDE 到 QODE 输入

最后一个案例来自 `examples/general_qham.py`：把一个带强迫的 Burgers 方程自动推导成量子线性演化系统的输入。依据是量子适配线性化的 QHAM 方法——"二次线性化"指第二次线性化，不是只能处理二次多项式。

```python
from functools import partial

from pyqecclang.algorithms.qode.ode import linear_qode
from pyqecclang.algorithms.common.hamiltonian import taylor_hamiltonian
from pyqecclang.applications.qham import Discretization, Field, Grid, Known, PolynomialPDE, QHAMPlan, qham_input_model, structured_fd_bindings

u = Field("u")
pde = PolynomialPDE.from_equations(
    {"u": 0.1 * u.d("x", 2) - u * u.d("x") + Known("forcing")},
    label="forced_burgers",
)
plan = QHAMPlan(pde, order=2)
grid = Grid(("x",), (4,), (1.0,), boundary="periodic")
space = Discretization(pde, grid, {"forcing": [0.05, 0, -0.05, 0]})
initial = space.encode_fields({"u": [0.1, 0.2, 0, -0.1]})
bindings = structured_fd_bindings(space, initial)
qode_input = qham_input_model(plan, bindings, eta=-0.4)

solver = linear_qode("schrodingerization",
                     hamiltonian_function=partial(taylor_hamiltonian, degree=1))
solution = qode_input.solve(solver, 0.01)

shifted = qode_input.dissipative_shift()
cbmd_solution = shifted.solve(
    linear_qode("cbmd", hamiltonian_function=partial(taylor_hamiltonian, degree=1)), 0.01
)
```

这段代码里有四个层次值得逐个指认。`PolynomialPDE.from_equations` 把 Python 表达式变成可序列化的 PDE 0.1 表示；乘积、数除、`d("x", 2)` 外导数都是表达式代数，非线性项 −u·u_x 不会被偷换成 2u·D_x(u)——有限差分一般不满足精确的 Leibniz 恒等式，外导数作为算子作用在乘积之后。`QHAMPlan(pde, order=2)` 是惰性对象：它只记住 PDE 与阶数，按需生成 HAM 递推与有序张量字；对二次无强迫问题，m = 20 的计划能描述 2^21 个函数块而不写一行 JSON。`structured_fd_bindings` 生成结构化差分端口：导数是可逆移位的 LCU，同点乘积是坐标 XOR 加矩形收缩，全程不物化任何 N^r 规模的矩阵。`qham_input_model` 把这些组装成 QODE 输入模型，`solve` 调用任意 QODE 生成器并取出物理输出通道（截断 HAM 各阶之和）。

脚本末尾的 `dissipative_shift()` 是一个诚实的设计点：Schrödingerization 不要求生成元耗散，而 LCHS/CBMD 要求，所以移位必须显式施加（G 减去 μI，μ ≥ ‖G‖），并且要记录对数幅值恢复因子——它不是免费的性能优化。

命令行入口做同样的事并落盘全套推导产物：

```bash
PYTHONPATH=src python3 -m pyqecclang.applications.qham --example burgers --order 2 --row physical
```

manifest 摘录（真实输出）：PDE 最高次 2、张量秩上限 3、函数块 9 个、算子端口 L/F/B_2、状态维数 4、提升后原始维数 **129**、生成元目标宽度 8。输出目录包含 pde.json、qcl-plan.json、rows.json、qode-input.json 与 derivation.md（同伦递推公式和耦合表）。内置案例除 Burgers 外还有 KdV（m=3，16 块）、强迫三次反应、双分量耦合与二维向量 Burgers（提升维数 35,968）。验证状态：链式法则见证残差约 1e-16（检验代数推导正确），真实 PySparQ 已执行过多分量非线性端口与完整小型 QHAM→QODE→物理输出链；HAM 收敛性与求解精度仍在待核验清单。

## 13. 案例目录全景

三个工具各生成一组可再生的案例产物（均在 gitignored 的 `out/` 下）：

`tools/build_catalog.py` 产出范式目录 `out/catalog/`，33 个案例，每个含开放、部分绑定、闭合三种 RIR、绑定清单、内存输入、OriginIR 文本与报告。代表性案例：`dj_gate/dj_qram`（本文案例三）、`grover_gate/grover_qram`（案例五）、`costa_gate/costa_qram/costa_sparse_qram`（Costa 三种输入绑定）、`qfvm_gate/qfvm_qram`（九槽 QFVM 玩具物理核，案例八的入口）、`qham_qode/qham_qpde`（m=1 特例，已被一般 QHAM 取代为阶段历史）、`be_algebra`、`arithmetic`、`banked_qram`（96 位逻辑字拆双 bank）、`measure_reset`（宿主读出计划）。加 `--native-parse` 用真实 uniqc 解析器验证全部产物。

`tools/build_stage2.py` 产出 `out/stage2/` 的 28 个案例：14 类算术、Roe face、完整 Roe 矩阵元、QRAM QFVM 块编码、QFVM+Costa+filtering 组装、Carleman/Schrödingerisation/LCHS/CBMD 四类 QODE 与 QPDE，以及两条 QHAM→QPDE 路径；`tools/verify_stage2_backend.py` 补真实后端验证。

`tools/build_general_qham.py` 产出 `out/qham-general/` 的五个一般 QHAM 案例（含 Burgers 的 Schrödingerization/CBMD 替代求解）；`tools/build_math_functions.py` 与 `tools/verify_math_functions.py` 产出并验证案例七的五函数与 Roe face 编译。

## 14. 验证命令与结果边界

日常验证三条命令（与 CI 一致）：

```bash
uv run python -m unittest discover -s tests/core -v      # 语言核心与全部生成层
uv run python -m unittest discover -s tests/schema -v    # JSON Schema 结构校验
uv run ruff check src tests examples tools
```

真实后端集成（解释器需装 uniqc 与 pysparq，不跳过）：

```bash
PYTHONPATH=src /path/to/backend/python -m unittest discover -s tests/integration -v
```

集成测试把同一个程序分别交给参考执行器、PySparQ 与真实 OriginIR-ext 模拟器并比较复幅度，这一层能抓住位序、相对相位、受控调用、伴随顺序与 QRAM 非零目标语义的错误。

最后重申每个案例都适用的边界。本指南验证的是表达与组装能力：程序能生成、结构合法、序列化确定、后端能消费、小实例与经典参照一致。算法数学正确性——块编码近似误差、态制备切线、定点溢出边界、求解器精度、HAM 收敛、范数恢复——都带着 `correctness: pending` 标记，完整清单见《规范》第 19 章与各阶段 validation JSON。这些待办不阻塞任何组装或导出，它们是下一阶段的作业单。
