# 算法自己的约定：从一个 gate 开始

**约定属于算法库，不属于 RIR 的语言类型系统。** 一个对象可以满足多个 Python 协议；上层算法检查自己需要的方法、参数和调用能力。新算法可以在自己的文件中定义新协议，无需修改 pyqecclang 的语法、RIR、序列化器或全局类型目录。

本章说明算法如何声明、检查和适配输入。完整可运行示例：[examples/algorithm_contracts.py](../../examples/algorithm_contracts.py)。较长的数学组装说明仍见 [QPDE/QODE 指南](differential-equations.md)。

## 1. 协议与算法生成器

Python 的 `typing.Protocol` 表示“一个对象满足什么接口”；库中的 QLSS/QODE protocol 表示“可替换的算法生成器”。二者相关但不同：

```text
输入对象满足某个 Python 接口
    → 算法生成器检查输入、选择实现
    → 生成新的 oracle/操作
    → 生成结果又满足某些 Python 接口
```

这些步骤发生在 Python 生成阶段，不是量子线路执行时的类型分派。RIR 仍保存寄存器、门、Call、Repeat 和开放声明，不保存 Python 回调。

## 2. 一个 gate 同时是 unitary、state preparation 和 LCU 输入

下面是完整代码片段：

```python
from pyqecclang import Bits, Builder, identity, requires
from pyqecclang.algorithms.input_model.interfaces import UnitaryProtocol, StatePreparationProtocol, BlockEncodingProtocol
from pyqecclang.algorithms.input_model.block_encoding import lcu

b = Builder("XGate", {"q": Bits(1)})
b.x(b["q"])
gate = b.finish()

requires(gate, UnitaryProtocol)
requires(gate, StatePreparationProtocol)
requires(gate, BlockEncodingProtocol)

initial = gate.state_preparation()   # U|0> = |1>
encoded_u = gate.block_encoding()  # U 自身：alpha=1、零信号位
sum_encoding = lcu([(1, gate), (1, identity(1))])  # 编码 U+I，alpha=2
```

这里没有给 `gate` 挂三个字符串标签。`Operation` 提供了三个实际方法，因此结构协议成立：`unitary()`、`state_preparation()`、`block_encoding()`。

**默认把 Operation 的全部公开寄存器视为完整 unitary 的目标空间。** 多个寄存器按签名顺序从低位到高位拼接。如果一个操作公开了 2 位 target 和 3 位 work，把整个 Operation 当作 unitary 态制备时，得到的是 5 位态，不会擅自认为那 3 位可以消失。

如果你知道 work 在零输入制备后复净，可以显式缩小态的目标解释：

```python
from pyqecclang.algorithms.input_model.oracles import StatePreparation

# operation 的接口必须恰好是 q / tmp。
# prep = StatePreparation.from_unitary(operation, target="q", work="tmp", clean_work=True)
```

`clean_work=True` 是你对这个实现的算法承诺，语言不证明它。默认的全公开寄存器适配无需这个承诺；当前 BE/state-prep 便捷包装器的单个 target 仍受 64 位限制。

这也澄清“任意 unitary 能否作为 b”：它一定定义了某个 `b=U|0>`。是否是你要解的问题中的那个 b，由应用定义；是否支持 QLSS 所需的 inverse/controlled，由具体 QLSS 检查。

## 3. QLSS 直接接收这个 gate 作为 b

```python
from pyqecclang import BlockSystem, LinearSystem, SpectralPromise, identity
from pyqecclang.algorithms.qlss.qlss import CostaConfig, make_costa_qlss

problem = LinearSystem(
    block=BlockSystem(identity(1), gate, SpectralPromise(1.0, 1.0)),
    rhs_norm=1.0,
)
solver = make_costa_qlss(CostaConfig(steps=1))
report = solver.check(problem)
report.require()
result = solver(problem)
```

`BlockSystem` 通过 `state_preparation()` 把 gate 取得为初态接口。Costa 的需求可直接查看：

```python
print(solver.contract.to_dict())
print(report.to_dict())
print([p.__name__ for p in solver.provides])  # ['StateOracleProtocol']
```

当前 Costa 组合入口要求 A 能提供 BE，或能提供可显式转换的 CKS 稀疏输入；b 能提供零输入、干净工作区的态制备；输入必须同宽，并支持组合中实际需要的 adjoint/controlled。谱界依然是问题的数学声明。

`result.state_oracle()` 得到输出态 oracle。输出 state oracle 包含成功信号，不会自动变成“没有后选择的干净初态制备”；若下一步只需要它的完整物理 unitary，可显式取 `result.operation`，此时其全部信号寄存器也属于完整空间。

`check` 不运行 QLSS 内核或量子模拟。它可以调用输入对象的方法取得/生成具体访问视图，因此自定义适配方法应当确定、无外部副作用；构造昂贵时可由对象自己缓存结果。

## 4. A 的参数从哪里读

已有 BE 包装类保留原 API，并增加只读别名：

| 读取 | 含义 |
|---|---|
| `A.width` / `A.main_qubit` | 目标寄存器位数 |
| `A.signal_qubits` / `A.anc_qubit` | BE 公开信号位数，不包含模块私有 locals 的资源峰值 |
| `A.alpha` | 当前编码的归一化常数 |
| `A.capabilities.adjoint` / `.controlled` | 从整个依赖图推导的有效调用能力 |
| `A.spec` | 可转 JSON 的描述快照，含寄存器、资源和 open/closed 状态 |
| `A.type` | 描述性角色名，例如 `block_encoding`；不作为唯一类型或分派依据 |

如果 A 是自己的算子对象，先通过 `A.block_encoding()` 取得该算法所需的具体视图，再读取这些参数。算法根据接口方法工作，不需要向一个中央枚举追加类型名。

CKS 输入是位置操作和元素操作的集合：`SparseAccess.sparse_access()` 返回自身；其 `spec.components` 分别描述 position/entry，参数中记录 sparsity 和 value_width。整个集合的 `anc_qubit=None`，避免把两个不同查询操作伪装成同一个 unitary 的总辅助位数。

## 5. 不继承 pyqecclang，也可以提供新输入

```python
from pyqecclang import identity, requires
from pyqecclang.algorithms.input_model.interfaces import BlockEncodingProtocol

class MyMatrix:
    def block_encoding(self):
        return identity(2)

A = MyMatrix()
requires(A, BlockEncodingProtocol)
```

算法库使用 `@runtime_checkable typing.Protocol` 做结构检查；`MyMatrix` 没有继承基类，也没有注册。

`requires` 检查接口存在，不证明方法签名或返回值的数学含义。具体算法调用 `block_encoding()` 后仍检查返回的 BE 类型、宽度、alpha 和所需能力。例如返回一个字符串会得到 `INPUT_ADAPTER` 报告，不会因为“恰好有同名方法”就被接受。

新约定可以完全放在应用里：

```python
from typing import Protocol, runtime_checkable
from pyqecclang import requires

@runtime_checkable
class HasDiagonal(Protocol):
    def diagonal_values(self) -> tuple[float, ...]: ...

def my_diagonal_algorithm(operator):
    requires(operator, HasDiagonal, path="my_diagonal_algorithm.operator")
    values = operator.diagonal_values()
    # 本算法继续检查长度、值域，并生成自己的操作。
    return values
```

这里的 `HasDiagonal` 不是新的语言类型，导出器从来不需要知道它。

## 6. Hermitian、Trotterizable 与 QSP 分别表示什么

数学算符与实现它的物理 unitary 必须分清：非 Hermitian 的 A 也可以有 unitary 的 block encoding `U_A`，但不能据此把 A 当作 Hermitian Hamiltonian。

当前 [hamiltonian.py](../api/algorithms/common/hamiltonian.rst) 定义了这一类算法自己的接口：

| 协议 | 提供的内容 | 由谁判断 |
|---|---|---|
| `HermitianProtocol` | `.hermitian` 声明 | HamSim 要求其为 True，数学真实性由实现方负责 |
| `BlockEncodingProtocol` | `.block_encoding()` | BE/QSP 型实现取得访问模型与 alpha |
| `TrotterizableProtocol` | `.trotter_list()` | Trotter 实现取得有序的带系数分解 |
| `EvolvableProtocol` | `.evolution(t)` | 单项提供 `exp(-it H_j)` 的具体酉演化 |

一个对象可以同时满足其中几个协议。下面的对象只有 Trotter 访问，没有 BE：

```python
from pyqecclang.algorithms.common.hamiltonian import PauliOperator, TrotterTerm, hamiltonian_simulation

class MyHamiltonian:
    hermitian = True

    def trotter_list(self):
        return (
            TrotterTerm(0.3, PauliOperator("I")),
            TrotterTerm(0.7, PauliOperator("X")),
        )

evolution = hamiltonian_simulation(MyHamiltonian(), 0.4, steps=3)
```

此处 `H=Σ c_j H_j`；生成器按列表顺序调用各项 `evolution(c_j*t/steps)`，再用 RIR `Repeat(steps)` 保留重复结构。返回 `BlockEncoding`，alpha=1。一般非对易情况下这是乘积公式近似，精度由配置与算法分析负责。

`TrotterTerm` 不只是一个矩阵名字。当前实现要求单项演化返回具体 `Operation`，其非零公开寄存器只有 `target`，且所有项同宽。需要辅助寄存器的精确 HamSim 可以后续扩展这个算法的契约；不能把有后选择信号的一般 Taylor BE 冒充精确单项酉演化。

`hamiltonian_simulation(method="auto")` 当前优先选择可用的 Trotter 分解，否则尝试 BE/QSP 路径。`method="qsp"` 要求注入实际的 `qsp(BE,time)->BlockEncoding` 生成器；库内目前没有通用 QSP-HamSim 内核，缺失时明确报错。这是局部、可替换的选择策略，不是“所有 Hermitian 输入都由语言自动实现 QSP”。具体 QSP 的额外前提仍由被注入实现负责检查。

`EncodedOperator(encoding, hermitian=False)` 可以表示非 Hermitian 算子并提供 BE；HamSim 会拒绝它，QODE 则可以消费其 BE，并按自身方法检查条件。

## 7. LCHS、Schrödingerization 各自负责约定

```python
from pyqecclang import QODEProblem, identity, scale
from pyqecclang.algorithms.qode.ode import linear_qode

problem = QODEProblem(
    generator=scale(-1, identity(1)),
    initial=gate,
    dissipative=True,
    initial_norm=1.0,
)
lchs = linear_qode("lchs")
schrodinger = linear_qode("schrodingerization")
lchs.check(problem, time=0.1).require()
state = lchs.solve(problem, 0.1)
other = schrodinger.solve(problem, 0.1)
```

`QODEProblem` 与 `QODEProtocol` 是普通算法库对象。LCHS/CBMD 的问题级 `solve` 要求显式 `dissipative=True`；Schrödingerization 不要求这个声明，但仍需应用保证辅助网格与恢复区有效。`check().ok` 的含义是结构与声明满足需求，绝不是已经证明 PDE 的性质。

旧 `(G, initial, time)` 调用仍兼容，照旧由调用者承担数学前提；新应用推荐 `.solve(QODEProblem(...),time)`，这样未声明耗散和明确非耗散不会被无声地接受。没有给 RIR 增加一个“耗散矩阵”类型。

Carleman 自己要求多线性系数端口 `F_p` 的布局与初态范数，再把生成的线性系统交给所选 QODE protocol。提升系统不一定耗散，需要时显式移位；这些数学边界仍见 [QPDE/QODE 指南](differential-equations.md)。

## 8. requires、报告、绑定各做一件事

`requires(value, Interface)` 是最小工具，用来写算法自己的检查。需要聚合错误和展示契约时，可用 `InputRequirement` 与 `ProtocolContract`，把本算法接受的 Python 协议类及适配函数传进去。它们不依赖闭合的字符串类型集合。

报告中的错误含 `code/path/expected/actual/message`。常见错误有 `INPUT_PROTOCOL`（缺接口）、`INPUT_ADAPTER`（取得视图失败）、`INPUT_WIDTH`、`INPUT_CAPABILITY`、`INPUT_PROMISE`、`OUTPUT_LAYOUT`。`report.require()` 会一次抛出全部已发现的问题，异常仍是 `ValidationError` 的子类，方便旧代码兼容。

`bind` 继续只负责已有 RIR 槽位的 ABI、alpha、调用能力和资源连接。它不负责认识 `HermitianProtocol`、`TrotterizableProtocol` 或你的新数学性质；声明的数学真实性仍由算法及实现方负责。同签名/alpha 可晚绑定，否则重新运行宿主生成器。

`A.spec` / `contract.to_dict()` / `report.to_dict()` 是检查报告，可存 JSON，但不是另一套可执行 IR。运行实现仍通过 `dumps/loads` 保存 RIR，恢复入口操作可用 `Operation.from_program(program)`。协议对象、缓存和任意宿主方法不会被偷偷序列化。

## 9. 运行与工程边界

```bash
PYTHONPATH=src .venv/bin/python examples/algorithm_contracts.py
PYTHONPATH=src .venv/bin/python -m pyqecclang inspect out/algorithm-contracts/qlss.closed.rir.json
```

示例输出包括接受/拒绝报告、开放/闭合 QLSS、同一 unitary 的 LCU、两个 QODE 描述和 Trotter 的模块化 OriginIR。

算法约定在 Python 层扩展，RIR 格式仍为 0.3。生产使用边界与完整验收方式见 [工程成熟度](limits.md)：语言与组装检查接受工程测试；Costa/CKS、LCHS、Schrödingerization、Carleman、QHAM 的完整数值正确性仍是实验状态。尚不能把整套量子求解算法称为已认证的生产求解器。
