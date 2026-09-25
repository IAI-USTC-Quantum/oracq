# 操作、寄存器与生成过程

<a href="../../manual/concepts.html">English</a> · **简体中文**

一个 oracq 程序经历三个不同的阶段。

1. Python 执行生成函数，选择算法、确定参数并构造模块。
2. RIR 保存模块定义、寄存器签名、调用关系和未完成的 oracle。
3. 后端导出或执行这份描述。

{obj}`Builder <oracq.infrastructure.builder.Builder>` 用来构造模块，{obj}`Operation <oracq.infrastructure.builder.Operation>` 表示一个模块及其依赖，{obj}`Program <oracq.infrastructure.ir.Program>` 则包含入口可用的完整模块集合。调用 `operation.program()` 会检查结构，不会运行量子算法。

## 寄存器与视图

寄存器按名字和位宽定义。下标 0 表示最低位，切片使用 Python 的左闭右开规则。{obj}`fuse <oracq.infrastructure.ir.fuse>` 可以把不同寄存器的片段拼成一个逻辑视图；`reinterpret` 只改变位模式的数值解释（完整语义见 [RIR 规范](../reference/rir.md#23-寄存器引用与视图)）。

一个寄存器或合并视图最多 64 位，这是与 PySparQ 存储模型一致的限制。程序可以包含多个寄存器，总量子位数不受 64 限制。

切片、拼接和重解释不生成门。调用模块时，实参必须与形参同宽、同存储类型，而且不得互相重叠。

## 模块与未完成的实现

RIR 中的 {obj}`Call <oracq.infrastructure.ir.Call>` 引用共享模块定义；{obj}`Repeat <oracq.infrastructure.ir.Repeat>` 保存重复次数和主体。它们在保存为 YAML/JSON 文本时不会被展开。

开放 oracle 的 `body=None` 表示尚未提供实现。空的指令列表则表示已经实现的恒等操作。这个区别贯穿验证、绑定和后端导出。

{obj}`bind <oracq.infrastructure.linking.bind>` 根据槽位签名绑定实现，并沿调用链传递新增的 QRAM 资源。寄存器布局或 block encoding 的 alpha 改变时，应重新运行生成函数。绑定与替换的完整流程见[教程：替换 oracle](../tutorials/oracle-binding.md)。

## 精度与数学前提

语言检查结构与可用调用能力。数值格式、积分节点、Taylor 阶数和线性化截断由算法生成器决定。矩阵是否 Hermitian、某个态是否满足物理模型，以及近似误差是否足够小，都需要算法层的声明和验证。

这一职责划分使新算法可以逐步扩展自己的 Python 协议，无需改变 RIR。详细约定见[算法输入与输出](contracts.md)。

入门教程见[第一个寄存器程序](../tutorials/first-program.md)；核心对象的完整清单见 [RIR 对象](../api/infrastructure/ir.rst)、[模块构造器](../api/infrastructure/builder.rst) 与 [绑定与能力分析](../api/infrastructure/linking.rst)。
