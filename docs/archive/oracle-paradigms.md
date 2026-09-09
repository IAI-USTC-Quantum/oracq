# Oracle 范式与普通实现

0.7 将这些范式作为可同时满足的算法库协议，而非排他的语言类型。先看 [算法自己的约定](../manual/contracts.md)。

在微分方程算法中如何使用这些范式，见 [QPDE/QODE 实作指南](../manual/differential-equations.md)：包含直接 BE、Hermitian parts、XOR/QRAM、CKS sparse、结构化算子和多线性系数输入，并说明每种适配的前提。

QFVM 的当前输入模型与 QLSS 替换契约已更新，见 [独立审查与修正记录](../manual/qfvm.md)。本文中旧的 QFVM 组装描述保留为阶段历史，不代表当前问题级入口。

本轮目标是确定接口、开放实现、可替换实现和组装方式。下表列出已提供的范式与构造。编译器检查结构，数学契约由库文档和下一阶段的验证支持。

| 范式 | 调用约定 | 已提供实现 |
|---|---|---|
| Block encoding | target、signal，以及 be_alpha | 门级 Pauli/LCU、查询角表构造的对角 BE、稀疏访问候选 BE、乘积/张量/直和/嵌入。 |
| XOR database | address 保留，data 按位 XOR 更新 | 显式 gate truth table、原生 QRAM、宽字 banked QRAM。 |
| State-prep isometry | 从零输入子空间制备 target，work 显式保留 | 基态、均匀态、复幅度门级制备、QRAM 角表制备。 |
| Sparse position | column 保留，index 原地变为位置索引，work 显式保留 | 完整置换的 gate 合成、正反向 QRAM 双表。 |
| Sparse entry | row/column 保留，data XOR 写入矩阵条目字 | gate truth table 与 QRAM 查询适配。 |
| Phase oracle | 对谓词匹配的基态施加相位 | 显式标记、由 XOR database 计算和反算。 |
| Reversible function | 输入与输出都是明确的寄存器参数 | 小型乘法/物理条目查表、查询字到旋转的适配。 |
| Algorithm stage | 自定义但固定的量子签名 | 未实现的物理核、嵌入或其他领域阶段可以保留为声明。 |

## CKS 稀疏访问

CKS 的位置接口是在给定列 j 后，将第 l 个非零位置原地写到索引寄存器。它与保留 l 的 XOR 查询不同。范式采用完整 n 位索引寄存器，稀疏度另用编译期参数记录。

gate 实现要求每列给出完整置换扩张。QRAM 实现依次计算正向位置到 work，交换 index/work，再用逆表清理 work。数据查询使用 row/column 两个坐标和 XOR 输出字。该区分依据 [CKS 论文 §1.1](https://arxiv.org/html/1511.02306#S1.SS1)。

用户可以先声明 SparseAccess 的两个 oracle，待查询方式确定后分别绑定。不能把一个 XOR 位置表直接当成原地置换实现。

## 等距角色与可逆扩张

state prep 的数学对象是从零输入子空间到目标态的映射。普通门实现给出了定义在完整寄存器空间上的 unitary 扩张，因此能够提供 inverse 和 controlled。

仅声明制备而没有声明 inverse 能力的接口，不能放进 Costa/Grover 的反射结构。零输入和 work 复净在本轮是明确契约，不是编译器已经证明的事实。

## 块编码

alpha 参与组合，保存在 IR 属性 be_alpha 中，不是可随意修改的成本标签。目标矩阵误差、相位求解误差和滤波精度不进入语言核心。

查询角表构造的对角 BE 明确编码 alpha*cos(encoded_angle/2)。稀疏 T†SWAP T 和 QFVM T_L†SWAP T_R 目前是结构候选，需要下一阶段验证矩阵解释、归一化和有效子空间。

门级矩阵 Pauli 展开是小实例实现，默认最多 5 个目标位。它用于使案例到达后端描述层级，不声称数据加载或经典展开具有量子加速。

## 使用示例

```python
from pyqecclang import Binding, bind, dumps, export_originir, unresolved
from pyqecclang.algorithms.oracles import abstract_database, qram_database
from pyqecclang.algorithms.oracle_algorithms import deutsch_jozsa

f = abstract_database("BooleanFunction", 2, 1)
open_program = deutsch_jozsa(f).program()
assert unresolved(open_program)

closed = bind(open_program, {
    "BooleanFunction": Binding(qram_database(2, 1).operation, {"table": "truth"})
})
assert not unresolved(closed)
print(export_originir(closed).text)
```

绑定前后都可以 dumps/loads。实现代码和数据表可以在后续工作中替换；如果接口宽度或 alpha 改变，应重新运行生成器。
