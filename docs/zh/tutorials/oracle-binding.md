# 给算法替换 oracle

[English](../../tutorials/oracle-binding.html) · **简体中文**

我们先编写一个不知道函数实现的 Bernstein–Vazirani 程序，再给它绑定门实现。这样可以看清“算法已经完整”和“输入 oracle 尚未完成”之间的区别。

## 先声明输入

```{testcode}
from oracq.algorithms.input_model.oracles import abstract_database
from oracq.algorithms.basics.oracle_algorithms import bernstein_vazirani
from oracq import unresolved

given = abstract_database("BooleanFunction", 3, 1)
opened = bernstein_vazirani(given).program()
print([item.name for item in unresolved(opened)])
assert [item.name for item in unresolved(opened)] == ["BooleanFunction"]
```

```{testoutput}
['BooleanFunction']
```

{obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>` 声明一个尚无实现的抽象数据库，{obj}`bernstein_vazirani <oracq.algorithms.basics.oracle_algorithms.bernstein_vazirani>` 在它上面组装算法。打印出的列表就是 {obj}`unresolved <oracq.infrastructure.linking.unresolved>` 找到的尚未绑定开放槽名字：算法本体已经完整，缺的只是名为 `BooleanFunction` 的输入实现。该声明提供三位地址和一位 XOR 结果。算法假设函数具有 `f(x)=s·x XOR c` 的形式；声明本身不证明这个前提。

## 绑定门实现

```{testcode}
from oracq.algorithms.basics.oracle_algorithms import affine_boolean_oracle
from oracq import bind, simulate

implementation = affine_boolean_oracle(3, secret=5, bias=1)
closed = bind(opened, {"BooleanFunction": implementation.operation})
assert not unresolved(closed)

state = simulate(closed)
probability = sum(abs(a)**2 for key, a in state.amplitudes.items() if key[0] == 5)
print(state.amplitudes)
print(probability)
assert abs(probability - 1) < 1e-12
```

```{testoutput}
{(5, 0): (-0.7071067811865471+0j), (5, 1): (0.7071067811865471+0j)}
0.9999999999999989
```

门实现由 {obj}`affine_boolean_oracle <oracq.algorithms.basics.oracle_algorithms.affine_boolean_oracle>` 用普通可逆门给出；打印出的振幅只在 input 读数为 `5` 的分支上非零；value 位上的 `±1/√2` 相位差来自仿射偏置，第二行的总概率约等于 `1`。读取 input 得到 `5`，即低位在前解释的秘密位串。仿射偏置改变了相位，但不影响该结果。

## 换成 QRAM

同一开放槽可以绑定 {obj}`qram_database(3,1) <oracq.algorithms.input_model.oracles.qram_database>`。使用 {obj}`Binding(..., {"table": "truth"}) <oracq.infrastructure.linking.Binding>` 将其资源映射到入口，运行时再提供 `truth` 表。

门实现和 QRAM 实现必须兑现同一 XOR 语义。{obj}`bind <oracq.infrastructure.linking.bind>` 检查接口与能力，函数的数学形式仍由应用负责。

## 相关页面

- 手册：[输入与算子](../manual/operators.md)（oracle 范式与视图）
- 规范：[开放 IR](../reference/open-ir.md)（开放声明、分批绑定与资源捕获）
- 算法页：[Bernstein–Vazirani](../manual/algorithms/bernstein-vazirani.md)、[XOR 数据库视图](../manual/algorithms/xor-database.md)
- API 参考：[Oracle 声明与实现](../api/algorithms/input_model/oracles.rst)、[绑定与能力分析](../api/infrastructure/linking.rst)
- 继续教程：[从论文访问模型到实现比较](algorithm-research.md)
