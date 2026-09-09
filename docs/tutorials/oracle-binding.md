# 给算法替换 oracle

我们先编写一个不知道函数实现的 Bernstein–Vazirani 程序，再给它绑定门实现。这样可以看清“算法已经完整”和“输入 oracle 尚未完成”之间的区别。

## 先声明输入

```{testcode}
from pyqecclang.algorithms.oracles import abstract_database
from pyqecclang.algorithms.oracle_algorithms import bernstein_vazirani
from pyqecclang import unresolved

given = abstract_database("BooleanFunction", 3, 1)
opened = bernstein_vazirani(given).program()
assert [item.name for item in unresolved(opened)] == ["BooleanFunction"]
```

该声明提供三位地址和一位 XOR 结果。算法假设函数具有 `f(x)=s·x XOR c` 的形式；声明本身不证明这个前提。

## 绑定门实现

```{testcode}
from pyqecclang.algorithms.oracle_algorithms import affine_boolean_oracle
from pyqecclang import bind, simulate

implementation = affine_boolean_oracle(3, secret=5, bias=1)
closed = bind(opened, {"BooleanFunction": implementation.operation})
assert not unresolved(closed)

state = simulate(closed)
probability = sum(abs(a)**2 for key, a in state.amplitudes.items() if key[0] == 5)
assert abs(probability - 1) < 1e-12
```

读取 input 得到 `5`，即低位在前解释的秘密位串。仿射偏置改变了相位，但不影响该结果。

## 换成 QRAM

同一开放槽可以绑定 `qram_database(3,1)`。使用 `Binding(..., {"table": "truth"})` 将其资源映射到入口，运行时再提供 `truth` 表。

门实现和 QRAM 实现必须兑现同一 XOR 语义。`bind` 检查接口与能力，函数的数学形式仍由应用负责。
