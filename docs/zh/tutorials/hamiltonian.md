# 提供自己的 Hamiltonian 分解

<a href="../../en/tutorials/hamiltonian.html">English</a> · **简体中文**

假设你已经知道 `H=0.3I+0.7X` 的项分解。可以通过普通 Python 方法将它交给 [Trotter 实现](../manual/algorithms/trotter.md)，不必先构造一个稠密矩阵，也不必注册新的语言类型。分解直接写成 {obj}`TrotterTerm <oracq.algorithms.common.hamiltonian.TrotterTerm>` 的元组，每个因子是一个 {obj}`PauliOperator <oracq.algorithms.common.hamiltonian.PauliOperator>`；{obj}`hamiltonian_simulation <oracq.algorithms.common.hamiltonian.hamiltonian_simulation>` 接收任何能给出这些项的对象。

```{testcode}
import cmath
import math
from oracq.algorithms.common.hamiltonian import PauliOperator, TrotterTerm, hamiltonian_simulation
from oracq import simulate

class MyHamiltonian:
    hermitian = True

    def trotter_list(self):
        return (TrotterTerm(0.3, PauliOperator("I")),
                TrotterTerm(0.7, PauliOperator("X")))

evolution = hamiltonian_simulation(MyHamiltonian(), time=0.4, steps=3)
state = simulate(evolution.operation.program())
expected = cmath.exp(-0.12j) * math.cos(0.28)
print(state.amplitudes)
print(evolution.alpha)
assert abs(state.amplitudes[(0, 0)] - expected) < 1e-12
```

```{testoutput}
{(0, 0): (0.9541441386892553-0.11505006784720223j), (1, 0): (-0.03308314468637069-0.274368274461448j)}
1.0
```

打印出的 `(0, 0)` 幅度就是解析值 `exp(-0.12j)·cos(0.28)` 的浮点形式，`X` 项把剩余概率转到 `(1, 0)` 上。这里两个项对易，因此可以直接使用解析结果检查幅度。`evolution.alpha=1`，因为生成的是完整酉演化。

对于非对易项，Trotter 的步数控制乘积公式近似。只有 Hermitian 声明还不够：每一项需要提供可调用的演化实现。若选择 QSP 路径，则应提供它需要的 BE 访问和实际 QSP 内核。

这一例子体现了算法约定的作用：Hamiltonian 对象提供 {obj}`trotter_list() <oracq.algorithms.common.hamiltonian.TrotterizableProtocol.trotter_list>`，Trotter 算法负责检查并调用它；RIR 只接收最后生成的模块与门。

## 相关页面

- 手册：[输入与算子](../manual/operators.md)（Hamiltonian 的算子视图与协议）
- 算法页：[Hamiltonian 演化](../manual/algorithms/hamiltonian-simulation.md)、[Trotter 乘积公式模拟](../manual/algorithms/trotter.md)、[截断 Taylor 块编码](../manual/algorithms/taylor-block-encoding.md)
- API 参考：[Hamiltonian 演化](../api/algorithms/common/hamiltonian.rst)
