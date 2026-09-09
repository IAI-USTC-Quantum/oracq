# 提供自己的 Hamiltonian 分解

假设你已经知道 `H=0.3I+0.7X` 的项分解。可以通过普通 Python 方法将它交给 Trotter 实现，不必先构造一个稠密矩阵，也不必注册新的语言类型。

```{testcode}
import cmath
import math
from pyqecclang.algorithms.hamiltonian import PauliOperator, TrotterTerm, hamiltonian_simulation
from pyqecclang import simulate

class MyHamiltonian:
    hermitian = True

    def trotter_list(self):
        return (TrotterTerm(0.3, PauliOperator("I")),
                TrotterTerm(0.7, PauliOperator("X")))

evolution = hamiltonian_simulation(MyHamiltonian(), time=0.4, steps=3)
state = simulate(evolution.operation.program())
expected = cmath.exp(-0.12j) * math.cos(0.28)
assert abs(state.amplitudes[(0, 0)] - expected) < 1e-12
```

这里两个项对易，因此可以直接使用解析结果检查幅度。`evolution.alpha=1`，因为生成的是完整酉演化。

对于非对易项，Trotter 的步数控制乘积公式近似。只有 Hermitian 声明还不够：每一项需要提供可调用的演化实现。若选择 QSP 路径，则应提供它需要的 BE 访问和实际 QSP 内核。

这一例子体现了算法约定的作用：Hamiltonian 对象提供 `trotter_list()`，Trotter 算法负责检查并调用它；RIR 只接收最后生成的模块与门。
