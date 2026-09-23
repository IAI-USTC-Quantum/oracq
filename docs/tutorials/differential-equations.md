# 为同一个线性问题替换 QODE 方法

本例考虑 `u'=-u`。我们用 `-I` 的 BE 表示生成元，用一个 X 门制备初态，再分别选择 LCHS 和 Schrödingerization。

```{testcode}
from pyqecclang import Bits, Builder, QODEProblem, identity, scale
from pyqecclang.algorithms.qode.ode import linear_qode
from pyqecclang.algorithms.common.hamiltonian import taylor_hamiltonian
from functools import partial

b = Builder("initial_one", {"q": Bits(1)})
b.x(b["q"])
problem = QODEProblem(scale(-1, identity(1)), b.finish(),
                      dissipative=True, initial_norm=1.0)
kernel = partial(taylor_hamiltonian, degree=1)
lchs = linear_qode("lchs", hamiltonian_function=kernel)
schrodinger = linear_qode("schrodingerization", hamiltonian_function=kernel)

assert lchs.check(problem, time=0.05).ok
assert schrodinger.check(problem, time=0.05).ok
first = lchs.solve(problem, 0.05)
second = schrodinger.solve(problem, 0.05)
assert first.width == second.width == 1
```

两个生成器都返回目标宽度为一位的态 oracle，但内部辅助寄存器、近似方式和恢复条件不同。这个例子检查的是组装接口；有限一阶 Taylor 并不代表已经得到任意精度的微分方程解。

`dissipative=True` 是问题声明。LCHS 的问题级入口要求它明确存在；Schrödingerization 使用其他数学条件，尤其需要选择合适的辅助窗口和恢复通道。

## 从 PDE 开始

线性 PDE 先经过空间离散化，形成 `DiscretePDE(generator, initial)`，再交给 `make_qpde(qode)`。非线性多项式 PDE 可以形成 `PolynomialODE`，由 Carleman 生成提升系统，再调用同一三参数线性求解协议。

完整的热方程、Burgers 方程与多种 given-oracle 示例位于：

```bash
uv run python examples/ode_input_models.py
```

各方法的符号、尺度和适用前提见[微分方程完整文档](../manual/differential-equations.md)。

从开放系数、分批绑定到不同数据路径的完整讲解见[科学计算工作流](scientific-workflows.md)。
