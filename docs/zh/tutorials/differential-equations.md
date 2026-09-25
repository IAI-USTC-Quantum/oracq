# 为同一个线性问题替换 QODE 方法

<a href="../../tutorials/differential-equations.html">English</a> · **简体中文**

本例考虑 `u'=-u`。我们用 `-I` 的 BE 表示生成元，用一个 X 门制备初态，再分别选择 [LCHS](../manual/algorithms/lchs.md) 和 [Schrödingerization](../manual/algorithms/schrodingerization.md)。生成元写成 {obj}`scale <oracq.algorithms.input_model.operators.scale>` 缩放 {obj}`identity <oracq.algorithms.input_model.operators.identity>`，问题连同初态一起装进 {obj}`QODEProblem <oracq.algorithms.qode.ode.QODEProblem>`。

```{testcode}
from oracq import Bits, Builder, QODEProblem, identity, scale
from oracq.algorithms.qode.ode import linear_qode
from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
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
print(first.width, second.width)
print(len(first.operation.program().modules), len(second.operation.program().modules))
assert first.width == second.width == 1
```

```{testoutput}
1 1
31 28
```

两个入口都由 {obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>` 按方法名生成，共用 {obj}`taylor_hamiltonian <oracq.algorithms.common.hamiltonian.taylor_hamiltonian>` 充当近似模拟核。打印结果显示：两个生成器给出同样的 target 宽度 `1`，但内部模块数不同（LCHS 为 31，Schrödingerization 为 28），对应不同的辅助寄存器、近似方式和恢复条件。这个例子检查的是组装接口；有限一阶 Taylor 并不代表已经得到任意精度的微分方程解。

`dissipative=True` 是问题声明。LCHS 的问题级入口要求它明确存在；Schrödingerization 使用其他数学条件，尤其需要选择合适的辅助窗口和恢复通道。

## 从 PDE 开始

线性 PDE 先经过空间离散化，形成 {obj}`DiscretePDE(generator, initial) <oracq.algorithms.qpde.pde.DiscretePDE>`，再交给 {obj}`make_qpde(qode) <oracq.algorithms.qpde.pde.make_qpde>`。非线性多项式 PDE 可以形成 {obj}`PolynomialODE <oracq.algorithms.qnlss.carleman.PolynomialODE>`，由 [Carleman](../manual/algorithms/carleman.md) 生成提升系统，再调用同一三参数线性求解协议。

完整的热方程、Burgers 方程与多种 given-oracle 示例位于：

```bash
uv run python examples/ode_input_models.py
```

## 相关页面

- 手册：[微分方程完整文档](../manual/differential-equations.md)（各方法的符号、尺度和适用前提）
- 算法页：[QODE 问题对象与协议](../manual/algorithms/qode-problem.md)、[LCHS](../manual/algorithms/lchs.md)、[Schrödingerization](../manual/algorithms/schrodingerization.md)、[Carleman 线性化](../manual/algorithms/carleman.md)
- API 参考：[QODE 组装接口](../api/algorithms/qode/ode.rst)、[PDE 模型与适配](../api/algorithms/qpde/pde.rst)、[Carleman 线性化](../api/algorithms/qnlss/carleman.rst)
- 继续教程：[科学计算工作流](scientific-workflows.md)（从开放系数、分批绑定到不同数据路径的完整讲解）
