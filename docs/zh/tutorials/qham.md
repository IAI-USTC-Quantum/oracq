# 从 PDE 表达式生成 QHAM 输入

<a href="../../en/tutorials/qham.html">English</a> · **简体中文**

这一例子先定义 Burgers 方程，再建立一个二阶 HAM 计划。创建计划不会展开所有张量块，也不会构造完整矩阵。未知场是 {obj}`Field <oracq.applications.qham.pde.Field>`，方程冻结成 {obj}`PolynomialPDE <oracq.applications.qham.pde.PolynomialPDE>`，二阶计划保存在 {obj}`QHAMPlan <oracq.applications.qham.linearization.QHAMPlan>` 里。

```{testcode}
from oracq.applications.qham import Field, PolynomialPDE, QHAMPlan

u = Field("u")
pde = PolynomialPDE.from_equations({"u": 0.1*u.d("x", 2) - u*u.d("x")})
plan = QHAMPlan(pde, order=2)
restored = QHAMPlan.loads(plan.dumps())
print([(port.name, port.arity) for port in pde.ports])
print(plan.block_count, plan.max_rank, plan.order)
assert restored == plan
assert any(port.arity == 2 for port in pde.ports)
```

```{testoutput}
[('L', 1), ('B_1', 2)]
8 3 2
```

打印结果显示：PDE 分解成一个线性端口 `L` 和一个二阶多线性端口 `B_1`；计划只记录 `8` 个张量块、最大秩 `3` 和截断阶 `2`，没有展开任何块，也没有构造完整矩阵。PDE 的非线性项被表示为多线性端口。选定网格和边界后，{obj}`structured_fd_bindings <oracq.applications.qham.stencils.structured_fd_bindings>` 可以通过移位、系数乘子和同点收缩生成这些端口的 BE；也可以先声明开放端口，等待数据访问实现。

量子组装入口 {obj}`qham_input_model <oracq.algorithms.input_model.qham.qham_input_model>` 生成 QODE 所需的提升算子和初态。此后可以替换线性求解方法，并选择表示 HAM 各阶之和的物理输出通道。

运行完整示例：

```bash
uv run python examples/general_qham.py
uv run python -m oracq.applications.qham --example burgers --order 2 --eta=-0.4
```

若要检查推导公式，阅读[QHAM 数学推导](../reference/qham-derivation.md)。若要绑定具体 oracle，阅读[QHAM 实现说明](../manual/qham.md)。截断 HAM 的代数表示和原始 PDE 的收敛问题需要分别验证。

完整脚本及其逐步解释见[科学计算工作流](scientific-workflows.md)，其中也对照了
Carleman、LCHS 和 CBMD 的输入、求解器选择与绑定边界。

## 相关页面

- 手册：[一般 QHAM 自动生成](../manual/qham.md)
- 规范：[QHAM 数学推导](../reference/qham-derivation.md)
- 算法页：[QHAM](../manual/algorithms/qham.md)、[QODE 问题对象与协议](../manual/algorithms/qode-problem.md)
- API 参考：[QHAM](../api/algorithms/input_model/qham.rst)、[PDE 模型与适配](../api/applications/qham/pde.rst)、[QHAM 有限闭包](../api/applications/qham/linearization.rst)、[结构化差分端口](../api/applications/qham/stencils.rst)
