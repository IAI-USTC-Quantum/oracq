# 从论文访问模型到实现比较

<a href="../../en/tutorials/algorithm-research.html">English</a> · **简体中文**

本教程面向实现新量子算法的研究者：先保存含开放 oracle 的算法，再选择实现，
最后对照独立数学参考并比较成本。核心不需要量子后端。

## 保存开放程序并分析调用

```{doctest}
>>> from oracq import bind_with_report, dumps, loads, estimate_resources
>>> from oracq.applications.oracle_study import oracle_study
>>> opened, implementations = oracle_study(width=2, repetitions=3)
>>> restored = loads(dumps(opened))
>>> cost = estimate_resources(restored, require_closed=False)
>>> cost.complete, cost.qubits
(False, None)
>>> sorted((call.adjoint, count) for call, count in cost.oracle_calls.items())
[(False, 3), (True, 3)]
```

{obj}`oracle_study <oracq.applications.oracle_study.oracle_study>` 一次给出开放程序和全部候选实现；{obj}`dumps <oracq.infrastructure.serialization.dumps>` 与 {obj}`loads <oracq.infrastructure.serialization.loads>` 的往返说明开放描述可以照常保存，{obj}`estimate_resources <oracq.infrastructure.estimate.estimate_resources>` 在未闭合时也返回成本台账。访问模型是 `data ^= (address + 1) mod 2**width`。同一个角字数据库用于构造
对角块编码，每次调用都会计算角字、旋转、反算。重复三次只保存 Repeat 和
模块调用；上面的台账说明每种实现要提供三次正向和三次伴随调用。

## 选择实现并保留绑定报告

```{doctest}
>>> implementation, memory = implementations["qram_table"]
>>> linked = bind_with_report(restored, {"AngleWord": implementation})
>>> linked.report.ok
True
>>> closed = linked.require()
>>> estimate_resources(closed).qram_queries
Counter({'angle_words': 6})
```

{obj}`bind_with_report <oracq.infrastructure.linking.bind_with_report>` 返回程序与 {obj}`BindingReport <oracq.infrastructure.linking.BindingReport>` 的组合，`require()` 只在报告无问题时交出闭合程序。将选择改为 `gate_table` 或 `arithmetic` 即可比较另外两种实现。
门表枚举同一整数函数，QRAM 表在运行时提供，可逆算术直接计算并清除私有
工作字。三者公开签名一致；私有辅助位和资源消耗可以不同。改变 width 需要
重新生成程序；这不属于保持接口的实现绑定。

## 验证与资源比较

```bash
python examples/research_workflow.py
PYTHONPATH=src /path/to/backend/python examples/research_workflow.py \
  --native --widths 2 3 -o out/research-workflow-native
```

独立参考使用三角公式：每个地址的角字为 k，重复 r 次后的成功与失败幅度
分别为 `cos(r*pi*k/2**width)/sqrt(2**width)` 和对应的 sin。
这是本例旋转扩张的性质，不能推广为任意块编码直接取幂的规则。
报告包含完整输出分支的幅度误差、成功概率、程序与数据指纹、绑定报告以及
开放和闭合成本。原生验证同时调用 PySparQ 两条路径与 OriginIR 后端。

本例验证均匀输入的全部地址分支；不把单一输入态的通过推广为任意酉算子的
等价证明。可逆整数算术的接口与复净另由核心回归检查。

更完整的科学计算链复用 `tools/build_qlss_comparison.py`：原始流场 QRAM、
编译后的 Roe 公式、稀疏访问和两种 QLSS 生成器。该脚本提供开放/闭合成本、
数据与程序指纹；数值验证继续使用 `tests/verification/verify_qham_qfvm.py`
和 `verify_mathfunc.py`，分别报告实现误差与方法误差。

## 相关页面

- 手册：[资源估计](../manual/resource-estimation.md)（成本台账口径）、[契约](../manual/contracts.md)
- 规范：[开放 IR](../reference/open-ir.md)（开放声明与分批绑定）
- 算法页：[XOR 数据库视图](../manual/algorithms/xor-database.md)（门表/QRAM 两种实现的视图基础）
- API 参考：[Oracle 实现比较](../api/applications/oracle_study.rst)、[资源估计](../api/infrastructure/estimate.rst)、[绑定与能力分析](../api/infrastructure/linking.rst)
- 继续教程：[科学计算工作流](scientific-workflows.md)（QHAM/Carleman/LCHS/CBMD 的完整链路）
