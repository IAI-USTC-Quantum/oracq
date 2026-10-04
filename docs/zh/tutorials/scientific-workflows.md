# 科学计算工作流：输入模型、求解器与复现

<a href="../../tutorials/scientific-workflows.html">English</a> · **简体中文**

本章串联 QHAM、Carleman、LCHS 和 CBMD 的完整示例，解释每一步产生的对象、
仍然开放的接口，以及结果能够支持的结论。论文保留机制与代表性实验；完整
用法、输入变体和逐步讲解在这里继续维护。

## 先区分三个层次

| 层次 | 决定什么 | 典型对象 |
|---|---|---|
| 数学输入与近似 | PDE、离散化、截断阶、初态范数、物理输出窗口 | {obj}`PolynomialPDE <oracq.applications.qham.pde.PolynomialPDE>`、{obj}`QHAMPlan <oracq.applications.qham.linearization.QHAMPlan>`、{obj}`PolynomialODE <oracq.algorithms.qnlss.carleman.PolynomialODE>` |
| 量子算法生成 | 选哪个求解器及内部模拟核，生成哪些模块调用 | {obj}`QODESolver <oracq.algorithms.qode.ode.QODESolver>`、{obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>`、{obj}`carleman_qode <oracq.algorithms.qnlss.carleman.carleman_qode>` |
| 实现绑定与执行 | 用门表、QRAM 或算术实现开放槽，提供运行期数据 | RIR {obj}`Program <oracq.infrastructure.ir.Program>`、{obj}`Binding <oracq.infrastructure.linking.Binding>`、内存快照 |

替换求解器通常需要重新生成程序；给既有槽换一个兼容实现可以使用 {obj}`bind <oracq.infrastructure.linking.bind>`。
改变公开位宽或已用于组装的 alpha 不能只替换标签。数学条件相同与接口兼容
也是两件事：一个生成元能提供 BE，并不意味着它满足所有 QODE 方法的条件。

## QHAM：从方程到提升后的线性输入

先阅读[最小 QHAM 教程](qham.md)。下面直接包含实际示例源码，因此导入、
方程、数据与求解器配置始终和仓库脚本一致。

```{literalinclude} ../../../examples/general_qham.py
:language: python
:start-at: from functools import partial
:caption: 强迫 Burgers 方程的 QHAM 组装与两种求解器选择。
```

按对象的生命周期理解这段代码：

1. {obj}`Field <oracq.applications.qham.pde.Field>` 与 {obj}`Known <oracq.applications.qham.pde.Known>` 构造未知场及按名字引用的强迫项。表达式描述待求场，
   不是把一个普通数值函数直接作用到量子振幅上。
2. {obj}`QHAMPlan(pde, order=2) <oracq.applications.qham.linearization.QHAMPlan>` 确定有限 HAM 截断与量子适配线性化规则。
   计划可保存、恢复和按行查询，不会立即物化完整提升矩阵。
3. {obj}`Grid <oracq.applications.qham.reference.Grid>` 与 {obj}`Discretization <oracq.applications.qham.reference.Discretization>` 决定周期网格、导数和已知数据。边界条件属于
   这一层，不由求解器猜测。
4. {obj}`structured_fd_bindings <oracq.applications.qham.stencils.structured_fd_bindings>` 用移位、局部系数和多线性收缩构造端口 BE；
   {obj}`qham_input_model <oracq.algorithms.input_model.qham.qham_input_model>` 再组装提升生成元和初态。提升初态各张量块的相对范数
   必须保留，不能逐块独立归一化。
5. {obj}`linear_qode <oracq.algorithms.qode.ode.linear_qode>` 选择线性演化方法；{obj}`partial(taylor_hamiltonian, degree=1) <oracq.algorithms.common.hamiltonian.taylor_hamiltonian>`
   显式指定本例的近似模拟核。`solve` 生成带信号和恢复信息的态 oracle，
   并不直接返回一个已认证收敛的经典 PDE 解。
6. `dissipative_shift` 显式改变生成元以满足另一种方法的前提。相应增长
   恢复因子不能丢弃，也不能把这个步骤当成免费提高成功率的方法。

原始方法见 [Xue 等的 QHAM 论文](https://doi.org/10.1007/s11433-024-2584-2)。
有限闭包、张量块权重与物理输出的推导见[QHAM 数学说明](../reference/qham-derivation.md)，
一般多项式、强迫与离散化限制见[完整实现文档](../manual/qham.md)。

## Carleman：开放系数端口与可替换线性求解器

Carleman 路径从多项式 ODE 的系数端口出发，按所选 cutoff 构造张量提升。
它与 QHAM 的截断对象不同，不能把两个截断阶当作同一精度参数。

```{literalinclude} ../../../examples/ode_input_models.py
:language: python
:start-at: "    # 6. Ordinary PDE"
:end-before: "    # 7."
:dedent: 4
:caption: Carleman 示例片段；导入、网格和线性求解器配置见完整脚本。
```

这里 `concrete.coefficients` 按非线性次数保存系数 BE。循环为每个系数建立
一个同宽、同信号布局、同 alpha 的开放声明，同时在 `bindings` 中保留对应
实现。`BurgersInitial` 同样是开放态制备槽，原始初态范数由宿主单独携带。

{obj}`PolynomialODE <oracq.algorithms.qnlss.carleman.PolynomialODE>` 因此描述数学输入；{obj}`carleman_qode <oracq.algorithms.qnlss.carleman.carleman_qode>` 负责提升；传入的
`linear_solver` 负责提升后的线性问题。脚本分别使用 Schrödingerization
与显式移位后的 LCHS，并把恢复因子写入记录。运行两次生成器不意味着两种
方法的成功概率、截断误差或线路成本相同。

数学约定与接口见 [Carleman 页面](../manual/algorithms/carleman.md)。
下载包含所有上下文的 {download}`完整输入模型脚本 <../../../examples/ode_input_models.py>`。

## LCHS：同一个开放图的门表与 QRAM 绑定

LCHS 的配置包含求积计划与 Hermitian 分支的模拟核。示例中的 Cauchy 节点
和一阶 Taylor 核是具体的近似选择，不是语言默认保证的目标精度。

```{literalinclude} ../../../examples/ode_input_models.py
:language: python
:start-at: "    # 3. Same open graph"
:end-before: "    # 4."
:dedent: 4
:caption: 不改调用方，只改变角数据库和初态制备的绑定。
```

`DiagonalAngles` 存储的是旋转角字，{obj}`diagonal_block_encoding <oracq.algorithms.input_model.oracles.diagonal_block_encoding>` 把它解释为
矩阵系数。`Initial` 是独立的初态制备槽。两个 `write` 调用接收同一个开放
`state`，但分别选择门表与 QRAM 实现；后者还必须提供 `memory` 中的两个表。

脚本的 `save_case` 将开放、第一次部分绑定和最终闭合的 RIR 分别保存，并
生成模块化后端文本。中间资源捕获仍保留为显式参数，不会把内存内容塞进 RIR。
更小的独立示例见[开放程序实现比较](algorithm-research.md)，方法配置与适用
条件见 [LCHS 页面](../manual/algorithms/lchs.md)。

## CBMD：复用 QHAM 输入，但显式处理耗散前提

```{literalinclude} ../../../examples/input_models.py
:language: python
:start-at: "    # ---- CBMD: QHAM"
:end-before: "    (args.output"
:dedent: 4
:caption: 已构造的 QHAM 端口通过显式适配交给 CBMD。
```

CBMD 的方法来源是 [Wang 等的轮廓矩阵分解](https://doi.org/10.1088/2058-9565/ae7b7e)。
这段代码复用 QHAM 的计划和开放系数输入，而不是在 CBMD 内重新解释 PDE。
`model.dissipative_shift()` 产生的移位与恢复信息属于数学适配；最终的
`qram_dict` 和 `qram_memory` 才是实现及运行期数据绑定。

节点、留数、有限截断及未包含项见 [CBMD 页面](../manual/algorithms/cbmd.md)。
不要根据 `check().ok` 推断轮廓截断误差或物理恢复误差已经得到认证。

## 输入模型变体与复现

{download}`输入模型变体脚本 <../../../examples/input_models.py>` 对四个家族提供
结构化端口、谱构造和 QRAM 数据路径。它是构造示例，不是跨方法准确率排行榜。

```bash
PYTHONPATH=src python examples/general_qham.py
PYTHONPATH=src python examples/ode_input_models.py
PYTHONPATH=src python examples/input_models.py
PYTHONPATH=src python examples/research_workflow.py
PYTHONPATH=src python tools/build_qlss_comparison.py
```

真实后端检查必须使用同时安装 `pysparq` 和 `uniqc` 的解释器：

```bash
PYTHONPATH=src /path/to/backend/python examples/research_workflow.py --native
PYTHONPATH=src /path/to/backend/python tools/run_verification.py \
  --group oracles --group blockencoding --group arithmetic \
  --group mathfunc --group qham_qfvm --group ode
```

数值报告应分别解释：线路相对独立有限精度参考的实现误差、有限精度或截断
相对原问题的方法误差，以及后选择成功率和物理范数恢复。程序可导出、后端
相互一致、方法对原 PDE 收敛是不同结论。完整边界见[适用范围](../manual/limits.md)。

## 格式细节、资源数据与扩展功能

| 要查的细节 | 维护位置 |
|---|---|
| RIR 节点、完整序列化文本与 OriginIR-ext 文法 | [RIR 规范](../reference/rir.md) |
| 开放声明、分批绑定、资源捕获 | [开放 IR](../reference/open-ir.md) |
| 数学函数的中间表示与序列化 | [数学 IR](../reference/math-ir.md) |
| Toffoli/旋转/QRAM 计数、通用规模扫描 | [资源估计](../manual/resource-estimation.md) |
| QFVM 的输入、Roe 算术和恢复链 | [QFVM 文档](../manual/qfvm.md) |
| 更广的算法清单与验证状态 | [算法手册](../manual/algorithms/index.md)、[验证覆盖](../development/validation-coverage.md) |
| 类、函数、参数名与返回对象 | [API 参考](../api/index.rst) |

这些页面承接完整实现与教程细节。非核心算法仍在库中维护，但它们的数量不被
用作论文中框架优越性的证据。跨框架代码量实验的脚本仍保留在
`tools/expressiveness/`，不属于论文保留的实验依据。
