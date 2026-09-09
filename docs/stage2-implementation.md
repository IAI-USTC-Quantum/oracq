# 第二阶段：寄存器级模拟、算术与微分方程组装

一般阶数与规则化 PDE 的 QHAM 新实现见 [数学推导](qham-general-derivation.md) 和 [实现说明](qham-general-implementation.md)；下文 m=1 部分为早期特例记录。

QFVM 的当前输入模型与 QLSS 替换契约已更新，见 [独立审查与修正记录](qfvm-qlss-input-model-review.md)。本文中旧的 QFVM 组装描述保留为阶段历史，不代表当前问题级入口。

本阶段已经按 [计划](stage2-plan.md) 实现一条可审阅的路径。判定对象是接口是否能表达算法、模块是否能组装以及后端是否能消费描述。**所有新算法的数学正确性、精度、成功概率和优势均待核验。** 当前代码不会把这些结论写成语言保证。

## 架构决策

Oracle 仍是带寄存器签名的操作，protocol 仍是普通 Python 生成函数。固定点字长、小数位、舍入、Taylor 阶数、积分节点、CBMD 极点和 Carleman 截断阶都属于生成参数。生成后进入 IR 的是具体寄存器、模块调用、门、QRAM 声明与描述属性。BE 的 alpha 继续随组合传播；不新增语言级 eps。

RIR 0.3 增加 Module.locals：私有寄存器以零态进入，在模块返回时应复净，不进入公开参数签名。公开的 signal/work 与私有 locals 用途不同：前者可以携带后选择条件，后者不能携带返回后的信息。结构验证不证明复净，执行器可检查，具体算术采用 compute/XOR/uncompute 构造。单个寄存器仍不超过 64 位，模块可以有多个寄存器。

一个模块可以同时具有：

- 可序列化的门级主体，用于普通执行与 OriginIR 导出。
- 宿主注册表中的 PySparQ 原生实现，用于跳过内部算术门和工作区。
- 或只有签名，等待后续绑定。native-only 实现可以模拟该模块，但不会使其获得门级导出能力。

| 层次 | 当前表示 | 保留边界 |
|---|---|---|
| 宿主生成 | Python 函数、配置记录 | 输入模型和 protocol 可替换 |
| RIR | Module、Call、Repeat、Ref、locals | 调用不自动展开 |
| OriginIR-ext | DEF、QRAMDECL、带工作参数的调用 | 重复和模块保留 |
| 严格目标门集 | Toffoli / U3 / CZ，另有 QRAM 指令 | 在 DEF 内分解门，保留 DEF 调用 |
| PySparQ | 寄存器事件执行、NativeRegistry | 命中原生模块时直接执行一个算子 |

## 算术实现

[arithmetic.py](../src/pyqecclang/arithmetic.py) 提供 14 类生成器：add、sub、neg、abs、mul、div、reciprocal、sqrt、lt、eq、select、and、or、xor。既有 add_const 保留寄存器模加；乘除和开方中的移位直接由位视图布线实现。

生成器先建立具有公共子表达式复用的 Boolean DAG。加法用进位网络，乘法用移位累加，除法用恢复除法，平方根用逐双位恢复算法。每个中间布尔量被计算到私有工作位，结果 XOR 到输出，随后反算。复杂度随位宽多项式增长，没有用整个函数的真值表替代分解。目前工作空间较多，尚未做可逆 pebbling 或最优算术电路优化。

FixedFormat 使用显式二补码和小数位。乘除采用幅值运算后恢复符号，截去低位，溢出按字长回绕。status[0] 表示除零或负数开方等定义域失效，status[1] 表示越出当前实现的可表示范围；这两个位也不是误差界。符号边界、溢出定义的一致性和全字长行为待下一阶段系统核验。高层将各节点状态按 OR 汇总。

BooleanCppFactory 从同一 Boolean 图生成 C++ 布尔求值和输出 XOR，支持切片与跨寄存器视图。宿主只序列化网络数据，原生代码不进入 RIR。导入网络会检查 DAG 和端口布局。

```python
from pyqecclang import FixedFormat, fixed_arithmetic, export_toffoli_u3_cz
from pyqecclang.arithmetic import arithmetic_native_registry
from pyqecclang.backends.pysparq import run_pysparq

operation = fixed_arithmetic("div", FixedFormat(width=12, fraction=6))
program = operation.program()
description = export_toffoli_u3_cz(program).text
registry = arithmetic_native_registry(program, cache_dir="out/native-cache")
# 实际应用先用 Builder 准备输入，再调用 operation。
state = run_pysparq(program, native_registry=registry)
```

严格门集 pass 处理多重控制、伴随和相位，包括 RZ 对应的受控相位修正。高阶控制使用干净的合取工作位并反算；共享的目标门集工作区作为附加 DEF 参数传递。QRAM 保留独立资源指令，尚未分解成物理 QRAM 网络。

## PySparQ 审阅与运行边界

实际调用的是相邻 QRAM-Simulator 的 pysparq.dynamic_operator.compile_operator。它需要 C++17 编译器和与已安装核心匹配的头文件。当前 CPU 版本的 System.registers 是可增长的存储；“64 位”限制针对单个寄存器。

已确认两点接口差异：

1. 动态算子提供调用和 dag，但不会自动拥有 conditioned_by 接口。适配器使用真实 split_systems，将满足量子控制条件的基态分区交给原生算子，再 combine_systems 合并；零控制在两端翻转。
2. 当前动态共享库中的 System.get() 会访问自身的静态寄存器表，可能与 Python 核心持有的表分离。生成的 C++ 算子使用由适配器提供的寄存器 ID 访问 s.registers.at(id)，避免在共享库内重新查表。真实受控叠加态、伴随、视图和 Roe 算术冒烟已执行。

这条路径依赖已审阅的 CPU ABI，用户自写 DynamicCppFactory 也需要遵守该存储契约。核心后端不在导入 pyqecclang 时加载。不修改 QRAM-Simulator 或 UnifiedQuantum。

## QFVM 的具体路径

以 [QFVM 原文](https://arxiv.org/abs/2102.03557) 的 QRAM/量子制备构造为依据，本轮实现专门化为周期一维 Euler、三守恒量、frozen-Roe Jacobian。它没有宣称完成原文所有空间维度、网格和边界条件。

[roe.py](../src/pyqecclang/roe.py) 从两侧 rho、动量、总能量计算速度、压力、总焓、Roe 平均、声速、三个特征值、特征向量及逆变换，并形成左右通量 Jacobian。熵修正采用显式 delta 参数。行/列分量是量子寄存器，以算术选择网络选择矩阵元。

[qfvm.py](../src/pyqecclang/qfvm.py) 的矩阵元模块查询西/中/东三个单元的原始守恒量，对两个界面调用 Roe 算术，得到西、中心、东块，加入中心质量项并考虑 dx。数据库不含按行列预计算的 A 元素。

组装包含：

- 4 位寄存器上的九槽位制备，三分量补齐到四分量。
- geometry QRAM 存稀疏位置、反向槽位和源单元/分量索引，支持 Hermitian dilation 的转置一侧。
- T_R 只做结构准备；T_L 调用量子矩阵元计算，再查询只以数值字为地址的 PTheta。
- PTheta 根据绝对值给旋转角，符号另作相位；矩阵元和角度工作区均反算。
- S 交换行列索引以及正反槽位，组成 T_L† S T_R，记录 alpha=9*Amax。
- 输出可直接注入已有 Costa QLSS，并保留 Dolph–Chebyshev filtering。

Amax 是调用者提供的矩阵元幅值上界声明，目前不证明它成立。PTheta 表截断/角度精度以及整条 BE 的角块等价性待核验。

[flow_data.py](../src/pyqecclang/flow_data.py) 在经典侧维护原始流场、界面通量、残差符号与残差平方范数树。修改一个单元只重算相邻两个界面、受影响的三个残差单元，以及叶子到根的树节点。QRAMStore 返回显式地址 patch。**当前 PySparQ QRAM 对象没有可用的局部写接口，因此原生物化时只重建变更的 bank；逻辑局部更新不能被表述为模拟器已支持物理局部写入。**

## 四类 QODE/QPDE

主接口位于 [differential.py](../src/pyqecclang/algorithms/differential.py)。

| 方法 | 已生成的结构 | 开放接口 / 当前专门化 |
|---|---|---|
| Carleman | F_p 的张量放置、各阶间转移、截断提升、初态张量幂、一级通道选择 | F_p 为 d×d^p 补齐 BE；常量项可用 p=0；线性求解器可注入 |
| Schrödingerization | warped 初态、QFT、P⊗H1−I⊗H2、演化函数、逆 QFT、物理 p 通道 | G=H1+iH2；周期辅助网格与恢复通道由生成参数确定 |
| LCHS | H+kL 分支、离散核权重、演化函数、LCU、初态 | A=L+iH，当前自治齐次；L≥0 是输入假设 |
| CBMD | QST Eq.12 的实极点级数权重、辅助极点记录、H+(k/a)L、LCU | 同上；有限截断和辅助极点贡献显式列为遗漏项 |

输入可以从一个 A 的 BE 通过 A/A† 组合得到 HermitianParts，也可直接提供两个 Hermitian BE；它们还可以来自稀疏 oracle、QRAM 算术或其他生成器。没有要求输入以稠密矩阵给出。

演化子接口 hamiltonian_function(H,t) 返回 BlockEncoding，保留其 alpha。当前普通实现 taylor_hamiltonian 是有限阶 Taylor/LCU 的可导出候选，不冒称精确的单位归一化 HamSim。后续可以注入 QSP、其他 HamSim 实现或在模块边界绑定原生模拟。有限阶参数是生成器参数。

linear_qode 返回既有的 (G, initial, time) 接口，能直接接入 QHAM。PDE 适配器将空间离散化结果保持为开放 oracle 或 PolynomialODE；本轮用该方式组装四种 QPDE，以及 QHAM→QPDE→CBMD/Schrödingerization。空间离散化本身没有被同名占位的“完整 PDE 算法”替代，仍是显式可替换的输入环节。

Carleman 原始数据维度按二次幂寄存器承载，系数的矩形零填充由输入模型声明。初态张量幂通过多次调用制备操作构造，并非复制未知量子态。Schrödingerization 的选定 p 是否处于有效恢复区域、时间依赖和非齐次扩展、Carleman 尾项及输入耗散条件均待后续处理。

## CBMD 对应文献

用户指定的论文是 [Quantum Simulation of Non-Hermitian Special Functions and Dynamics via Contour-based Matrix Decomposition](https://arxiv.org/abs/2511.10267)，QST 11, 035027 (2026)，DOI [10.1088/2058-9565/ae7b7e](https://doi.org/10.1088/2058-9565/ae7b7e)。实现采用 v3 的一阶动力学 Eq.12–13，未将 CBMD 简化成任意积分节点的 LCHS 别名。ContourPlan 保存实节点、非实辅助极点、两类系数与省略项；通用 cbmd_function 另外提供 f(sH+L) 的 Hermitian function 注入点。

其他方法参考：[LCHS](https://arxiv.org/abs/2303.01029)、[Schrödingerization](https://arxiv.org/abs/2212.14703)、[Carleman 量子算法](https://arxiv.org/abs/2011.03185)。这些资料为构造依据，不意味着当前实现已经复现论文的全部假设或结论。

## 重建和审阅

```bash
PYTHONPATH=src .venv/bin/python tools/build_stage2.py
.venv/bin/pytest tests/core tests/schema -q
.venv/bin/ruff check src tests examples tools
# 以下使用安装了真实 pysparq/uniqc 且可找到 g++ 的 Python 环境。
PYTHONPATH=src python tools/verify_stage2_backend.py
PYTHONPATH=src python -m unittest discover -s tests/integration -v
```

本工作区可使用 ../QECC.Lang/.venv/bin/python 执行真实后端；已有 GNU 编译器是 ../.tools/envs/devenv/bin/x86_64-conda-linux-gnu-g++。本轮只在 out/toolchain/g++ 建了本地符号链接供动态编译器发现，没有改动 .tools。

out/stage2/ 包含 28 个案例的 open.rir.json、closed.rir.json、modular.originir、toffoli_u3_cz.originir、report.json，QFVM 另有 memory.json 和增量 patch。开放版本和闭合版本分别保存；模块及 Repeat 不在 JSON 中展开。后端实际执行/解析证据见 [阶段验收记录](stage2-validation.json)，算法正确性待核验项见该记录及上文。
