# 应用组装与本轮边界

一般阶数与规则化 PDE 的 QHAM 新实现见 [数学推导](qham-general-derivation.md) 和 [实现说明](qham-general-implementation.md)；下文 m=1 部分为早期特例记录。

QFVM 的当前输入模型与 QLSS 替换契约已更新，见 [独立审查与修正记录](qfvm-qlss-input-model-review.md)。本文中旧的 QFVM 组装描述保留为阶段历史，不代表当前问题级入口。

所有例子由 pyqecclang.workloads.build_case 生成。tools/build_catalog.py 保存开放、部分绑定、闭合 IR、绑定清单、内存输入和 OriginIR-ext。输出位于被 Git 忽略的 out/catalog，源码和计划文件保留在仓库中。

## Costa QLSS 与 filtering

general walk 按插值 H(s) 块编码、指定辅助位反射和全局相位组装。H(s) 的 R(s) 矩阵用 RY 与 Z 分解，并保留 A、A†、B、B† 的调用。主体结构参考本地 QRAM-Simulator 的 general 分支和 [Costa 等论文](https://arxiv.org/html/2111.08152v1)。

filtering 是独立模块。它制备 unary 权重态，按 unary 控制位调用 walk 的幂，随后反制备。Dolph–Chebyshev 计划通过 Laurent 多项式递推给出权重、stride 和 offset，负幂通过 Adjoint/Repeat 表达。所有滤波与 walk 信号位均保留。

没有迁移 PySparQ Filtering 中尚未完成且直接返回成功概率 1 的占位逻辑。本实现确实生成上述量子结构，但尚未认证其解态、调度误差或成功概率。

## QFVM

qfvm_inputs 给出位置、反向槽、列向位置、流场、边界、物理条目、幅度转导和残差制备九个开放槽。访问制备模块明确查询数据、调用物理核、转导幅度，并反算临时数据。BE 组装保留左右制备和交换，最后注入 Costa/filter 生成器。

qfvm_gate/qfvm_qram 分别提供门级或 QRAM 绑定。物理核绑定成明确的玩具查表函数，用于判定输入、输出、工作区和资源传播是否可用，不冒充 Roe 通量或完整流体计算。

可以只绑定数据访问，继续把 QfvmPhysicalEntry 保留为 null 主体。外层算法、控制、过滤、寄存器布局和已经绑定的资源都不会丢失。这样的 IR 是本语言的正常工作状态。

## QHAM、QODE 与 QPDE

QHAM 实例采用 m=1、两个空间点。提升生成元由三块线性项、两个 Kronecker 项和两个折叠耦合项构造。全局维数 10 补齐到 16，物理输出为第一个通道。初值按 u、u、0、u⊗u 构造，保留通道相对幅度。

QODE 将有限时间步组成历史线性系统，使用注入的 QLSS 后选择终点通道。QPDE 提供离散化接口，再委托 QODE。qham_qode 和 qham_qpde 分别验证两条生成器链。

本轮开放槽是 QhamLinear、QhamFold、QhamInitial，具体小型实现使用 Pauli/LCU 和门级态制备。没有改用 synthetic tridiagonal 夹具冒充该提升系统。

此外提供 LCHS、Schrodingerisation 的阶段接口，以及一个标量 Carleman 截断和两点热方程、Poisson 实例。Schrodingerisation 的 Hamiltonian lift 仍是可替换 oracle；其小绑定用于范式检查，不是特定 PDE 的精度认证。

## 旧文档覆盖

coverage.md 对每个旧用例给出新例子或设计变更说明。它保证每类范式有实现路径和可生成证据，不声称已经重写并执行全部旧 .qec，也不以旧 CLIR 逐字相同为目标。

旧文档的六组参考负载分别对应：
- QFVM 的访问和求解链。
- BE 算术、多项与异构 LCU。
- 可逆算术开放槽和普通查表实现。
- 稀疏访问、Costa walk 与 filtering。
- oracle 目录与可调用 QLSS。
- QHAM 的提升和 QODE/QPDE 注入。

## 下一阶段

下一阶段按算子矩阵、状态、成功分支、读出与外层行为分层验证。重点包括 CKS 稀疏归一化、一般矩阵 BE、Costa 的初态与反射、Dolph 系数约定、相位求解、QFVM 物理核与 QHAM 高阶构造。

本轮不把这些尚未认证的结论作为施工前置条件。所有候选构造和具体小绑定均在报告中标记为 paradigm，correctness 为 not_assessed。
