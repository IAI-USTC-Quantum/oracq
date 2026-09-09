# 一般 QHAM 自动推导与输入模型：实施计划

依据 arXiv:2411.06759v2（量子适配线性化）及其 v1 secondary linearization 名称。二次表示第二次线性化，不等于只能处理二次多项式。

| 编号 | 工作项 | 状态 |
|---|---|---|
| H1 | 从同伦方程推导高阶递推、张量积规则和有限闭包证明 | done |
| H2 | 定义可序列化的多项式演化 PDE 规则与惰性 QCL 生成器 | done |
| H3 | 矩形算子、张量位置和强迫注入的模块化 BE 降低 | done |
| H4 | 保留初值范数的结构化制备、开放 QODE 输入与输出通道 | done |
| H5 | Burgers/KdV/三次反应/耦合系统的数学与后端案例 | done |

首批规则：一阶时间演化；任意有限空间维、有限分量；未知场及其空间导数的有限多项式，已知空间系数和强迫；自治离散算子，固定 eta=hH。数学公式也指出时变 eta/算子的位置，但当前量子组装不把时变问题静默冻结成自治问题。

核心推导：U0'=L U0+f；Ui'=L Ui-eta*sum(l=0..i-1)(1+eta)^(i-1-l)*C_l，C_l 是各多线性非线性项的 q^l 系数。对最高次数 D>=2，令 p=D-1，张量字 a=(a0,...,ak-1) 的权重为 p*sum(a)+k，取上界 p*m+1。非线性替换让阶数和严格下降、权重不增，强迫移除一个动态因子；因此形成有限的精确 HAM 截断系统闭包。

必须检验链式法则恒等式 dLift(U)/dt = G*Lift(U)，覆盖 m>1、eta!= -1、非零强迫与三次非线性。HAM 收敛、量子求解精度及复杂度优势单独列为待核验，不把小规模稠密参考当作量子实现。

## 施工结果

H1–H5 已完成本轮范围。通用递推和闭包见 [数学推导](qham-general-derivation.md)，完整生成入口和输入契约见 [实现说明](qham-general-implementation.md)，验证证据见 [记录](qham-general-validation.json)。

- pde.py / linearization.py：规则化 PDE、惰性 QCL plan 与逐行耦合。
- quantum.py / stencils.py：矩形 BE、移位/收缩端口、初态权重与 QODE 接入。
- reference.py：独立 HAM 求值、链式法则和逐块线性作用。
- python -m pyqecclang.qham：从 PDE JSON 自动导出推导和 QODE manifest。
- 五组案例包含 m=3 KdV 和 35,968 维二维向量系统；全部生成开放/闭合 RIR。
- 数学与结构测试、真实后端均通过；收敛与量子求解精度不在完成声明中。
