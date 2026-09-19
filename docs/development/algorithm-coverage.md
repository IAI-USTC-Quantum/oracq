# 量子算法覆盖工作板

本文档跟踪 pyqecclang 相对开源生态的量子算法实现覆盖面，按"input model 自由度优先"的原则安排实现顺序。核心原则：算法与输入模型解耦，一个算法定义通过开放声明与分批绑定服务多种 input model（abstract / gate / qram 三层）。

## Input model 词汇表

| 缩写 | Input model | 仓库实现位置 |
|---|---|---|
| SP | 态制备 / 幅值编码 | `algorithms/oracles.py` StatePreparation（gate/qram） |
| QRAM | QRAM/QROM 数据查找 | `infrastructure/ir.py` Load、`oracles.py` XorDatabase/banked |
| SO | 稀疏矩阵 oracle（位置+元素分离） | `oracles.py` SparseAccess |
| BE | Block encoding / LCU | `algorithms/operators.py`、`block_encoding.py` |
| UO | 酉算子 oracle（受控 U） | `interfaces.py` UnitaryProtocol |
| HAM | 哈密顿量表示（Pauli/稀疏/低秩） | `algorithms/hamiltonian.py` |
| FO | 函数 oracle（相位/幅度/布尔） | `infrastructure/mathfunc/`、`oracles.py` |
| CP | 经典数据直接参数化（QUBO、变分参数） | `algorithms/variational.py` 等 |
| ODE | 微分方程输入（A 的访问模型 + 初态 + b） | `algorithms/ode_models.py`、`ode.py`、`qham.py` |
| DM | 密度矩阵 / Gibbs 态 | 未实现 |
| EST | 采样/测量访问（shadow、层析） | 未实现 |

## 工作板

状态：✅ 已实现；🚧 进行中；📄 候选论文（无可靠开源实现或值得以本仓库抽象重写）。

### A 类：算法与 input model 显式解耦的框架（最高优先级）

| 条目 | 算法 | Input model | 状态 | 说明 |
|---|---|---|---|---|
| A1 | QSVT 标准变换库：相位合成（奇偶分解）、定点搜索、矩阵求逆/符号函数、特征态过滤、哈密顿模拟相位 | BE | ✅ | Gilyén et al. 2019；实现于 `algorithms/qsvt.py`，在 `transforms.py` 序列骨架上组装 |
| A2 | PREPARE–SELECT 标准分解：prepare/select 独立 oracle、系数经 gate/QRAM/alias 绑定、与 qubitization walk 打通 | BE/HAM | ✅ | Low–Chuang 2019、Babbush et al. 2018；实现于 `algorithms/prepare_select.py`，含 alias sampling |
| A3 | 量子行走框架：图 oracle input model、Szegedy/MNRS marked-vertex 搜索框架 | FO（邻接 oracle） | ✅ | MNRS 2011；实现于 `algorithms/graph_walks.py`，含邻接 oracle 三层与 hitting time 工具 |
| — | Element distinctness / 碰撞查找 | FO | 📄 | Ambainis；依赖 A3 框架 + Johnson 图行走（需非均匀 setup 与附加数据寄存器） |
| — | 非 Abelian HSP | FO | 📄 | 数学门槛高，暂缓 |

### B 类：新 input model 打开新算法族

| 条目 | 算法 | Input model | 状态 | 说明 |
|---|---|---|---|---|
| B1 | Gibbs/热态制备、量子 Metropolis | DM + HAM | ✅ | 实现于 `algorithms/density.py`：纯化访问三层范式 + QSVT 纯化路线 Gibbs 制备（DM input model 已建立） |
| B2 | 量子 SDP / 凸优化（Brandão–Svore、van Apeldoorn–Gilyén） | DM/BE | ✅ | 实现于 `algorithms/qsdp.py`：迹估计探针电路 + MMW 驱动 + 逐轮量子子程序生成 |
| B3 | Heinrich 量子求和/积分 | FO + QRAM | ✅ | 实现于 `algorithms/integration.py`：比较器线性读出 + QAE，三层绑定对拍 |
| B4 | Jordan 梯度估计 | FO（概率 oracle） | ✅ | 实现于 `algorithms/gradient.py`：相位 oracle 三层 + 单查询 d 分量读出 |
| B5 | Alias sampling 态制备、Select-Swap QROM | CP→SP / QRAM | ✅ | alias 见 `prepare_select.py`；Select-Swap 实现于 `algorithms/data_loading.py`（含 `qrom_cost` 资源对拍） |

### C 类：抽象层次能更正确地重写的算法

| 条目 | 算法 | Input model | 状态 | 说明 |
|---|---|---|---|---|
| C1 | 量子化学低秩分解（DF/THC/SF → LCU → BE） | HAM（低秩张量）→ BE | ✅ | 实现于 `algorithms/lowrank.py`：DF（两能级合成旋转 + 受控旋转对角编码）与 THC，产物可进 qubitization_walk |
| C2 | 非线性/随机 ODE 新进展（Fokker–Planck、SDE） | ODE | ✅ | 实现于 `algorithms/sde.py`：FP 守恒形式离散化接 QODEProblem 契约，OU 矩解析对拍 |
| C3 | DQI（解码量子干涉优化） | CP（编码约束） | ✅ | 实现于 `algorithms/dqi.py`：Dicke 态 + 译码器开放声明，分布与 Krawtchouk 闭式逐点对拍 |
| C4 | QRAM 系 QML（QPCA、推荐系统） | QRAM + SP | ✅ | QPCA 实现于 `algorithms/qpca.py`：LMR 密度矩阵指数化 + QPE，一阶收敛率见证；推荐系统暂缓 |

### 已实现基线（对拍目录见 `applications/catalog.py`）

| 领域 | 算法 |
|---|---|
| Oracle/搜索/估计 | Deutsch–Jozsa、Bernstein–Vazirani、Simon、Grover、振幅放大、QAE/量子计数、QPE、Hadamard/SWAP test |
| Fourier/数论 | QFT、QFT 加法、模乘、order finding |
| 模拟 | Trotter、Taylor 化、qubitization walk、QSVT 序列骨架、OAA |
| 线性系统 | CKS（Chebyshev）、Costa（离散绝热 walk）、CKS §5 VTAA 变时求解器 |
| 微分方程 | LCHS、CBMD、Schrödingerization、Carleman、Euler 历史态、结构化 FD、QHAM、QFVM |
| 其他 | 变分（QAOA-MaxCut/VQE）、重复码纠错、定点算术、mathfunc 前端 |

## 维护规则

- 新增算法归入相应类别文件，不写旧导入兼容层；寄存器名、位宽、视图保留到后端降低。
- 修改本工作板时同步更新状态列；算法实现完成后在对应行注明实现文件。
- 选择新论文的标准：论文明确定义输入 oracle；可被协议契约体系机械转写；正确性可由 contracts 检查背书。
- 逐算法的实现与验证详情见[算法页面](../manual/algorithms/index.md)，修改工作板时同步对应页面。
