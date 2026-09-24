# 算法目录

算法库按用途组织为十个子包：`input_model`（输入模型与数据访问）、`common`（通用原语）、`qlss`（线性系统）、`qnlss`（非线性系统）、`qode`（常微分方程）、`qpde`（偏微分方程）、`qml`（量子机器学习）、`optimization`（量子优化与变分方法）、`basics`（基础示例算法）与 `qec`（量子纠错）。下表给出入口文件、已实现的内容以及使用时需要留意的边界。API 参考列出了完整签名。

| 子包 | 文件 | 实现与范围 |
|---|---|---|
| input_model | [`contracts.py`](../../api/algorithms/input_model/contracts.rst) | 可检查输入契约：Oracle 能力/规格、协议契约与验收报告 |
| input_model | [`operators.py`](../../api/algorithms/input_model/operators.rst) | BlockEncoding 类与 identity/product/scale/LCU 等基本组合 |
| input_model | [`oracles.py`](../../api/algorithms/input_model/oracles.rst) | XorDatabase、StatePreparation、StateOracle、SparseAccess 四大 Oracle 范式与 abstract/gate/qram 工厂 |
| input_model | [`interfaces.py`](../../api/algorithms/input_model/interfaces.rst) | StatePreparation/Unitary/BlockEncoding 等共享协议与适配函数 |
| input_model | [`block_encoding.py`](../../api/algorithms/input_model/block_encoding.rst) | 块编码代数组合：tensor、direct_sum、projector、lcu、pauli_word 等 |
| input_model | [`sparse.py`](../../api/algorithms/input_model/sparse.rst) | 稀疏访问到块编码的适配，含 Chebyshev 块与可逆查表 |
| input_model | [`spectral.py`](../../api/algorithms/input_model/spectral.rst) | 谱对角块编码、稀疏谱块编码与谱态制备 |
| input_model | [`lowrank.py`](../../api/algorithms/input_model/lowrank.rst) | 量子化学低秩分解（DF/THC）哈密顿量 LCU/块编码 |
| input_model | [`data_loading.py`](../../api/algorithms/input_model/data_loading.rst) | Select-Swap QROM 数据加载与代价模型 |
| input_model | [`qdata.py`](../../api/algorithms/input_model/qdata.rst) | QVector/QMatrix 量子数据结构（平方范数树与 sample-and-query） |
| input_model | [`density.py`](../../api/algorithms/input_model/density.rst) | 密度矩阵纯化访问、Gibbs 态制备与迹距离等经典工具 |
| input_model | [`qham.py`](../../api/algorithms/input_model/qham.rst) | 从有限 HAM 闭包构造 QODE 输入及物理输出通道 |
| input_model | [`graph_walks.py`](../../api/algorithms/input_model/graph_walks.rst) | 图 oracle 输入模型与 Szegedy/MNRS 行走搜索框架 |
| common | [`prepare_select.py`](../../api/algorithms/common/prepare_select.rst) | LCU 的 PREPARE–SELECT 标准分解与 alias 采样 |
| common | [`qsvt.py`](../../api/algorithms/common/qsvt.rst) | QSVT 标准变换：QSP 相位、矩阵求逆、特征滤波、定点搜索 |
| common | [`transforms.py`](../../api/algorithms/common/transforms.rst) | qubitization、显式相位序列 QSVT、oblivious amplification 的组装 |
| common | [`state_preparation.py`](../../api/algorithms/common/state_preparation.rst) | 初态组合：扩展初态、物理子空间选择、BE 作用到态 |
| common | [`hamiltonian.py`](../../api/algorithms/common/hamiltonian.rst) | Pauli 项演化、Trotter 组合、Taylor BE、可注入的 QSP 接口 |
| common | [`fourier.py`](../../api/algorithms/common/fourier.rst) | 正/逆 QFT、零宽 work 适配、模 2^n 的 Fourier 加法 |
| common | [`arithmetic.py`](../../api/algorithms/common/arithmetic.rst) | 可逆定点算术：Boolean SSA 网络、compute/uncompute、原生注册 |
| common | [`estimation.py`](../../api/algorithms/common/estimation.rst) | QPE、标准振幅估计、Hadamard test、Swap test |
| common | [`search.py`](../../api/algorithms/common/search.rst) | Grover、Grover iterate、成功子空间振幅放大 |
| common | [`walks.py`](../../api/algorithms/common/walks.rst) | 周期格点上的 Hadamard coined walk |
| common | [`spectral_synthesis.py`](../../api/algorithms/common/spectral_synthesis.rst) | 谱线路的算子级合成优化（均匀受控制备、扇出合并） |
| common | [`integration.py`](../../api/algorithms/common/integration.rst) | Heinrich 量子求和与数值积分 |
| qlss | [`qlss.py`](../../api/algorithms/qlss/qlss.rst) | 问题契约、Costa walk/filter、CKS 基础 Chebyshev/LCU 路线 |
| qlss | [`vtaa_cks.py`](../../api/algorithms/qlss/vtaa_cks.rst) | CKS §5 变时幅度放大：QSP 判决时钟、分频带逆 LCU、Ambainis 嵌套放大与 A' 反计算 |
| qnlss | [`newton.py`](../../api/algorithms/qnlss/newton.rst) | 量子牛顿法：M_F 数据结构、差分 Jacobian oracle |
| qnlss | [`carleman.py`](../../api/algorithms/qnlss/carleman.rst) | 多项式 ODE 的 Carleman 有限阶张量提升 |
| qode | [`ode.py`](../../api/algorithms/qode/ode.rst) | QODE 协议和 Euler history 组装；具体方法在独立文件中 |
| qode | [`ode_models.py`](../../api/algorithms/qode/ode_models.rst) | Hermitian 反厄米分解与 LinearODE 共享输入模型 |
| qode | [`cbmd.py`](../../api/algorithms/qode/cbmd.rst) | 轮廓分解矩阵函数与演化组装 |
| qode | [`lchs.py`](../../api/algorithms/qode/lchs.rst) | 耗散线性演化的 Hermitian 分支有限加权和 |
| qode | [`schrodingerization.py`](../../api/algorithms/qode/schrodingerization.rst) | 辅助坐标 + Fourier 变换的非酉演化量子表示 |
| qode | [`sde.py`](../../api/algorithms/qode/sde.rst) | Fokker–Planck/SDE 到线性 ODE 的输入模型与经典见证 |
| qode | `legacy.py` | 早期演化工厂兼容层（make_lchs_qode 等） |
| qpde | [`pde.py`](../../api/algorithms/qpde/pde.rst) | QODE 协议之上的 QPDE 输入与求解薄封装 |
| qml | [`recommendation.py`](../../api/algorithms/qml/recommendation.rst) | Kerenidis–Prakash 量子推荐系统 |
| qml | [`qpca.py`](../../api/algorithms/qml/qpca.rst) | 密度矩阵指数化 + 相位估计的量子主成分分析 |
| qml | [`qsdp.py`](../../api/algorithms/qml/qsdp.rst) | Gibbs 采样 + 迹估计 + 矩阵乘权的量子半定规划框架 |
| qml | [`qcnn.py`](../../api/algorithms/qml/qcnn.rst) | 量子卷积神经网络的经典侧张量/卷积/池化逻辑 |
| qml | [`qcnn_layer.py`](../../api/algorithms/qml/qcnn_layer.rst) | 量子卷积神经网络的量子构件 |
| optimization | [`dqi.py`](../../api/algorithms/optimization/dqi.rst) | DQI 解码量子干涉：GF(2) max-XORSAT |
| optimization | [`variational.py`](../../api/algorithms/optimization/variational.rst) | 参数化 ansatz、MaxCut QAOA、Pauli 测量和 VQE 测量电路集合 |
| optimization | [`gradient.py`](../../api/algorithms/optimization/gradient.rst) | Jordan 量子梯度估计 |
| basics | [`oracle_algorithms.py`](../../api/algorithms/basics/oracle_algorithms.rst) | D-J、Bernstein–Vazirani、Simon 采样；Simon 的 GF(2) 消元在经典侧进行 |
| basics | [`number_theory.py`](../../api/algorithms/basics/number_theory.rst) | 有限规模模乘置换、QPE 求阶、连分数因子候选后处理 |
| qec | [`error_correction.py`](../../api/algorithms/qec/error_correction.rst) | 三位 bit/phase flip 重复码的编码与相干恢复 |

## 如何选择起点

如果输入是布尔函数，先看查询与搜索类别。如果需要估计概率或期望值，先看 `estimation`。已有 Hamiltonian 项分解时可以选择 Trotter；只有算子 BE 时，需要选择能够消费该访问模型的实现。

PDE 算法先构造空间离散化的输入，再选择 QODE 或线性化路线。数据在 QRAM 中并不意味着矩阵已经被 block encoded，输入适配步骤仍需明确。

## 运行展示目录

```bash
uv run python examples/algorithm_gallery.py
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

展示目录包含 22 个小实例，也包含同一算法的不同应用，例如振幅估计与量子计数。它们用于说明接口和读出方法，不能视为 22 种互不相关的算法。

## 实现依据

Fourier 加法采用 [Draper 的 QFT 加法构造](https://arxiv.org/abs/quant-ph/0008033)。振幅放大和估计采用 [Brassard 等人的框架](https://arxiv.org/abs/quant-ph/0005055)。QAOA 的 cost/mixer 分层依据 [Farhi 等人的原始算法](https://arxiv.org/abs/1411.4028)。求阶和经典因子后处理依据 [Shor 的构造](https://arxiv.org/abs/quant-ph/9508027)。

当前模乘使用有上限的置换合成，默认最多 8 位；它用于检查求阶接口与线路，不代表已实现可扩展的 Shor 模算术。VQE 和 QAOA 提供量子电路，经典优化器由应用选择。QSVT 接收调用方提供的相位序列，不包含通用相位求解器。

## 算法页面

每个算法一页，说明接口、输入模型、实现要点与验证证据的位置。页面按文件名排序；每页首行下方注明类别（C1–C6 或应用层）与所属模块，类别定义见[验证计划](../../development/validation-plan.md)。

```{toctree}
:maxdepth: 1
:caption: 算法页面
:glob:

*
```

## 数值验证

量子半定规划（`qsdp.py`，暂无专页）及其 `density.py` 内层原语（纯化访问、Gibbs 态制备）的论文级数值实验见 `tests/verification/verify_misc_algorithms.py`（misc_algorithms 组），全部在真实后端上执行，经典预言机全部独立（numpy 特征分解 / `eigvalsh` / 闭式公式）。

**实验设计**：
(a) 纯化见证——2×2 复密度矩阵（非对角元 0.1±0.05j）经 {obj}`gate_purification <oracq.algorithms.input_model.density.gate_purification>` 制备纯化态，偏迹环境后与原矩阵逐矩阵元对拍；
(b) Gibbs 纯化——对角 $H = \mathrm{diag}(1, -0.5)$（$\beta = 0.6$）与稠密非对角 $H$（$\beta = 0.8$）两个实例，`error = 0.05`，后选 signal == 0 分支的归一化约化态对照 numpy Gibbs 态 $e^{-\beta H}/Z$，另做 error = 0.4/0.2/0.1 的收敛扫描；
(c) 迹估计——{obj}`trace_estimate_circuit <oracq.algorithms.qml.qsdp.trace_estimate_circuit>` 对 Z/X/Y 三个 Pauli 观测量在复密度矩阵上的探针读数对照 numpy 迹 $\mathrm{Tr}(P\rho)$，以及 Gibbs 近似纯化路径的 {obj}`trace_from_joint <oracq.algorithms.qml.qsdp.trace_from_joint>` 联合解码；
(d) MMW 驱动——四约束可行性实例（等式 $r_z = 0.2$、$r_x = 0.1$ 拆成双边不等式，$\varepsilon = 0.08$），收敛后对返回的平均迭代 $\bar\rho$ 做独立的 PSD/迹 1/违反量核算；
(e) 单轮量子迭代——{obj}`iteration_circuits <oracq.algorithms.qml.qsdp.iteration_circuits>` 生成的 Gibbs 纯化 + 逐约束迹估计电路在参考执行器上给出的估计与经典估计器逐约束对照。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| purification-partial-trace | 2 量子位 | 四路径全振幅对拍 | 矩阵元最大误差 | 1.1e-16 |
| gibbs-purification-diagonal | 9 量子位，qsp_degree 3 | reference + rir-pysparq + originir-ext | 迹距离 / 后选成功率 | 3.0e-5 / 0.0965 |
| gibbs-purification-nondiagonal | 13 量子位，qsp_degree 3 | 同上 | 迹距离 / 后选成功率 | 1.9e-6 / 0.1130 |
| gibbs-error-scaling | error = 0.4/0.2/0.1 | reference | 迹距离 | 5.6e-4 / 8.7e-5 / 8.7e-5（均 ≤ error 且单调不增） |
| trace-estimate-pauli-xyz | ≤5 量子位 | 三路径 | Tr(Zρ)/Tr(Xρ)/Tr(Yρ) 估计 | 0.4 / 0.2 / 0.1（max 误差 5.3e-16） |
| trace-estimate-joint-gibbs | 10 量子位 | reference + rir-pysparq | Tr(Zρ_gibbs) 联合解码 | −0.42196 vs −0.42190（误差 6.0e-5） |
| qsdp-mmw-driver-feasibility | 4 约束，256 迭代 | 经典驱动 | 收敛 / 最大违反 / min 特征值 | ✓ / 0.0738 ≤ 0.08 / 0.4295 |
| qsdp-mmw-quantum-round | 2 约束单轮 | reference + rir-pysparq | 量子 vs 经典估计最大差 | 2.3e-6 |

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_misc_algorithms.py
```

产物：`out/verification/misc_algorithms.json`（24 个案例全过，本节对应 `purification-*`、`gibbs-*`、`trace-estimate-*`、`qsdp-*` 八个案例）。组内其余算法的数值验证见各自页面（密度矩阵指数化、QPCA、DQI、变分拟设、VQE、QAOA、重复码）。

## 数值验证（应用目录与结构层模块）

`applications/catalog.py` 与 `applications/gallery.py` 的全部条目由 nt_qlss_sde 组的 `tests/verification/verify_nt_qlss_sde.py` 覆盖：每个条目都在真实后端（rir-pysparq、adapter-pysparq，部分含 reference 执行器）跑通并与参考执行器逐振幅对拍，oracle 层面的数值正确性由底层算法各自页面的验证节承担（数论、CKS 等见本目录对应页）。

**实验设计**：gallery 22 例三后端（reference / rir-pysparq / adapter-pysparq）全振幅对拍；catalog 33 例中 31 例做 reference vs rir-pysparq 对拍（`stateprep_qram` 兼作 pysparq.rir 工作位复净的回归哨兵，另对 adapter 全态对拍），`qham_qode` 与 `qham_qpde` 两例因展开步数超 $10^6$、单后端运行约 60–90 秒，采用 spawn 子进程并发完成 reference vs adapter 对拍。这些条目是**组装示例**：其数值内容（模乘置换、求阶分布、CKS 解态等）以各算法页的独立 oracle 验证为准，目录层只钉死组装确定性，不再重复物理 oracle——这是目录/画廊条目以"结构测试 + 跨后端对拍"而非独立数值预言机覆盖的原因。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| gallery-*（22 例） | 小实例 | 三后端 | 跨后端最大振幅偏差 | 0.0（全部） |
| catalog-*（29 例） | 小实例 | reference + rir-pysparq | 跨后端最大振幅偏差 | ≤ 1.4e-16 |
| catalog-stateprep_qram | QRAM 角表 | 三后端 | 全态偏差 / 工作位残留 | 0.0 / 0.0 |
| catalog-qham_qode、qham_qpde | 39296 态 | reference + adapter | 跨后端最大振幅偏差 | 1.8e-9 |

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

产物：`out/verification/nt_qlss_sde.json`（`gallery-*` 22 例、`catalog-*` 33 例）。
