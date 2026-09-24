# 量子算法实现验证计划（审阅稿）

本文档系统规划 oracq 的算法验证体系：先按"正确的判定方式"对全部算法分类，再为每一类规定标准见证技术，最后给出全算法验证矩阵与实施路线。审阅重点：第 2 节的分类是否完整、第 4 节的横切机制是否值得作为基础设施投入、第 5 节的阶段划分与优先级。

## 1. 验证总原则

- **验证对象是 RIR 程序的语义**。算法是生成器，见证必须作用在生成的程序上（经参考模拟器或真实后端执行），而不是作用在生成器的内部数据结构上。
- **每个算法至少三层证据**：结构合法（IR 校验、签名、契约属性）、小规模数值见证（参考模拟器对拍解析解或经典参考）、可绑定性（抽象声明经 bind 绑定后语义不变）。
- **input model 三层一致性**是本仓库的特有不变量：同一算法的 abstract / gate / qram 绑定必须产生一致的可观察语义（允许 QRAM 字宽量化误差，需给出误差界）。
- **确定性优先**：概率性断言必须落在分布层面（闭式期望、支撑集、矩），禁止单一样本阈值断言；模拟必须使用固定种子或精确态矢量。
- 现有层级保持：L0 静态检查（ruff、布局、文档 warning-as-error）→ L1 结构 → L2 数值见证（tests/core）→ L3 对拍目录（applications/catalog.py）→ L4 真实后端（tests/integration，uniqc/pysparq 真机对拍，不用 mock）→ L5 文档 doctest。

## 2. 算法分类（按"正确"的判定方式）

| 类别 | 判定准则 | 标准见证技术 | 主要风险 |
|---|---|---|---|
| C1 精确离散语义 | 输出分布/作用酉与真值表或精确构造**逐点相等** | 穷举小实例全输入对拍；酉性恒等式 W†W=I；零空间/秘密串精确恢复 | 基本无 flaky |
| C2 近似连续语义 | 作用算子 ≈ 目标连续函数，误差有显式上界 | 与经典参考（解析解/矩阵指数/高精度数值）在容差内比较；**收敛性扫描**：误差 vs 度数/κ/Δt 拟合理论率 | 容差选取掩盖系统误差；数值病态（如 QSP 高阶合成） |
| C3 概率分布语义 | 输出**分布**等于闭式期望（均值、方差、支撑集） | 精确态矢量求分布，逐结果概率对拍闭式（Krawtchouk/二项等）；QAE 类验证估计值落在理论置信区间 | 采样抖动 → 禁止样本断言，只用精确分布 |
| C4 启发式/优化语义 | 严格优于随机基线，且小实例达到已知最优/理论分数 | 与随机猜测期望的分布级比较；小实例穷举最优解对拍；理论分数公式对拍 | 边界断言抖动（C3 类 DQI 已发生一次：恰好等于基线时 assertGreater 失败） |
| C5 数据访问层 | 查询语义逐点正确 + 三层绑定一致 + 复净 | 逐地址真值表对拍；abstract/gate/qram 三绑定一致性；LocalExit 强制复净检查；叠加态查询边缘分布检查 | QRAM 量化误差需显式误差界 |
| C6 组合骨架 | 组件契约满足 + 端到端小实例语义正确 + bind 不变 | 端到端模拟 + 绑定前后结果一致 + 结构属性（寄存器布局、Repeat 不展开）检查 | 组件各自正确但组合语义错位（寄存器映射错误） |

## 3. 全算法验证矩阵

状态：✅ 已有充分见证；🟡 有见证但有缺口；❌ 缺见证。模块按分层列出。

### infrastructure / 数据访问层（C5）

| 模块 | 内容 | 状态 | 缺口与计划 |
|---|---|---|---|
| oracles.py | XorDatabase、StatePreparation、SparseAccess、BlockEncoding 三层 | ✅ | 补"三绑定一致性"参数化测试（同一声明逐绑定对拍），目前各绑定独立见证 |
| data_loading.py | QROM / Select-Swap | ✅ | 资源公式已对拍；建议把 {obj}`qrom_cost <oracq.algorithms.input_model.data_loading.qrom_cost>` 纳入目录报告属性 |
| mathfunc/ | 经典函数→可逆线路 | ✅ | 补定点量化误差界的显式见证 |
| prepare_select.py | PREPARE/SELECT、alias | ✅ | alias 的 clean_work=False 已诚实标注；补 alias 三绑定一致性 |

### 精确离散算法（C1）

| 模块 | 算法 | 状态 | 缺口与计划 |
|---|---|---|---|
| oracle_algorithms.py | DJ/BV/Simon | ✅ | — |
| fourier.py | QFT、QFT 加法 | ✅ | — |
| arithmetic.py | 定点加减乘、布尔网络 | ✅ | — |
| number_theory.py | 模乘、order finding | ✅ | 小规模见证已注明不可外推为完整 Shor |
| error_correction.py | 重复码 | ✅ | — |
| walks.py | coined cycle walk | ✅ | — |
| graph_walks.py | 邻接 oracle、Szegedy、MNRS | ✅ | 环的失败已证实为理论结果；补 Johnson 图后需新见证（element distinctness） |

### 近似连续算法（C2）

| 模块 | 算法 | 状态 | 缺口与计划 |
|---|---|---|---|
| transforms.py | qubitization walk、QSVT 序列、OAA | ✅ | 序列约定已由 qsvt.py 钉死（逐点一致 1e-10） |
| qsvt.py | 相位合成、求逆、过滤、HS 模拟、定点搜索 | ✅ | **收敛性扫描缺失**：误差 vs 度数/κ 的批量曲线未自动化；合成度数病态边界（≳16）已有负例覆盖：tests/core/test_qsvt.py:PhaseSynthesisTests.test_input_validation 与 TransformWitnessTests.test_transform_input_validation（含奇偶性违例、|coeff|>1、端点未饱和、虚部奇偶、κ 病态、度数 0、t=0、Δ≥1、奇偶度数违例等场景） |
| hamiltonian.py | Trotter、Taylor、hamsim | ✅ | Trotter 阶数误差率未做扫描 |
| sparse.py / block_encoding.py | 稀疏 BE、BE 代数 | ✅ | — |
| qlss.py | CKS、Costa | 🟡 | 有范式见证；缺与 HHL 论文参考值的端到端对拍案例（可入目录） |
| ode.py / lchs.py / cbmd.py / schrodingerization.py / carleman.py | ODE/PDE 求解器 | 🟡 | 各自有见证；缺统一的"解析可解 ODE 族"收敛性基准（同一问题过全部求解器） |
| sde.py | Fokker–Planck | ✅ | OU 矩对拍；补 SDE 求解器与 LCHS 的端到端目录案例 |
| density.py | 纯化、Gibbs | ✅ | 收敛扫描已补：tests/core/test_density.py:GibbsTests.test_error_convergence_decreases（error ∈ {0.4,0.2,0.1} 单调不增+每档≤error，实测 5.6e-4/8.7e-5/8.7e-5）与 test_error_bound_uniform_in_beta（β ∈ {0.2,0.5,1.0} 网格每点 ≤0.1） |
| lowrank.py | DF/THC → BE | ✅ | α 上界紧性已补：tests/core/test_lowrank.py:DoubleFactorizationTests.test_alpha_matches_closed_form_eigenvalues（2×2 闭式独立对拍 places=10）、test_df_alpha_tighter_than_pauli（DF α=1.4 ≤ Pauli α=1.5）；ThcTests.test_thc_alpha_matches_hand_computed_bound（THC α=1.996 独立手算） |

### 概率分布算法（C3）

| 模块 | 算法 | 状态 | 缺口与计划 |
|---|---|---|---|
| estimation.py | QPE、QAE、量子计数、Hadamard/SWAP test | ✅ | QAE 置信区间声明未见证（只验证了点估计） |
| search.py | Grover、振幅放大 | ✅ | — |
| integration.py | Heinrich 求和/积分 | ✅ | 闭式均值对拍 + heinrich_rate 参数校验已有见证（tests/core/test_integration.py:SumPreparationTests / QuantumSumTests / RateTests） |
| gradient.py | Jordan 梯度 | ✅ | 线性精确 + 扰动收敛已有见证（tests/core/test_gradient.py:GradientTests.test_perturbed_linear_concentrates_with_grid_bits 含失败概率衰减率 q_{m+1} ≤ 0.34·q_m） |
| qpca.py | QPCA、密度矩阵指数化 | ✅ | 特征值读出峰对拍 + Δt 一阶误差率已有见证（tests/core/test_qpca.py:DensityMatrixExponentiationTests / QpcaTests） |

### 启发式/优化算法（C4）

| 模块 | 算法 | 状态 | 缺口与计划 |
|---|---|---|---|
| variational.py | QAOA-MaxCut、VQE、ansatz | ✅ | 单边见证；补"达到小实例最优割"的强见证 |
| dqi.py | DQI | ✅ | 分布逐点对拍闭式（Krawtchouk）；identity 译码器基线断言已修为允许边界等值 |

### 应用层

| 模块 | 内容 | 状态 | 缺口与计划 |
|---|---|---|---|
| qfvm.py / qham/ | QFVM、QHAM | ✅ | 目录对拍已有 |
| catalog.py / gallery.py | 33 案例对拍目录 | ✅ | 新算法（qsvt 求逆、graph_walks 搜索、data_loading、integration、qpca）未登记 |

## 4. 横切验证机制（基础设施投入项）

1. **不变量测试库**（新建议，tests/core/witness.py 或类似）：
   - `assert_unitary(program)`：参考模拟器上 W†W=I（抽样基态）。
   - `assert_uncomputation(program, work_regs)`：工作位复原（复用 LocalExit 机制）。
   - `assert_bind_invariant(abstract_program, bindings)`：同一开放声明的各绑定产生一致语义（误差界参数化）。
   - `assert_block_equals(be, matrix, alpha)`：BE 的 (0,0) 块对拍。
   收益：新算法见证从"各写一套"变为声明式组合，且回归口径统一。
2. **确定性策略**（规则化）：所有概率断言必须是对精确模拟分布的闭式对拍；理论值恰好落在断言边界时用 ≥/≤ 加数值容差，不允许严格不等（DQI 边界案例已按此修复）。
3. **收敛性扫描框架**：对 C2 类算法提供参数扫描工具（误差 vs 度数/κ/Δt），断言拟合斜率 ≥ 理论率的下界；结果存 JSON 供 validation.json 汇总。
4. **资源估算对拍**：凡提供成本函数（qrom_cost、QSP 度数公式、Trotter 步数）的模块，生成程序的结构属性必须与成本公式一致（双向检查）。
5. **目录登记规范**：新算法完成后必须登记 catalog 案例（gate/qram 两版输入），并接入 tests/integration 的真实后端对拍（不用 mock）。

## 5. 实施路线

| 阶段 | 内容 | 前置 |
|---|---|---|
| V1 | 不变量测试库（4 个断言原语）+ 在途模块（density/gradient/lowrank/integration/qpca）的见证按矩阵要求补齐（已完成，0.x 迭代） | 无 |
| V2 | 收敛性扫描框架 + qsvt/hamsim/ODE 求解器的收敛基准接入 validation.json | V1 |
| V3 | catalog 登记新算法（qsvt 求逆、MNRS 搜索、Select-Swap、Heinrich 积分、QPCA）+ tests/integration 真实后端对拍扩展 | V1 |
| V4 | 三绑定一致性参数化测试铺开（oracles/prepare_select/graph_walks/data_loading） | V1 |

## 6. 与现有验收的关系

`validation.md` / `validation.json` 继续作为版本验收记录；本计划落地后，V1–V4 的产物按版本写入其中。算法覆盖工作板（algorithm-coverage.md）记录"实现了什么"，本文档记录"如何证明正确"，两者行级对应。
