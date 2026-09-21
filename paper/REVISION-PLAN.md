# pyqecclang 论文修改计划

2026-09-21：当前工作以算法研究者的实现与分析框架为主线，见
`docs/development/research-workflow-plan.md`。以下内容保留为前一轮审阅历史；
其中定位、能力判断和待办状态应与当前主稿及新验收记录区分。

日期：2026-09-16。来源：外部审稿意见（5 Major + 5 Minor）+ 编辑部前对本轮"对标软件调研"的扩充要求（商业平台与主流开源框架，见 C 节）。

状态图例：`[已改 2026-09-16]` 本轮已修改并重新编译验证；`[待办]` 已列入计划尚未执行；`[需作者决策]` 涉及发布、署名等只有作者能定的事项。

---

## A. 审稿 Major issues

### A1 语言本体不可获取（最严重）`[需作者决策 → 待办]`
- 事实：正文无 repository URL / DOI / availability statement；GitHub 搜 `pyqecclang` 0 结果、PyPI 404。473 实验 campaign 与全部 doctest 清单目前无法核查。
- 行动：
  1. 作者发布仓库（含 tag 对应论文版本 0.8.0）与 PyPI 包；
  2. 论文增加 Data/Code availability 小节：repo URL、版本号/commit pin、PyPI 链接、后端 revision ledger（`backend-revisions.json`）的公开位置；
  3. Abstract 或 Introduction 首段加一句 artifact 链接。
- 注：后端可获取性没问题（PySparQ 在 PyPI、uniqc 在 GitHub），仅需语言本体发布。

### A2 两个已知缺陷未修复即投稿 `[已修复 2026-09-17]`
- 缺陷 1：OAA（oblivious amplitude amplification）迭代器缺收尾 U。
  - 修复：`src/pyqecclang/algorithms/transforms.py::oblivious_amplification` 由 `[R U† R U]^it` 改为 `U·[R U† R U]^it`（标准三查询形式的迭代推广）；零信号块满足切比雪夫恒等式 `ΠWΠ = B(4B†B − 3I)`，V/2 夹具幅值 0.5 → 1.0 精确（sin 3θ）。
  - 回归：`verify_fourier.py`（unitary ≤4.6e-16；半块三档迭代幅值 sin 3/5/7θ 精确）、`verify_hamiltonian.py`（一般块恒等式 3.6e-16；库算子 = 脚本组装标准序列，逐振幅偏差 0.0）。
- 缺陷 2：Schrödingerization 动量符号错误（恢复时间反演流）。
  - 修复：`src/pyqecclang/algorithms/schrodingerization.py` 生成元动量项系数 `+1 → −1`（K' = −P⊗H1 − I⊗H2）；正向流恢复：精确演化网格误差 1.3e-16（反演流失配 0.43 作方向性证据），端到端残余 = Nyquist 模 Taylor 余项 3.3e-2（degree 4，随阶数下降）。
  - 回归：`verify_ode.py` decay 案例改为通过性断言；`schrodinger_sign_flipped_qode` 独立重组装作为符号回归钉。
- 重跑：全 13 组 campaign 重跑，475 例全 PASS（473→475：OAA 半块拆为 m=1,2,3 三档）；out/catalog 重生成；tests/core 275 绿；ruff 绿；docs 构建通过。
- 论文同步：摘要恢复强表述并改为"两缺陷已修复 + 回归钉"；verify-matrix Fourier/ODE 行与 status-matrix 脚注 c、findings 段、discussion limitations、附录计数同步；docs/manual 三页（schrodingerization / oblivious-amplification / differential-equations）同步修复叙事与新数字。

### A3 漏引最接近的先行工作 `[已改 2026-09-16]`
- LIQUi|>（arXiv:1402.4467）：最早的语言+编译+资源估计+线性系统演示的垂直栈。已加入 §related 的语言谱系段并引用。
- Qualtran：bloq 抽象是最接近的结构先例，此前只在资源估计语境被引用。已在 §related 新增"容错算法分析库"段落，把 bloq 与我们的 open module 做正面对位。
- pyqsp（ichuang/pyqsp，2021 起）：QSP 相位角数值生成包。已在 QSP 引擎与 §related 处引用，作为相位计算先行工作；我们的差异点（自验证 + 电路内集成）相应写明。
- qsp4pde（2026，QSP-PDE 求解 Python 包 + 配套论文）：已加入 §related 覆盖表脚注/正文。
- 同时新增：OpenFermion、PennyLane（含其 QSVT/FABLE/GQSP 子例程模板的准确表述）、pyLIQTR、isQ、CUDA-Q、Amazon Braket SDK、Classiq、UnitaryLab、Azure Quantum RE、InQuanto、Phasecraft（详见 C 节）。

### A4 覆盖对比表计分规则不对称 `[已改 2026-09-16]`
- 问题：Table（ode-coverage）只把"框架官方标准库"计入对手、排除其生态（Qualtran/pyLIQTR/Classiq library/PennyLane 模板等），而 pyqecclang 全库计入。
- 修改：
  1. 表题与正文明确"官方标准库口径"，另加一段"生态口径"说明：按生态计，QSVT/block-encoding 在多个生态存在（Qualtran GQSP、pyLIQTR、PennyLane 模板、Classiq library），我们的差异是语言级集成与六族 ODE/PDE 求解器；
  2. 新增平台对比表（tab:platform-comparison）覆盖上述生态，逐项给出 ✓/partial/– 并附引文；
  3. 结论 "to our knowledge, unmatched" 已软化为"就我们所知没有系统同时做到……"的枚举式表述。

### A5 表达力证据全部自写，缺外部锚点 `[已完成 2026-09-17]`
- 生态部署：`~/projects/pyqecclang-dev/`（9 个独立 venv：qiskit/qualtran/pyliqtr/pennylane/qrisp/pyqsp/qsp4pde/unitarylab/classiq，全部 pip freeze 存档；classiq 综合需账号，按协议降级为不做综合）。
- 基准套件：6 个已发表锚点任务 T1–T6，预注册协议（LOC=非空非注释行、各框架文档化最高层路径、infeasible 不记 0、每行有正确性见证与冻结版本），pyqecclang 侧脚本入仓 `tools/expressiveness/`，多框架 bundle 在 pyqecclang-dev/benchmarks/（含每任务 SPEC.md 与 results.json）。
- 结果（LOC / 状态）：T1 pyqsp 30·ok vs pyqecclang 57·ok（但相位合成度数上限暴露）；T2 PennyLane 65·ok / Qualtran 78·ok / pyqecclang 70·**failed**（κ=8 规格需 b=293，库上限 b=9——真实能力缺口，已写入论文与 future work）；T3 Qiskit-from-scratch 77 / Qrisp 41（Qrisp 自带 HHL）/ pyqecclang 107·ok（4.7e-3）；T4 pyqecclang 87·ok 原生 LCHS vs Qiskit/Qualtran **infeasible**（Qiskit 从零 84 行可做最小变体）；T5 UnitaryLab 算法库 59·ok / Qiskit 从零 59 / pyqecclang 68·ok；T6 qsp4pde 30·ok / pyqecclang 98·ok。
- 论文落点：附录新增 \section{A reimplementation study}（表 tab:expressiveness + 三条解读 + 效度威胁）；related.tex 算法库段落加指向；discussion future work 加 Laurent 相位合成。审稿人要的"外部锚点 + 重实现对比代码量"完成。

---

## B. 审稿 Minor issues

1. `wang2026cbmd` 作者不全（Crossref/期刊页核实共 7 位）`[已改 2026-09-16]`：已补全为 Wang, Liu, Xue, Zhuang, Dou, Chen, Guo（详见 refs.bib）。
2. Table 1（qpl-comparison）引用了未印刷的脚注 e `[已改 2026-09-16]`：脚注 e 已补定义。
3. catalog 计数 32 vs 33 不一致 `[已改 2026-09-16]`：仓库实际 33 个 case（`out/catalog/` 33 个目录 + index.json），`casestudies.tex` 的 32 已改 33。Appendix B "four families" 与正文 "six" 的口径冲突 `[已改 2026-09-16]`：已改为"六个 ODE 族与 PDE 路径中的四个完整清单（Carleman、LCHS、CBMD、QHAM）"。
4. pyqecclang/QECC.Lang 命名易被读成量子纠错 `[需作者决策]`：选项 (a) Introduction 加一句命名说明（成本低，建议）；(b) 改名（成本高，影响引用连续性）。等作者定夺后执行。
5. 作者占位符与空 acknowledgments `[需作者决策]`；提交包内 `paper/qecc-lang/`（Overleaf 同步副本，含重复的 main.tex/sections/）与 `paper/` 并存 `[待办]`：提交前清理或同步，避免审稿人看到重复文件。

---

## C. 本轮新增：对标软件调研（已写入 §related）`[已改 2026-09-16]`

调研时间 2026-09-16，事实均经官网/GitHub/arXiv 核实；已写入论文的部分：

### 商业/工业平台
| 平台 | 定位 | 与本文的关系（写入 related work 的要点） |
|---|---|---|
| UnitaryLab（酉术量子，上交大背景，2025-11 发布 1.0） | "量子科学计算平台"：Agent 驱动 Web Studio + 闭源模拟器 SDK + MIT 算法库 | 最接近的商业同类：Schrödingerization PDE 求解器、QSVT-QLSA、HHL、多种 Ham-sim；无 QRAM/资源估计、无开放程序工件概念 |
| Classiq（以色列） | Qmod 高级 DSL + 约束驱动电路综合（闭源平台、开算法示例库） | 商业侧 FT 科学算法覆盖最全（QSVT 求逆、HHL、ODE 目录、GQSP、QPE）；无开放 oracle/范式契约，综合即展开 |
| Microsoft Azure Quantum Resource Estimator | 语言无关的容错资源估计（Q#/QIR/logical counts 输入，QDK 开源） | 资源估计轴的最强对手，论文已在 reserved 轴呼应；引用官方文档 + Beverland et al. |
| Quantinuum InQuanto | 专用量子化学平台（闭源 Python） | 领域垂直栈代表（QPE 含 QEC-QPE 叙事）；非语言、非通用科学计算 |
| Phasecraft（英国） | 算法公司，材料模拟论文交付 | 无语言/平台产品；作为"算法-资源缩减"研究路线对照 |
| NVIDIA CUDA-Q、Amazon Braket | 量子-HPC 执行基础设施（CUDA-Q 开源 Apache-2.0） | 一句话对照：执行与异构层，把算法/语言留给上层 |

### 主流开源框架（GitHub）
| 框架 | stars（2026-09-16） | 写入 related work 的要点 |
|---|---|---|
| Qualtran（Google） | ≈389，活跃 | bloq 抽象 + block-encoding/GQSP/qubitization/Ham-sim 构件 + 资源估计；无打包 QLSS/ODE 求解器；非语言（构造库） |
| pyLIQTR（MIT-LL/USC ISI） | ≈47 | BlockEncodings + QSP/QSVT 电路 + Clifford+T 资源估计；clam 仅经典 ODE 积分；引用 Zenodo DOI |
| OpenFermion（Google） | ≈1736 | 电子结构专用；无 QLSS/QSVT/PDE；作为"领域输入模型库"对照 |
| PennyLane（Xanadu） | ≈3466 | 有 QSVT/FABLE/GQSP/Qubitization 子例程模板与实验性 ftqc 模块；定位可微分混合编程，无打包 FT 科学算法——注意论文表述已按此校准，勿写成"PennyLane 无 QSVT" |
| Qiskit 生态 | 主仓 ≈7801 | 主库确认无 QSVT/block-encoding（2026-09 main 分支核实）；HHL 随旧 algorithms 模块退役且未进 qiskit-algorithms；qiskit-nature 仍在维护但非 FT |
| isQ（中科院系） | 编译器仓 ≈15（2023-10 后沉寂） | 独立 QPL 软件栈，IEEE TQE 2023；作为独立语言存续性讨论点 |
| pyQPanda（本源） / MindQuantum | ≈1210 / Gitee ≈3124 | NISQ 变分为主，HHL 仅演示级；无 QSVT/block-encoding/资源估计 |
| pyqsp（ichuang） | ≈141 | QSP/QSVT 相位角数值生成先行工作（审稿点 A3） |
| qsp4pde | — | 2026 QSP-PDE 求解包 + 配套论文（审稿点 A3） |

另：Mitiq（GPL-3.0，误差缓解）、TFQ、QuTiP、cuQuantum、bloqade、Catalyst、Qibo、Tangelo 等经逐一判定不属于"量子科学计算编程框架"范畴，不写入。

### 论文内落点
- §related 新增两段（容错算法分析库；商业与工业平台）+ 平台对比表；
- §related 的 PennyLane/Qiskit 表述按调研校准；
- refs.bib 新增 16 条（LIQUi|>、isQ、OpenFermion、PennyLane、pyLIQTR、pyqsp、qsp4pde、UnitaryLab、Classiq、Azure RE、CUDA-Q、Braket、QPanda、MindQuantum、InQuanto、Phasecraft），全部经 arXiv/Crossref/Zenodo/官方页核实（2026-09-16）。

---

---

## E. T2 性能强化设计（QSP 高度数相位合成；已设计未执行，2026-09-17）

### E0'. 分层原则（2026-09-17 修订）

本论文是语言特性论文；T2 相关事项必须先分层，避免把"捆绑数值例程的质量"当成"语言能力"来比：

| 事项 | 层面 | 论文应有的处理 |
|---|---|---|
| 相位以普通值进入 QSVT 构造器；协议生成器接受注入函数 | **语言层（已存在）** | 这是论点本身：接缝已留好，落后例程可整体替换而语言/IR/验证不动 |
| DK 求根条件、度数护栏 40、自检网格 | **实现层** | 如实披露为 bundled routine 的属性（primitives.tex 已改为此口径） |
| Laurent 分治 / GQSP / Newton 精修等算法 | **实现层（别人已解决）** | **复制成熟实现**（vendoring，首选 pyqsp 的 Laurent，MIT，注明来源），经我们的往返验证管线收编；不以"自研数值更强"为叙事 |
| C2 Newton–Schulz 分段组合 | 算法演示层 | 可选的语言组合性展示，不作为 T2 修复 |
| C1 κ 配额 | 基准卫生 | 独立记录为 SPEC 修订 |

### E0. 问题精确定位

T2 失败链（κ=8, ε=1e-2 规格，求逆目标 b=293 → 度 585）：

1. **度数护栏**：`qsvt.py::_MAX_DEGREE = 40`（`_GRID=4096` 自检网格，每行 3-4 次矩阵乘，度 585 的自检本身可承受；护栏是保守值不是硬极限）。
2. **真正的瓶颈是数值条件**：补多项式路线（Haah 型 complement + 求根 + 剥离）在 b≥10 失败——求逆 J_b 多项式的补多项式根呈**簇状**（oscillatory 目标的典型病态），纯 Python Durand–Kerer（`_roots`，float64）在簇内相对间距 ~1e-7 处收敛不到 `tol=1e-30`，`_cluster`（rel 2e-5）随后配对失败。实测最大可合成 b=9（度 17）。**b=9 时 J 多项式在谱边相对误差 0.75⁹≈7.5% → 方向误差 3.2e-2**，超 1e-2 阈值。
3. **参照系**：pyqsp（Laurent 精确分治，Chao et al. arXiv:2003.02831）合成到度 853；Qualtran（GQSP 闭式相位 + Chebyshev 系数截断，有效度 199）直接吃 b=293；PennyLane 单项式基在 b≳60 就 float64 灾难消去（b=59 为其上限）——**病根是"基/求根"选择而非 QSP 本身**，各家解法不同。
4. **一个被忽略的快速修正**：T2 的实际矩阵 A=[[1,−1/3],[−1/3,1]] 的 κ=2.25，按库规则 b=ceil(ln(1/ε)/−ln(1−1/κ²))=**21（度 41）**；基准 SPEC 把 T1 的 κ=8 一并套给了 T2。κ 正确配额后所需度从 585 降到 41——只比护栏高 1，**唯一障碍仍是 b≥10 的求根失败**。

目标分两档：**T2'（κ 配额修正后）**= 度 ~41 稳定合成；**T2（原规格 κ=8）**= 度 ~585 机精度合成（对齐 pyqsp/Qualtran）。

### E1. 方案选项

**方案 A：修现有管线的数值条件（小改，覆盖 T2'）**
- A1 求根器升级：Durand–Kerer → **Aberth–Ehrlich 同时迭代**（簇状根收敛显著更好）或 companion 矩阵特征值（`numpy` 已是执行环境依赖？注意：核心纯 Python，此处可退化为纯 Python Jacobi？不可——建议 A1 只在可选 numpy 存在时启用，保持依赖策略），外加**Newton 单根抛光**（DK 输出做 3-5 步 Newton）。
- A2 精度提升：`mpmath`（可选依赖，try-import，float64 回退）下构造补多项式 + 求根，50 位精度，相位最后舍入 float64。自检（`qsp_response` 往返）仍在 float64 上做——**验证门槛不变**，只放宽合成路径的内部精度。
- A3 护栏：`_MAX_DEGREE` 40 → 64（A1/A2 落地后），护栏语义改为"自检通过即接受"。
- 验收：b=21..30 合成往返残差 ≤1e-9（现 b=9 已达 8.7e-7）；T2' 通过（度 41，方向误差预计 <1e-2 的数倍裕量——度 41 时 J 多项式谱边误差 (1−1/κ²)^21≈0.8%…0.01 级，需实测）。
- 工作量：1-2 天。风险：低；`_q_from_roots` 的簇配对逻辑（rel=2e-5）在高精度根上可能要重新调参。

**方案 B：换合成算法（中改，覆盖 T2 原规格）**
- B1 **Laurent 精确分治**（pyqsp 同款，arXiv:2003.02831）：把 P(x) 写成 z=e^{iθ} 的 Laurent 多项式，补多项式 Q 由对称结构逐因子分解——不对病态实多项式求根，而对结构化对象做**精确**低秩分解。pyqsp 实测度 853 机精度，是最被验证的路线。落点：`qsvt.py` 新增高阶路径（度 >64 走 B1，低度保现有路径与行为兼容）；可移植 pyqsp 的 `angle_sequence.py` 核心逻辑（MIT），**必须去依赖化重写并保留我们的往返自检**（这是我们相对 pyqsp 的差异点：合成进验证管线，不合格相位进不了 IR）。
- B2 **相位 Newton 精修**（QSPPACK，arXiv:2002.11649）：任意初值（A 方案的粗相位）上对相位变量直接 Newton——需要 QSP 乘积对相位的导数（layer-stripping 递推的伴随反传，~100 行）。作为 B1 的兜底/交叉验证，或独立小改方案。
- B3 **GQSP 闭式构造**（Motlagh–Lin）：广义 QSP 的相位有闭式公式（每步两个复数的辐角公式），**完全无需求根**；Qualtran 正是这样吃到 b=293 的。代价：需要新增 GQSP 装配器（三相位约定 θ,φ,λ 的交错序列，builder 层面与 `qsvt_sequence` 同构）+ 对应响应验证公式；理论上限最高、实现面也最大。
- 验收（任一）：b=293（度 585）合成往返 ≤1e-9；T1/T2 原规格通过；新增合成路径的单元测试 + 与现有低度路径的交叉一致（同目标多项式两条路径相位一致到相位等价类）。
- 工作量：B1 移植+重写 3-5 天；B2 2-3 天；B3 5-8 天（含装配器与验证公式）。

**方案 C：算法级绕行（不改合成器，改"怎么到达 1/x"）**
- C1 **κ 配额修正**（零代码）：T2 基准 SPEC 按实际矩阵 κ=2.25 配额 b=21。这是基准协议修正而非性能改进，单独立项记录（诚实注明原 SPEC 过度配额）。
- C2 **Newton–Schulz 分段求逆**（语言层亮点方案）：迭代 X←X(2I−AX̃)，每段只是度 3 的多项式模式（相位合成毫无压力），段间用**刚修复的 OAA** 放大成功幅值——合成度数恒为个位数，调用次数 2^k·O(1) 增长。这是最能体现 pyqecclang 组合性主张的方案（"协议组合替代单块高度数"），也是唯一把 Major 2 修复直接转化为能力展示的路线；数值上 κ=8 约 3-4 段到 1e-2~1e-3。缺点：不是"单序列 QSVT 求逆"，与 T2 锚点（Martyn 单序列演示）形式不同，只能作为 T2 的补充行（T2-NS）而非替换。
- C3 任务级替代：同一线性系统已由 QLSS/CKS 路径在 T3 以 4.7e-3 通过——库"能解这个任务"，缺的只是 T2 特定的相位数值层。作为论文表述的上下文，不作为修复方案。

**方案 D：外部相位注入 + 我方验证（架构层，1 天）**
给 `qsvt_matrix_inversion`/QSVT 构造器加 `phases=` 注入点（协议组件替换，R8 叙事的现成实例）：外部相位（如 pyqsp laurent 输出）经我们的 `qsp_response` 往返验证后才进入 `qsvt_sequence` 装配。基准侧以"pyqecclang(装配+验证) × pyqsp(相位)"复测 T2，LOC 与相位来源如实标注。论文叙事从"引擎度数不足"转为"**合成可插拔、验证不可绕过**"——这与 §primitives 的协议论述完全一致。

### E2. 推荐组合（决策建议，按分层原则修订）

1. **主路线（2-3 天）：vendoring**——把 pyqsp 的 Laurent 精确分治（`pyqsp.angle_sequence`，MIT）去依赖化移植为 `qsvt.py` 的 `synthesize_phases_laurent(...)` 高度数后端（度 >64 走此路径，低度保现有补根路径与逐位兼容），来源与许可在模块 docstring 与 docs 注明；**往返自检照旧**（我们验证、不信任来源）。这不是"补短板的权宜"，而是论文论点的实例化：合成可插拔、验证不可绕过。
2. **语言层收尾（半天）**：给 `qsvt_matrix_inversion` 等构造器加 `phases=` 直通参数（纯 API 补充，接缝本就存在）；基准侧以 "pyqecclang(装配+验证) × pyqsp(相位)" 复测 T2 一行（LOC 与来源如实标注）。
3. **基准卫生（零代码）**：C1——T2 SPEC 按实际矩阵 κ=2.25 配额 b=21（度 41）另测一行，修订注明。
4. **可选展示（2-3 天）**：C2（Newton–Schulz + OAA 分段求逆）作为语言组合性案例（T2-NS 行），顺带强化 Major 2 修复的价值叙事。
5. **不做**：自研 Aberth/高精度求根来"提升我们数值性能"（A 方案降级为 vendoring 期间的临时缓解，仅当移植受阻才启用）；B3 GQSP 自研（与 vendoring 收益重叠）。

### E2'. 论文口径（已同步 2026-09-17）

primitives.tex（"引擎是可能来源之一，度数上限是本例程属性而非设计属性"）、related.tex 指针句、附录 findings 第二条（"分层报告：失败属于捆绑例程，语言层接缝已在"）、discussion future work（"vendoring 属实现升级，语言/IR/验证零改动"）已全部改为分层口径。附录对照表的 Err 列定位为**通过性见证**而非性能比较；三条解读以第一条（各系统能让你说什么）为主。

### E3. 验证与同步（任一阶段落地后）

- `tests/core`：新增合成路径单元测试（簇状根 fixture、b 扫描、双路径交叉一致）；
- `tests/verification`：hamiltonian 组 QSP 案例扩展高度数用例（现有度 ≤40 行为必须逐位兼容——estimation/eigenstate-filter 等消费方不回归，重跑 hamiltonian+estimation 组）；
- 基准：T1/T2 按原预注册协议复测，results.json 更新（C1 修正另记 SPEC 修订）；
- 论文：primitives.tex 度数表述、附录 T2 行与 findings、future work 删"Laurent 合成"或改为"已落地"；docs（qsp 页）同步。

### E4. 风险

- B1 移植的边界情形（非对称目标、奇偶性混装）需要 pyqsp 测试语料对拍；
- 高度数自检成本：度 585 × 4096 网格的 `qsp_response` 往返 ~分钟级（可稀疏化网格或分批验证）；
- 与 `_cluster`/`_q_from_roots` 的耦合：A 方案调参可能牵动低度路径，需回归（低度相位逐位不变作为验收项）。

## D. 验证清单
- [x] `pdflatex + bibtex` 全文重编译，无未定义引用/新 warning（本轮）
- [x] A2 修复后（2026-09-17）：全 13 组 campaign 重跑 475/475 PASS（out/verification/ 刷新；OAA 半块拆 m=1,2,3 三档，计数 473→475）；out/catalog 重生成；tests/core 275 PASS；tests/integration 30 PASS（真实后端）；ruff 绿；docs 构建（HTML+doctest，warning 即错）绿
- [ ] A1 发布后：availability 小节 + 链接校验
- [x] 运行手册：campaign/integration 解释器 = `~/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python`（同时可 import pysparq/uniqc/numpy/scipy）+ `PYTHONPATH=src`；PySparQ JIT 需要 PATH 上有 `g++`（本机无系统 g++，用 `ln -s ~/projects/qcfd-dev/.tools/envs/devenv/bin/x86_64-conda-linux-gnu-g++ /tmp/gxx-shim/g++` 后 `PATH=/tmp/gxx-shim:$PATH`）
- [x] A5 表达力锚点实验完成（2026-09-17）：T1–T6 全部落地，结果见 ~/projects/pyqecclang-dev/benchmarks/t*/results.json 与论文附录 reimplementation study；T2 我方失败（相位合成度数上限）——强化方案见 E 节（已设计未执行）
- [ ] 提交前清理 `paper/qecc-lang/` Overleaf 副本与占位符
