# 框架调研笔记（第一阶段）

<a href="../../development/framework-survey.html">English</a> · **简体中文**

为审稿驱动的修订收集的工作笔记。每节末尾给出一句话结论，论文修改时不再需要回头查证。以下全部为 2026-10-02 收集的文档与论文层面证据，未做任何基准实测。"文档中未找到"指截至该日在公开文档（Qualtran 另含仓库目录结构）中未见，不涉及未发布代码。

## Qrisp BlockEncoding 论文（arXiv:2604.18276）

书目信息（据 arXiv 摘要页）：Matic Petrič 与 René Zander，《Block-encodings as programming abstractions: The Eclipse Qrisp BlockEncoding Interface》，2026 年，arXiv:2604.18276 [quant-ph]，v1 于 2026-04-20，DOI 10.48550/arXiv.2604.18276。配套文档：[Qrisp BlockEncoding 参考](https://qrisp.eu/reference/Block%20Encodings/BlockEncoding.html)。

逐点对照 oracq：

- 构造接口：`BlockEncoding(alpha, ancillas, unitary, num_ops, is_hermitian)`，工厂方法 `from_array`、`from_operator`（QubitOperator/FermionicOperator）、`from_lcu`（带平衡二叉树 `q_switch` 的 PREP/SELECT）、`from_projector`、`from_eye`。辅助位以模板声明，算子形状在作用时动态确定。oracq 最接近的构造器是显式 Pauli 展开（`matrix_pauli_encoding`，≤ 5 比特）、`projector` 以及稀疏/低秩输入模型——oracq 没有稠密矩阵或物理算子构造器。
- 代数：运算符重载 `+`/`−`（LCU）、标量 `*`（α → |c|·α）、`@`（乘积，α 可乘）、`kron`、`dagger`，另有 `qubitization()` 与 `chebyshev(k)`。这与 oracq 的块编码代数（`linear_combination`、`product`、`scale`、`tensor`、`adjoint_be`、`qubitization_walk`、`chebyshev_block`）直接对应，包括被跟踪的 α。
- 矩阵反演：`inv(eps, kappa, method={"QET","QSVT","GQSVT"})`，声明复杂度 𝒪(κ² log(κ/ε)) = 𝒪(κ log(κ/ε)) 次多项式 × 𝒪(κ) 后选择重复；相位合成委托给内部求解器，**未文档化度数上限或精度极限**。oracq 的 `qsvt_matrix_inversion` 文档化度数守卫（40），且每次合成经往返自检。Qrisp 另有成品求解器（`CKS`、`dalzell_inversion`、`lanczos_alg`）以及 `pseudo_inv`、`svt`、`poly`、`sim(t, N)`。
- 资源估计：`resources(*operands)` 返回一次编码幺正执行的门计数、深度与比特数，构建在 Jasp 追踪上；整程序计数用 `@count_ops`。oracq 的 `estimate_resources` 不执行任何内容，输出容错资源台账（Toffoli/Clifford+T 原子、待合成旋转、QRAM 查询/写入、逐 oracle 调用计数、开放槽位的不完整标记）。
- 执行：Qrisp 的 BE 可直接运行（`apply`、`apply_rus`、`expectation_value`）并集成 JAX；oracq 把生成与执行分离，经独立执行器路径读出。

oracq 仍独有的差异：带"已声明未绑定"oracle 与 checked late binding 的持久化开放程序；具名运行时 QRAM 资源；带结构化 issue 的静态契约检查；Python 公式到可逆线路的编译器。一句话结论：Qrisp BlockEncoding 接口是组件与分析层的直接对照——同样跟踪 α 的代数、反演与资源查询——差异在持久化开放程序工作流及其周围的科学约定与 QRAM 捕获。

## Qualtran 与 Catalyst 证据核对

用于收紧相关工作措辞。来源：[Qualtran Bloq 参考](https://qualtran.readthedocs.io/en/latest/reference/qualtran/Bloq.html)、[Register](https://qualtran.readthedocs.io/en/latest/reference/qualtran/Register.html) 与 [Signature](https://qualtran.readthedocs.io/en/latest/reference/qualtran/Signature.html) 页、[资源计数模块](https://qualtran.readthedocs.io/en/latest/reference/qualtran/resource_counting.html)、[序列化模块页](https://qualtran.readthedocs.io/en/latest/reference/qualtran/serialization.html)及仓库内 `qualtran/serialization/bloq.py` 与 `qualtran/protos/bloq.proto`、Qualtran 论文（arXiv:2409.04643），以及 Catalyst 的[首页](https://docs.pennylane.ai/projects/catalyst/en/stable/index.html)、[架构](https://docs.pennylane.ai/projects/catalyst/en/stable/dev/architecture.html)与 [dialects](https://docs.pennylane.ai/projects/catalyst/en/stable/dev/dialects.html) 页。

- Qualtran **已有**：寄存器签名（名称、dtype、形状、side）；层次化分解（`decompose_bloq`、`CompositeBloq`、`BloqBuilder`，构造期报寄存器不匹配）；无完整分解的符号调用计数（`build_call_graph`，SymPy 类型计数）；protobuf 序列化。未分解的 bloq 仅以"名称+签名+属性"的浅层条目存活序列化，反序列化须解析到已注册的具体 Python 类。
- Qualtran **未见文档化**的 late binding（在程序对象建成后替换子 bloq 实现并复查接口；CompositeBloq 文档明确不可变），也**未见文档化**的具名运行时 QRAM 资源（QROM/QROAM 是线路 bloq）。
- Catalyst **已有**高层结构化混合 IR（MLIR dialects：Quantum、Gradient、Catalyst 等），但它是编译器内部的瞬态产物：流水线降至 LLVM/QIR 并产出二进制，中间阶段只能经调试工具查看。**未见文档化**的"声明但不实现、事后在接口检查下解析"的量子 oracle 槽位（文档索引无 "oracle" 条目；`quantum.custom` 是完整指定的具名门）。

一句话结论：这些机制存在但停在不同位置——Qualtran 有签名、层次、符号计数与抽象 bloq 的浅层序列化，Catalyst 在编译器内部有高层混合 IR；两者均未文档化"持久化开放程序、绑定时科学约定检查、具名运行时 QRAM 资源捕获、重复部分绑定"这一完整流程。

## QSP 相位生成接口

来源：[pyqsp README](https://github.com/ichuang/pyqsp/blob/master/README.md) 与源码 docstring（`pyqsp/angle_sequence.py`、`pyqsp/sym_qsp_opt.py`、`pyqsp/response.py`），[pyLIQTR Features.md](https://github.com/isi-usc-edu/pyLIQTR/blob/main/docs/Features.md) 及其 `phase_factors` 包。

- pyqsp v0.2.0（MIT；依赖 numpy/scipy/matplotlib；`pip install pyqsp`）：输入为 **Chebyshev 基**系数；主入口 `QuantumSignalProcessingPhases(poly, method=…, chebyshev_basis=True)` 对度数 d 返回 d+1 个相位。相位约定：Wx 信号算子，U = e^{iφ₀Z} ∏ W·e^{iφ_k Z}，φ₀ 为最后作用的旋转，无 π/4 偏移；`method="sym_qsp"` 时目标落在 Im⟨0|U|0⟩。方法：`laurent`（求根补全，中度数机器精度，高度数不稳）与 `sym_qsp`（简约相位上带 FFT Jacobian 的 Newton 迭代，源自 arXiv:2002.11649 与 arXiv:2307.12468；文档称可稳定到"数千个相位"）。响应自评：`pyqsp.response.ComputeQSPResponse` / `PlotQSPResponse`。
- pyLIQTR 不包装 pyqsp：其 `Angler_opt` 是"非常直接地"取自 QSPPACK 的 L-BFGS 对称 QSP 生成器，目标是 **Re⟨0|U|0⟩** 且带显式 π/4 端偏移——其相位不经转换不能与 pyqsp 互换。旧的 scipy/mpsolve 求根路线正在弃用。
- PennyLane 更正："相位角外部提供"的说法已过时——`qml.poly_to_angles`（求根路线约到度数 1000，迭代 L-BFGS-B 更高）已存在，文档记载对 P(x) 约 1e-10 的线路级验证。

oracq 本阶段落地的契约选择：PhaseSynthesizer 接收实目标 f 的升幂实系数（单项式基）与可选的钉死虚部补全 h，返回反射约定时间正序相位；钉死模式要求精确复现 p = f + i·h，自由模式只钉死 Re p = f。pyqsp 适配器（仅自由模式）把目标转到 Chebyshev 基并把 Wx 相位映射到反射约定（逆转列表、端相位减 π/4、中间减 π/2、按 d mod 4 给首作用相位加偏移，d ≡ 3 mod 4 时目标取反）；每次调用由 `qsp_response` 网格校验独立核验。pyLIQTR 的 Re 目标 π/4 偏移约定未采用。一句话结论：相位约定映射是精确且机器校验的；适配器只越过合成度数守卫——单项式系数管线自身在度数约 40 处有可靠性极限，高度数目标必须在 Chebyshev 基下构造与校验（κ=8、ε=1e-2 反演目标即在度数 585 以此方式达到，往返残差 ≈ 1e-13）。

## 成熟度分级素材

四级：**L1** 有接口 · **L2** 能生成线路 · **L3** 能闭合执行（端到端模拟/运行）· **L4** 达到指定精度且有验证。文档层面分级，来源随框架附链接；"—" 表示所查文档中未见。

| 框架 | L1 接口 | L2 线路 | L3 执行 | L4 精度验证 | QLS | ODE/PDE |
|---|---|---|---|---|---|---|
| Qualtran | 是（[BE bloq](https://qualtran.readthedocs.io/en/latest/bloqs/block_encoding/block_encoding.html)） | 是（[GeneralizedQSP](https://qualtran.readthedocs.io/en/latest/bloqs/qsp/generalized_qsp.html)） | 部分（bloq 经典模拟；目的在于资源估计） | — | — | — |
| pyLIQTR | 是（[Features.md](https://github.com/isi-usc-edu/pyLIQTR/blob/main/docs/Features.md)） | 是（QSVT 线路类，OpenQASM 导出） | 部分（内部用 cirq；交付物是 Clifford+T 估计） | — | — | 仅经典（`clam`） |
| PennyLane | 是（[qml.QSVT](https://docs.pennylane.ai/en/stable/code/api/pennylane.QSVT.html)、`poly_to_angles`、FABLE、Qubitization） | 是（BlockEncode 仅模拟器） | 是（default.qubit 文档示例） | 是（`poly_to_angles` ~1e-10 验证；QSVT 演示对拍经典解） | QSVT 反演演示 | — |
| Classiq | 是（[qsvt_inversion 库](https://docs.classiq.io/latest/explore/algorithms/quantum_linear_solvers/qsvt_matrix_inversion/qsvt_matrix_inversion)、`qsp_approximate`） | 是（封闭综合引擎） | 是（模拟器与云后端） | 是（QSVT/HHL 示例对拍经典解） | HHL、QSVT、VQLS、绝热 | Poisson、LCHS、time-marching、CFD 示例 |
| UnitaryLab | 是（[算法手册](https://docs.unitarylab.com/en/docs/unitarylab-algorithms-user-manual/)） | 是（`.run()` 建线路） | 是（numpy/torch 后端） | 部分——Schrödingerization 的 `method='block'` 文档明确回退经典计算（[手册页](https://docs.unitarylab.com/en/docs/unitarylab-algorithms-user-manual/schrodingerization/)） | HHL、QSVT-QLSA、VQLS、AQC | Schrödingerization 热/平流（trotter + 经典） |
| Qrisp | 是（[BlockEncoding](https://qrisp.eu/reference/Block%20Encodings/BlockEncoding.html)） | 是 | 是（`apply_rus`，模拟+硬件） | 是（[inv() QSLP 示例](https://www.qrisp.eu/reference/Block%20Encodings/methods/inv.html)对拍经典解） | CKS、`inv()` | — |

oracq 自校准行：L1 是（类型化 oracle 接口与契约）、L2 是、L3 是（参考执行器加后端路径）、L4 对已见证 fixture 是（自校验相位合成；κ=8/ε=1e-2 反演经可替换合成器在相位级与 2×2 块编码级达到），QLS 是（CKS、VTAA-CKS），ODE/PDE 是（LCHS、Schrödingerization、CBMD、Carleman、QHAM/QFVM 族）。

一句话结论：四级分级能干净地区分模板、演示与已验证求解器；所查文档显示 PennyLane（模拟器）、Classiq（封闭引擎）、Qrisp 与 oracq 在 QSVT/QLS 上有完整 L1–L4 链，Qualtran 与 pyLIQTR 为部分链（面向资源估计），UnitaryLab 的 Schrödingerization `block` 选项有文档化的经典回退；ODE/PDE 求解器族仍是 oracq 的差异化覆盖。
