# 一般量子同伦分析（General QHAM）

> 类别 C6 · 模块 `pyqecclang.applications.qham` · 阶段 V1

## 概述

QHAM 把有限多项式演化 PDE 变成量子可解的线性系统。输入是自治的一阶时间演化系统

$$
u' = f + Lu + \sum_r B_r(u,\dots,u), \qquad u(0)=u_{in},
$$

允许任意有限分量、有限空间维、未知场及其空间导数的有限多项式、已知系数/强迫及其导数、以及完整单项式上的外层导数。生成器按 HAM（同伦分析）截断阶 m 递推出各阶 $U_i$，再对截断后的 HAM 系统做量子适配线性化（QCL），得到提升空间上的 $Y'=GY$，交给任意 QODE 求解器；输出选取物理块 $u_{sum}=\sum_i U_i$。构造对照 QHAM 论文（[arXiv:2411.06759](https://arxiv.org/html/2411.06759v2)），数学规则见推导参考。

张量字 $Y_{(a_0,\dots,a_{k-1})}$ 的有限闭包取 $\mathrm{grade}\cdot\sum a_j + \mathrm{rank} \le \mathrm{grade}\cdot m + 1$（grade = max(1, D−1)，D 为 PDE 最高多项式次数），非线性替换使 HAM 阶数和严格下降；闭包精确表示截断后的 HAM 系统，不再做一次 Carleman 式高阶截断。

## 接口与输入模型

入口分四段（`from pyqecclang.applications.qham import ...`）：

```python
PolynomialPDE.from_equations(equations, *, axes=("x",), label="polynomial_pde")   # Field/Known 表达式
QHAMPlan(pde, order)                                    # 惰性闭包：blocks / row_terms / offset / locate
Discretization(pde, Grid(axes, shape, spacing), known)  # 空间离散与经典参考
qham_input_model(plan, bindings, *, eta=-1.0, max_blocks=256, max_terms=4096)
```

- PDE 用 `Field`/`Known` 表达式代数构造（加减乘、非负整幂、`d(axis, order)` 导数；未知场分母与非整幂抛 `ValidationError`），`dumps`/`loads` 与 JSON Schema 往返。
- `bindings` 由 `gate_bindings(discretization, initial)`（小规模物化端口做对照，单端口 target ≤ 5 位）或 `structured_fd_bindings(discretization, initial)`（移位 + 分量选择 + 同点收缩，不物化 $N^r\times N^r$ 的端口矩阵）构造；也可用 `QHAMBindings.declare(plan, state_width, port_specs, *, initial_norm, ...)` 声明开放端口。
- input model 为 ODE（生成元 + 初态 + 物理输出窗口）：生成元以 BE 承载、初值以 SP 承载，开放槽位经 `bind` 晚绑定。

`qham_input_model` 返回 `QHAMInputModel`：

| 成员 | 含义 |
|---|---|
| `generator` / `initial` | 提升生成元的 BE 与提升初态制备 |
| `state_width` / `eta` | 物理态位宽 / 同伦参数 |
| `log_initial_norm` | 初始范数的对数（生成期用缩放权重避免直接对大幂求和） |
| `solve(qode, time)` | 调用 QODE 协议并选取物理块，返回 `StateOracle` |
| `dissipative_shift(shift=None)` | 显式整体移位 G − μI（μ ≥ α_G），适配 LCHS/CBMD，累计 `growth_shift` |

另有 `open_qham_input(plan, bindings, *, generator_alpha, generator_signal, eta=-1.0, name="QhamGenerator")`：整个生成元保留为未解析的抽象 BE（不给它伪造空主体），初态仍结构化生成；`taylor_qode(generator, initial, time, *, degree=2)` 提供不要求 Hermitian/耗散前提的有限 Taylor 候选，用于打通门级执行与对照。生成元模块属性含 `algorithm="qham_generated_qcl"`、`ham_order`、`pde_degree`、`tensor_rank_limit`、`raw_dimension`、`explicit_couplings`、`qham_plan` 与 `correctness="linearization_witnessed; solver pending"`。

## 实现要点

生成链为 PDE → HAM 递推（修正项权重 $-\eta(1+\eta)^{p}$、物理项权重 $1-(1+\eta)^{p}$）→ 有序张量字闭包 → 矩形端口放置与嵌入 → LCU 组合成 G。每个端口按多线性映射解释（L: V→V，B_τ: V 的 r 重张量积→V，F: 常量→V，端口系数已含 PDE 的符号与常数）：`place_port` 把端口放在指定张量位置并置换其余坐标，`embed_rectangular` 同时约束输入/输出两个窗口，避免矩形零填充污染相邻块。`QHAMPlan` 是惰性对象——二次无强迫、m = 20 时可描述 $2^{21}$ 个函数块而不物化，`row_terms`/`offset`/`locate` 可单独查询；显式生成超过 `max_blocks`/`max_terms` 预算时抛错，计划本身仍有效。

提升初态保留各张量字的相对范数 $r, r, r^2, \dots, r^K$（r = ‖u_in‖，有强迫再加 one 分支）：按缩放权重制备分支标签，对相应目标区间多次调用可重复的初值 oracle（不是克隆未知量子态），再从地址区间经布尔网络反算标签；每次对已知零输入的制备都复净 work，同一段工作寄存器顺序复用。u_in = 0 且有强迫时初态只剩 one 分量；无强迫的零初值对应零解，应在经典侧直接返回。

适用边界：单个 BE 的 target 包装受 64 位限制（`plan.max_rank * n > 64` 报错，惰性 plan 与逐行算法不受此限）；`correctness="linearization_witnessed; solver pending"` 表明已见证的是线性化代数，HAM/PDE 收敛认证与 QODE 求解精度属于求解器一侧；自动选择收敛的 h/m、时间依赖绑定与 IQHAM 外层重启未实现（完整边界见用户指南第 9 节）。

## 验证方案

类别 C6（组合骨架，判定准则见 `../../development/validation-plan.md` §2：组件契约满足 + 端到端小实例语义正确 + bind 不变）。三层证据位于 `tests/core/test_qham_general.py:QhamGeneralTests`：

- 结构：类的构造与属性断言；`test_rules_roundtrip_and_unsupported_nonlinearity`（PDE JSON 往返相等、`1/u` 与 `u**0.5` 拒绝）。
- 数值：`test_lazy_closure_and_rank_unrank`——度 2–4 × 阶 0–3 的闭包自洽（非线性替换使阶数和严格下降、offset/locate 互逆）；$2^{21}$ 块计划的序列化 < 3000 字符；显式块数超预算报错。`test_quadratic_dimension_and_nontrivial_eta_weight`——u′ = u² 时 block_count = 2^(m+1)、raw_dimension(4) = 5^(m+1)+3；η = −0.4 的同伦权重 0.24 / 0.4 / 0.4。`test_multidimensional_coupled_identity`——二维耦合 PDE 上 G·lift 与独立链式法则求值的残差 < 1e-10。`test_m1_reduces_to_previous_special_case`——m = 1 无强迫的生成元与旧特例 `qham_lift_m1` 的 BE 逐列对拍（16×16 幅度乘 α，places = 10）。另有 `test_chain_rule_identity_general_cases`（5 组 PDE 家族 × η ∈ {−1, −0.4, 0.2}，同一残差界）、`test_lazy_stencil_rows_match_explicit_linearization`（惰性行规则与显式矩阵 places = 11）、`test_structured_initial_preserves_relative_tensor_norms`、`test_rectangular_forcing_and_contraction_do_not_spill`、`test_dissipative_adaptation_is_explicit`、`test_zero_initial_forcing_and_reused_work` 补充初值、矩形嵌入、移位与工作区复用细节。
- 绑定：`test_open_qode_inputs_keep_operator_and_initial_oracles`——`QHAMBindings.declare` 声明的开放端口（`Qham_L`、`Qham_B_*`、`Qham_initial`）在 `model.solve` 之后保持未解析，经 `bind` 换成 gate 实现后无未解析槽且 IR 序列化往返不变；`open_qham_input` 保留的槽位恰为 `{QhamGenerator, Qham_initial}`。

## 已知缺口与计划阶段

与 `validation-coverage.md` 的应用层 `qfvm.py`/`qham/` 行一致：无登记缺口，阶段 V1。求解精度与收敛认证按 `correctness` 属性标注为 pending，属 QODE 求解器一侧的实现边界，不在本模块缺口列。

## 相关链接

- 源码：`src/pyqecclang/applications/qham/`（pde / linearization / reference / stencils）与 `src/pyqecclang/algorithms/qham.py`（量子组装）
- 用户指南：[一般 QHAM 自动生成](../qham.md)；数学推导：[QHAM 推导](../../reference/qham-derivation.md)
- API 参考：[QHAM](../../api/algorithms/qham.rst)、[PDE 模型与适配](../../api/applications/qham/pde.rst)、[结构化差分端口](../../api/applications/qham/stencils.rst)
- 同组页面：[量子有限体积法](qfvm.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
