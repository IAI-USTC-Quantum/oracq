# 量子有限体积法（Quantum Finite Volume Method, QFVM）

> 类别 C6 · 模块 `pyqecclang.applications.qfvm` · 阶段 V1

## 概述

QFVM 从经典流场数据构造量子线性系统求解（QLSS）所需的全部输入：矩阵元素 oracle、稀疏位置 oracle 与右端残差态制备，再把同一个问题交给可替换的 QLSS（声明稀疏输入的 CKS 基础路线，或声明块编码输入的 Costa 路线）。接口划分对照 QFVM 论文（[arXiv:2102.03557](https://arxiv.org/html/2102.03557v1)，见输入模型审查）。当前实现针对一维周期网格上的 Euler 方程：三个守恒量（密度、动量、能量）加 frozen-Roe Jacobian，不覆盖原论文的全部网格、边界和物理模型。

线性系统取物理坐标上的 Hermitian 扩张 `D = [[0, M], [M.T, 0]] + padding_value·I_pad`（M 为 Roe 块），求解后在输出侧选取物理坐标块；这个扩张服务于线性求解，不能当作原生成元的等价 QODE 演化。

## 接口与输入模型

```python
roe_qfvm_inputs(*, cell_width=2, fmt=None, angle_width=8, prefix="RoeQfvm")
qfvm_sparse_access(inputs, *, padding_value=1.0, **entry_options)
roe_qfvm_block_encoding(inputs, *, amax=8.0, padding_value=1.0, **entry_options)
roe_qfvm_problem(inputs, *, spectrum, rhs_norm=None, amax=8.0, padding_value=1.0, **entry_options)
roe_qfvm_step(inputs, qlss, *, spectrum=None, rhs_norm=None, **options)
bind_qfvm(program, inputs)
qfvm_memories(inputs, flow, *, amax=8.0)
```

- `roe_qfvm_inputs` 返回 `RoeQfvmInputs`：六个抽象数据库（rho / momentum / energy / geometry / theta / 残差值）加一个抽象右端态制备。`cell_width` 范围 2..25（至少四个单元），`angle_width` 范围 1..64；矩阵维数为 `2**width`，其中 `width = cell_width + 3`（3 个守恒量 + 1 个补齐分量）。
- input model 分层：流场、几何与角度数据是 QRAM（`qram_database` 绑定）；矩阵访问产物是 SO（`qfvm_sparse_access` 返回 `SparseAccess`：位置 + 元素 oracle，稀疏度 9）；右端是 SP（`rhs_qram_preparation` 的符号残差树）。
- `roe_qfvm_problem` 组装 `LinearSystem`：`spectrum` 必须是显式 `SpectralPromise`（与 padding_value 合并取更宽的界），`amax`、谱界与矩阵性质是调用者声明，须适用于实际量化后的矩阵；数据假设记录在 `data_assumptions`（相干 XOR QRAM、快照固定、原始场量而非矩阵表等）。
- `roe_qfvm_step` 要求声明 input_model 的 `QLSSProtocol`（裸 callable 拒绝），返回 `qlss.solve(problem)`；`bind_qfvm` 把抽象槽替换为 QRAM 数据库绑定，`qfvm_memories` 从 `RoeFlowData` 快照组装内存表。

`RoeQfvmInputs` 的属性：`width`（= cell_width + 3）与 `geometry_width`（= width + 4 + cell_width + 7，打包 neighbor / reverse / source / rowvar / colvar / band / valid 七段）。

## 实现要点

矩阵元素不预存成表：`roe_entry` 查询 source 单元及其左右邻居的守恒量数据库，对两个界面调用可逆编译的 `roe_face`（`frozen_roe_face` 纯函数经 mathfunc `compile_function` 编译为定点算术，含 sqrt / div / mul / select），按 band 选 west / center / east 通量差，质量项经 `component_equal` 布尔网络只写入分量对角；算术失败（status ≠ 0）时元素 totalize 为零。位置 oracle 只编码邻居、槽位与分量索引（`geometry_cells`，不含任何矩阵值）：每列 9 个结构槽位（3 邻居 × 3 分量）经 9 次相干值转置（`value_transposition`）给出 CKS 需要的原地索引置换；补齐分量（变量 3）为对角 `padding_value`，其余无效槽位返回零元素，周期边界。右端用 `ptheta_cells` 把残差值量化为旋转角 `2·arccos(min(1, |v|/amax))`，`rhs_qram_preparation` 以 QRAM 角度表制备幅度、以符号数据库的 Z 反冲写入符号。

经典侧的 Riemann 计算与局部流场更新由 `RoeFlowData`（`applications/flow_data.py`）管理：单点更新只重算相邻界面与受影响的树节点；原生后端的 bank 修改仍需重新物化，不等同于论文假设的常数时间物理写。模块属性 `correctness="pending"`：当前见证的是输入适配与两条 QLSS 路线的可替换性，数值求解正确性未认证；完整假设与边界见 [QFVM/QLSS 输入模型审查](../../reference/qfvm-input-models.md)。

## 验证方案

类别 C6（组合骨架，判定准则见 `../../development/validation-plan.md` §2：组件契约满足 + 端到端小实例语义正确 + bind 不变）。`validation-coverage.md` 将 `qfvm.py`/`qham/` 登记于应用层行（C6、V1、无缺口）；QFVM 的具体见证位于 ODE/PDE 求解器行的证据列：

- 结构与数值：`tests/core/test_differential.py:DifferentialStructureTests.test_qfvm_uses_raw_flow_and_modular_arithmetic`——`FixedFormat(4, 1)`、`amax=4` 的输入下 `roe_qfvm_block_encoding` 的 α = 36；BE 程序恰含 4 个未解析槽，闭合后资源为 rho / momentum / energy / geometry 四个数据库；程序包含 sqrt / div / mul / select 算术模块；geometry 表大小恰为 `2**(width+4)`；稀疏编码标注稀疏度 9 与自伴扩展。
- 经典侧：同类 `test_local_riemann_updates_only_neighbor_faces_and_tree`——16 单元流场单点更新只重算界面 (4, 5) 与单元 (4, 5, 6)，rhs 角度树更新数少于全量，且不物化任何 matrix bank。
- 绑定：`test_qfvm_uses_raw_flow_and_modular_arithmetic` 内 `bind_qfvm` 替换全部抽象槽后无未解析声明。

## 已知缺口与计划阶段

与 `validation-coverage.md` 的应用层行一致：无登记缺口，阶段 V1。模块属性中的 `correctness="pending"` 是实现边界记录（数值求解正确性认证不在当前验证范围），不是覆盖矩阵上的缺口。

## 数值验证

2026-09-16 由 `tests/verification/verify_qham_qfvm.py` 执行的论文级数值验证（真实后端：PySparQ 原生 RIR 解释器、UniQC 全振幅态向量；无 mock、无 skip）。本机环境无 C++ 编译器，`arithmetic_native_registry` 不可用，Roe 算术经 PySparQ RIR 门级展开执行（单条目电路约 1.6×10⁷ 展开门），故条目级验证采用叠加采样而非全矩阵穷举。

### 实验设计

- **定点格式选择**：`roe_face` 的 compile_function 要求 `Index(2)` 的行列输入可表示，最小可用格式为 `FixedFormat(5,2)`；`entropy_delta=0.5` 使 2δ、δ²、δ 在该格式下全部精确可表示。注意：若小数位不足以表示 2δ（如 fmt=(4,1) 配默认 δ=0.125，2δ=0.25 截断为 raw 0），熵修正分支会发生除零，条目按文档化的 totalize 行为静默归零——集成测试此前只覆盖 padding 对角，未覆盖非零 Roe 条目。
- **Roe 面通量**：编译后的 `roe_face` 在 16+8+8=32 个叠加分支（全 row/col 对、密度扫描、动量扫描）下与脚本内独立定点仿真（按 `fixed_arithmetic` 文档语义 toward_zero/modular_wrap 重实现）做**位级**对拍；方法误差（对照 float64 `roe_formulas`）单独报告。
- **稀疏条目/位置**：条目 oracle 在列 5 的 8 分支叠加（三带非零元 + 零结构 + 补齐零）下逐位对照定点仿真矩阵；padding 对角独立案例；位置 oracle 全 32 列叠加验证完整置换与 9 槽位映射（对照独立几何语义）。
- **RHS 与经典恒等式**：RHS 残差态制备（QRAM 角度树 + 符号反冲）对照独立残差计算；flow_data 的 F*(L,R)=left·U_L+right·U_R 对照矩阵求逆实现 `riemann_flux`；M·u−mass·u=−residual（隐式 FVM 符号约定）；局部更新 patch 与全量重算一致；ptheta 角表逐点真值。统一四单元流场 ρ=(1,1,1.25,1.25)、m=(0,0.25,0,−0.25)、E=(1,1,1,1)。

### 关键指标

| 案例 | 规模 | 后端路径 | 指标值 |
|---|---|---|---|
| roe_face 位级对拍 | 32 分支，fmt=(5,2) | rir-pysparq | raw/status 逐位一致（32/32）；方法误差 0.427 |
| 条目矩阵采样 | 列 5，8 分支 | rir-pysparq | 实现误差 0.0（逐位一致）；方法误差 0.481；零结构正确 |
| 补齐对角 padding | 坐标 (7,7)/(23,7) | rir-pysparq | 对角 raw=4（1.0），扩张块外为 0 |
| 位置 oracle | 32 列 × 32 输入 | rir-pysparq | 每列完整置换；288 槽位映射 0 失配；几何表逐点 0 失配；work 复净 |
| RHS 制备 | 16 地址，angle_width 8 | rir-pysparq | 振幅误差 1.11e-16，符号 0 失配，总概率 1.0；角度树/符号 bank 逐点 0 失配 |
| 通量恒等式 | 4 组采样态 | 经典独立 | F* 块形式 vs 求逆实现 3.33e-16；F*(U,U)=F(U) 2.22e-16 |
| 单步更新 | 4 单元 | 经典独立 | M·u−mass·u=−residual 2.22e-16；残差 bank 0 失配；局部 patch 与全量一致 |

方法误差是 5 位定点流水线的固有量化误差（约 1–2 个量子），不是实现缺陷：实现侧（对照定点语义仿真）为逐位精确。

### 复现

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_qham_qfvm.py
```

产物：`out/verification/qham_qfvm.json`（15 个案例的全部指标与判据）。

## 相关链接

- 源码：`src/pyqecclang/applications/qfvm.py`（Roe 算术 `applications/roe.py`、公式 `applications/roe_formulas.py`、经典数据 `applications/flow_data.py`）
- 用户指南：[QFVM 输入模型与求解器替换](../qfvm.md)
- API 参考：[QFVM 应用](../../api/applications/qfvm.rst)
- 同组页面：[一般量子同伦分析](qham.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
