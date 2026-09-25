# QFVM 输入模型与求解器替换

[English](../../manual/qfvm.html) · **简体中文**

QFVM 应用从流场数据构造线性系统，并把它交给可替换的 QLSS。当前实现针对周期一维 Euler 方程，使用三个守恒量和 frozen-Roe Jacobian。它不覆盖原论文的全部网格、边界和物理模型。端到端工作流见[教程：科学计算工作流](../tutorials/scientific-workflows.md)。

## 数据与量子访问

流场的密度、动量和能量保存在 QRAM 数据中。几何表记录邻居、槽位与分量索引，不预存完整 Jacobian 矩阵。

量子矩阵元 oracle 查询相关单元，调用可逆 Roe 算术，并根据行列索引选择元素。位置 oracle 给出每列的结构位置，其接口是 [CKS](algorithms/cks.md) 使用的原地索引置换。RHS 制备使用残差数据结构提供的角表。

经典 Riemann 计算和局部流场更新由 {obj}`RoeFlowData <oracq.applications.flow_data.RoeFlowData>` 管理。逻辑 QRAM patch 可以局部更新；当前原生后端物化仍可能重建被修改的 bank，不能把两者等同。

## 构造问题

```python
from oracq import SpectralPromise
from oracq.applications.qfvm import roe_qfvm_inputs, roe_qfvm_problem
from oracq.algorithms.qlss.qlss import CKSConfig, CostaConfig, make_cks_qlss, make_costa_qlss

inputs = roe_qfvm_inputs()
problem = roe_qfvm_problem(
    inputs,
    spectrum=SpectralPromise(norm_upper=10.0, sigma_min_lower=0.5),
    amax=4.0,
)
cks = make_cks_qlss(CKSConfig(order=2))
costa = make_costa_qlss(CostaConfig(steps=1))
cks.check(problem).require()
costa.check(problem).require()
```

`amax`、谱界和矩阵性质是调用者声明，应适用于实际量化后的矩阵。示例中的数值用于说明接口，并不代替某个真实流场的谱分析。

## 替换 QLSS

[CKS](algorithms/cks.md) 入口消费稀疏问题，[Costa](algorithms/costa-walk.md) 入口请求 BE。当前实对称稀疏适配通过 `T†ST` 构造相应 BE，并记录 alpha。没有从任意 BE 反向恢复高效稀疏访问的通用适配。

QFVM 的非对称物理矩阵使用明确的 Hermitian 扩张，并在输出时选择物理坐标。这个扩张服务于线性求解；不能直接将它当作原生成元的等价 QODE 演化。

两条 QLSS 路线可以保持相同的问题输入与物理输出解释，但辅助位宽度和内部模块不同。更换实现后应重新生成上层布局。只有在 ABI 和 alpha 不变时，才适合对已有开放槽晚绑定。

## 输出与范数

{obj}`SolveResult <oracq.algorithms.qlss.qlss.SolveResult>` 返回态 oracle、矩阵范数探针、输入 alpha、谱声明和适配记录。它将求解成功概率与独立矩阵探针的条件概率分开，避免误用过滤成功率作为解向量范数。

解态方向、成功通道、数值误差以及 CFD 外循环仍需按具体案例验证。更完整的数学假设和原始论文对照见[QFVM/QLSS 输入模型审查](../reference/qfvm-input-models.md)。

## 实现逐行讲解（含 QRAM 数据结构）

本节把 QFVM 的实现完整拆开讲：先列全部 QRAM bank 的契约（谁写、谁读、存什么），再逐字段讲几何表与残差树，然后逐段讲三段量子电路（矩阵元 oracle、稀疏位置 oracle、RHS 制备），最后给端到端装配的逐行清单。源码对应 `src/oracq/applications/qfvm.py` 与 `src/oracq/applications/flow_data.py`。

### QR1. QRAM bank 总览

{obj}`roe_qfvm_inputs <oracq.applications.qfvm.roe_qfvm_inputs>` 声明的是**抽象槽位**（开放模块），数据到绑定与执行期才出现。绑定后每个 bank 的完整契约：

| bank | 地址宽 → 数据宽 | 经典侧谁写 | 量子侧谁查 | 内容 |
|---|---|---|---|---|
| `rho` / `momentum` / `energy` | cell_width → fmt.width | `RoeFlowData.update` | {obj}`roe_entry <oracq.applications.qfvm.roe_entry>`：每元素查 3 邻居 × 3 场 = 9 次 | 每单元守恒量的定点编码 |
| `geometry` | width+4 → geometry_width | {obj}`geometry_cells() <oracq.applications.qfvm.geometry_cells>` 静态生成（与流场值无关） | 位置 oracle 9 次；条目 oracle 9 次 | 七段打包的纯几何字（QR2） |
| `rhs_values` | cell_width+2 → fmt.width | `RoeFlowData.update` | —（经典范数树的数据面） | 每 (单元, 分量) 的量化残差值 |
| `rhs_sign` | cell_width+2 → 1 | `RoeFlowData.update` | RHS 制备查 2 次（写入 + 反查询） | 残差符号位 |
| `rhs_angles` | cell_width+2 → angle_width | `RoeFlowData.update`（树缓存） | {obj}`qram_state_prep <oracq.algorithms.input_model.oracles.qram_state_prep>` 每层 2 次、共 2·(cell_width+2) 次 | 平方范数树内部节点的旋转角字 |
| `theta` | fmt.width → angle_width | {obj}`ptheta_cells() <oracq.applications.qfvm.ptheta_cells>` 静态生成 | —（预留数据面，见下） | θ=2·acos(min(1,|v|/amax)) 的角字 |

三个要点：

1. **矩阵元从不预存**。QRAM 里只有原始场量与几何；矩阵元素由可逆 Roe 算术在查询时现场算出（QR4）。
2. **`theta` 是预留的**。它是"残差值定点字 → 旋转角字"的换算表，随 {obj}`qfvm_memories <oracq.applications.qfvm.qfvm_memories>` 一并物化；当前 CKS/Costa 两条路线没有电路查询它，文档如实记录这一点，不把它说成已接入的 oracle。
3. **内存快照必须固定**。同一段相干计算中的 P、P†、受控调用与范数探针必须读同一版本的数据；这写在 {obj}`roe_qfvm_problem <oracq.applications.qfvm.roe_qfvm_problem>` 的 `data_assumptions` 属性里，语言不替硬件证明它。

### QR2. 几何表 `geometry_cells` 逐字段

先定坐标布局：完整矩阵坐标是 `half·2^width + 4·cell + var`——最低 2 位是分量 `var`（0..3，3 为补齐分量），中间 cell_width 位是单元号，最高位 `half` 选择 Hermitian 扩张 `D=[[0,M],[Mᵀ,0]]` 的上/下块。QRAM 地址 = `row + slot·2^width`，其中 row 打包 (half, cell, var)，slot 打包 (band, other_var)。

逐行读 `geometry_cells`：

```python
cell, var, half = (row >> 2) % n, row % 4, row >> (cw + 2)
```

拆出行坐标：低 2 位 `var`、中段 `cell`、最高位 `half`（上块 = M，下块 = Mᵀ）。

```python
band, other_var = divmod(slot, 3)
valid = int(slot < 9 and var < 3)
```

16 个 slot 值中只有前 9 个是结构槽位（band∈{西,中,东} × other_var∈{0,1,2}）；补齐分量（var=3）与无效槽位标 `valid=0`。

```python
neighbor_cell = (cell + band - 1) % n
neighbor = other_var + 4 * neighbor_cell + ((1 - half) << (cw + 2))
```

对手坐标：本行经槽位 (band, other_var) 耦合到**另一半空间**的 (neighbor_cell, other_var)。上块行指下块列、下块行指上块列——这正是 D 的转置耦合结构；`% n` 实现周期边界。

```python
reverse = (2 - band) * 3 + var if valid else 0
source = cell if half == 0 else neighbor_cell
rowvar, colvar, source_band = (
    (var, other_var, band) if half == 0 else (other_var, var, 2 - band)
)
```

`reverse` 是同一矩阵元的对称槽位；`source` 是计算该元素时该用的源单元（下块行取邻居）；最后一行把下块条目折算成上块的同一物理元素（行列互换、band 取 2−band）。

```python
values = ((neighbor, w), (reverse, 4), (source, cw), (rowvar, 2),
          (colvar, 2), (source_band, 2), (valid, 1))
```

七段按低位到高位打包成 `geometry_width` 位的一个字：对手坐标 | 对称槽位 | 源单元 | 行分量 | 列分量 | 源带 | 有效位。位置 oracle 与条目 oracle 查到的都是这个字，再按 `_geometry_refs` 的同一段宽切开。

### QR3. RHS 的 QRAM 符号残差树

**经典侧**（`RoeFlowData.update`）维护 1-based 堆数组树（`8n` 项），叶存残差分量的平方：

```python
node = 4 * self.n + address          # 叶：堆 1..8n，叶子从 4n+1 起（address=4*cell+j）
self.tree[node] = value * value      # 残差平方
while node > 1:                      # 沿祖先向上收集受影响节点
    node //= 2
    ancestors.add(node)
```

```python
self.tree[node] = self.tree[2*node] + self.tree[2*node + 1]      # 自底向上重算子树和
angle = 2 * math.acos(math.sqrt(self.tree[2*node] / self.tree[node]))
changes["rhs_angles"][node - 1] = round(angle * (1 << aw) / (2 * math.pi)) % (1 << aw)
```

角度公式来自 RY 半角约定：`RY(θ)|0> = cos(θ/2)|0> + sin(θ/2)|1>`，要让左右子树幅度比为 `sqrt(S_left/S_node)`，就取 `θ = 2·acos(sqrt(S_left/S_node))`。角度量化成 `angle_width` 位角字写入 `rhs_angles` bank（地址 = node−1）。单点流场更新只触碰相邻界面与 O(log) 个树节点——这就是"局部更新"声明的来源。

**量子侧幅度**（`qram_state_prep`，按层制备）：

```python
for depth in range(width):
    bit = width - depth - 1                     # 本层决定的目标位（从最高位往下）
    if depth:
        b.xor(b["target"][bit + 1:], addr[:depth])   # 已确定前缀拷入地址
    offset = (1 << depth) - 1
    b.add_const(addr.reinterpret("uint"), offset)    # 地址 = 前缀 + 树节点偏移
    b.qram("angles", addr, angle)                    # 查该内部节点的角字
    for k in range(angle_width):
        with b.control(angle[k]):
            b.ry(b["target"][bit], 2 * math.pi * (1 << k) / (1 << angle_width))
    b.qram("angles", addr, angle)                    # 反查询复净 angle
    b.add_const(addr.reinterpret("uint"), (-offset) % (1 << address_width))
```

逐行要点：第 depth 层的树节点编址是 `(1<<depth)−1+前缀`（与 {obj}`qram_state_angles <oracq.algorithms.input_model.oracles.qram_state_angles>` 的表编址一致）；一个 `angle_width` 位角字不必编译成 `2^angle_width` 种旋转门——按位加权 `2π·2^k/2^aw` 用 `angle_width` 次受控 RY 合成任意角；每层一查一反查，总 QRAM 查询 `2·width` 次（对照逐叶制备的 O(2^n) 次查询）。

**量子侧符号**（{obj}`rhs_qram_preparation <oracq.applications.qfvm.rhs_qram_preparation>`）：幅度树只产生非负振幅，符号用独立的 1 位库经**相位反冲**写入：

```python
invoke(b, prep.operation, "prep", target=b["target"][:n], work=b["work"][: n + aw])  # 幅度
flag = b["work"][n + aw :]
invoke(b, sign.operation, "sign", address=b["target"][:n], data=flag)  # 符号 XOR 进 flag
b.z(flag)                                                              # flag=1 的分支乘 −1
invoke(b, sign.operation, "sign", address=b["target"][:n], data=flag)  # 反查询复净 flag
```

{obj}`Load <oracq.infrastructure.ir.Load>` 自逆，所以"写入符号 → Z → 还原符号"的净效果只剩分支相位：负残差分量得到 π 相位，flag 与 work 全部复净（`zero_input=True`）。

### QR4. 矩阵元 oracle `roe_entry` 逐段

寄存器：`source`（单元）、`row`/`col`（分量）、`band`（西/中/东）、`value`（定点输出）、`status`（算术状态）。资源是三个守恒量库。

```python
for offset in (-1, 0, 1):
    addr = g.local(cw)
    b.xor(b["source"], addr)                                   # addr = source
    b.add_const(addr.reinterpret("uint"), offset % (1 << cw))  # addr += 偏移（周期回绕）
    addresses.append(addr)
    for name, op in fields:
        word = g.local()
        invoke(b, op, name, address=addr, data=word)           # 查 rho/momentum/energy
        state.append(word)
```

对 source 及其左右邻居各查三个守恒量库——元素计算的输入是 3 个单元 × 3 个场，共 9 次 QRAM 查询。

```python
face = roe_face(fmt=fmt, gamma=gamma, entropy_delta=entropy_delta)
for i in range(2):
    b.call(face, rho_l=..., m_l=..., e_l=..., rho_r=..., m_r=..., e_r=...,
           row=b["row"], col=b["col"], left=left, right=right, status=flag)
```

{obj}`roe_face <oracq.applications.roe.roe_face>` 是 {obj}`frozen_roe_face <oracq.applications.roe_formulas.frozen_roe_face>` 纯函数经 {obj}`compile_function <oracq.infrastructure.mathfunc.compile_function>` 编译出的可逆定点电路（含 sqrt/div/mul/select）；调用两次得到左、右界面的通量左右特征分量。`status≠0`（溢出/除零）时元素按文档约定 totalize 为零。

```python
west = -faces[0][0] / dx
center = (faces[1][0] - faces[0][1]) / dx + g.choose(equal, mass, 0)
east = faces[1][1] / dx
value = g.choose3(b["band"], [west, center, east])
```

按有限体积通量差组装西/中/东三个候选；质量项只在 `row == col`（`component_equal` 布尔网络）时加到中心块对角；最后按 `band` 三选一。

### QR5. 稀疏位置 oracle（CKS 原地置换）逐段

CKS 要求"给定列与稀疏序号，原地得到行索引"。{obj}`qfvm_sparse_access <oracq.applications.qfvm.qfvm_sparse_access>` 的 locator 用 9 次相干值转置构造完整置换：

```python
for rank in range(9):
    ...  # slot 置为 rank；查 geometry(column, slot) 得七段字；neighbor = data[:width]
    padded = locator.local(...); locator.xor(locator["column"], padded)
    locator.add_const(padded.reinterpret("uint"), rank)
    with locator.control(locator["column"][:2], 3):      # 补齐分量（var=3）
        locator.xor(data[:n], neighbor)
        locator.xor(padded, neighbor)                    # 对手坐标改成 column+rank（padding 对角）
    for bit in range(n):
        if (rank >> bit) & 1:
            locator.x(left[bit])                         # left 初值 = rank
    for previous in range(rank):
        locator.call(value_transposition(n), index=left, a=lefts[previous], b=neighbors[previous])
```

补齐分量没有真实的邻居结构：对手坐标被改写为 `column+rank`，使 9 个序号恰好映射到 9 个不同的补齐坐标（padding 对角）。`left` 从 rank 出发，对之前每个 `(left_j, neighbor_j)` 做值转置——若 rank 已被占用就被换到尚未使用的槽位，保证 9 个槽位映射构成**完整置换**。

```python
setup = tuple(locator._frames[0])
for left, right in zip(lefts, neighbors, strict=True):
    locator.call(value_transposition(n), index=locator["index"], a=left, b=right)
locator.emit(Adjoint(setup))
```

对 `index` 依次执行 9 次值转置（`index` 中等于 `left_j` 的字被换成 `neighbor_j`），再用 {obj}`Adjoint(setup) <oracq.infrastructure.ir.Adjoint>` 复净全部工作位。逆映射 = 同一电路的 adjoint，这正是 CKS 接口需要的原地语义；代价是稀疏度次可逆比较/换位，不是免费单位门。

**条目 oracle**同用几何表：对 9 个槽位各做"查询 → {obj}`compare_words <oracq.algorithms.input_model.sparse.compare_words>` 匹配 row → XOR 累加出 selected 七段字 → 反比较 → 反查询"的 compute/uncompute 对；拆段后调用 QR4 的 `roe_entry` 得到 `(value, status)`：

```python
with b.control(valid):
    with b.control(status, 0):
        b.xor(value, b["data"])                    # 有效且算术成功的元素 XOR 写入
with b.control(fuse(same, b["column"][:2]), 7):    # 对角且列分量为补齐（3）
    ...                                            # 写 padding_value
```

前后各包一层 `Adjoint(forward)` 保持可逆。最终 {obj}`SparseAccess(location, entry, width, fmt.width, 9) <oracq.algorithms.input_model.oracles.SparseAccess>`。

### QR6. 端到端装配逐行

```python
from oracq import FixedFormat, SpectralPromise
from oracq.applications.flow_data import RoeFlowData
from oracq.applications.qfvm import (
    bind_qfvm, qfvm_memories, roe_qfvm_inputs, roe_qfvm_problem,
)
from oracq.algorithms.qlss.qlss import CKSConfig, CostaConfig, make_cks_qlss, make_costa_qlss

# 定点格式：数值验证表明 (5,2) 是 roe_face 常数全部可精确表示的最小格式。
fmt = FixedFormat(5, 2)
# 1. 声明输入模型：六个抽象数据库 + 一个抽象 RHS 制备（全是开放槽，尚无数据）。
inputs = roe_qfvm_inputs(fmt=fmt, angle_width=6)
# 2. 经典侧流场：RoeFlowData 写 rho/momentum/energy/rhs_* 六个 bank 并维护范数树。
flow = RoeFlowData([(1, 0, 1), (1, 0.25, 1), (1.25, 0, 1), (1.25, -0.25, 1)],
                   fmt=fmt, angle_width=6)
# 3. 组装问题：Hermitian 扩张 + 补齐对角合并进谱声明；谱界是调用者的义务。
problem = roe_qfvm_problem(
    inputs, spectrum=SpectralPromise(norm_upper=10.0, sigma_min_lower=0.5), amax=4.0,
)
# 4. 同一个问题交两条路线：CKS 消费 SparseSystem；Costa 触发稀疏→BE 适配。
cks = make_cks_qlss(CKSConfig(order=2))(problem)
costa = make_costa_qlss(CostaConfig(steps=1))(problem)
# 5. 闭合开放槽：五个数据库绑定 QRAM 实际库，RHS 槽绑定符号残差树制备。
cks_rir = bind_qfvm(cks.operation.program(), inputs)
# 6. 物化运行时内存表：流场快照六个 bank + geometry 表 + theta 表。
memories = qfvm_memories(inputs, flow, amax=4.0)
# 7. 执行时把 memories 按名交给后端（run_pysparq(cks_rir, memory=memories)）；
#    替换 Costa 只需换第 4 步的 protocol，其余各行不变。
```

## 数值验证

2026-09-16 的论文级数值验证（`tests/verification/verify_qham_qfvm.py`，真实后端无替身）从数值上确认了本章的数据与访问构造：

- **数值通量**：编译后的 Roe 面电路与独立定点仿真逐位一致（32/32 分支），对照 float64 Roe 公式的方法误差 0.427（5 位定点固有量化）；稀疏条目 oracle 在采样列上逐位一致，补齐对角精确给出 padding_value，结构域外为零。
- **单步更新**：F*(L,R)=left·U_L+right·U_R 与矩阵求逆实现 {obj}`riemann_flux <oracq.applications.flow_data.riemann_flux>` 相差 3.33e-16；矩阵恒等式 M·u−mass·u=−residual（隐式 FVM 符号约定）残差 2.22e-16；`RoeFlowData` 局部更新只重算相邻界面与受影响树节点，且与全量重算逐 bank 一致。
- **位置访问**：全 32 列叠加下每列都是完整置换，9 个结构槽位映射与独立几何语义 0 失配，几何 QRAM 表逐点真值。
- **RHS 制备**：残差态振幅与独立残差/范数计算一致（误差 1.11e-16），符号经 Z 反冲精确写入，角度树与符号 bank 逐点真值。

验证用 {obj}`FixedFormat(5,2) <oracq.algorithms.common.arithmetic.FixedFormat>` 与 `entropy_delta=0.5`（全部常数精确可表示）；更低精度格式下 2δ 可能截断为零导致熵修正分支除零、条目按文档行为归零（详见[算法页数值验证](algorithms/qfvm.md#数值验证)）。QLSS 两条路线的数值求解精度属求解器一侧，不在本页验证范围。产物：`out/verification/qham_qfvm.json`。

## QMem 直连并行路径

`applications/qfvm_qmem.py` 提供同一套 QFVM 数据访问的并行实现，数据面全部改经[指针式 QRAM 访问](qdata.md)：

- 三个守恒量库合并为一张 `(场, 单元)` 状态表 {obj}`QMem(b, "state", shape=(3, n)) <oracq.infrastructure.qmem.QMem>`；周期邻居在 cell_width 位 scratch 上做模加（与旧路径相同），再经二维指针 `[field, addr]` 查询。
- 几何表按 `(槽位, 列)` 二维寻址（旧路径的 {obj}`fuse(column, slot) <oracq.infrastructure.ir.fuse>` 编址逐点一致）。
- 残差态由 {obj}`QVector <oracq.algorithms.input_model.qdata.QVector>` 平方范数树制备（取代 `qram_state_prep` + 符号库组合）。
- 模块直接声明 QRAM 形式资源，不再经过抽象数据库槽位与 {obj}`bind_qfvm <oracq.applications.qfvm.bind_qfvm>`。

可逆 Roe 算术、九槽位几何选择、补齐对角与原地位置置换的电路结构与本路径完全相同。等价性证据（`tests/core/test_qfvm_qmem.py`）：几何/状态/周期邻居的叠加探针逐字一致、位置 oracle 全列置换逐振幅一致、残差态振幅一致、含编译 Roe 算术的物理层在单点上逐位一致。本页上文描述的槽位-绑定路径仍是 QLSS 集成的主路径；QMem 路径当前定位是数据访问层的等价重写与后续演进的基线。
