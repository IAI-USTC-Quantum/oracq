# QFVM 输入模型与求解器替换

QFVM 应用从流场数据构造线性系统，并把它交给可替换的 QLSS。当前实现针对周期一维 Euler 方程，使用三个守恒量和 frozen-Roe Jacobian。它不覆盖原论文的全部网格、边界和物理模型。

## 数据与量子访问

流场的密度、动量和能量保存在 QRAM 数据中。几何表记录邻居、槽位与分量索引，不预存完整 Jacobian 矩阵。

量子矩阵元 oracle 查询相关单元，调用可逆 Roe 算术，并根据行列索引选择元素。位置 oracle 给出每列的结构位置，其接口是 CKS 使用的原地索引置换。RHS 制备使用残差数据结构提供的角表。

经典 Riemann 计算和局部流场更新由 `RoeFlowData` 管理。逻辑 QRAM patch 可以局部更新；当前原生后端物化仍可能重建被修改的 bank，不能把两者等同。

## 构造问题

```python
from pyqecclang import SpectralPromise
from pyqecclang.applications.qfvm import roe_qfvm_inputs, roe_qfvm_problem
from pyqecclang.algorithms.qlss import CKSConfig, CostaConfig, make_cks_qlss, make_costa_qlss

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

CKS 入口消费稀疏问题，Costa 入口请求 BE。当前实对称稀疏适配通过 `T†ST` 构造相应 BE，并记录 alpha。没有从任意 BE 反向恢复高效稀疏访问的通用适配。

QFVM 的非对称物理矩阵使用明确的 Hermitian 扩张，并在输出时选择物理坐标。这个扩张服务于线性求解；不能直接将它当作原生成元的等价 QODE 演化。

两条 QLSS 路线可以保持相同的问题输入与物理输出解释，但辅助位宽度和内部模块不同。更换实现后应重新生成上层布局。只有在 ABI 和 alpha 不变时，才适合对已有开放槽晚绑定。

## 输出与范数

`SolveResult` 返回态 oracle、矩阵范数探针、输入 alpha、谱声明和适配记录。它将求解成功概率与独立矩阵探针的条件概率分开，避免误用过滤成功率作为解向量范数。

解态方向、成功通道、数值误差以及 CFD 外循环仍需按具体案例验证。更完整的数学假设和原始论文对照见[QFVM/QLSS 输入模型审查](../reference/qfvm-input-models.md)。

## 数值验证

2026-09-16 的论文级数值验证（`tests/verification/verify_qham_qfvm.py`，真实后端无替身）从数值上确认了本章的数据与访问构造：

- **数值通量**：编译后的 Roe 面电路与独立定点仿真逐位一致（32/32 分支），对照 float64 Roe 公式的方法误差 0.427（5 位定点固有量化）；稀疏条目 oracle 在采样列上逐位一致，补齐对角精确给出 padding_value，结构域外为零。
- **单步更新**：F*(L,R)=left·U_L+right·U_R 与矩阵求逆实现 `riemann_flux` 相差 3.33e-16；矩阵恒等式 M·u−mass·u=−residual（隐式 FVM 符号约定）残差 2.22e-16；`RoeFlowData` 局部更新只重算相邻界面与受影响树节点，且与全量重算逐 bank 一致。
- **位置访问**：全 32 列叠加下每列都是完整置换，9 个结构槽位映射与独立几何语义 0 失配，几何 QRAM 表逐点真值。
- **RHS 制备**：残差态振幅与独立残差/范数计算一致（误差 1.11e-16），符号经 Z 反冲精确写入，角度树与符号 bank 逐点真值。

验证用 `FixedFormat(5,2)` 与 `entropy_delta=0.5`（全部常数精确可表示）；更低精度格式下 2δ 可能截断为零导致熵修正分支除零、条目按文档行为归零（详见[算法页数值验证](algorithms/qfvm.md#数值验证)）。QLSS 两条路线的数值求解精度属求解器一侧，不在本页验证范围。产物：`out/verification/qham_qfvm.json`。
