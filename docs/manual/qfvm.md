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
