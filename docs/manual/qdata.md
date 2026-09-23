# 量子数据结构与 KP 推荐系统

{obj}`QVector <pyqecclang.algorithms.input_model.qdata.QVector>` 与 {obj}`QMatrix <pyqecclang.algorithms.input_model.qdata.QMatrix>` 封装 QFVM（arXiv:2102.03557 式 15–22 残差平方和树）与 Kerenidis–Prakash 推荐系统（arXiv:1603.08675 Thm 5.1 + 附录 A.1）共用的量子数据结构，数据面全部经由 [QRAM 指针访问](qmem.md)。API 见 [量子数据结构（qsample 与 sample-and-query）](../api/algorithms/input_model/qdata.rst)。

## QVector：qsample 向量

平方范数二叉树：叶存分量平方，内部节点缓存 RY 旋转角字（θ = 2·acos(√(S_left/S_node))）。量子侧按层制备归一化态，每层两次 QRAM 查询；符号经 1 位 bank 的相位反冲写入。经典侧单点更新只重算叶到根的路径。

```{doctest}
>>> from pyqecclang import QVector, simulate
>>> from pyqecclang.algorithms.common.arithmetic import FixedFormat
>>> vector = QVector((3.0, 4.0), fmt=FixedFormat(8, 4), angle_width=12)
>>> state = simulate(vector.preparation().operation.program(), vector.snapshot())
>>> [round(a.real, 3) for _, a in sorted(state.amplitudes.items())]
[0.6, 0.8]
```

## QMatrix：sample-and-query 矩阵

条目 bank 支持任意叠加查询（地址按 (行, 列) 二维指针寻址）；每行一棵 QVector 树给出 Ũ: |i⟩|0⟩→|i⟩|Ā_i⟩；行范数根树给出 Ṽ: |0⟩|j⟩→|Ã⟩|j⟩。条目限非负定点值。

```{doctest}
>>> from pyqecclang import QMatrix, simulate
>>> from pyqecclang.algorithms.common.arithmetic import FixedFormat
>>> matrix = QMatrix([[0.5, 0.25], [0.5, 0.25]], fmt=FixedFormat(8, 4), angle_width=12)
>>> state = simulate(matrix.query().operation.program(),
...                  {"entries": matrix.snapshot()["entries"]}, initial={"address": 0b10})
>>> sorted(state.amplitudes.items())
[((2, 8), (1+0j))]
```

## KP 量子推荐系统

`kp_recommendation` 按论文 Lemma 5.3 构造 W = Ũ R₁ Ũ⁻¹ · Ṽ R₀ Ṽ⁻¹（R₀/R₁ 为零基矢反射），对其做相位估计：cos(θᵢ/2) = σᵢ/‖A‖_F。逐相位字估计 σ̂ 并按阈值翻 flag（Alg 2 的确定性投影），逆相位估计后测 (flag, item)；flag=1 分支的条目分布即推荐采样分布。镜像相位 2^p−t 与 t 对应同一奇异值（W 特征值成对 e^{±iθ}）。

```{doctest}
>>> from pyqecclang import KPRecommendationConfig, QMatrix, simulate, kp_recommendation
>>> from pyqecclang.algorithms.common.arithmetic import FixedFormat
>>> matrix = QMatrix([[0.5, 0.25], [0.5, 0.25]], fmt=FixedFormat(8, 4), angle_width=12)
>>> result = kp_recommendation(matrix, 0, KPRecommendationConfig(precision=4, sigma=0.5))
>>> success, distribution = result.readout(simulate(result.operation.program(), result.memories()))
>>> round(success, 3), {k: round(v, 3) for k, v in sorted(distribution.items())}
(1.0, {0: 0.8, 1: 0.2})
```

秩一矩阵下推荐分布精确落在主奇异向量 v₁²=(0.8, 0.2) 上（残差为角量化误差）。带小扰动的矩阵上，阈值 σ 滤掉小奇异值分量，成功概率低于 1，条件分布趋近 v₁²——数值验证见 `tests/core/test_recommendation.py`。已声明的假设：非负定点条目、角字量化、无幅度放大（KP §6 的变时放大记为后续工作）。

## QFVM 的 QMem 直连路径

`applications/qfvm_qmem.py` 把 QFVM 的数据访问重写到同一套语法糖上：三守恒量合并为 (场, 单元) 状态表 `QMem(b, "state", shape=(3, n))`，邻居 cell±1 用模 2^cell_width 的指针算术实现周期边界；几何表按 (槽位, 列) 二维寻址；残差态由 QVector 树制备。模块直接声明 QRAM 形式资源，不再经过抽象槽位与 bind。电路语义（可逆 Roe 算术、九槽位几何、补齐对角、原地位置置换）与既有 `applications/qfvm.py` 路径一致，等价性由 `tests/core/test_qfvm_qmem.py` 在 simulate 上逐振幅对拍。
