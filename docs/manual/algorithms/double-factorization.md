# 双因子分解块编码（Double Factorization）

> 类别 C2 · 模块 `pyqecclang.algorithms.lowrank` · 阶段 V1

## 概述

量子化学电子哈密顿量的双因子分解（DF）把二体积分张量压缩为低秩形式

$$
H = \mathrm{scalar}\cdot I + \sum_r U_r\, \mathrm{diag}(g_r)\, U_r^\dagger ,
$$

每个秩项是对角谱 $g_r$ 经旋转 $U_r$ 共轭的厄米矩阵。实现依据 Berry et al. 2019 与 von Burg et al. 2021 的"低秩张量 → LCU → BE"管道（文献引用见模块 docstring）：积分张量作为经典输入数据给出，本模块不做真实的化学积分计算，产出的块编码可直接进入 qubitization。

## 接口与输入模型

```python
DoubleFactorization(scalar, rotations, spectra)
DoubleFactorization.from_symmetric(scalar, rotations, factors)
diagonalize_symmetric(matrix, *, tolerance=1e-12, max_sweeps=100)
double_factorized_encoding(df, *, name=None)
```

- `DoubleFactorization`：DF 输入模型（input model 为 CP——张量数据作为经典参数直接给出，即 algorithm-coverage 中"低秩张量（HAM）→ BE"管道）。`rotations` 为显式小酉矩阵（维度 2..32 内二的幂，按列正交校验酉性，容差 1e-9），`spectra` 为对应实谱；数据类的属性 `width` 给出目标量子位数，`rank` 给出秩项数。
- `from_symmetric`：物理 DF 的经典预处理入口——接收显式酉 $U_r$ 与实对称 $G_r$，把 $G_r = V_r\,\mathrm{diag}(g_r)\,V_r^{\mathsf T}$ 的特征向量矩阵折叠进旋转（$U_r V_r$）。
- `diagonalize_symmetric`：实对称矩阵的 Jacobi 特征分解，返回 `(特征值, 特征向量矩阵)`；非对称输入在生成期抛 `ValidationError`。
- `double_factorized_encoding`：组装 LCU 块编码（`name` 形参保留但当前不参与命名）。

返回 `BlockEncoding`，模块属性：

| 属性 | 含义 |
|---|---|
| `be_alpha` | $\lvert\mathrm{scalar}\rvert + \sum_r \lVert g_r\rVert_1$（各谱 1-范数之和） |
| `be_form` | `"double_factorization"` |
| `df_rank` / `lcu_terms` | 秩项数 / 外层 LCU 的非零项数 |
| `df_lambda` | 与 `be_alpha` 同值的 λ 口径报告 |

## 实现要点

外层 LCU 的 PREPARE 在秩指标 $r$ 上，权重 $\propto \lVert g_r\rVert_1$；`scalar` 非零时以 `(scalar, identity)` 并入同一 LCU，全部系数为零（含全零谱）时抛 `ValidationError`。每个秩项的项块编码是 $U_r^\dagger\cdot\mathrm{diag}\cdot U_r$：对角部分在单比特 signal 上对每个基态 $\lvert t\rangle$ 施加受控 $R_y(\theta_t)$，$\cos(\theta_t/2) = g_t/\alpha_r$、$\alpha_r = \sum_p\lvert g_p\rvert$ 取谱 1-范数，(0,0) 块恰为 $\mathrm{diag}(g_r)/\alpha_r$；两侧共轭施加合成的 $U_r$。

$U_r$ 由两能级分解合成：逐列消元为对角相位后按逆序回放，每个两能级酉经 ZYZ 分解（$e^{i\varphi}R_z(\alpha)R_y(\beta)R_z(\gamma)$）加沿 Gray 路径的多控 X 转置实现。寄存器布局：项块编码为 `target(n) | signal(1)`，外层 LCU 的 signal 为 selector（项数所需位数）拼接 `work(1)`。

适用边界：显式矩阵路径仅支持 2..32 维（不超过 5 个量子位）的验证小实例；大规模时 $U_r$/$G_r$ 的谱与旋转角度表应改走 QRAM 数据绑定（见 `prepare_select` 的三层范式）。归一化口径上，单比特正定 $g$ 的 DF α 等于 $\mathrm{tr}(g)$，不引入 Pauli 展开的非对角冗余（紧性见证见验证节）。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）：块编码的零信号块须与目标哈密顿量在容差内一致。三层证据：

- 结构：`tests/core/test_lowrank.py:DiagonalizeSymmetricTests.test_reconstructs_factor`（$V\,\mathrm{diag}(\lambda)\,V^{\mathsf T}$ 重构与特征向量正交性，places = 10）与 `test_rejects_non_symmetric`；`DoubleFactorizationTests` 的构造属性断言（`be_form`、`df.rank == 1` / `df.width == 1`）；`ThcTests.test_invalid_inputs_fail_at_generation` 覆盖 DF 负例（空秩项、谱长度不符、非酉旋转、`scalar` 非数值、`from_symmetric` 数量不一致）。
- 数值：`DoubleFactorizationTests.test_single_rank_block_equals_hamiltonian`——单秩 DF 经 `assert_block_equals` 与 $G$ 逐列对拍（places = 9）；`test_scalar_and_rotated_ranks`——scalar 与 Hadamard 旋转的第二秩，期望矩阵在测试内独立组装后对拍，α 与 $\lvert\mathrm{scalar}\rvert + \sum_r\lVert g_r\rVert_1$ 对拍（places = 10）；`test_feeds_qubitization_walk`——产物进入 `transforms.qubitization_walk` 后无未绑定槽，寄存器为 `["target", "signal"]`。
- 绑定：`DoubleFactorizationTests.test_matches_pauli_encoding_block`——同一哈密顿量的 DF 编码与 Pauli LCU 编码经 `block_column` helper 逐列跨 BE 对拍（places = 9），覆盖两条独立编码路径的一致性。

α 紧性见证（V1 新增）：`test_alpha_matches_closed_form_eigenvalues` 用 2×2 闭式特征值 $\lambda_\pm = (t \pm \sqrt{t^2-4d})/2$（$t$ 为迹、$d$ 为行列式）独立计算 $\sum\lvert\lambda\rvert$ 与 `be.alpha` 对拍（places = 10，实测 2.0）；`test_df_alpha_tighter_than_pauli` 取非对角主导的正定 $g$，实测 DF α = 1.4 ≤ Pauli LCU α = 1.5（理论条件：单比特 PSD 矩阵上 Pauli α − DF α = $\lvert g_{01}\rvert + \lvert\Delta/2\rvert - \mathrm{tr}/2$，非对角主导时为正）。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（块对拍 + 闭式特征值与 Pauli 对照的 α 紧性 + 跨编码绑定一致性 + qubitization 衔接）。

## 相关链接

- 源码：`src/pyqecclang/algorithms/lowrank.py`
- 同模块算法：[THC 块编码](thc.md)
- API 参考：[化学低秩分解块编码](../../api/algorithms/lowrank.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
