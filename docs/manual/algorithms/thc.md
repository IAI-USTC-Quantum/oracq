# 张量超收缩块编码（Tensor Hypercontraction）

> 类别 C2 · 模块 `pyqecclang.algorithms.input_model.lowrank` · 阶段 V1

## 概述

张量超收缩（THC）把电子哈密顿量的二体积分张量压成叶算符形式

$$
H = \sum_{\mu\nu} \zeta_{\mu\nu}\, L_\mu L_\nu^\dagger ,
$$

其中 $\zeta$ 为实对称系数矩阵，$L_\mu$ 为显式小矩阵叶算符。实现依据 Lee et al. 2021 的 THC 表示接入本仓库的"低秩张量 → LCU → BE"管道（文献引用见模块 docstring）：积分张量作为经典输入数据给出，本模块不做真实的化学积分计算，产出的块编码可直接进入 qubitization。

## 接口与输入模型

```python
THCDecomposition(coefficients, leaves)
thc_encoding(thc)
```

- `THCDecomposition`：THC 输入模型（input model 为 CP——ζ 与叶矩阵作为经典参数直接给出，即 algorithm-coverage 中"低秩张量（HAM）→ BE"管道）。`coefficients` 必须是实对称方阵（对称性容差 1e-12），`leaves` 为维度一致（2..32 内二的幂）的显式小矩阵，无需酉或 Hermitian。数据类的属性 `width` 给出目标量子位数，`leaf_count` 给出叶算符个数。
- `thc_encoding`：组装 THC 哈密顿量的 LCU 块编码。

返回 `BlockEncoding`，模块属性：

| 属性 | 含义 |
|---|---|
| `be_alpha` | $\sum_{\mu\nu}\lvert\zeta_{\mu\nu}\rvert\,\alpha_\mu\alpha_\nu$（$\alpha_\mu$ 为叶算符的 Pauli l1 上界） |
| `be_form` | `"thc"` |
| `thc_leaves` / `lcu_terms` | 叶算符个数 / 外层 LCU 的非零项数 |
| `thc_lambda` | 与 `be_alpha` 同值的 λ 口径报告 |

## 实现要点

每个叶算符先经 `matrix_pauli_encoding` 展开为 Pauli LCU 块编码（Pauli 系数低于 1e-12 的字丢弃；显式展开仅支持不超过 5 位的小实例），$L_\nu^\dagger$ 由 `adjoint_be` 取共轭方向，$(\mu,\nu)$ 项取两者的 BE 乘积——signal 拼接、α 相乘。外层 LCU 的 PREPARE 在 $(\mu,\nu)$ 对上，权重 $\propto\lvert\zeta_{\mu\nu}\rvert\alpha_\mu\alpha_\nu$，signal 为 selector（项数所需位数）拼接各叶编码的 signal；ζ 全零时抛 `ValidationError`。

适用边界：叶矩阵 2..32 维；THC 分解相对真实哈密顿量的近似误差属于经典预处理，本模块对给定的 ζ 与叶做精确组装。大规模 ζ 与叶数据应改走 QRAM 数据绑定，产物可直接交给 `transforms.qubitization_walk`。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）：块编码的零信号块须与目标哈密顿量在容差内一致。三层证据：

- 结构：`tests/core/test_lowrank.py:ThcTests` 的构造属性断言（`be_form == "thc"`、`thc_leaves`、`leaf_count`）；`ThcTests.test_invalid_inputs_fail_at_generation` 覆盖 THC 负例（ζ 非对称、ζ 非方阵、全零系数）。
- 数值：`ThcTests.test_thc_block_equals_hamiltonian`——两个叶矩阵与 2×2 ζ 组装，期望矩阵按 $\sum_{\mu\nu}\zeta_{\mu\nu}L_\mu L_\nu^\dagger$ 在测试内独立组装后经 `assert_block_equals` 逐列对拍（places = 9）。
- 绑定：本模块的跨编码绑定见证是同文件 `DoubleFactorizationTests.test_matches_pauli_encoding_block`（同一哈密顿量的 DF 编码与 Pauli LCU 编码逐列对拍）；THC 的叶算符本身即显式 Pauli LCU 编码路径，不单独设 abstract 声明。

α 紧性见证（V1 新增）：`ThcTests.test_thc_alpha_matches_hand_computed_bound` 在测试内用独立的 `pauli_l1` 闭式（$M = c_I I + xX + yY + zZ$ 的系数由矩阵元线性表出，不经 `matrix_pauli_encoding`）手算叶上界 $\alpha_0 = 1.2$、$\alpha_1 = 1.3$，则 α = 0.7·1.44 + 0.1·1.56 + 0.1·1.56 + 0.4·1.69 = 1.996，与 `be.alpha` 对拍（places = 10）。

## 已知缺口与计划阶段

无已知缺口，阶段 V1 见证已齐（块对拍 + 独立手算的 α 口径对拍）。

## 相关链接

- 源码：`src/pyqecclang/algorithms/lowrank.py`
- 同模块算法：[双因子分解块编码](double-factorization.md)
- API 参考：[化学低秩分解块编码](../../api/algorithms/input_model/lowrank.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值验证脚本：`tests/verification/verify_blockencoding.py`（真实后端执行，无 mock、无 skip；2026-09-16 共 73 个案例全部通过），产物 `out/verification/blockencoding.json`。本页对应 `thc-encoding-leaves2-d2` 与 `thc-encoding-leaves3-d4` 两个案例。

实验设计：2×2 稠密叶（2 叶，一般复矩阵、非酉非 Hermitian）经 OriginIR-ext + UniQC `to_matrix` 全幺正提取，(0,0) 块乘 α 后与 numpy 独立组装的 $\sum_{\mu\nu}\zeta_{\mu\nu}L_\mu L_\nu^\dagger$ 逐元对拍，α 对照不经 `matrix_pauli_encoding` 的独立 Pauli l1 手算口径 $\sum_{\mu\nu}|\zeta_{\mu\nu}|\alpha_\mu\alpha_\nu$；4×4 三叶（对角叶 + 全耦合 3×3 ζ，7 个非零 LCU 项）走 reference + rir-pysparq + adapter-pysparq 三后端逐列，检验多叶多规模组装。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `thc-encoding-leaves2-d2` | 2 叶 2×2，α=1.381 | originir-ext + to_matrix，三后端交叉 | max_error | 1.7e-16 |
| 同上 | — | — | α − 独立 l1 口径 | 0（恰相等） |
| `thc-encoding-leaves3-d4` | 3 叶 4×4（对角叶），α=1.782 | 三后端 | max_error / 后端偏差 | 5.6e-16 / 0 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_blockencoding.py
```

产物：`out/verification/blockencoding.json`。
