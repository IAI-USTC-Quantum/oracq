# CKS Chebyshev 线性系统求解器（CKS Chebyshev QLSS）

> 类别 C2 · 模块 `pyqecclang.algorithms.qlss` · 阶段 V3

## 概述

稀疏 Hermitian 线性系统 $Ax = b$ 的基础 Chebyshev/LCU 求解路线（Childs–Kothari–Somma, arXiv:1511.02306 §4 的基础构造）。第 5 节的 VTAA 变时层见 [VTAA-CKS 求解器](vtaa-cks.md)。把 $1/x$ 的截断 Chebyshev 展开作用到编码矩阵 $A/\alpha$ 的谱上：

$$
A^{-1} \approx \sum_{j=0}^{d-1} c_j\, T_{2j+1}(A/\alpha), \qquad
c_j = 4\,(-1)^j\, 2^{-2d} \sum_{i=j+1}^{d} \binom{2d}{d+i},
$$

经[稀疏矩阵块编码](sparse-block-encoding.md)得到 BE，`chebyshev_block` 生成奇次 Chebyshev 幂，LCU 以 $c_j$ 为权重组合逆算子 BE，最后作用到 RHS 制备上。矩阵性质（谱界、Hermitian）由调用者声明，语言不证明。

## 接口与输入模型

```python
make_cks_qlss(config=None)          # 返回 QLSSProtocol（input_model="sparse"）
cks_chebyshev(system, config=None)  # 底层内核，接收 SparseSystem
CKSConfig(order=2, terms=None)      # order 1..128；terms 为截断项数，≤ order
```

问题输入是 `LinearSystem(sparse=...)`，只能指定一种源输入模型：

- `SparseSystem(access, value_format, entry_bound, rhs, spectrum, diagonal_nonnegative, hermitian)`——input model 为 SO（位置 + 元素 oracle）加 SP（RHS 态制备）；`hermitian` 与 `diagonal_nonnegative` 必须显式声明 `True`，宽度与值格式须匹配。
- `SpectralPromise(norm_upper, sigma_min_lower, evidence)`——调用者声明的谱界（缺省 `evidence = "caller_declared_unverified"`），要求 $0 < \sigma_{\min}^{\text{lower}} \le$ 范数上界；`inverse_bound(alpha) = max(1, alpha / sigma_min_lower)`，与 alpha 冲突时拒绝。

`protocol(problem)` 返回 `SolveResult`：

| 属性 | 含义 |
|---|---|
| `state` | 物理子空间解态（`StateOracle`，`.operation` 为其 RIR） |
| `norm_probe` | 独立矩阵范数探针（`StateOracle`） |
| `input_alpha` / `encoded_inverse_bound` | 消费的 BE 归一化与 $\max(1, \alpha/\sigma_{\min})$ |
| `rhs_norm` / `adapter_trace` | 经典右端范数与适配轨迹 |
| `recover_norm(p_solver, p_joint)` | $\lVert x\rVert = \lVert r\rVert / (\alpha\sqrt{p_{joint}/p_{solver}})$ |

内核输出属性：`algorithm = "cks_chebyshev_basic"`、`input_model = "sparse_location_inplace_and_entry_xor"`、`polynomial_order` / `polynomial_terms`、`inverse_lcu_normalization`、`implementation_scope = "CKS section 4 basic LCU; no VTAA"`、`correctness = "pending"`。

## 实现要点

生成链：`real_symmetric_sparse_encoding`（CKS $T^\dagger S T$，$\alpha = s \cdot a_{\max}$）→ `chebyshev_block(a, 2j+1)` 奇次幂 → `lcu` 组合逆算子 BE → `apply_be_to_state` 作用到 RHS。协议层 `solve()` 统一做：契约检查（`check()` 不运行内核，按 `INPUT_TYPE` / `INPUT_ZERO_RHS` / `INPUT_ADAPTER` / `INPUT_PROMISE` / `INPUT_SPECTRUM` 出报告）、`select_subspace` 选取物理通道、组装独立范数探针——Costa filtering 的成功率不能套用 CKS 逆算子 LCU 的归一化因子，两者成功分支形成方式不同；探针属性 `probability_contract = "joint solver success and fresh matrix-ancilla success"`。输出要求 adjoint 与 controlled 能力（供 QFVM 组合）。

适用边界：仅论文 §4 基础 LCU，无 VTAA；解精度与成功通道标注 prototype/pending，不构成求解精度承诺。BE 输入不能自动恢复稀疏 oracle（反向适配不存在）；本 protocol 也没有两参数 legacy 调用形状，只接收 `LinearSystem`。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。与验证覆盖矩阵 `qlss.py` 行一致：

- 结构：`tests/core/test_qlss_input_models.py:QLSSInputTests`（构造与契约断言）。
- 数值：`test_signed_sparse_encoding_and_chebyshev`（负非对角 2×2 的 BE 角块与 alpha = 2 对拍、`chebyshev_block(be, 3)` 与 $4h^3 - 3h$ 对拍，places = 11）；`test_matrix_probe_recovers_scalar_system_norm`（标量系统经探针恢复 $\lVert x\rVert = 4$，places = 10）；`test_norm_recovery_uses_conditional_matrix_probe`（`recover_norm(0.5, 0.125) = 2.0`，概率次序非法时拒绝）；`test_protocols_consume_different_models_and_scale_kappa`（同一问题两种 protocol 输出同宽解态、相同 input_alpha，程序无未绑定槽位）。
- 绑定：`test_costa_rhs_reflection_is_independent_of_unitary_extension`（U 扩展独立性，登记于本行）。拒绝路径由 `test_no_implicit_be_to_sparse`（BE 输入不得走稀疏 protocol）与 `test_zero_rhs_is_not_a_state_preparation_problem`（零右端经典侧处理）钉死。

## 已知缺口与计划阶段

与 HHL 论文参考值的端到端对拍缺失（阶段 V3，可入 catalog 目录）：当前见证覆盖输入适配、kappa 换算与范数恢复的局部语义，求解精度本身仍是 prototype 声明。VTAA 变时层已由 [VTAA-CKS 求解器](vtaa-cks.md) 单独实现；本页保持 §4 基础路线的定位。与验证覆盖矩阵 `qlss.py` 行的缺口列一致。

## 相关链接

- 源码：`src/pyqecclang/algorithms/qlss.py`
- 同族页面：[Costa 行走求解器](costa-walk.md)、[VTAA-CKS 变时求解器](vtaa-cks.md)、[稀疏矩阵块编码](sparse-block-encoding.md)
- API 参考：[量子线性系统](../../api/algorithms/qlss.rst)
- 输入模型审查：[QFVM 输入模型审查](../../reference/qfvm-input-models.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_nt_qlss_sde.py`（nt_qlss_sde 组），全部在真实后端执行。经典预言机全部独立：Chebyshev 矩阵多项式 $P(M)=\sum_j c_j T_{2j+1}(M)$ 用 numpy 以 $T_2$ 递推直接求值（系数由 `math.comb` 闭式重算），真解由 `numpy.linalg.solve` 给出，均不经过被测实现的辅助函数。

**实验设计**：(a) 内核收敛扫描——$\kappa=3$ 的 $2\times2$ 有符号稀疏系统（$A=[[0.75,-0.25],[-0.25,0.75]]$，$\alpha=1.5$），`cks_chebyshev` 在 order 2/4/8/16 下于 reference 与 rir-pysparq 运行，条件解态（signal = 0 分支归一化）与成功概率对照多项式预言（实现误差），方法误差对照 numpy 真解；(b) 协议级——`make_cks_qlss(CKSConfig(order=8))` 求解 `LinearSystem`，$p_{\text{solver}}$ 与独立矩阵范数探针的 $p_{\text{joint}}$ 对照多项式预言，`recover_norm` 对照 $\lVert A^{-1}b\rVert$；(c) Costa 构件——Dolph–Chebyshev 权重对照闭式窗 $\gamma\,T_d(\beta\cos\theta)$（numpy chebval 4097 点采样）、`schedule` 闭式与单调性、`unary_weight_preparation` 概率分布；`costa_qlss` 装配程序做三后端对拍（其 `kernel_status` 为库内声明的 prototype，求解精度不作判据）。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| cks-kernel-order-2/4/8 | $2\times2$、$\kappa=3$ | 双路径 | 实现误差 / 成功概率误差 | ≤ 4.9e-16 / ≤ 1.7e-16 |
| cks-kernel-order-16 | 同上 | 双路径 | 实现误差 / 成功概率误差 / 跨后端 | 7.2e-12 / 8.7e-13 / 2.5e-7 |
| cks-method-convergence | order 2→16 | numpy oracle | 方法误差 | 0.5537 → 0.4084 → 0.2130 → 0.0662（严格递减，约 $e^{-4/3}$ 每阶） |
| cks-protocol-norm-recovery | order 8 | 双路径 | 探针概率误差 / recover_norm 相对误差 / 保真度 | 1.1e-16 / 0.1456 / 0.9531 |
| costa-plan-schedule-unary | degree ≤ 6 | numpy + 双路径 | DC 权重 / schedule / unary 制备 | 2.7e-15 / 2.2e-16 / 1.7e-16 |
| costa-assembly-cross-backend | steps=1 | 三路径 | 逐振幅偏差 | 0.0 |

recover_norm = 2.7018 对照 $\lVert A^{-1}b\rVert = 3.1623$（$b=(2,0)$）：0.1456 的相对误差与 order=8 的方法误差 0.213 自洽，随阶数按几何率收敛（kernel 扫描）。order=16 处 reference 与多项式预言保持 ~1e-16 一致，rir-pysparq 在约 $10^6$ 展开步后出现 2.5e-7 的单振幅浮点漂移（物理可观测量一致到 ≤ 8.7e-13），已在最终报告记录。

**结构层模块的覆盖说明**：`contracts.py` 与 `interfaces.py` 是契约/协议结构层，本身不含独立数值语义——其代码路径（`ProtocolContract.check`、协议适配器 `as_sparse_access` / `as_block_encoding` 等）在协议级案例中被真实执行，正确性由 `tests/core` 结构测试与本组的跨后端对拍联合覆盖；`algorithms/legacy.py`（LCHS/Schrödingerization 工厂）与 `applications/legacy.py`（QFVM/QHAM 组装）同理，其数值内容归属底层算法页，组装确定性由 catalog 案例（`lchs`、`schrodingerisation`、`qfvm_*`、`qham_*`）的跨后端对拍覆盖。

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_nt_qlss_sde.py
```

产物：`out/verification/nt_qlss_sde.json`（`cks-*`、`costa-*` 共 7 个案例）。
