# CKS Chebyshev 线性系统求解器（CKS Chebyshev QLSS）

> 类别 C2 · 模块 `pyqecclang.algorithms.qlss` · 阶段 V3

## 概述

稀疏 Hermitian 线性系统 $Ax = b$ 的基础 Chebyshev/LCU 求解路线（Childs–Kothari–Somma, arXiv:1511.02306 §4 的基础构造，无 VTAA）。把 $1/x$ 的截断 Chebyshev 展开作用到编码矩阵 $A/\alpha$ 的谱上：

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

与 HHL 论文参考值的端到端对拍缺失（阶段 V3，可入 catalog 目录）：当前见证覆盖输入适配、kappa 换算与范数恢复的局部语义，求解精度本身仍是 prototype 声明。与验证覆盖矩阵 `qlss.py` 行的缺口列一致。

## 相关链接

- 源码：`src/pyqecclang/algorithms/qlss.py`
- 同族页面：[Costa 行走求解器](costa-walk.md)、[稀疏矩阵块编码](sparse-block-encoding.md)
- API 参考：[量子线性系统](../../api/algorithms/qlss.rst)
- 输入模型审查：[QFVM 输入模型审查](../../reference/qfvm-input-models.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
