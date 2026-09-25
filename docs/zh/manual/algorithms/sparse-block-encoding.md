# 稀疏矩阵块编码（Sparse Matrix Block Encoding）

<a href="../../../index.html">English</a> · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.input_model.sparse`](../../api/algorithms/input_model/sparse.rst) · 阶段 V1

## 概述

把稀疏访问 oracle（位置与元素分离，input model 为 SO）适配为块编码。正式入口 {obj}`real_symmetric_sparse_encoding <oracq.algorithms.input_model.sparse.real_symmetric_sparse_encoding>` 面向实对称（Hermitian）、非负对角矩阵，构造 CKS 型 $T^\dagger S T$（Childs–Kothari–Somma, arXiv:1511.02306，接口约定见 [QFVM 输入模型审查](../../reference/qfvm-input-models.md)）：先做幅度转导制备 $T$——在 `neighbor` 上制备稀疏位置的前缀均匀叠加、查位置 oracle 得行索引、查元素 oracle 得条目字、按 $\sqrt{|v|/a_{\max}}$ 旋转成功旗标；再交换两侧坐标与两侧失败旗标得到 $S$。若每列枚举 $s$ 个结构位置且 $a_{\max}$ 确实约束元素幅值，则 $T^\dagger S T$ 的零信号块为 $A/(s \cdot a_{\max})$。

## 接口与输入模型

```python
real_symmetric_sparse_encoding(access, fmt, amax, *, diagonal_nonnegative=False, rotation=None)
```

API 入口：{obj}`real_symmetric_sparse_encoding <oracq.algorithms.input_model.sparse.real_symmetric_sparse_encoding>`

- `access`：{obj}`SparseAccess <oracq.algorithms.input_model.oracles.SparseAccess>`——`location`（column/index/work 原地置换）与 `entry`（row/column/data XOR）两个操作。
- `fmt`：{obj}`FixedFormat <oracq.algorithms.common.arithmetic.FixedFormat>` 条目定点格式，宽度须等于 `access.value_width`；`fmt.signed` 决定是否叠加符号相位。
- `amax`：元素幅值上界 $a_{\max}$，有限正数。
- `diagonal_nonnegative`：必须显式传 `True`，否则报"当前对称稀疏适配要求非负对角"。
- `rotation`：幅度转导操作，缺省 {obj}`magnitude_rotation(fmt, amax) <oracq.algorithms.input_model.sparse.magnitude_rotation>`。

返回 {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`，模块属性：

| 属性 | 含义 |
|---|---|
| `be_alpha` | $s \cdot a_{\max}$（$s$ 为 `access.sparsity`） |
| `construction` | `"CKS_Tdag_S_T"` |
| `self_adjoint_extension` | `True`（Chebyshev walk 的前置条件） |
| `correctness` | `"pending"` |

辅助入口：`magnitude_rotation(fmt, amax)`（值宽 ≤ 12 位给显式门实现，超宽返回带 `amplitude_contract` 属性的待绑定声明）、{obj}`prefix_state(width, count) <oracq.algorithms.input_model.sparse.prefix_state>`（前缀均匀叠加）、{obj}`compare_words(width, kind) <oracq.algorithms.input_model.sparse.compare_words>`（eq/lt 布尔网络，缓存复用）、{obj}`chebyshev_block(a, degree) <oracq.algorithms.input_model.sparse.chebyshev_block>`（见下）。旧入口 {obj}`sparse_block_encoding <oracq.algorithms.input_model.sparse.sparse_block_encoding>` 已标注 legacy（`legacy_input_model=True`、`matrix_contract` 声明"legacy transduction unspecified"），仅供旧目录描述，不是一般稀疏输入适配器。

## 实现要点

寄存器布局：`target = Bits(n)`、`signal = Bits(n+2)`；`signal[:n]` 为 neighbor，`signal[n]` / `signal[n+1]` 为两侧失败旗标，条目字与位置 work 是局部寄存器、复净后交还。制备序列为叠加 → 位置查询 → 元素查询 → 幅度旋转（→ 带符号格式的方向相位）；外层 `swap(target, signal[:n])` 与 `swap(signal[n], signal[n+1])` 实现 $S$，两侧旗标交换是必要部分，只换索引不换旗标会破坏自伴性。负非对角元素按 `target < neighbor` 的方向约定叠加 $\pi$ 相位（`sign_convention` 属性记录）。

`chebyshev_block(a, degree)`：要求输入 BE 显式声明 `self_adjoint_extension`（普通 BE 不足以保证 Chebyshev 语义），{obj}`Repeat(degree) <oracq.infrastructure.ir.Repeat>` 交替"信号零态正反射 + 调用 a"，返回 `be_alpha = 1.0`、`argument_scale = a.alpha` 的行走幂，其零信号块实现 $T_k(A/\alpha)$。

适用边界：仅实 Hermitian 加非负对角；一般矩阵须先做显式 Hermitian dilation。值宽超过 12 位时幅度转导是开放声明，绑定前程序不可执行。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。与验证覆盖矩阵 `sparse.py` / `block_encoding.py` 行一致：

- 结构：`tests/core/test_language.py:BlockEncodingTests`（BE 构造与属性断言）。
- 数值：本入口的直接见证登记于 `qlss.py` 行——`tests/core/test_qlss_input_models.py:QLSSInputTests.test_signed_sparse_encoding_and_chebyshev`：2×2 含负非对角矩阵（$s = 2$、$a_{\max} = 1$）逐列读出零信号块乘回 `alpha = 2` 与经典矩阵对拍（places = 11），且 `chebyshev_block(be, 3)` 的零信号块与 $4h^3 - 3h$（$h = A/\alpha$，即 $T_3$）逐元素对拍。
- 绑定：`tests/core/test_language.py:BlockEncodingTests.test_alpha_survives_ir_serialization`（JSON 往返保留 `be_alpha`）。

## 已知缺口与计划阶段

无已知缺口，阶段 V1（验证覆盖矩阵 `sparse.py` / `block_encoding.py` 行缺口列为空）。符号相位约定与 alpha 换算已有小实例见证；大字长幅度转导以显式待绑定声明的形式保留，其量化误差界未单列见证。

## 相关链接

- 源码：`src/oracq/algorithms/input_model/sparse.py`
- 同族页面：[块编码组合代数](block-encoding-algebra.md)、[CKS Chebyshev 求解器](cks.md)、[Costa 行走求解器](costa-walk.md)、[稀疏矩阵访问](sparse-access.md)、[Select-Swap QROM](select-swap.md)、[VTAA-CKS 变时求解器](vtaa-cks.md)
- API 参考：[稀疏访问适配](../../api/algorithms/input_model/sparse.rst)
- 概念：[Oracle 与算子表示](../operators.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值验证脚本：`tests/verification/verify_blockencoding.py`（真实后端执行，无 mock、无 skip；2026-09-16 共 73 个案例全部通过），产物 `out/verification/blockencoding.json`。本页对应 `sparse-be-*`、`chebyshev-walk-*`、`cross-tridiagonal-*`、`cross-qram-be-*` 共 14 个案例。

实验设计：

1. **稀疏 BE 多矩阵多规模**：三对角矩阵 $\alpha I + \beta T$（带符号格式 $\beta<0$ 检验方向相位约定、无符号格式检验无符号路径），dim ∈ {2, 4, 8, 16}，每列 $s$ 个结构位置、边界列以零条目位置补足。零信号块逐列提取（reference 全列 + rir-pysparq / adapter-pysparq 交叉），乘回 $\alpha = s\cdot a_{\max}$ 后与经典矩阵逐元对拍；dim = 2（15 qubits）与无符号 dim = 4（12 qubits）另走 OriginIR-ext 态向量逐列提取，带符号 dim = 4 达 24 qubits OriginIR 预算、dim = 8 为 34 qubits 超预算，均按预算只走 pysparq/reference 路径并在案例参数注明。
2. **QRAM 数据绑定稀疏 BE**（dim = 4）：位置正/反表与元素表全部经 memory 绑定（不进 IR），reference 与 rir 逐振幅一致。
3. **Chebyshev 行走**：2×2 稀疏 BE 的 `chebyshev_block` 阶数 $k=1..4$，零信号块对照 $T_k(A/\alpha)$。
4. **与 pysparq 自带块编码的独立交叉验证**：同一三对角矩阵分别经 oracq 稀疏 BE（及 dim ≤ 8 时的 Pauli LCU BE 第二路径）与 pysparq `BlockEncodingTridiagonal` 编码，各自提取有效块乘自身归一化后互相对拍（dim ∈ {2, 4, 8, 16}）；另以 pysparq `BlockEncodingViaQRAM`（C++ 判据配置 data_size=50、rational=51、exponent=15，矩阵 Frobenius 归一化入表）对拍三对角与非三对角（对匹配稀疏图 $s=2$）两个 dim = 4 实例。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `sparse-be-signed-tridiagonal-d2` | dim 2，s=2，α=3.0 | originir-ext + 三路径 | max_error | 4.4e-16 |
| `sparse-be-signed-tridiagonal-d4` | dim 4，s=3，α=4.5 | 三路径 | max_error | 1.3e-16 |
| `sparse-be-signed-tridiagonal-d8` | dim 8，s=3，α=4.5 | reference + rir | max_error | 1.3e-16 |
| `sparse-be-unsigned-unitary-d4` | dim 4，s=3，α=4.5 | originir-ext + 三路径 | max_error | 1.1e-16 |
| `sparse-be-unsigned-d16` | dim 16，s=3，α=4.5 | reference + rir | max_error | 1.1e-16 |
| `sparse-be-qram-access-d4` | dim 4，QRAM 绑定 | reference + rir | max_error / 后端偏差 | 1.3e-16 / 0 |
| `chebyshev-walk-k1..k4` | dim 2，$\alpha=3.0$ | 三路径 | max_error | 1.7e-16 / 2.5e-16 / 5.6e-16 / 6.7e-16 |
| `cross-tridiagonal-*`（5 组） | dim 2–16 | reference × pysparq | qecc / pysparq / 交叉误差 | ≤4.4e-16 / ≤4.4e-16 / ≤6.7e-16 |
| `cross-qram-be-tridiagonal-d4` | dim 4，s=3 | reference × pysparq QRAM | qecc 误差 / pysparq 量化误差 | 1.3e-16 / 2.8e-5 |
| `cross-qram-be-matched-pairs-d4` | dim 4，s=2 非三对角 | 同上 | 交叉误差 | 2.8e-15 |

交叉验证发现（如实记录，非断言库缺陷）：pysparq `BlockEncodingTridiagonal` 的 (0,0) 块归一化在主寄存器 1 位（dim = 2）且 $\beta \neq 0$ 时实测为 $|\alpha|+2|\beta|$，而非其构造文档的 Frobenius 范数（1 位寄存器上加一/减一都触发溢出分支，anc==0 角块失去 Frobenius 归一；$\beta=0$ 时无移位分支、仍等于 $A/\lVert A\rVert_F$）。其 C++ 正确性测试域 `randint(2,5)`（dim 4–16）不覆盖该规模；脚本按实证归一化对拍并在案例参数 `psparq_alpha_effective` 记录两种口径。

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_blockencoding.py
```

产物：`out/verification/blockencoding.json`。
