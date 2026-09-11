# 稀疏矩阵块编码（Sparse Matrix Block Encoding）

> 类别 C2 · 模块 `pyqecclang.algorithms.sparse` · 阶段 V1

## 概述

把稀疏访问 oracle（位置与元素分离，input model 为 SO）适配为块编码。正式入口 `real_symmetric_sparse_encoding` 面向实对称（Hermitian）、非负对角矩阵，构造 CKS 型 $T^\dagger S T$（Childs–Kothari–Somma, arXiv:1511.02306，接口约定见 [QFVM 输入模型审查](../../reference/qfvm-input-models.md)）：先做幅度转导制备 $T$——在 `neighbor` 上制备稀疏位置的前缀均匀叠加、查位置 oracle 得行索引、查元素 oracle 得条目字、按 $\sqrt{|v|/a_{\max}}$ 旋转成功旗标；再交换两侧坐标与两侧失败旗标得到 $S$。若每列枚举 $s$ 个结构位置且 $a_{\max}$ 确实约束元素幅值，则 $T^\dagger S T$ 的零信号块为 $A/(s \cdot a_{\max})$。

## 接口与输入模型

```python
real_symmetric_sparse_encoding(access, fmt, amax, *, diagonal_nonnegative=False, rotation=None)
```

- `access`：`SparseAccess`——`location`（column/index/work 原地置换）与 `entry`（row/column/data XOR）两个操作。
- `fmt`：`FixedFormat` 条目定点格式，宽度须等于 `access.value_width`；`fmt.signed` 决定是否叠加符号相位。
- `amax`：元素幅值上界 $a_{\max}$，有限正数。
- `diagonal_nonnegative`：必须显式传 `True`，否则报"当前对称稀疏适配要求非负对角"。
- `rotation`：幅度转导操作，缺省 `magnitude_rotation(fmt, amax)`。

返回 `BlockEncoding`，模块属性：

| 属性 | 含义 |
|---|---|
| `be_alpha` | $s \cdot a_{\max}$（$s$ 为 `access.sparsity`） |
| `construction` | `"CKS_Tdag_S_T"` |
| `self_adjoint_extension` | `True`（Chebyshev walk 的前置条件） |
| `correctness` | `"pending"` |

辅助入口：`magnitude_rotation(fmt, amax)`（值宽 ≤ 12 位给显式门实现，超宽返回带 `amplitude_contract` 属性的待绑定声明）、`prefix_state(width, count)`（前缀均匀叠加）、`compare_words(width, kind)`（eq/lt 布尔网络，缓存复用）、`chebyshev_block(a, degree)`（见下）。旧入口 `sparse_block_encoding` 已标注 legacy（`legacy_input_model=True`、`matrix_contract` 声明"legacy transduction unspecified"），仅供旧目录描述，不是一般稀疏输入适配器。

## 实现要点

寄存器布局：`target = Bits(n)`、`signal = Bits(n+2)`；`signal[:n]` 为 neighbor，`signal[n]` / `signal[n+1]` 为两侧失败旗标，条目字与位置 work 是局部寄存器、复净后交还。制备序列为叠加 → 位置查询 → 元素查询 → 幅度旋转（→ 带符号格式的方向相位）；外层 `swap(target, signal[:n])` 与 `swap(signal[n], signal[n+1])` 实现 $S$，两侧旗标交换是必要部分，只换索引不换旗标会破坏自伴性。负非对角元素按 `target < neighbor` 的方向约定叠加 $\pi$ 相位（`sign_convention` 属性记录）。

`chebyshev_block(a, degree)`：要求输入 BE 显式声明 `self_adjoint_extension`（普通 BE 不足以保证 Chebyshev 语义），`Repeat(degree)` 交替"信号零态正反射 + 调用 a"，返回 `be_alpha = 1.0`、`argument_scale = a.alpha` 的行走幂，其零信号块实现 $T_k(A/\alpha)$。

适用边界：仅实 Hermitian 加非负对角；一般矩阵须先做显式 Hermitian dilation。值宽超过 12 位时幅度转导是开放声明，绑定前程序不可执行。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。与验证覆盖矩阵 `sparse.py` / `block_encoding.py` 行一致：

- 结构：`tests/core/test_language.py:BlockEncodingTests`（BE 构造与属性断言）。
- 数值：本入口的直接见证登记于 `qlss.py` 行——`tests/core/test_qlss_input_models.py:QLSSInputTests.test_signed_sparse_encoding_and_chebyshev`：2×2 含负非对角矩阵（$s = 2$、$a_{\max} = 1$）逐列读出零信号块乘回 `alpha = 2` 与经典矩阵对拍（places = 11），且 `chebyshev_block(be, 3)` 的零信号块与 $4h^3 - 3h$（$h = A/\alpha$，即 $T_3$）逐元素对拍。
- 绑定：`tests/core/test_language.py:BlockEncodingTests.test_alpha_survives_ir_serialization`（JSON 往返保留 `be_alpha`）。

## 已知缺口与计划阶段

无已知缺口，阶段 V1（验证覆盖矩阵 `sparse.py` / `block_encoding.py` 行缺口列为空）。符号相位约定与 alpha 换算已有小实例见证；大字长幅度转导以显式待绑定声明的形式保留，其量化误差界未单列见证。

## 相关链接

- 源码：`src/pyqecclang/algorithms/sparse.py`
- 同族页面：[块编码组合代数](block-encoding-algebra.md)、[CKS Chebyshev 求解器](cks.md)、[Costa 行走求解器](costa-walk.md)
- API 参考：[稀疏访问适配](../../api/algorithms/sparse.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
