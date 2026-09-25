# 一般 QHAM：从 PDE 到量子适配线性系统的推导

<a href="../../reference/qham-derivation.html">English</a> · **简体中文**

本文先固定数学规则，再据此实现生成器。主要对照 [QHAM v2](https://arxiv.org/html/2411.06759v2) II.1–II.3；其 [v1](https://arxiv.org/html/2411.06759v1) 使用 secondary linearization 名称。这里“二次”表示第二次线性化，不是将所有 PDE 先降为二次多项式。工程描述见手册[一般 QHAM 自动生成：PDE → HAM → QCL → QODE](../manual/qham.md)，可运行示例见教程[从 PDE 表达式生成 QHAM 输入](../tutorials/qham.md)。

## 1. 输入范围与约定

考虑一阶时间演化系统，在空间离散后统一记为

```text
u' = f + L u + sum_{r=2}^D B_r(u,...,u),    u(0)=u_in.
```

u 可以包含多分量和多维网格。每个 B_r 是已知的 r 线性映射 V^tensor(r) → V，例如 Burgers 中的 -u*u_x，KdV 中的 -6*u*u_x，以及三次反应项。空间导数、已知空间系数、分量选择和同点乘积归入这些映射，不作为量子态的非线性门。

首批自动表示规则是自治、有限阶多项式演化 PDE，允许已知空间强迫与系数。未知场分母、任意超越非线性、隐式时间方程和未消去的约束不在此规则内。规则内的多项式次数、HAM 截断阶数、分量数和空间维数不固定为某个案例。时变 L、f、B_r 和 eta 可以沿用下述代数，但需要时间依赖 QODE 输入，不能把它们直接冻结后称为同一个问题。

选择 L_H = partial_t - L，U0'=L U0+f，U0(0)=u_in，Ui(0)=0（i>0）。令 eta=hH，当前量子组装中取常数。

## 2. 从同伦方程推导高阶递推

从 (1-q)L_H(U-U0)-q*eta*(partial_t U-N(U))=0 出发，写 U(q)=sum_i q^i Ui，并定义

```text
C_l = sum_{r=2}^D sum_{a1+...+ar=l} B_r(U_a1,...,U_ar),
D_i = (partial_t-L) Ui.
```

比较 q 的系数得到 D_1=-eta*C_0，D_i=(1+eta)D_(i-1)-eta*C_(i-1)（i>1）。因此

```text
Ui' = L Ui - eta * sum_{l=0}^{i-1} (1+eta)^(i-1-l) C_l.       (A)
```

这条式子是生成器的依据。特别地，m=2 时

```text
U0' = L U0 + f
U1' = L U1 - eta B_2(U0,U0)
U2' = L U2 - eta [B_2(U0,U1)+B_2(U1,U0)]
              - eta(1+eta) B_2(U0,U0)
```

最后一项在 eta=-1 时恰好消失；只用这个参数测试会漏掉一般 h 下的错误。不能仅将论文局部耦合示意式中的 -eta 复制到每条高阶边上。

强迫只显式出现在 U0 方程，其他 Ui 通过 U0 间接受到强迫影响。这依赖于上述初猜的选择。

## 3. 独立坐标上的张量变量

对有序非负整数字 a=(a0,...,a_(k-1))，定义

```text
Y_a(x0,...,x_(k-1)) = product_j U_aj(xj).
```

这些因子有独立空间坐标。Y_(0,1) 与 Y_(1,0) 不能未经坐标置换直接合并。同点乘法必须在对相应因子作用空间算子后，再进行对角收缩。离散情况下，B_r 就是这样的矩形线性映射。实现上，{obj}`structured_fd_bindings <oracq.applications.qham.stencils.structured_fd_bindings>` 从移位与收缩生成基本矩阵，不物化端口矩阵。

用乘积法则只求一次时间导数：

```text
Y_a' = sum_p I ⊗ ... ⊗ (U_ap') ⊗ ... ⊗ I.
```

代入 (A)，得到三类线性边：

1. 线性项：同一字到同一字，在第 p 个坐标作用 L。
2. a_p=0 的强迫项：从删除该因子的字出发，在第 p 个位置插入已知向量 f。
3. a_p>0 的非线性项：将该因子替换为有序 r 元组 b，sum(b)=l<a_p；对这 r 个输入坐标作用 B_r 收缩到一个输出坐标，系数为 -eta(1+eta)^(a_p-1-l)。

因此，不需要把 U_i 的量子态作为另一个求解器可反复复制的输入。整个系统只需要初始数据与固定线性算子/多线性映射的访问。

## 4. 有限闭包证明

令 p=max(1,D-1)，设字的权重为

```text
weight(a)=p*sum(a)+len(a),    W=p*m+1.
```

最初的单因子 Ui（i≤m）都满足 weight≤W。非线性替换后

```text
weight(new)-weight(old)
 = p*(l-a_p)+(r-1) <= -p+(D-1) <= 0.
```

并且阶数和严格下降。强迫边删除一个零阶因子，权重减少 1；线性边保持不变。因此所有边都落在有限集合 weight(a)≤W 中，最多需要 W 个动态因子。全零字只产生线性项和强迫降阶项，闭包终止。

这个闭包对给定的 HAM 截断系统是精确的，不是再做一次 Carleman 阶数截断。生成器侧的 {obj}`QHAMPlan <oracq.applications.qham.linearization.QHAMPlan>` 按同一权重判据惰性枚举闭包，不物化矩阵。与原始非线性 PDE 的差异仍来自 HAM 截断、空间离散和后续数值/量子求解。

对于二次 PDE，p=1，字数为 2^(m+1)-1。另加物理输出块后是论文的 2^(m+1) 个函数块。每个单坐标空间维数为 N 时，原始总维数为

```text
N + sum_{k=1}^{m+1} binom(m+1,k) N^k
  = (N+1)^(m+1)+N-1.
```

高次情形采用上述有限闭包的规范超集；不声称所有块都必要或达到论文的最优资源规模。

## 5. 物理解和强迫齐次化

额外引入物理块 Y_phys=sum_{i=0}^m Ui。对 (A) 求和得到

```text
Y_phys' = L Y_phys + f
        + sum_{l=0}^{m-1} [1-(1+eta)^(m-l)] C_l.              (B)
```

该关系在 eta=0 也成立，不需要除以 eta。初值为 Y_phys(0)=u_in。张量字初值只有全零字非零：Y_(0,...,0)(0)=u_in^tensor(k)，其余初值为零。

若有强迫，引入常量分量 c'=0、c(0)=1，把来自空字的 f 注入也写入线性生成元，得到一个齐次增广 QODE。强迫对张量块的作用是低秩字到高秩字的矩形注入，不能只给 U0 加一个右端向量而漏掉张量方程的强迫项。

## 6. 验证方法与边界

数学验证不依赖 QLSS 输出。对任意随机 Ui，分别计算：

```text
chain_rule = derivative_of_Lift(U) applied to HAM_rhs(U)
linear_rhs = G * Lift(U).
```

两者应逐分量一致。测试需要覆盖 m>1、eta 不等于 -1、非零 f、三次项、混合多项式和耦合分量。m=1 无强迫时还应还原已有 [u,U0,U1,U0⊗U0] 模型。

闭包证明不证明 HAM 收敛，也不保证任意生成元适合 LCHS/CBMD 的耗散性假设。显式生成所有块的实现具有组合规模；后续可利用惰性逐行规则设计更高效的稀疏 oracle。MIR/RIR 中保留模块并不自动带来论文复杂度优势。
