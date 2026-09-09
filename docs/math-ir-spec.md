# MIR 0.1：纯数学函数图

MIR 是 Python 纯函数前端与 RIR 之间的生成层表示。它可以独立 JSON 往返，完整格式见 [Schema](math-ir.schema.json)。它不替代 RIR 0.3，也不把 Python 回调写入量子 IR。

## 对象与类型

MathProgram 包含 version="0.1"、entry 和 functions。每个 MathFunction 包含唯一 name、源码 label、parameters、有序 nodes 和 returns。函数图无递归，helper 调用保留为 call 节点。

参数类型为 real、complex、bool 或 index。index 另带 1..64 的无符号位宽；其他参数 width=0，具体量子字长由 lowering 的 FixedFormat 决定。Index 在进入数学计算时转换为当前定点表示，格式必须容纳它的完整范围。值节点只有 real、complex、bool 三类；复数 lowering 到两个独立实数寄存器。

MathNode 包含 op、kind、args 和 data。节点编号为其在 nodes 中的位置，args 只能引用编号更小的节点。returns 是非空节点编号数组；返回 tuple 在这里表示多个独立结果，不允许递归嵌套 tuple 输出。函数参数、节点类型、元数、调用接口以及 DAG 均由 MathProgram.validate 核对。

| op | 参数 / data | 结果 |
|---|---|---|
| input | data=[参数名] | 参数的数学类型；index 转 real |
| const | data=[实数或布尔]；复数为 [实部,虚部] | 声明的有限数值类型 |
| add/sub/mul/div/pow | 两个数值节点 | real 或提升到 complex |
| neg/conj | 一个数值节点 | 与输入同类型 |
| real/imag/abs | 一个数值节点 | real |
| lt/eq | 两个可比较节点 | bool；complex 不支持 lt |
| and/or/not | 布尔节点 | bool |
| complex | 两个 real 节点 | complex |
| select | bool、true 值、false 值 | 两分支的公共类型 |
| intrinsic | data=[数学函数名] | 按内建签名 |
| call | data=[函数符号,返回索引]，args 按形参顺序 | 被调函数对应返回类型 |

同一 helper 的同参数调用与公共表达式可以复用。一个多结果调用在 MIR 中具有不同返回索引，lowering 将它们合并为一次 RIR Call。Python 的静态 range 循环在有界生成阶段展开成 SSA；这不改变 helper 模块或已有 RIR Repeat 的保留规则。

## 有限数值与生成配置

MIR 表达函数结构；生成配置包括 FixedFormat 和 MathConfig。前者指定二补码字长/小数位，后者指定数学核阶数与近似区间。这些参数及采样得到的系数记录在 RIR 属性中。相同数学图在不同配置下得到不同的 RIR 模块符号，可同时组装。

lowering 使用现有定点加减乘除/开方等 Boolean 电路。非多项式实函数采用 Chebyshev 采样系数与 Clenshaw 递推；采样数由阶数决定，独立于输入位模式总数，不生成整个函数的真值表。复数函数以实数核和代数关系分解。

具体量子电路精确实现其有限布尔函数 F_tilde 的可逆扩张；F_tilde 与理想数学函数 F 的接近程度不由本语言证明。整数常量幂采用平方-乘，绝对指数上限 128；一般实数幂采用 exp(y log(x))，实数路径要求正底数。复数路径使用选定的复对数分支。

## 输出和状态

公开接口为输入端口、输出端口和 status[2]：

```text
|inputs, outputs, status, 0_private>
  ↦ |inputs, outputs XOR F_tilde(inputs),
     status XOR flags(inputs), 0_private>
```

输出不是赋值覆盖；非零输出同样合法。编译器将所需中间量计算到私有 locals，复制结果和标志后反算，输入端口保持不变。结构检查不构成复净证明；本阶段用小规模参考与真实后端检查这条契约。

status[0] 记录定义域失效，例如除零、实数负数开方和非正实数取对数。status[1] 汇总底层溢出/常量越界以及数学核超出配置区间。它们不是精度界，也不表示 IEEE NaN/Inf。

select 的结果状态为条件状态 OR 被选中分支状态。未选中分支虽在可逆实现中计算和反算，其定义域标志不污染最终状态。helper 的返回状态是其所有返回值的状态合并，调用者再合并已求值实参的状态。无返回依赖的死计算可以消除；不模拟 Python 浮点异常时序。

没有 IEEE 有符号零，复函数在切线上的边界约定无法逐位复现 Python cmath；当前零虚部采用非负一侧的约定。复杂分支切线、近似精度和全域状态行为待核验。

## 前端边界

compile_function 读取普通 Python 函数源码并解释受限 AST；不执行待编译函数，也不使用动态 eval/exec。源码字符串可指定 entry；没有指定时使用最后一个 def。支持数值常量、局部赋值/解包、算术、比较、条件表达式、结构化 if、静态有界 range、tuple 返回、纯 helper 和白名单 math/cmath 调用。

拒绝 I/O、对象突变、任意对象方法、动态循环、递归、异常处理、生成器、lambda 及不能读取源码的可调用对象。闭包和显式 constants 只捕获有限数值。源码不可用时可传 def 字符串。math 和 cmath 别名以及 from 导入均可识别；数学 intrinsic 当前用位置参数。

数学函数名目录依据 [Python cmath 文档](https://docs.python.org/3/library/cmath.html)，AST 节点依据 [Python AST 文档](https://docs.python.org/3/library/ast.html)。支持的语法/数值子集和定点行为由本规范限定，不能视作任意 Python 或完整 cmath 兼容实现。
