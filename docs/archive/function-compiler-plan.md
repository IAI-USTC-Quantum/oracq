# 纯数学函数到可逆量子模块：实施计划

目标：受限、无副作用的普通 Python 数学函数 → 版本化数学计算图 → 模块化寄存器 RIR → gate / PySparQ。映射为 |x,y,s,0〉→|x,y XOR F(x),s XOR status(x),0〉，保持输入并复净私有工作区。

| 编号 | 工作项 | 状态 |
|---|---|---|
| F1 | 定义数学图、源码白名单前端、纯 helper 调用与诊断 | done |
| F2 | 补充多项式/Chebyshev 数学核、状态传播和复数分解 | done |
| F3 | 自动 out-of-place lowering、重复输入别名处理和反算 | done |
| F4 | 将普通 Roe 数学函数接入 QFVM oracle 组装 | done |
| F5 | 构造/往返/小规模门与真实 native 冒烟、说明与案例 | done |

数值采用显式定点配置。cmath 的基本函数名称是源语言识别目录，不调用 cmath 去模拟未知量子输入，也不枚举输入真值表。非多项式实函数以有配置区间和阶数的 Chebyshev 多项式合成；复数使用两条实数寄存器及明确的代数分解。第一阶段不实现 IEEE NaN/Inf/有符号零，分支切线和近似精度仍待核验。

接受赋值、算术、比较、条件表达式/结构化 if、静态有界 range 循环、tuple 返回和纯 helper。拒绝 I/O、对象突变、任意方法调用、动态循环、递归。读取函数源码并解释 AST，编译过程不执行用户函数。helper 保留独立数学图和 RIR Module；相同调用结果可复用。量子分支的状态只合并被选择路径。

RIR 本体维持 0.3；数学图是新增的生成层表示，必须能独立 JSON 往返。正确性/精度不作为本阶段的研究结论，接口、编译行为、可逆更新及后端消费需要真实验证。

## 施工结果

F1–F5 已完成当前阶段的范式交付。数学正确性继续标为 pending，不以这些测试推导全域精度。

- 数学图与前端：src/oracq/mathfunc/{graph,frontend}.py。
- 分解与自动量子模块：src/oracq/mathfunc/{numeric,lowering}.py。
- QFVM 普通公式：src/oracq/mathfunc/roe_formulas.py；原 roe_face 接口转接自动编译结果。
- 七组产物：out/math-functions/；三个描述通过真实 OriginIR 解析，自动 Roe 的 QRAM 矩阵元路径通过真实 PySparQ 执行。
- 验证：75 项核心/Schema、21 项真实后端测试通过；MIR 约束补充后的 11 项定向检查通过。

使用方式见 [纯函数编译](../manual/math-functions.md)，规范见 [MIR 0.1](../reference/math-ir.md)，证据见 [验证记录](function-compiler-validation.json)。
