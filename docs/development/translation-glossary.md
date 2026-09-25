# Translation Glossary (Chinese → English)

**English** · <a href="../zh/development/translation-glossary.html">简体中文</a>

Single source of truth for translating this repository from Chinese to English.
Every translation pass (code comments, docstrings, docs pages, exception
messages) MUST use the exact English term on the right for the Chinese term on
the left. Append new entries as they are discovered; never introduce a second
translation for a term that already has an entry.

Punctuation: replace 、with `,`; replace full-width ：，（） with ASCII
`:`, `,`, `()`; keep the content of formulas and identifiers untouched.

## Core architecture

| Chinese | English |
|---|---|
| 寄存器级中间表示（RIR） | register-level intermediate representation (RIR) |
| 中间表示 | intermediate representation (IR) |
| 生成阶段 / 生成期 | generation stage / generation-time |
| 降低 | lowering |
| 模块 | module |
| 模块调用 | module call |
| 入口模块 | entry module |
| 依赖模块 | dependency module |
| 构造器 | builder |
| 指令 | instruction |
| 原语 | primitive |
| 寄存器 | register |
| 寄存器规格 | register specification |
| 位宽 | bit width |
| 视图 | view |
| 重解释 | reinterpretation |
| 切片 | slice / slicing |
| 下标 0 是最低位 | index 0 is the least significant bit |
| 开放 oracle / 未实现的实现 | open oracle / oracle placeholder |
| 声明 | declaration |
| 绑定 | binding (verb: bind) |
| 能力 | capability |
| 能力合取 | capability conjunction |
| 能力规格 | capability specification |
| 候选实现 | candidate implementation |
| 块编码 | block encoding |
| 信号寄存器 | signal register |
| 缩放因子 | scaling factor |
| 结构验证 | structural validation |
| 语义测试 | semantic test |
| 序列化 / 序列化器 | serialization / serializer |
| 生成层边界 | generation-layer boundary |
| 兼容层 / 兼容路径 | compatibility layer / compatibility path |
| 严格网表 | strict netlist |
| 后端导出 | backend export |
| 后端降低 | backend lowering |
| 参考执行器 | reference executor |
| 宿主读出 | host readout |
| 读出 | readout |
| 布局 | layout |
| 资源估计 | resource estimation |
| 原生实现 | native implementation |
| 注册 | registration (of native implementations) |

## Input models and algorithms

| Chinese | English |
|---|---|
| 访问模型 | access model |
| 输入模型 | input model |
| 算法契约 | algorithm contract |
| 输入协议 | input protocol |
| 验收报告 | acceptance report |
| 可检查的 | checkable |
| 算子 | operator |
| 算子包装 | operator wrapper |
| 组合 / 可组合 | composition / composable |
| 态制备 | state preparation |
| 稀疏访问 | sparse access |
| 可逆算术 | reversible arithmetic |
| 可逆线路 | reversible circuit |
| 量子数据结构 | quantum data structure |
| 密度矩阵 | density matrix |
| Gibbs 态 | Gibbs state |
| 谱 / 谱分解 | spectral / spectral decomposition |
| 低秩分解 | low-rank decomposition |
| 振幅放大 | amplitude amplification |
| 振幅估计 | amplitude estimation |
| 相位估计 | phase estimation |
| 重叠估计 | overlap estimation |
| 量子计数 | quantum counting |
| 求阶 | order finding |
| 模乘 | modular multiplication |
| 量子行走 | quantum walk |
| 变时（VTAA） | variable-time (VTAA) |
| 线性系统求解器 | linear-system solver |
| 演化 | evolution |
| 变分 | variational |
| 误差恢复 | error recovery |
| 重复码 | repetition code |
| 数据加载 | data loading |
| 求和与积分 | summation and integration |
| 梯度估计 | gradient estimation |
| 推荐系统 | recommendation system |
| 主成分分析 | principal component analysis (PCA) |
| 卷积神经网络 | convolutional neural network |
| 半定规划 | semidefinite programming (SDP) |

## Applications (QFVM / QHAM)

| Chinese | English |
|---|---|
| 流场 | flow field |
| 矩阵元 | matrix element |
| 有限闭包 | finite closure |
| 线性化 | linearization |
| 空间离散 | spatial discretization |
| 经典参考 | classical reference |
| 差分 / 差分端口 | finite-difference / stencil port |
| 推导报告 | derivation report |
| 参考工作负载 | reference workload |
| 展示目录 | gallery |
| 直连数据路径 | direct data path |
| 数据路径 | data path |
| 输入模型审阅 | input-model review |

## Validation and process

| Chinese | English |
|---|---|
| 验证 | validation |
| 数值验证 | numerical validation |
| 对拍 | cross-check (against a real backend) |
| 真实后端 | real backend |
| 模拟替身 | mock substitute |
| 验证覆盖矩阵 | validation coverage matrix |
| 适用边界 | applicability boundary |
| 已知缺口 | known gap |
| 验证方案 | validation approach |
| 实现要点 | implementation notes |
| 手册 | manual |
| 教程 | tutorial |
| 规范 | specification |
| 历史资料 / 存档 | archive material |
| 分层 | layering |
| 迁移 | migration |
| 占位 | placeholder |
| 网表 | netlist |
| 逐步讲解 | step-by-step walkthrough |

## Exception message style

Exception messages and CLI help become plain English sentences:
`同名模块定义冲突：X` → `conflicting definitions for module name: X`.
Keep the interpolated value and its f-string expression unchanged; translate
only the surrounding prose. Start lowercase only if the original starts a
sentence fragment; prefer a complete capitalized sentence.
