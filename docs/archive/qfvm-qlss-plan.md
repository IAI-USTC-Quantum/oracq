# QFVM / QLSS 输入模型审查计划

| 编号 | 项目 | 状态 |
|---|---|---|
| Q1 | 对照 QFVM III/IV、CKS 1.1/4 和 Costa 输入及 general walk | done |
| Q2 | 分离 LinearSystem、SparseAccess、BEInput 与 QLSSProtocol，记录缩放/物理子空间/范数恢复 | done |
| Q3 | 为 QFVM 提供真正的原地位置 oracle 和任意 row/column 元素 oracle；明确填充空间 | done |
| Q4 | 稀疏输入适配、Costa 参数换算、CKS Chebyshev 基础路线与替换见证 | done |
| Q5 | 修正残差态制备的逐层 QRAM 查询和局部数据更新实现 | done |
| Q6 | 独立说明文档、针对性测试与后端描述 | done |

当前确认的问题：QFVM 强制 BE 输入；几何 XOR 查询不等于 CKS 原地位置 oracle；alpha 未进入 Costa 有效逆谱界；补齐的变量存在零模式；归一化解与物理更新幅值的契约缺失；旧 QRAM 角树制备按所有 prefix 枚举查询。

本轮修正接口与可审阅的实现链。CKS 路线限于论文的 Chebyshev/LCU 基础构造，不声称完成 VTAA 或获得最优复杂度；算法全域精度和 CFD 正确性继续单列待核验。

本轮 Q1–Q6 已完成输入模型审查与范式修正。独立结论见 [QFVM / QLSS 审查文档](../manual/qfvm.md)，验证范围见 [验收记录](qfvm-qlss-validation.json)。数值求解器等价、VTAA、Costa 内核认证和 CFD 闭环不在完成声明之内。
