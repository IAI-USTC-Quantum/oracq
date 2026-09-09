# 第二阶段实施计划

本阶段扩展实现路径和算法结构。正确性、精度、收敛性和资源优势在文档中列为待核验，不作为施工前置条件。

- CBMD 按用户指定 QST 11, 035027 (2026), DOI 10.1088/2058-9565/ae7b7e；对应 arXiv:2511.10267v3。
- QFVM 先实现 1D Euler 三守恒量的 frozen-Roe 专门化，并保留几何/数据访问接口。
- 固定点字长和舍入策略显式配置；本阶段只做接口、构造和后端运行冒烟。

| 编号 | 工作项 | 状态 |
|---|---|---|
| S0 | 审阅 PySparQ 自定义算子、QFVM/Roe 和 CBMD 输入模型 | done |
| S1 | RIR 局部工作区、PySparQ 模块级 native registry 与动态 C++ 适配 | done |
| S2 | Boolean 算术合成与 Toffoli/U3/CZ 目标门集降低 | done |
| S3 | 定点 Roe 物理量、特征结构和接口 Jacobian 的算术生成 | done |
| S4 | QRAM 数据结构、经典 Riemann 残差和局部增量更新 | done |
| S5 | 几何查询、量子矩阵元、P_Theta 与 T_L† S T_R/QFVM 组装 | done |
| S6 | Carleman、Schrodingerization、LCHS、CBMD 的开放输入模型与生成器 | done |
| S7 | 案例面板、gate/native 描述产物、实际后端冒烟和待核验清单 | done |

施工依赖为 S1、S2 支撑 S3；S3、S4 接入 S5；S6 复用现有 BE/QLSS 与新后端；S7 汇总可审阅证据。

验收要求：

- 原生自定义算子能在 PySparQ 执行，并明确不等于门级实现。
- 算术具备实际算法分解，降低后量子门只含 Toffoli、U3、CZ，QRAM 保持独立资源指令。
- QFVM 查询原始流场后计算矩阵元，经典 Riemann 只更新数据和残差结构，不预存矩阵元。
- 四类微分方程方法具有不同的结构生成器，输入 oracle 可保持开放。
- 所有新算法报告 correctness=pending，并记录当前专门化和未覆盖的参数范围。

施工结果见 [工作面板](stage2-board.md)，验收范围和待核验项见 [实施说明](stage2-implementation.md)。本计划在实现前建立，状态随施工更新。
