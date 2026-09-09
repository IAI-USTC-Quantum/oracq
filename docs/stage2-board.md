# 第二阶段工作面板

本面板的 done 表示当前阶段的范式实现与描述交付完成；算法正确性仍为 pending。范围、专门化和待核验项见 [实施说明](stage2-implementation.md)。

| 编号 | 工作项 | 状态 | 实现 / 证据 |
|---|---|---|---|
| S0 | 审阅 PySparQ 自定义算子、QFVM/Roe 和 CBMD 输入模型 | done / correctness pending | [stage2-implementation.md](../docs/stage2-implementation.md)、[backend-review.md](../docs/backend-review.md) |
| S1 | RIR 局部工作区、PySparQ 模块级 native registry 与动态 C++ 适配 | done / correctness pending | [native.py](../src/pyqecclang/native.py)、[layout.py](../src/pyqecclang/layout.py)、[test_stage2_native.py](../tests/integration/test_stage2_native.py) |
| S2 | Boolean 算术合成与 Toffoli/U3/CZ 目标门集降低 | done / correctness pending | [arithmetic.py](../src/pyqecclang/arithmetic.py)、[basis.py](../src/pyqecclang/backends/basis.py)、[toffoli_u3_cz.originir](../out/stage2/arithmetic_div/toffoli_u3_cz.originir) |
| S3 | 定点 Roe 物理量、特征结构和接口 Jacobian 的算术生成 | done / correctness pending | [roe.py](../src/pyqecclang/roe.py)、[open.rir.json](../out/stage2/roe_face/open.rir.json) |
| S4 | QRAM 数据结构、经典 Riemann 残差和局部增量更新 | done / correctness pending | [flow_data.py](../src/pyqecclang/flow_data.py)、[qfvm_patch.json](../out/stage2/qfvm_patch.json)、[test_differential.py](../tests/core/test_differential.py) |
| S5 | 几何查询、量子矩阵元、P_Theta 与 T_L† S T_R/QFVM 组装 | done / correctness pending | [qfvm.py](../src/pyqecclang/qfvm.py)、[toffoli_u3_cz.originir](../out/stage2/qfvm_costa_filter/toffoli_u3_cz.originir)、[qfvm-native-result.json](../out/stage2/qfvm-native-result.json) |
| S6 | Carleman、Schrodingerization、LCHS、CBMD 的开放输入模型与生成器 | done / correctness pending | [differential.py](../src/pyqecclang/algorithms/differential.py)、[open.rir.json](../out/stage2/carleman_qode/open.rir.json)、[toffoli_u3_cz.originir](../out/stage2/qham_cbmd_qpde/toffoli_u3_cz.originir) |
| S7 | 案例面板、gate/native 描述产物、实际后端冒烟和待核验清单 | done / correctness pending | [build_stage2.py](../tools/build_stage2.py)、[verify_stage2_backend.py](../tools/verify_stage2_backend.py)、[stage2-validation.json](../docs/stage2-validation.json) |

## 交付案例

14 类算术、Roe face、Roe entry、QRAM QFVM BE、QFVM + Costa + filtering、四类 QODE 与四类 QPDE、两条 QHAM→QPDE 路径，共 28 个。全部生成开放/闭合 RIR 和两种模块化 OriginIR 描述。

11 个代表性目标门集/QRAM 描述通过真实解析；17 个真实后端集成测试通过。完整 Roe entry 已从 QRAM 流场完成 PySparQ 模拟，使用 1211 次原生模块调用和 305 个普通事件；没有展开内部 Boolean 算术工作区。数值解精度不由这些结果推出。
