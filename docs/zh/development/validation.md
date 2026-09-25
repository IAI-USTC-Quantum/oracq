# 本版验收记录

[English](../../development/validation.html) · **简体中文**

## 算法研究工作流改进（2026-09-21）

单次输入适配、求解器命名兼容、部分绑定报告与资源来源、开放成本分析及
QRAM/算术实现比较已完成。核心与 Schema 391 项测试、278 个子测试，真实
integration 41 项、六组相关数值验证 193 个案例，以及新增 9 个实现比较案例
全部通过。静态、类型和 Sphinx 检查通过。
详细范围、环境与产物路径见[工作流实施与验收记录](research-workflow-plan.md)。

## 0.8 发布验收快照

本页记录 0.8 的目录整理、算法扩展和文档构建范围。最终测试数量与产物摘要由同目录的 `validation.json` 保存。

验证分为四类：原有核心与 Schema 回归；新算法的独立数学见证；真实后端的复幅度对拍；Sphinx HTML、教程和安装包检查。

数学见证覆盖正号 QFT、所有三位输入的模加、BV 秘密恢复、Simon 零空间、振幅放大和估计、Hadamard/Swap test、单边 QAOA、Pauli 测量、周期行走、模乘求阶以及重复码恢复。已有 QLSS/QODE/QHAM 的验证边界没有因此扩大为全面数值认证。

文档构建使用 warning-as-error。API 从规范源码路径生成；历史文档不混入构建。旧导入路径另有兼容性见证。

## 验收结果

- 核心与 Schema：144 项测试、135 个 subtest 通过。
- 真实后端：29 项测试通过，其中算法展示目录包含 22 个案例，三种执行结果的最大幅度差为 `3.89e-16`。
- 文档：56 个 API 模块页面、7 篇教程；HTML 构建零 warning，10 个教程测试与 2 个搜索测试通过。
- 浏览器：首页、教程、API 导航、中英文搜索、公式和 Mermaid 图均已检查。
- 安装包：源码包与 wheel 构建通过，隔离环境中的规范导入、旧导入兼容和基础算法执行通过。

机器可读结果见 [validation.json](validation.json)。

## V1 验证覆盖进展（2026-09-10）

以下进展按 [`validation-plan.md`](validation-plan.md) §5 的 V1 阶段落地；上一节仍是 0.8 发布时点的验收快照，本节反映 0.8 之后的迭代进展而不重写历史。完整的"每个算法的证据在哪里"见 [`validation-coverage.md`](validation-coverage.md)。

- 不变量测试库：`tests/core/witness.py` 提供 `assert_unitary` / `assert_uncomputation` / `assert_bind_invariant` / `assert_block_equals` 四个断言原语，自证测试位于 `tests/core/test_witness.py`（4 个测试类，11 个用例，覆盖正确/故意错误两种路径）。
- 在途模块见证按矩阵要求补齐：
  - density：Gibbs 误差扫描 (`GibbsTests.test_error_convergence_decreases` 实测 5.6e-4 / 8.7e-5 / 8.7e-5，单调不增 + 每档 ≤ error) 与 β 网格 (`test_error_bound_uniform_in_beta`)。
  - gradient：失败概率衰减率 (`GradientTests.test_perturbed_linear_concentrates_with_grid_bits` 断言 q_{m+1} ≤ 0.34·q_m，实测比率 0.292 / 0.268)。
  - lowrank：2×2 闭式独立对拍 places=10 (`DoubleFactorizationTests.test_alpha_matches_closed_form_eigenvalues`)、DF α=1.4 ≤ Pauli α=1.5 (`test_df_alpha_tighter_than_pauli`)、THC α=1.996 独立手算 (`ThcTests.test_thc_alpha_matches_hand_computed_bound`)。
  - integration：闭式均值对拍 + heinrich_rate 校验 (`SumPreparationTests` / `QuantumSumTests` / `RateTests`)。
  - qpca：特征值读出峰对拍 + Δt 一阶误差率 (`DensityMatrixExponentiationTests` / `QpcaTests`)。
- 核心测试总数：237 → 253（`python -m unittest discover -s tests/core` 全绿）。`src/` 本轮未改动。
