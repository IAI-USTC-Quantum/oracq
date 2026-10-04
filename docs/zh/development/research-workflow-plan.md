# 算法研究工作流改进

<a href="../../development/research-workflow-plan.html">English</a> · **简体中文**

核心受众是根据论文实现、验证和比较量子科学计算算法的研究者。
实现顺序按依赖推进，每阶段同步规范、示例与验证结果。

| 阶段 | 交付内容 | 状态 |
|---|---|---|
| 1 | 单次输入适配、结构化诊断、求解器与契约命名兼容迁移 | 已实现，回归通过 |
| 2 | 部分绑定契约、结构化绑定报告、资源与能力回归 | 已实现，回归通过 |
| 3 | 紧凑旋转计数、QRAM 实参分账、开放 oracle 调用台账 | 已实现，回归通过 |
| 4 | 实现选择案例、Roe/QFVM 研究报告、可追溯验证 | 9 个新案例四路对拍及 193 个领域案例通过 |
| 5 | 算法作者教程、论文主线与证据更新、完整验收 | 已完成 |

已复现的问题：一次生成重复调用提供方方法；开放估计仍报闭合错误；
旋转计数随 Repeat 物化；QRAM 查询按内部形参而非入口实参合并。

RIR 保留具体位宽与组装常量；开放的是模块实现。改变公开布局或 alpha
需要显式适配或重新生成。Python 提供方不进入序列化结果；分析报告与 RIR
分别保存。未知 oracle 成本与私有工作区不能作为零成本。

核心验收使用 unittest、Schema、ruff 与 mypy；语义执行使用真实 uniqc/pysparq。
Sphinx HTML 与 doctest 以 warning 为错误，论文用 LaTeX 构建核对引用和图表。
生成产物只写入 out；不修改相邻后端仓库，不提交或推送。

## 原生验证环境

本轮使用 `out/native-validation/bin/python`，通过 uv 创建和安装。
PySparQ 的 PyPI 0.1.1 尚无原生 RIR 接口，因此从 QRAM-Simulator 的
已提交版本 `4e4c9f16fe3e912a828915231ecb4a6663048bd2` 提取快照到本仓
`out/native-sources/`，用 uv 构建为 `0.0.0+g4e4c9f16`。
UniQC 对应发行包为 `unified-quantum==0.1.1`，C++ 模拟器为
`uniqc-cppsimulator==1.0.1`。不加载相邻仓库的未提交改动。

新工作流案例在 width=2,3,4、repetitions=3 下比较门表、QRAM 表和算术：
四条实际执行路径相对独立参考的最大幅度误差为 `2.5e-16`。
详细报告位于 `out/research-workflow-native/report.json`。

QMatrix 原生测试按执行预算分层：10 位角字行制备的公开寄存器占 14 位，
地址计算私有工作区占 16 位，总计 30 位；保留两个真实 PySparQ 路径对拍。
新增 4 位角字案例，总计 24 位，用参考执行器、两条 PySparQ 路径与 OriginIR
完成四路对拍。此分层经用户确认；没有使用模拟替身或 skip。

## 验收结果

- `tools/check_project.py --docs`：全部检查通过；核心与 Schema 共 391 项测试、278 个子测试通过。
- 真实后端 `tests/integration`：41 项通过；原 10 位 QMatrix 案例保留。
- `tools/run_verification.py`：oracles、blockencoding、arithmetic、mathfunc、qham_qfvm、ode 六组共 193 个案例通过。该结果不替代历史十三组报告的版本范围。
- 新实现比较案例：9 个配置、四条执行路径全部通过，最大幅度误差 `2.5e-16`。
- ruff、mypy、Sphinx HTML 与 55 项 doctest 通过；论文生成于 `out/paper/main.pdf`。

工作流检查记录在 `out/checks/research-workflow/report.json`；原生执行日志为
`out/research-native-integration.log` 与 `out/research-native-verification.log`。
后者生成的六组 JSON 报告包含相同源码指纹、验证脚本指纹和后端版本。
论文主稿为 `paper/main.tex`；旧的 Overleaf 副本未作同步，不作为本次交付入口。
