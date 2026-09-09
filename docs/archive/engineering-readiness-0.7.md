# 工程成熟度与发布边界

0.7 是算法约定与工程加固版本。它把输入接口、能力检查、错误报告、配置验证和可重复检查变成明确的开发流程。**不能把它描述为所有量子求解算法已经达到生产数值正确性的版本。** 本轮用户的阶段目标仍是可靠表达、组装与检查算法范式。

| 层次 | 本轮交付 | 当前边界 |
|---|---|---|
| RIR 与操作 | 版本化序列化、模块调用、开放主体、结构校验 | 仍是 RIR 0.3；不引入算法专用类型 |
| 算法输入 | Python 结构协议、requires、可聚合 check 报告、单位操作的多种视图 | 方法存在不证明数学含义；具体适配器验证返回接口 |
| 算法输出 | StateOracleProtocol / BlockEncodingProtocol，可检查 provides | 不将后选择态冒充干净制备 |
| 配置 | Costa/CKS、积分/轮廓、辅助网格、Carleman 的输入检查 | 阶数、耗散声明、零填充不等于误差证明 |
| 实现替换 | 模块 ABI/alpha/能力约束；保留 gate/QRAM 绑定 | 未知数据更新可能使数学声明失效，应用负责更新证据 |
| Hamiltonian | 实际 Pauli 项演化、Trotter 列表组合、可注入 QSP 接口 | 通用 QSP-HamSim 内核未提供，缺少时明确拒绝 |
| 发布工程 | 集中检查工具、Gitea CI 配置、源码/wheel 构建与隔离安装 | CI 尚需远端实际执行；本地通过不声称 CI 已绿 |
| 量子求解器 | Costa/CKS、LCHS、Schrödingerization、Carleman、QHAM 组装 | 完整数值精度、成功通道、收敛和优势仍实验性 |

## 可重复检查

```bash
uv sync --locked --extra dev
uv run python tools/check_project.py

PATH="$PWD/out/toolchain:$PATH" uv run python tools/check_project.py \
  --backend-python ../QECC.Lang/.venv/bin/python

uv build --out-dir out/release
```

检查工具失败即返回非零退出码，报告位于 `out/checks/latest/report.json`。不提供后端解释器只代表核心检查，不表示跳过了“完整验收所要求的真实后端”。原生测试不使用 mock/skip 代替真实 API。

本轮验证结果以 [contract-hardening-validation.json](contract-hardening-validation.json) 为记录；对应环境的后端修订仍在 [backend-revisions.json](../../backend-revisions.json)。新增的单项 Trotter 与 unitary-LCU 测试比较参考执行器、真实 PySparQ 与真实 OriginIR 的复幅度。

## 0.6 → 0.7 迁移

原有 BE/state-prep/sparse 包装类、`bind`、`linear_qode(...)(G,initial,t)` 和 QLSS 入口保留。新增 `block_encoding()`、`state_preparation()`、`sparse_access()` 等宿主方法；普通 Operation 可以直接用于支持相应协议的输入位置。

`A.type` 仅用于阅读描述，不能写成唯一类型分派。推荐检查 Python 协议；自定义对象不用继承库类，也不用修改语言。

新应用使用 `QODEProtocol.solve(QODEProblem(...),t)` 明确数学声明；旧三参数形式仍由调用方承担前提。QLSS/QODE 现在会更早拒绝缺接口、宽度不匹配、缺受控/逆能力或错误的内核输出。先前通过的非法 bool/NaN/非整数阶数配置可能因此报错。

`annotate` 不再无意把 `supports_adjoint=False` 或 `supports_controlled=False` 改回 True。能力仍由 RIR 依赖图保守推导。state-prep 库生成器明确标注零输入与干净工作区的承诺；这是算法库元数据，不是新增的 RIR 验证规则。

RIR 格式和 Schema 没有变更；原有 0.1/0.2/0.3 描述继续按既有读取规则工作。新的属性可能改变生成的模块名和描述字节，测试应验证语义/结构，不把哈希名字当作稳定 API。

## 下一阶段的验收对象

选具体问题后，逐项验证其真实 input model、实际矩阵/态、算法输入前提、成功子空间、读出尺度以及截断误差。一个对象满足 Python Protocol 并不能替代这些检查。推荐从小型非对易 Hamiltonian、带符号/复数矩阵、真实 QFVM 条目和可解析 PDE 开始，使用独立数学参考验证，而不是让参考结果回填量子输出。
