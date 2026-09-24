# 开发与验收

## 放置新代码

RIR、验证、序列化、执行和后端代码放在 `infrastructure`。量子生成器放在 `algorithms` 的对应类别文件；共享的算法接口和报告工具也属于算法库。物理模型、数据准备与应用级组合放在 `applications`。

不要把新的算法分支继续加入 `elementary` 或 `differential` 兼容入口。它们只用于旧导入，不是新的实现位置。

## 添加算法

公开生成器需要说明输入模型、寄存器布局、生成参数、输出和适用范围。输出必须是普通 {obj}`Operation <oracq.infrastructure.builder.Operation>` 或明确的 oracle 包装类；Python 回调只能在生成阶段执行。

先给出可以独立计算的最小见证，再验证后端对同一线路的解释。不要用经典参考结果回填量子输出。状态制备、成功概率和物理尺度应分别检查。

## 运行检查

```bash
uv sync --locked --extra dev --extra docs
uv run python tools/check_project.py --docs
PATH="$PWD/out/toolchain:$PATH" uv run python tools/check_project.py --docs \
  --backend-python ../QECC.Lang/.venv/bin/python
uv build --out-dir out/release
```

第一条检查命令不需要外部量子模拟器。第二条增加真实后端测试，后端缺失时会失败，不以 skip 代替验收。CI 配置执行核心、文档与构建检查；完整原生验收需要相应环境。

生成结果、环境和构建文件都写入被忽略的 `out/` 或现有构建目录。提交使用 Conventional Commits；没有用户授权不执行 commit 或 push。
