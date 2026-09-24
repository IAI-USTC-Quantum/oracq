# 参与开发

使用 `uv sync --locked --extra dev --extra docs` 安装开发与文档环境。代码标识符使用英文，文档与注释使用中文。

- 基础设施放在 `src/oracq/infrastructure/`。
- 量子算法放在 `src/oracq/algorithms/` 的相应类别文件。
- 领域模型、数据准备和应用级组合放在 `src/oracq/applications/`。
- 文档使用 Sphinx；完整文档与教程分开，API 从源码生成。历史记录放入 `docs/archive/`。

运行核心、案例、类型与文档检查（含 `ruff`、`mypy src examples` 与 Sphinx 构建）：

```bash
uv run python tools/check_project.py --docs
```

完整原生验收需要真实后端解释器：

```bash
PATH="$PWD/out/toolchain:$PATH" uv run python tools/check_project.py --docs \
  --backend-python ../QECC.Lang/.venv/bin/python
```

修改 RIR 必须同步规范、Schema、序列化和语义测试。算法自己的 Python 协议不需要修改 RIR。详情见[开发与验收](docs/development/contributing.md)和[文档写作](docs/development/writing-docs.md)。

生成结果写入被忽略的 `out/`。没有明确授权不 commit/push；提交采用 Conventional Commits，不提交虚拟环境、运行数据或二进制。
