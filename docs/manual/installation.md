# 安装与环境

语言核心要求 Python 3.11 或更高版本，没有第三方运行时依赖。量子模拟器和文档工具单独安装。

在仓库目录中建立开发环境：

```bash
uv sync --locked --extra dev --extra docs
uv run python examples/algorithm_gallery.py
```

第二条命令会生成一组小型算法的 RIR、OriginIR-ext 和验证摘要，输出位于 `out/algorithm-gallery/`。

如果只需要从源码调用核心接口，可以将 `src` 加入 Python 路径：

```bash
PYTHONPATH=src python examples/algorithm_gallery.py
```

## 可选后端

执行真实后端测试需要另一个已安装 `pysparq` 和 `uniqc` 的 Python 环境。本仓库不会在安装核心包时下载或编译它们。执行算术原生算子还需要可用的 C++17 编译器。

```bash
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

`--native` 会比较参考执行器、PySparQ 和 OriginIR 后端的完整复幅度。它不只是检查导出文本能否解析。

## 构建文档

```bash
uv run sphinx-build -W --keep-going -b html docs out/docs/html
uv run sphinx-build -W --keep-going -b doctest docs out/docs/doctest
```

打开 `out/docs/html/index.html` 即可浏览本站。HTML 构建失败或教程断言失败都会返回非零退出码。
