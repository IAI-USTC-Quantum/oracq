# 安装与环境

<a href="../../manual/installation.html">English</a> · **简体中文**

语言核心要求 Python 3.11 或更高版本，运行期仅依赖 PyYAML 做文本序列化。量子模拟器和文档工具单独安装。可选 `pyqsp` extra（`uv sync --extra pyqsp`）安装 pyqsp，用于可替换的 QSP 相位合成适配器。

在仓库目录中建立开发环境：

```bash
uv sync --locked --extra dev --extra docs
uv run python examples/algorithm_gallery.py
```

第二条命令会生成一组小型算法的 RIR、OriginIR-ext 和验证摘要，输出位于 `out/algorithm-gallery/`。展示案例的使用与修改见[运行与修改算法展示目录](../tutorials/gallery.md)。

如果只需要从源码调用核心接口，可以将 `src` 加入 Python 路径：

```bash
PYTHONPATH=src python examples/algorithm_gallery.py
```

## 可选后端

执行真实后端测试需要另一个已安装 `pysparq` 和 `uniqc` 的 Python 环境。本仓库不会在安装核心包时下载或编译它们。执行算术原生算子还需要可用的 C++17 编译器。

```bash
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

`--native` 会比较参考执行器、PySparQ 和 OriginIR 后端的完整复幅度。它不只是检查导出文本能否解析。环境就绪后的常用命令行用法见[命令行](cli.md)。

## 构建文档

```bash
uv run python tools/build_docs.py --lang all
```

该命令以 warning 视为错误的方式运行两种语言树的 HTML 与 doctest 构建：英文 HTML 输出到 `out/docs`（即发布站点根目录），中文 HTML 输出到 `out/docs/zh`，doctest 产物在 `out/docs-doctest/{en,zh}`。打开 `out/docs/index.html` 即可浏览英文站点，中文站点从 `out/docs/zh/index.html` 进入。HTML 构建失败或教程断言失败都会返回非零退出码。构建配置与写作规范（含交叉链接约定）见[编写文档](../development/writing-docs.md)。
