# 编写文档

文档分为完整文档和教程。完整文档解释稳定的规则、接口和限制；教程围绕一个具体任务，给出输入、代码、结果解释和下一步。API 参考由源码生成。

## 写作要求

先说明读者要解决的问题，再解释为什么采用这些步骤。使用完整句子和具体对象名称，避免让读者从开发日志里推断当前行为。术语第一次出现时给出含义；保留必要的数学定义，但不要用缩写代替解释。

当前行为直接写成当前规则。历史设计、阶段计划和已放弃的路线放入 `archive`。对未完成的算法，准确说明缺少哪一部分，例如相位求解、统计读出或数值误差验证。

## 可执行教程

MyST 的 `testcode` 块会由 Sphinx doctest builder 执行。优先用断言检查稳定的数学结果，避免依赖随机采样或格式不稳定的打印输出。

```bash
uv run sphinx-build -W --keep-going -b doctest docs out/docs/doctest
```

## API 页面

公开模块按分类列在 `docs/api/`。新增模块后运行：

```bash
uv run python tools/generate_api_docs.py
uv run sphinx-build -W --keep-going -b html docs out/docs/html
```

API 通过 autodoc 导入实际源码，不使用 mock 导入。可选后端必须继续在执行入口导入，以便在核心环境中构建文档。完整使用方式参见 [Sphinx autodoc](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html)。

## 中文与 API 名称搜索

本站在 `docs/_ext/search_support.py` 中维护一致的索引与查询分词：中文使用字和二元片段，Python 标识符保留完整名称及下划线分段。该扩展避免所用 Sphinx 版本的中文 stemmer 脚本不匹配问题，不修改第三方安装目录。搜索规则改动后应同时运行 `tests/docs` 和浏览器查询检查。

## 算法页面

`docs/manual/algorithms/` 下每个算法一页，由 `index.md` 的 glob toctree 自动收录；写新页不需要改 `index.md`。文件名用 kebab-case（如 `qsvt-matrix-inversion.md`），与入口函数或算法英文名对应。

页面首行下方用一行引用注明类别与所属模块（`> 类别 Cn · 模块 pyqecclang.algorithms.<module> · 阶段 Vn`），类别与阶段取值必须与 `validation-coverage.md` 一致。正文固定六节：

1. **概述**：问题陈述、数学定义（可用 MyST dollarmath）、文献依据。文献只写 `docs/manual` 现有内容或源码 docstring 明确引用的，不许编造。
2. **接口与输入模型**：入口函数签名（以源码为准）、input model 类型（词汇见 `algorithm-coverage.md`）、返回对象属性表。
3. **实现要点**：生成策略、寄存器布局、设计决策与适用边界；未实现的部分准确说明缺什么。
4. **验证方案**：类别与判定准则（引 `validation-plan.md` §2）、三层证据位置（`tests/core/<file>:<TestClass.test_method>`）、见证技术与实测口径（容差、实测值）。本节事实必须与 `validation-coverage.md` 一致。
5. **已知缺口与计划阶段**：与 `validation-coverage.md` 缺口列一致。
6. **相关链接**：源码模块、API 参考页（`docs/api/algorithms/` 下，先确认实际文件名再链接）、`../../development/validation-coverage.md`。

维护约定：新增算法必须同时新增本页，并同步 `validation-coverage.md` 与 `algorithm-coverage.md`；三者描述的类别、阶段、缺口与证据位置必须一致。
