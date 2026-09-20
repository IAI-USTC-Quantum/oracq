"""Sphinx 配置；仅导入纯 Python 源码，不伪造可选后端。"""

import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "docs/_ext"))
project = "pyqecclang"
author = "pyqecclang developers"
copyright = "2026, pyqecclang developers"
release = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
language = "zh_CN"
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.doctest",
    "sphinxcontrib.mermaid",
    "search_support",
]
html_search_language = "pyqecclang"
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["archive/**", "_build/**", "README.md"]
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist", "fieldlist"]
myst_fence_as_directive = ["mermaid"]
myst_heading_anchors = 4
html_theme = "furo"
html_title = "pyqecclang 文档"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "light_css_variables": {"color-brand-primary": "#2458a6", "color-brand-content": "#2458a6"},
    "dark_css_variables": {"color-brand-primary": "#89b8ff", "color-brand-content": "#89b8ff"},
}
autodoc_typehints = "none"
autodoc_member_order = "bysource"
autodoc_default_options = {"members": True, "undoc-members": True, "show-inheritance": True}
autosummary_generate = True
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_use_ivar = True
# 注解中的内建 type（如 ``tuple[type, ...]``）会与多个 ``*.type`` 属性目标产生
# 歧义引用提示；Sphinx 仅在多个目标间做展示选择，这里抑制该类别提示。
suppress_warnings = ["ref.python"]
doctest_global_setup = "from pyqecclang import *"
