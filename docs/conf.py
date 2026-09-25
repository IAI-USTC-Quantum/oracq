"""Sphinx configuration for the English documentation tree.

Only pure-Python sources are imported here; optional backends are never faked.
The Chinese mirror tree lives under ``zh/`` with its own ``zh/conf.py`` and is
excluded from this build.
"""

import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "docs/_ext"))
project = "oracq"
author = "oracq developers"
copyright = "2026, oracq developers"
release = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
language = "en"
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
html_search_language = "en"
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["archive/**", "_build/**", "README.md", "zh/**"]
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist", "fieldlist"]
myst_fence_as_directive = ["mermaid"]
myst_heading_anchors = 4
html_theme = "furo"
html_title = "oracq documentation"
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
# Built-in ``type`` in annotations (e.g. ``tuple[type, ...]``) collides with
# several ``*.type`` attribute targets; Sphinx only picks among duplicates for
# display, so this hint class is suppressed here.
suppress_warnings = ["ref.python", "myst.xref_missing"]
doctest_global_setup = "from oracq import *"
