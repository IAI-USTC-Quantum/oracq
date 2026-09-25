"""Sphinx configuration for the Chinese documentation tree.

Shares every setting with the English root configuration and overrides only
the language-specific knobs. Build with the repository root ``docs`` directory
as this tree's parent, e.g. ``sphinx-build -b html docs/zh out/docs/zh``.
"""

from pathlib import Path

_en_conf = Path(__file__).resolve().parents[1] / "conf.py"
_namespace = {"__file__": str(_en_conf)}
exec(compile(_en_conf.read_text(encoding="utf-8"), str(_en_conf), "exec"), _namespace)
globals().update(_namespace)

# Chinese overrides: Chinese UI strings, the custom zh+identifier search
# tokenizer provided by docs/_ext/search_support.py, and paths that are
# relative to this zh/ source directory.
language = "zh_CN"
html_search_language = "oracq"
html_title = "oracq 文档"
html_static_path = ["../_static"]
exclude_patterns = ["_build/**"]
