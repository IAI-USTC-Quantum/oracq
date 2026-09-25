"""Prepend language-switcher links to Chinese docs pages (run from repo root).

For every docs/zh/**/*.md page, insert directly under the H1 title:
    <a href="...">English</a> · **简体中文**
linking to the English twin when one exists at the mirrored path under docs/,
otherwise to the English landing page. Href depths follow the DEPLOYED site
layout — a built page sits at <site>/<lang>/<relpath>.html, so from
/site/zh/manual/x.html the English twin is ../../en/manual/x.html (one level
more than source-tree intuition). Existing switcher lines are skipped.
"""

import posixpath
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ZH = REPO / "docs/zh"
DOCS = REPO / "docs"

changed = 0
skipped = 0
for path in sorted(ZH.rglob("*.md")):
    rel = path.relative_to(ZH)
    dirname = posixpath.dirname(rel.as_posix())
    page_dir = f"/oracq/zh/{dirname}" if dirname else "/oracq/zh"
    en_twin = DOCS / rel
    if en_twin.exists():
        deployed = "/oracq/en/" + rel.as_posix()[: -len(".md")] + ".html"
    else:
        deployed = "/oracq/en/index.html"
    target = posixpath.relpath(deployed, page_dir)
    line = f'<a href="{target}">English</a> · **简体中文**'
    text = path.read_text(encoding="utf-8")
    if text.startswith('<a href='):
        skipped += 1
        continue
    m = re.match(r"^(#\s+.*?\n)", text)
    if m:
        text = m.group(1) + "\n" + line + "\n" + text[m.end():]
    else:
        text = line + "\n\n" + text
    path.write_text(text, encoding="utf-8")
    changed += 1

print("switcher added to", changed, "pages;", skipped, "already had one")
