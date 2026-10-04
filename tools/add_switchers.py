"""Maintain language-switcher links on Chinese docs pages (run from repo root).

For every docs/zh/**/*.md page, place directly under the H1 title:
    <a href="...">English</a> · **简体中文**
linking to the English twin when one exists at the mirrored path under docs/,
otherwise to the English landing page. Href depths follow the DEPLOYED site
layout — the English tree builds to the site root and the Chinese tree to
<site>/zh/, so from /site/zh/manual/x.html the English twin is
../../manual/x.html. An existing switcher line is rewritten in place, making
the script idempotent.
"""

import posixpath
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ZH = REPO / "docs/zh"
DOCS = REPO / "docs"

SWITCHER = re.compile(r'(?m)^<a href="[^"]*">English</a> · \*\*简体中文\*\*[ \t]*$')

changed = 0
for path in sorted(ZH.rglob("*.md")):
    rel = path.relative_to(ZH)
    dirname = posixpath.dirname(rel.as_posix())
    page_dir = f"/oracq/zh/{dirname}" if dirname else "/oracq/zh"
    en_twin = DOCS / rel
    if en_twin.exists():
        deployed = "/oracq/" + rel.as_posix()[: -len(".md")] + ".html"
    else:
        deployed = "/oracq/index.html"
    target = posixpath.relpath(deployed, page_dir)
    line = f'<a href="{target}">English</a> · **简体中文**'
    text = path.read_text(encoding="utf-8")
    if SWITCHER.search(text):
        new = SWITCHER.sub(line, text, count=1)
    else:
        m = re.match(r"^(#\s+.*?\n)", text)
        if m:
            new = m.group(1) + "\n" + line + "\n" + text[m.end():]
        else:
            new = line + "\n\n" + text
    if new != text:
        path.write_text(new, encoding="utf-8")
        changed += 1

print("switcher written on", changed, "pages")
