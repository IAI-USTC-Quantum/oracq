"""Prepend language-switcher links to Chinese docs pages (run from repo root).

For every docs/zh/**/*.md page, insert directly under the H1 title:
    [English](<relative .html path to the EN twin>) · **简体中文**
when an English twin exists at the mirrored path under docs/, otherwise link
to the English landing page. Existing switcher lines are skipped.
"""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ZH = REPO / "docs/zh"
DOCS = REPO / "docs"

changed = 0
skipped = 0
for path in sorted(ZH.rglob("*.md")):
    rel = path.relative_to(ZH)
    depth = len(rel.parts) - 1
    up = "../" * (depth + 1)
    en_twin = DOCS / rel
    if en_twin.exists():
        target = f"{up}{rel.with_suffix('.html').as_posix()}"
    else:
        target = f"{up}index.html"
    line = f'<a href="{target}">English</a> · **简体中文**'
    text = path.read_text(encoding="utf-8")
    if text.startswith("[English]"):
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
