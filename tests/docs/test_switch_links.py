"""Cross-language raw anchors must resolve under the deployed site layout.

Built HTML sits at <site>/<lang>/<relpath>.html with both language trees as
siblings, so a raw <a href> from one tree to the other must climb out of the
language root first. This test recomputes the deployed-relative path for every
raw anchor that targets a docs-tree page and asserts the href matches —
catching depth bugs that only show up on the published site.
"""

import posixpath
import re
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DOCS = REPO / "docs"
ANCHOR = re.compile(r'<a href="([^"]+)">')
FENCE = re.compile(r"```.*?```", re.DOTALL)


def strip_code(text: str) -> str:
    """Example anchors inside fences or inline code are documentation, not
    live links; skip them."""
    text = FENCE.sub("", text)
    text = re.sub(r"``[^`]+``", "", text)
    return re.sub(r"`[^`\n]+`", "", text)


def lang_and_rel(path: Path) -> tuple[str, str]:
    parts = path.relative_to(DOCS).parts
    if parts[0] == "zh":
        return "zh", "/".join(parts[1:])
    return "en", "/".join(parts)


def iter_pages():
    for path in sorted(DOCS.rglob("*.md")):
        rel = path.relative_to(REPO).as_posix()
        if rel.startswith(("docs/archive/", "docs/_build/")) or path.name == "README.md":
            continue
        yield path


class SwitchLinkTests(unittest.TestCase):
    def test_raw_anchors_use_deployed_relative_depths(self):
        errors = []
        for path in iter_pages():
            lang, rel = lang_and_rel(path)
            dirname = posixpath.dirname(rel)
            page_dir = f"/oracq/{lang}/{dirname}" if dirname else f"/oracq/{lang}"
            for href in ANCHOR.findall(strip_code(path.read_text(encoding="utf-8"))):
                if href.startswith(("http://", "https://", "mailto:", "#")):
                    continue
                file_part, _, fragment = href.partition("#")
                deployed_target = posixpath.normpath(posixpath.join(page_dir, file_part))
                m = re.fullmatch(r"/oracq/(en|zh)/(.+)", deployed_target)
                if not m or m.group(2).startswith("../"):
                    errors.append(f"{path.relative_to(REPO)}: anchor leaves the site layout: {href}")
                    continue
                tlang, trel = m.group(1), m.group(2)
                root = DOCS / "zh" if tlang == "zh" else DOCS
                if trel.endswith(".html"):
                    trel = trel[: -len(".html")] + ".md"
                if not (root / trel).is_file():
                    errors.append(f"{path.relative_to(REPO)}: deployed target has no source twin: {href}")
                    continue
                expected = posixpath.relpath(f"/oracq/{tlang}/{m.group(2)}", page_dir) + (
                    f"#{fragment}" if fragment else ""
                )
                if href != expected:
                    errors.append(f"{path.relative_to(REPO)}: {href} (expected {expected})")
        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main()
