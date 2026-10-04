"""Cross-language raw anchors must resolve under the deployed site layout.

The English tree builds to the site root (<site>/<relpath>.html) and the
Chinese tree nests under it (<site>/zh/<relpath>.html). This test recomputes
the deployed-relative path for every raw anchor that targets a docs-tree page
and asserts the href matches — catching depth bugs that only show up on the
published site.
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


def page_dir_of(lang: str, rel: str) -> str:
    dirname = posixpath.dirname(rel)
    if lang == "zh":
        return f"/oracq/zh/{dirname}" if dirname else "/oracq/zh"
    return f"/oracq/{dirname}" if dirname else "/oracq"


def deployed_target_lang(target: str) -> tuple[str, str] | None:
    """Map a normalized site-absolute target to (lang, relpath), or None when
    it does not land inside either published tree."""
    if target.startswith("/oracq/zh/"):
        return "zh", target[len("/oracq/zh/"):]
    if target.startswith("/oracq/"):
        return "en", target[len("/oracq/"):]
    return None


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
            page_dir = page_dir_of(lang, rel)
            for href in ANCHOR.findall(strip_code(path.read_text(encoding="utf-8"))):
                if href.startswith(("http://", "https://", "mailto:", "#")):
                    continue
                file_part, _, fragment = href.partition("#")
                target = posixpath.normpath(posixpath.join(page_dir, file_part))
                located = deployed_target_lang(target)
                if located is None or located[1].startswith("../"):
                    errors.append(f"{path.relative_to(REPO)}: anchor leaves the site layout: {href}")
                    continue
                tlang, trel = located
                root = DOCS / "zh" if tlang == "zh" else DOCS
                if trel.endswith(".html"):
                    trel = trel[: -len(".html")] + ".md"
                if not (root / trel).is_file():
                    errors.append(f"{path.relative_to(REPO)}: deployed target has no source twin: {href}")
                    continue
                if trel.endswith(".md"):
                    trel = trel[: -len(".md")] + ".html"
                deployed = f"/oracq/zh/{trel}" if tlang == "zh" else f"/oracq/{trel}"
                expected = posixpath.relpath(deployed, page_dir) + (
                    f"#{fragment}" if fragment else ""
                )
                if href != expected:
                    errors.append(f"{path.relative_to(REPO)}: {href} (expected {expected})")
        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main()
