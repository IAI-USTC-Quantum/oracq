"""README files sit outside the Sphinx build, so nothing else catches a
relative link or anchor that rots. This test resolves every relative link in
README.md and README.zh-CN.md against the repository root and, for Markdown
targets, checks the anchor against the GitHub-style slugs of the file's
headings."""

import re
import unittest
from pathlib import Path
from urllib.parse import unquote

REPO = Path(__file__).resolve().parents[2]
LINK = re.compile(r"\[[^\]]+\]\(([^)\s]+)\)")
FENCE = re.compile(r"```.*?```", re.DOTALL)
HEADING = re.compile(r"^#+\s+(.*?)\s*#*\s*$", re.MULTILINE)


def slug(heading: str) -> str:
    """GitHub-style anchor: lowercase, drop non-word characters except
    hyphens, spaces become hyphens (CJK characters are kept)."""
    text = heading.strip().lower()
    text = re.sub(r"[^\w\s-]", "", text, flags=re.UNICODE)
    return re.sub(r"\s", "-", text)


def headings_of(text: str) -> set[str]:
    body = FENCE.sub("", text)
    return {slug(match) for match in HEADING.findall(body)}


class ReadmeLinkTests(unittest.TestCase):
    def test_relative_links_and_anchors_resolve(self) -> None:
        problems: list[str] = []
        for name in ("README.md", "README.zh-CN.md"):
            readme = REPO / name
            text = FENCE.sub("", readme.read_text(encoding="utf-8"))
            own_anchors = headings_of(readme.read_text(encoding="utf-8"))
            for target in LINK.findall(text):
                if target.startswith(("http://", "https://", "mailto:")):
                    continue
                fragment = ""
                if "#" in target:
                    target, fragment = target.split("#", 1)
                    fragment = unquote(fragment)
                if not target:
                    if fragment and fragment not in own_anchors:
                        problems.append(f"{name}: anchor #{fragment} has no heading")
                    continue
                resolved = REPO / target
                if not resolved.exists():
                    problems.append(f"{name}: {target} does not exist")
                    continue
                if fragment and resolved.suffix == ".md":
                    anchors = headings_of(resolved.read_text(encoding="utf-8"))
                    if fragment not in anchors:
                        problems.append(f"{name}: {target}#{fragment} has no matching heading")
        self.assertEqual(problems, [])


if __name__ == "__main__":
    unittest.main()
