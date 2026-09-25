"""Structural validation of documentation cross-links.

Because ``suppress_warnings = ["ref.python", "myst.xref_missing"]`` in conf.py
lets broken py-domain references and cross-language links degrade silently to
plain text under ``-W`` builds, this suite scans the Markdown sources
directly: ``{obj}`` targets must be importable, relative links and
literalinclude directives must point at existing files, relative-link anchors
must exist in the target page's headings (MyST does not warn on missing
anchors), algorithm entry pages must link symmetrically, and each page's
"source" line must point at an existing path. Both language trees are
validated; ``.html`` targets are cross-language switcher links checked against
their ``.md`` twins.
"""

import importlib
import re
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DOCS = REPO / "docs"

ROLE_OBJ = re.compile(
    r"\{(?:obj|class|func|mod|data|exc|meth|attr)\}`(?:[^<`]*?<)?([A-Za-z_][\w.]*)`"
)
DOWNLOAD = re.compile(r"\{download\}`[^<`]*?<([^<`\s]+)>`")
MD_LINK = re.compile(r"\]\(([^)\s]+)\)")
LITERALINCLUDE = re.compile(r"^(`{0,3}\{literalinclude\}\s+)(\S+)", re.MULTILINE)
FENCE = re.compile(r"```.*?```", re.DOTALL)
SRC_PATH = re.compile(r"^- (?:源码：|Source:)\s.*?`((?:src|tools|examples|tests)/[^`]+)`", re.MULTILINE)
HEADING = re.compile(r"^(#{1,4})\s+(.+?)\s*#*\s*$", re.MULTILINE)
SLUG_CLEAN = re.compile(r"[^\w\u4e00-\u9fff\-]")
EXTERNAL = ("http://", "https://", "mailto:")


def iter_doc_files():
    for path in sorted(DOCS.rglob("*.md")):
        rel = path.relative_to(REPO).as_posix()
        if rel.startswith(("docs/archive/", "docs/_build/")) or path.name == "README.md":
            continue
        yield path


def strip_fences(text: str) -> str:
    return FENCE.sub("", text)


def strip_inline_code(text: str) -> str:
    """Link syntax inside inline code is not rendered as a link; skip it when
    checking markdown links."""
    text = re.sub(r"``[^`]+``", "", text)
    return re.sub(r"`[^`\n]+`", "", text)


def resolve_object(target: str) -> str | None:
    """Return an error message; None when the target resolves."""
    parts = target.split(".")
    for split in range(len(parts), 0, -1):
        try:
            obj = importlib.import_module(".".join(parts[:split]))
        except ImportError:
            continue
        try:
            for attr in parts[split:]:
                obj = getattr(obj, attr)
        except AttributeError:
            return f"{target}: {'.'.join(parts[:split])} has no {target}"
        return None
    return f"{target}: no importable module prefix"


def heading_slug(title: str) -> str:
    """Reproduce MyST default_slugify: lowercase, spaces to hyphens, drop
    everything outside word characters, hyphens, and CJK."""
    title = re.sub(r"\[[^\]]*\]\([^)]*\)", "", title)  # link text stays out of the slug
    title = re.sub(r"\$[^$]+\$", "", title)  # math stays out of the slug
    title = re.sub(r"\*\*([^*]+)\*\*", r"\1", title)
    title = re.sub(r"\*([^*]+)\*", r"\1", title)
    title = re.sub(r"`([^`]+)`", r"\1", title)
    return SLUG_CLEAN.sub("", title.lower().replace(" ", "-"))


def page_slugs(path: Path) -> set[str]:
    """Anchors of the target page's level 1-4 headings; duplicates get -N
    appended per MyST rules."""
    slugs: list[str] = []
    for match in HEADING.finditer(path.read_text(encoding="utf-8")):
        slug = heading_slug(match.group(2))
        if slug in slugs:
            base, i = slug, 1
            while f"{base}-{i}" in slugs:
                i += 1
            slug = f"{base}-{i}"
        slugs.append(slug)
    return set(slugs)


def disk_target(link_file: str) -> str:
    """Cross-language switcher links use .html targets in the built site;
    validate them against their .md source twins."""
    return link_file[:-5] + ".md" if link_file.endswith(".html") else link_file


class XrefTargetTests(unittest.TestCase):
    def test_obj_role_targets_resolve(self):
        errors = []
        for path in iter_doc_files():
            text = strip_fences(path.read_text(encoding="utf-8"))
            for target in ROLE_OBJ.findall(text):
                problem = resolve_object(target)
                if problem:
                    errors.append(f"{path.relative_to(REPO)}: {problem}")
        self.assertEqual(errors, [])

    def test_relative_link_targets_exist(self):
        errors = []
        for path in iter_doc_files():
            raw = path.read_text(encoding="utf-8")
            targets = list(DOWNLOAD.findall(raw))
            targets += [m.group(2) for m in LITERALINCLUDE.finditer(raw)]
            targets += MD_LINK.findall(strip_inline_code(strip_fences(raw)))
            for target in targets:
                if target.startswith(EXTERNAL + ("#",)) or not target:
                    continue
                location = (path.parent / disk_target(target.split("#")[0])).resolve()
                if not location.exists():
                    errors.append(f"{path.relative_to(REPO)}: {target}")
        self.assertEqual(errors, [])

    def test_anchor_targets_exist(self):
        """Anchors of relative links must exist in the target .md page's
        heading set (rst targets and cross-language .html links excluded —
        headings differ across languages)."""
        errors = []
        cache: dict[str, set[str]] = {}
        for path in iter_doc_files():
            text = strip_inline_code(strip_fences(path.read_text(encoding="utf-8")))
            for target in MD_LINK.findall(text):
                if "#" not in target or target.startswith(EXTERNAL):
                    continue
                file_part, _, anchor = target.partition("#")
                if not anchor or file_part.endswith((".rst", ".json", ".html")):
                    continue
                target_path = path.parent / disk_target(file_part) if file_part else path
                if not target_path.exists():
                    continue  # file existence is covered by test_relative_link_targets_exist
                key = str(target_path)
                if key not in cache:
                    cache[key] = page_slugs(target_path)
                if anchor not in cache[key]:
                    errors.append(f"{path.relative_to(REPO)}: {target}")
        self.assertEqual(errors, [])

    def test_algorithm_page_links_are_symmetric(self):
        for pages in (DOCS / "manual" / "algorithms", DOCS / "zh" / "manual" / "algorithms"):
            edges = {}
            for path in sorted(pages.glob("*.md")):
                text = strip_fences(path.read_text(encoding="utf-8"))
                edges[path.stem] = {
                    match[:-3]
                    for match in re.findall(r"\]\(([a-z0-9-]+\.md)(?:#[^)]*)?\)", text)
                }
            errors = []
            for source, targets in edges.items():
                for target in targets:
                    if target == source or not (pages / f"{target}.md").exists():
                        continue
                    if source not in edges.get(target, {}):
                        errors.append(f"{pages.relative_to(REPO)}: {target}.md lacks a link back to {source}.md")
            self.assertEqual(errors, [])

    def test_algorithm_page_source_lines_exist(self):
        errors = []
        for path in iter_doc_files():
            for ref in SRC_PATH.findall(strip_fences(path.read_text(encoding="utf-8"))):
                if not (REPO / ref).exists():
                    errors.append(f"{path.relative_to(REPO)}: {ref}")
        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main()
