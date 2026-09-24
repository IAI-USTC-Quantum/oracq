"""文档交叉链接的结构校验。

conf.py 的 ``suppress_warnings = ["ref.python"]`` 会让失效的 py-domain 引用在
``-W`` 构建中静默降级为纯文本，因此这里直接扫描 Markdown 源文件校验：
``{obj}`` 目标必须可导入，相对链接与 literalinclude 必须指向存在的文件，
相对链接的标题锚点必须存在于目标页（MyST 不对缺失锚点报警），算法词条页
之间的链接必须对称，算法页「源码」行必须指向存在的路径。
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
LITERALINCLUDE = re.compile(r"^\{literalinclude\}\s+(\S+)", re.MULTILINE)
FENCE = re.compile(r"```.*?```", re.DOTALL)
SRC_PATH = re.compile(r"^- 源码：.*?`((?:src|tools|examples|tests)/[^`]+)`", re.MULTILINE)
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
    """内联代码中的链接语法不会被渲染成链接，检查 md 链接时应忽略。"""
    text = re.sub(r"``[^`]+``", "", text)
    return re.sub(r"`[^`\n]+`", "", text)


def resolve_object(target: str) -> str | None:
    """返回错误信息；目标可解析时返回 None。"""
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
            return f"{target}: {'.'.join(parts[:split])} 没有 {target}"
        return None
    return f"{target}: 找不到可导入的模块前缀"


def heading_slug(title: str) -> str:
    """复现 MyST default_slugify：小写、空格转连字符、剔除字母数字连字符以外字符。"""
    title = re.sub(r"\[[^\]]*\]\([^)]*\)", "", title)  # 链接文本不进入 slug
    title = re.sub(r"\$[^$]+\$", "", title)  # 数学式不进入 slug
    title = re.sub(r"\*\*([^*]+)\*\*", r"\1", title)
    title = re.sub(r"\*([^*]+)\*", r"\1", title)
    title = re.sub(r"`([^`]+)`", r"\1", title)
    return SLUG_CLEAN.sub("", title.lower().replace(" ", "-"))


def page_slugs(path: Path) -> set[str]:
    """目标页 1-4 级标题的锚点集合；重复标题按 MyST 规则追加 -N。"""
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
            targets += LITERALINCLUDE.findall(raw)
            targets += MD_LINK.findall(strip_inline_code(strip_fences(raw)))
            for target in targets:
                if target.startswith(EXTERNAL + ("#",)) or not target:
                    continue
                location = (path.parent / target.split("#")[0]).resolve()
                if not location.exists():
                    errors.append(f"{path.relative_to(REPO)}: {target}")
        self.assertEqual(errors, [])

    def test_anchor_targets_exist(self):
        """相对链接的 #锚点 必须存在于目标 .md 页的标题集合（rst 目标除外）。"""
        errors = []
        cache: dict[str, set[str]] = {}
        for path in iter_doc_files():
            text = strip_inline_code(strip_fences(path.read_text(encoding="utf-8")))
            for target in MD_LINK.findall(text):
                if "#" not in target or target.startswith(EXTERNAL):
                    continue
                file_part, _, anchor = target.partition("#")
                if not anchor or file_part.endswith((".rst", ".json")):
                    continue
                target_path = path.parent / file_part if file_part else path
                if not target_path.exists():
                    continue  # 文件存在性由 test_relative_link_targets_exist 负责
                key = str(target_path)
                if key not in cache:
                    cache[key] = page_slugs(target_path)
                if anchor not in cache[key]:
                    errors.append(f"{path.relative_to(REPO)}: {target}")
        self.assertEqual(errors, [])

    def test_algorithm_page_links_are_symmetric(self):
        pages = DOCS / "manual" / "algorithms"
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
                if source not in edges[target]:
                    errors.append(f"{target}.md 缺少指回 {source}.md 的链接")
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
