"""源码分类与 API 入口的维护约束。"""

import ast
import importlib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


class RepositoryLayoutTests(unittest.TestCase):
    def test_root_only_has_public_facade_cli_and_compatibility(self):
        files = {p.name for p in (ROOT / "src/pyqecclang").glob("*.py")}
        self.assertEqual(files, {"__init__.py", "__main__.py", "_compat.py"})

    def test_public_algorithm_modules_have_api_pages(self):
        for path in (ROOT / "src/pyqecclang/algorithms").glob("*.py"):
            if path.stem.startswith("_") or path.stem == "legacy":
                continue
            page = ROOT / "docs/api/algorithms" / (path.stem + ".rst")
            self.assertTrue(page.exists(), str(page))
            self.assertIn("automodule:: pyqecclang.algorithms." + path.stem, page.read_text())

    def test_internal_imports_do_not_use_compatibility_paths(self):
        from pyqecclang._compat import ALIASES, SPLIT_EXPORTS

        old = set(ALIASES) | set(SPLIT_EXPORTS)
        for path in (ROOT / "src/pyqecclang").rglob("*.py"):
            if path.name == "_compat.py":
                continue
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node, ast.ImportFrom) and node.module:
                    self.assertNotIn(node.module, old, str(path))

    def test_old_package_attributes_reference_canonical_modules(self):
        package = importlib.import_module("pyqecclang")
        self.assertIs(package.ir, importlib.import_module("pyqecclang.infrastructure.ir"))
        qham = importlib.import_module("pyqecclang.qham")
        self.assertIs(qham.quantum, importlib.import_module("pyqecclang.algorithms.qham"))
