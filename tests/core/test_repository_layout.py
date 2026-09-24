"""源码分类与 API 入口的维护约束。"""

import ast
import importlib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


class RepositoryLayoutTests(unittest.TestCase):
    def test_root_only_has_public_facade_cli_and_compatibility(self):
        files = {p.name for p in (ROOT / "src/oracq").glob("*.py")}
        self.assertEqual(files, {"__init__.py", "__main__.py", "_compat.py"})

    SUBPACKAGES = {
        "input_model",
        "common",
        "qlss",
        "qnlss",
        "qode",
        "qpde",
        "qml",
        "optimization",
        "basics",
        "qec",
    }

    def test_algorithms_flat_level_only_has_categories(self):
        base = ROOT / "src/oracq/algorithms"
        files = {p.name for p in base.glob("*.py")}
        self.assertEqual(files, {"__init__.py"})
        dirs = {p.name for p in base.iterdir() if p.is_dir() and p.name != "__pycache__"}
        self.assertEqual(dirs, self.SUBPACKAGES)

    def test_public_algorithm_modules_have_api_pages(self):
        base = ROOT / "src/oracq/algorithms"
        for path in sorted(base.rglob("*.py")):
            if path.name == "__init__.py" or path.stem.startswith("_") or path.stem == "legacy":
                continue
            module = "oracq.algorithms." + ".".join(path.relative_to(base).with_suffix("").parts)
            page = ROOT / "docs/api/algorithms" / path.relative_to(base).with_suffix("").parent / (path.stem + ".rst")
            self.assertTrue(page.exists(), str(page))
            self.assertIn("automodule:: " + module, page.read_text())

    def test_internal_imports_do_not_use_compatibility_paths(self):
        from oracq._compat import ALIASES, SPLIT_EXPORTS

        old = set(ALIASES) | set(SPLIT_EXPORTS)
        for path in (ROOT / "src/oracq").rglob("*.py"):
            if path.name == "_compat.py":
                continue
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node, ast.ImportFrom) and node.module:
                    self.assertNotIn(node.module, old, str(path))

    def test_old_package_attributes_reference_canonical_modules(self):
        package = importlib.import_module("oracq")
        self.assertIs(package.ir, importlib.import_module("oracq.infrastructure.ir"))
        qham = importlib.import_module("oracq.qham")
        self.assertIs(qham.quantum, importlib.import_module("oracq.algorithms.input_model.qham"))
