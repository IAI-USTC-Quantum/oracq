"""搜索索引的中文词片段与 Python 标识符保持可检索。"""

import importlib.util
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "search_support", Path(__file__).resolve().parents[2] / "docs/_ext/search_support.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SearchTests(unittest.TestCase):
    def test_chinese_phrase_tokens_are_in_sentence_index(self):
        query = set(module.split_terms("振幅估计"))
        document = set(module.split_terms("使用振幅估计算法"))
        self.assertTrue(query <= document)

    def test_python_identifier_and_parts(self):
        terms = module.split_terms("amplitude_estimation")
        self.assertIn("amplitude_estimation", terms)
        self.assertIn("amplitude", terms)
        self.assertIn("estimation", terms)
        self.assertEqual(module.split_terms(""), [])
