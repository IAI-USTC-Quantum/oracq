"""中英文与 Python 标识符的搜索分词，两端使用相同规则。"""

import re

from sphinx.search import SearchLanguage

TOKEN = re.compile(r"[A-Za-z0-9_]+|[\u3400-\u9fff]+")


def split_terms(text):
    terms = []
    for word in TOKEN.findall(text.lower()):
        if "\u3400" <= word[0] <= "\u9fff":
            terms.extend(word)
            terms.extend(word[i : i + 2] for i in range(len(word) - 1))
        else:
            terms.append(word)
            if "_" in word:
                terms.extend(part for part in word.split("_") if part)
    return list(dict.fromkeys(terms))


class MixedSearchLanguage(SearchLanguage):
    lang = "oracq"
    language_name = "Chinese and Python identifiers"
    js_stemmer_code = (
        "var Stemmer = function () { this.stemWord = function (w) { return w.toLowerCase(); }; };"
    )
    js_splitter_code = r"""
    function splitQuery(text) {
      const terms = [];
      for (const word of (text.toLowerCase().match(/[A-Za-z0-9_]+|[\u3400-\u9fff]+/g) || [])) {
        if (word[0] >= '\u3400' && word[0] <= '\u9fff') {
          terms.push(...word);
          for (let i = 0; i < word.length - 1; i++) terms.push(word.slice(i, i + 2));
        } else {
          terms.push(word);
          if (word.includes('_')) terms.push(...word.split('_').filter(Boolean));
        }
      }
      return [...new Set(terms)];
    }
    """

    def split(self, input):
        return split_terms(input)

    def stem(self, word):
        return word.lower()


def setup(app):
    app.add_search_language(MixedSearchLanguage)
    return {"version": "1", "parallel_read_safe": True, "parallel_write_safe": True}
