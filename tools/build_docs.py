"""Build the bilingual Sphinx documentation.

English tree:  docs    -> out/docs/en
Chinese tree:  docs/zh -> out/docs/zh

Warnings are errors (-W); both HTML and doctest builders run per language.
The two HTML output roots are siblings so that relative ../en|zh cross-language
links resolve in the built HTML. HTML builds also drop a landing redirect and
a .nojekyll marker at out/docs so the directory can be published as a GitHub
Pages artifact as-is; doctest output goes to out/docs-doctest to keep the
publishable tree clean.
"""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TREES = {"en": ROOT / "docs", "zh": ROOT / "docs/zh"}

LANDING = """<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>oracq documentation</title>
    <meta http-equiv="refresh" content="0; url=en/index.html">
    <link rel="canonical" href="en/index.html">
  </head>
  <body>
    <p>Redirecting to the <a href="en/index.html">English documentation</a>.
    The Chinese mirror is at <a href="zh/index.html">zh/index.html</a>.</p>
  </body>
</html>
"""


def build(lang: str, builder: str) -> None:
    if builder == "html":
        out = ROOT / "out/docs" / lang
    else:
        out = ROOT / "out/docs-doctest" / lang
    command = [
        sys.executable,
        "-m",
        "sphinx",
        "-W",
        "--keep-going",
        "-b",
        builder,
        str(TREES[lang]),
        str(out),
    ]
    print("Building", lang, builder, "->", out.relative_to(ROOT), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def add_landing() -> None:
    site = ROOT / "out/docs"
    (site / "index.html").write_text(LANDING, encoding="utf-8")
    (site / ".nojekyll").write_text("", encoding="utf-8")
    print("Landing redirect + .nojekyll written to", site.relative_to(ROOT))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lang", choices=["en", "zh", "all"], default="all")
    parser.add_argument("--builder", choices=["html", "doctest", "all"], default="all")
    args = parser.parse_args()
    langs = ["en", "zh"] if args.lang == "all" else [args.lang]
    builders = ["html", "doctest"] if args.builder == "all" else [args.builder]
    for lang in langs:
        for builder in builders:
            build(lang, builder)
    if "html" in builders:
        add_landing()


if __name__ == "__main__":
    main()
