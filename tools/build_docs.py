"""Build the bilingual Sphinx documentation.

English tree:  docs    -> out/docs     (the published site root)
Chinese tree:  docs/zh -> out/docs/zh

Warnings are errors (-W); both HTML and doctest builders run per language.
The Chinese tree is nested inside the English site root so that relative
../zh cross-language links resolve in the built HTML. HTML builds also drop
a .nojekyll marker and a 404 page at out/docs so the directory can be
published as a GitHub Pages artifact as-is; doctest output goes to
out/docs-doctest to keep the publishable tree clean.

The publishable tree is wiped on every HTML run so stale pages from older
layouts can never leak into a deployment; doctest output stays incremental.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TREES = {"en": ROOT / "docs", "zh": ROOT / "docs/zh"}
SITE = ROOT / "out/docs"

# Before the English tree moved to the site root it lived under /en/, so
# external pages may still deep-link there. GitHub Pages serves this file for
# any unknown path; map the legacy /en/ prefix onto the root, otherwise show
# plain not-found links anchored at the deepest known directory.
NOT_FOUND = """<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>oracq documentation — page not found</title>
    <script>
      (function () {
        var m = location.pathname.match(/^(.*?\\/en\\/?)(.*)$/);
        if (m) location.replace(
          m[1].slice(0, -3) + m[2] + location.search + location.hash);
      })();
    </script>
  </head>
  <body>
    <p>Page not found.</p>
    <p><a id="home" href="./">English documentation</a> ·
       <a id="zh" href="zh/">Chinese mirror (zh/)</a></p>
    <script>
      (function () {
        var m = location.pathname.match(/^(.*?\\/en\\/?)(.*)$/);
        var base = m ? m[1].slice(0, -3)
                     : location.pathname.replace(/[^/]*$/, "");
        document.getElementById("home").href = base;
        document.getElementById("zh").href = base + "zh/";
      })();
    </script>
  </body>
</html>
"""


def build(lang: str, builder: str) -> None:
    if builder == "html":
        out = SITE if lang == "en" else SITE / lang
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


def add_site_markers() -> None:
    (SITE / ".nojekyll").write_text("", encoding="utf-8")
    (SITE / "404.html").write_text(NOT_FOUND, encoding="utf-8")
    print(".nojekyll + 404.html written to", SITE.relative_to(ROOT))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lang", choices=["en", "zh", "all"], default="all")
    parser.add_argument("--builder", choices=["html", "doctest", "all"], default="all")
    args = parser.parse_args()
    langs = ["en", "zh"] if args.lang == "all" else [args.lang]
    builders = ["html", "doctest"] if args.builder == "all" else [args.builder]
    if "html" in builders:
        shutil.rmtree(SITE, ignore_errors=True)
    for lang in langs:
        for builder in builders:
            build(lang, builder)
    if "html" in builders:
        add_site_markers()


if __name__ == "__main__":
    main()
