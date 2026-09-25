"""Build the bilingual Sphinx documentation.

English tree:  docs    -> out/docs/en
Chinese tree:  docs/zh -> out/docs/zh

Warnings are errors (-W); both HTML and doctest builders run per language.
The two output roots are siblings so that relative ../en|../.. cross-language
links resolve in the built HTML.
"""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TREES = {"en": ROOT / "docs", "zh": ROOT / "docs/zh"}


def build(lang: str, builder: str) -> None:
    if builder == "html":
        out = ROOT / "out/docs" / lang
    else:
        out = ROOT / "out/docs" / f"{lang}-{builder}"
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


if __name__ == "__main__":
    main()
