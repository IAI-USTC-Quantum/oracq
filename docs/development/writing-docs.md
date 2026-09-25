# Writing documentation

**English** · [简体中文](../zh/development/writing-docs.html)

The documentation is bilingual. The English pages under
`docs/{manual,tutorials,reference,development}/` are the primary source, and a
Chinese mirror with the identical layout lives under `docs/zh/`: every page
keeps the same relative path in both trees, and a page edited on one side must
be synchronized on the other. The API pages are generated into both trees
(`docs/api` and `docs/zh/api`); the algorithm manual pages are currently
Chinese-only (68 pages under `docs/zh/manual/algorithms/`) with an English
index placeholder at `manual/algorithms/index.md`.

Documentation splits into manuals and tutorials. Manuals explain stable rules,
interfaces, and limitations; a tutorial centers on one concrete task and gives
the input, the code, an interpretation of the results, and the next step. The
API reference is generated from the source.

## Writing requirements

State the problem the reader is trying to solve first, then explain why these
steps are taken. Use complete sentences and concrete object names; do not make
the reader infer current behavior from development logs. Give the meaning of a
term at its first appearance; keep the necessary mathematical definitions, but
do not use abbreviations in place of explanations.

Write current behavior directly as current rules. Historical designs, stage
plans, and abandoned routes go into `archive`. For unfinished algorithms,
state exactly which part is missing, for example phase solving, statistical
readout, or numerical error validation.

## Language switcher

Every page carries a language-switcher line directly under its H1: English
pages use `**English** · [简体中文](<relative path>)`, Chinese pages use
`[English](<relative path>) · **简体中文**` (the Chinese side is prepended
automatically by `tools/add_switchers.py`). Cross-tree targets always use the
`.html` form — for example
`[简体中文](../zh/development/contributing.html)` from this directory — so
the link resolves in the built HTML of the other tree.

## Executable tutorials

MyST `testcode` blocks are executed by the Sphinx doctest builder. Prefer
assertions that check stable mathematical results and avoid depending on
random sampling or print output with unstable formatting. When the output is
deterministic, follow the `testcode` block with a `testoutput` block showing
the executed result; the build verifies it verbatim, and long output can be
elided with `...` (ELLIPSIS). Tutorial examples should be paired with their
executed output wherever possible, so readers see the expected result before
running anything.

```bash
uv run python tools/build_docs.py --lang all --builder doctest
```

Both language trees build through `tools/build_docs.py --lang all`: the
English tree `docs` builds to `out/docs/en` and the Chinese tree `docs/zh`
builds to `out/docs/zh`, each with warnings treated as errors.

## Cross-linking

When running text first mentions a public API object, link it to the API
reference with the MyST `{obj}` role, writing the target as the
fully-qualified path (consistent with the `:obj:` targets in
`docs/api/toplevel.rst`), for example
`` {obj}`Builder <oracq.infrastructure.builder.Builder>` ``. Do not use short
root-package paths (such as `oracq.Builder`); link only the first occurrence
of an object; never add links inside code blocks.

Links between pages keep using relative markdown links: tutorials and manuals
point at API pages with real paths such as
`[quantum linear systems](../api/algorithms/qlss/qlss.rst)`, and at manual
pages with `[core concepts](../manual/concepts.md)`; `#heading-anchor` deep
links work (`myst_heading_anchors = 4`). The same-family / same-group links
in an algorithm page's "Related links" section must be symmetric: if A links
to B, B must also link to A. When a tutorial demonstrates the algorithm, add a
tutorial back-link line to "Related links" in this format:

```markdown
- Tutorial: [swap the QODE method for the same linear problem](../../tutorials/differential-equations.md)
```

Tutorial pages end with a fixed "Related pages" section listing the related
manual chapters, reference pages, algorithm pages, and API pages.

Fine-grained link density: core concepts and API names must carry a link at
their first mention in running text; leaving them bare over the long term is
not allowed. API objects use `{obj}`; concepts, chapters, and specification
content use relative links, adding a heading anchor when pointing at a
specific subsection (e.g.
`[registers and views](../manual/concepts.md#registers-and-views)`). On an
algorithm page, the citation line links the module path to the corresponding
API page; under the signature block in "Interface and input model", an
"API entries:" line lists the `{obj}` entry points; the "Concepts:" line in
"Related links" points at the manual concept pages the algorithm critically
depends on. Hub pages (concepts, operators, qdata, qmem, contracts) and the
specification pages (reference) must be back-linked by their consumer pages,
so they do not become orphans. New pages add links under the same standard —
do not link only to the section index.

All `{obj}` targets and relative links are validated by
`tests/docs/test_xrefs.py` (target importable / file exists / anchor present
among the target page's headings); `suppress_warnings = ["ref.python"]` in the
build configuration keeps dead py-domain references from failing the `-W`
build, so that test is the safety net.

## API pages

Public modules are listed by category under `docs/api/` (Chinese titles under
`docs/zh/api/`). After adding a module, run:

```bash
uv run python tools/generate_api_docs.py --lang all
uv run python tools/build_docs.py --lang all
```

The generator also produces `docs/api/toplevel.rst` (and its Chinese twin
`docs/zh/api/toplevel.rst`) from the root package's `oracq.__all__` — the
package-overview page that links to each module page grouped by defining
module; after adding a root export name or adjusting `__all__`, rerun the
generator the same way. The `TITLES` and `TITLES_ZH` tables inside the
generator maintain each module's page title per language; remember to add an
entry for a new module. Legacy import paths are kept for compatibility only
and never appear in the API documentation.

The API pages import the actual source through autodoc, with no mock imports.
Optional backends must keep being imported at the execution entry points, so
the documentation can be built in the core environment. See
[Sphinx autodoc](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html)
for full usage.

## Chinese and API-name search

Both language trees maintain consistent indexing and query tokenization in
`docs/_ext/search_support.py`: Chinese uses characters and bigram segments,
while Python identifiers keep full names plus underscore-separated segments.
The extension works around the Chinese stemmer-script mismatch of the Sphinx
version in use and does not modify third-party installation directories.
After changing the search rules, run both `tests/docs` and browser query
checks.

## Algorithm pages

Each algorithm gets one page under `docs/zh/manual/algorithms/`, collected
automatically by the glob toctree in that directory's `index.md`; writing a
new page requires no change to `index.md`. File names are kebab-case (e.g.
`qsvt-matrix-inversion.md`), matching the entry function or the algorithm's
English name. These pages are currently Chinese-only (68 pages); the English
tree carries a placeholder index at `docs/manual/algorithms/index.md` that
links into the Chinese catalog until the translations land.

Directly under the page's first line, a one-line blockquote records the
category and owning module
(`> 类别 Cn · 模块 oracq.algorithms.<module> · 阶段 Vn`; in English
`> Category Cn · module oracq.algorithms.<module> · stage Vn`); the category
and stage values must agree with `validation-coverage.md`. The body has six
fixed sections:

1. **Overview**: problem statement, mathematical definitions (MyST dollarmath
   allowed), literature basis. Cite only literature already present in
   `docs/manual` or explicitly referenced by source docstrings — never invent
   references.
2. **Interface and input model**: entry-function signature (as in the source),
   the input model type (vocabulary in `algorithm-coverage.md`), and a table
   of the returned object's attributes.
3. **Implementation notes**: generation strategy, register layout, design
   decisions and applicability boundary; for unimplemented parts, state
   exactly what is missing.
4. **Validation approach**: category and acceptance criteria (citing
   `validation-plan.md` §2), the locations of the three evidence layers
   (`tests/core/<file>:<TestClass.test_method>`), the witness technique and
   the measured figures (tolerances, measured values). The facts in this
   section must agree with `validation-coverage.md`.
5. **Known gaps and planned stage**: consistent with the gap column of
   `validation-coverage.md`.
6. **Related links**: source module, API reference page (under
   `docs/api/algorithms/`, confirm the actual file name before linking), and
   `../../development/validation-coverage.md`.

Maintenance agreement: adding an algorithm must add its page at the same time
and synchronize `validation-coverage.md` and `algorithm-coverage.md`; the
category, stage, gaps, and evidence locations described by the three must
agree.
