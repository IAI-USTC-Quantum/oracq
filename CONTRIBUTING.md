# Contributing

**English** · [简体中文](CONTRIBUTING.zh-CN.md)

Install the development and documentation environments with
`uv sync --locked --extra dev --extra docs`. Write code identifiers, comments,
and user-facing strings in English; documentation is English-first with a
Chinese mirror under `docs/zh/`.

- Infrastructure lives in `src/oracq/infrastructure/`.
- Quantum algorithms live in the matching category files under `src/oracq/algorithms/`.
- Domain models, data preparation, and application-level composition live in `src/oracq/applications/`.
- Documentation uses Sphinx; full manuals and tutorials are separate, and API pages are generated from source. Historical records go into `docs/archive/`.

Run the core, case, typing, and documentation checks (including `ruff`,
`mypy src examples`, and the Sphinx builds):

```bash
uv run python tools/check_project.py --docs
```

Full native acceptance requires a real-backend interpreter:

```bash
PATH="$PWD/out/toolchain:$PATH" uv run python tools/check_project.py --docs \
  --backend-python ../QECC.Lang/.venv/bin/python
```

Changing RIR requires updating the specification, the JSON Schema, the
serializers, and the semantic tests in lockstep. Python protocols owned by an
algorithm do not require RIR changes. See
[development and acceptance](docs/development/contributing.md) and
[writing documentation](docs/development/writing-docs.md) for details.

Generated output goes into the ignored `out/` directory. Do not commit or push
without explicit authorization; commits follow Conventional Commits in English,
and virtual environments, run data, or binaries are never committed. GitHub is
the primary publishing remote; day-to-day development is synced to the Gitea
intermediate remote.
