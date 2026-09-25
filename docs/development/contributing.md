# Development and acceptance

**English** · [简体中文](../zh/development/contributing.html)

## Placing new code

RIR, validation, serialization, execution, and backend code lives in
`infrastructure`. Quantum generators go into the matching category file under
`algorithms`; shared algorithm interfaces and reporting tools also belong to
the algorithm library. Physical models, data preparation, and application-level
composition live in `applications`.

Do not keep adding new algorithm branches to the `elementary` or
`differential` compatibility entry points. They exist only for legacy imports
and are not locations for new implementations.

## Adding an algorithm

A public generator must document its input model, register layout, generation
parameters, outputs, and scope of applicability. The output must be a plain
{obj}`Operation <oracq.infrastructure.builder.Operation>` or an explicit
oracle wrapper class; Python callbacks may only execute during the generation
stage.

Provide the smallest independently computable witness first, then verify how
the backend interprets the same circuit. Do not backfill quantum outputs with
classical reference results. State preparation, success probability, and
physical scales should each be checked separately.

## Running the checks

```bash
uv sync --locked --extra dev --extra docs
uv run python tools/check_project.py --docs
PATH="$PWD/out/toolchain:$PATH" uv run python tools/check_project.py --docs \
  --backend-python ../QECC.Lang/.venv/bin/python
uv build --out-dir out/release
```

The first check command needs no external quantum simulator. The second adds
real-backend tests; it fails when the backend is missing — skips are not
accepted in place of acceptance. The CI configuration runs the core,
documentation, and build checks; full native acceptance requires the
corresponding environment.

Generated results, environments, and build files all go into the ignored
`out/` directory or existing build directories. Commits follow Conventional
Commits; never commit or push without user authorization.
