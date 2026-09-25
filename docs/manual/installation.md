# Installation and environment

**English** · <a href="../zh/manual/installation.html">简体中文</a>

The language core requires Python 3.11 or newer and has no third-party runtime
dependencies. Quantum simulators and documentation tooling are installed
separately.

Set up a development environment in the repository directory:

```bash
uv sync --locked --extra dev --extra docs
uv run python examples/algorithm_gallery.py
```

The second command generates the RIR, OriginIR-ext, and validation summaries
of a set of small algorithms; the output lands in `out/algorithm-gallery/`.
For using and modifying the showcase examples see
[Running and modifying the algorithm gallery](../tutorials/gallery.md).

If you only need to call the core interfaces from source, add `src` to your
Python path:

```bash
PYTHONPATH=src python examples/algorithm_gallery.py
```

## Optional backends

Running real-backend tests requires a separate Python environment with
`pysparq` and `uniqc` installed. This repository does not download or compile
them when the core package is installed. Executing arithmetic native operators
additionally requires a working C++17 compiler.

```bash
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

`--native` compares the full complex amplitudes of the reference executor,
PySparQ, and the OriginIR backend. It does more than check that the exported
text can be parsed. For everyday command-line usage once the environment is
ready see [Command line](cli.md).

## Building the documentation

```bash
uv run sphinx-build -W --keep-going -b html docs out/docs/html
uv run sphinx-build -W --keep-going -b doctest docs out/docs/doctest
```

Open `out/docs/html/index.html` to browse this site. A failed HTML build or a
failed tutorial assertion both return a non-zero exit code. For build
configuration and writing conventions (including the cross-linking rules) see
[Writing documentation](../development/writing-docs.md).
