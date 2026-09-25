# oracq Development Conventions

This repository is a standalone Python package. All work happens here; do not modify the neighboring UnifiedQuantum, QRAM-Simulator, or QECC.Lang projects for compatibility testing.

RIR is the architectural center. When changing the IR, update docs/reference/rir.md, the JSON Schema, the serializers, and the semantic tests in lockstep. Register names, bit widths, and views must survive until backend lowering. Module calls and Repeat structures must not be unconditionally expanded during generation or text serialization.

The core does not depend on quantum backends. Backend export and execution are separated; native dependencies are imported only at execution entry points. Code identifiers, comments, docstrings, and user-facing strings (exception messages, CLI help, printed output) are written in English. Documentation is English-first with a maintained Chinese mirror under docs/zh/. Commit messages are written in English.

The core validation command is `python -m unittest discover -s tests/core -v`. Real-backend validation runs tests/integration with an interpreter that has uniqc and pysparq installed; never substitute mocks or skips for a real cross-check. Static checks use `ruff check src tests examples`; type checking uses `mypy src examples`.

Do not commit virtual environments, caches, build artifacts, or simulated output. Do not run git commit or push without an explicit request.

The source tree is layered into infrastructure, algorithms, and applications. New algorithms go into the matching category file, not into legacy import-compatibility layers. Documentation uses Sphinx with manual, tutorials, reference, and api sections per language tree; historical material lives in archive. Documentation changes run both the HTML and doctest builders, and warnings are treated as errors.
