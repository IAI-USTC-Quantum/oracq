# oracq

**English** · [简体中文](zh/index.html)

oracq organizes quantum algorithms in Python and compiles them to a saveable,
composable register-level intermediate representation (RIR). Algorithms may
depend on oracle placeholders that are not implemented yet; once concrete
implementations are supplied, the same description exports to OriginIR-ext or
runs on PySparQ (see [backends and export](manual/backends.md) for details).

If this is your first encounter with the project, start with the
[core concepts](manual/concepts.md) and then follow a tutorial. If you are
implementing or reviewing algorithms, the full manuals explain interfaces,
input preconditions, composition rules, and backend limitations. The API
reference is generated directly from the current source.

```{toctree}
:maxdepth: 1
:caption: Manuals

manual/index
reference/index
api/index
```

```{toctree}
:maxdepth: 2
:caption: Tutorials

tutorials/index
```

```{toctree}
:maxdepth: 2
:caption: Development

development/index
```

The project currently provides assembly paths ranging from oracle queries,
search, estimation, and Hamiltonian evolution to QLSS, QODE, QFVM, and QHAM.
Maturity varies per implementation; check
[applicability and validation status](manual/limits.md) before choosing an
algorithm. The Chinese mirror of this site lives under
[zh/](zh/index.html); algorithm manual pages are currently Chinese-only there
while their English translations are in progress.
