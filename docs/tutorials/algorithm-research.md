# From a paper's access model to implementation comparison

**English** · <a href="../../zh/tutorials/algorithm-research.html">简体中文</a>

This tutorial is for researchers implementing new quantum algorithms: first save the algorithm with its open oracles, then choose an implementation, and finally compare against an independent mathematical reference and compare costs. The core of this needs no quantum backend.

## Save the open program and analyze the calls

```{doctest}
>>> from oracq import bind_with_report, dumps, loads, estimate_resources
>>> from oracq.applications.oracle_study import oracle_study
>>> opened, implementations = oracle_study(width=2, repetitions=3)
>>> restored = loads(dumps(opened))
>>> cost = estimate_resources(restored, require_closed=False)
>>> cost.complete, cost.qubits
(False, None)
>>> sorted((call.adjoint, count) for call, count in cost.oracle_calls.items())
[(False, 3), (True, 3)]
```

{obj}`oracle_study <oracq.applications.oracle_study.oracle_study>` returns the open program and all candidate implementations at once; the round trip through {obj}`dumps <oracq.infrastructure.serialization.dumps>` and {obj}`loads <oracq.infrastructure.serialization.loads>` shows that open descriptions save as usual, and {obj}`estimate_resources <oracq.infrastructure.estimate.estimate_resources>` returns a cost ledger even while the program is unclosed. The access model is `data ^= (address + 1) mod 2**width`. The same angle-word database is used to build the diagonal block encoding; every call computes the angle word, rotates, and uncomputes. Repeating three times stores only a Repeat and module calls; the ledger above shows that each implementation must provide three forward and three adjoint calls.

## Choose an implementation and keep the binding report

```{doctest}
>>> implementation, memory = implementations["qram_table"]
>>> linked = bind_with_report(restored, {"AngleWord": implementation})
>>> linked.report.ok
True
>>> closed = linked.require()
>>> estimate_resources(closed).qram_queries
Counter({'angle_words': 6})
```

{obj}`bind_with_report <oracq.infrastructure.linking.bind_with_report>` returns the program together with a {obj}`BindingReport <oracq.infrastructure.linking.BindingReport>`; `require()` hands over the closed program only when the report shows no problems. Change the selection to `gate_table` or `arithmetic` to compare the other two implementations.
The gate table enumerates the same integer function, the QRAM table supplies it at run time, and reversible arithmetic computes it directly and clears the private
work word. All three expose the same signature; their private ancillas and resource consumption may differ. Changing width requires
regenerating the program; that is not an interface-preserving implementation binding.

## Validation and resource comparison

```bash
python examples/research_workflow.py
PYTHONPATH=src /path/to/backend/python examples/research_workflow.py \
  --native --widths 2 3 -o out/research-workflow-native
```

The independent reference uses trigonometric formulas: with angle word k at each address, after r repetitions the success and failure amplitudes
are `cos(r*pi*k/2**width)/sqrt(2**width)` and the corresponding sin.
This is a property of this example's rotation dilation and must not be generalized into a rule for taking arbitrary block encodings to a power directly.
The report covers the amplitude error over all output branches, the success probability, program and data fingerprints, the binding report, and
both the open and closed costs. Native validation additionally calls both PySparQ paths and the OriginIR backend.

This example validates all address branches of the uniform input; a pass on a single input state is not promoted to a proof of equivalence
for arbitrary unitary operators. The interface and workspace cleanup of reversible integer arithmetic are covered separately by core
regression tests.

A more complete scientific-computing chain reuses `tools/build_qlss_comparison.py`: the raw flow-field QRAM,
the compiled Roe formulas, sparse access, and two QLSS generators. That script provides open/closed costs and
data and program fingerprints; numerical validation continues to use `tests/verification/verify_qham_qfvm.py`
and `verify_mathfunc.py`, reporting implementation error and method error respectively.

## Related pages

- Manual: [Resource estimation (Toffoli+Clifford+T+QRAM)](../manual/resource-estimation.md) (cost-ledger accounting), [Conventions owned by algorithms](../manual/contracts.md)
- Specification: [Open IR](../reference/open-ir.md) (open declarations and batched binding)
- Algorithm pages: <a href="../../zh/manual/algorithms/xor-database.html">XOR database views</a> (the view basis shared by the gate-table and QRAM implementations)
- API reference: [Oracle implementation comparison](../api/applications/oracle_study.rst), [Resource estimation](../api/infrastructure/estimate.rst), [Binding and capability analysis](../api/infrastructure/linking.rst)
- Continue with: [Scientific computing workflows](scientific-workflows.md) (the complete chains of QHAM/Carleman/LCHS/CBMD)
