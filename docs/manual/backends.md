# Export and execution backends

**English** · [简体中文](../zh/manual/backends.html)

oracq separates export from execution.
{obj}`export_originir <oracq.infrastructure.backends.originir.export_originir>`
and
{obj}`export_toffoli_u3_cz <oracq.infrastructure.backends.basis.export_toffoli_u3_cz>`
only produce descriptions and do not require a quantum simulator to be
installed. API reference pages for each backend: [OriginIR-ext](../api/infrastructure/backends/originir.rst),
[strict gate set basis](../api/infrastructure/backends/basis.rst),
[PySparQ](../api/infrastructure/backends/pysparq.rst), [reference executor](../api/infrastructure/execution.rst),
and [readout](../api/infrastructure/readout.rst).

## OriginIR-ext

The exporter preserves `DEF`, module calls, and `QRAMDECL`. Register
operations are lowered to concrete gates at export time;
{obj}`Repeat <oracq.infrastructure.ir.Repeat>` is expressed through a reusable
auxiliary definition. QRAM memory is passed in as an external resource and is
not embedded in the RIR.

`export_toffoli_u3_cz` lowers ordinary gates to Toffoli, U3, and CZ; QRAM
remains a standalone resource instruction. It does not expand the device
network of a physical QRAM.

{obj}`run_originir <oracq.infrastructure.backends.originir.run_originir>`
executes with UnifiedQuantum. That backend's parser expands `DEF`s, so there
is an expansion budget before execution; this does not change the saved
modular RIR.

## PySparQ

{obj}`run_pysparq <oracq.infrastructure.backends.pysparq.run_pysparq>`
executes modules as register events.
{obj}`NativeRegistry <oracq.infrastructure.native.NativeRegistry>` can invoke
custom native operators at module boundaries, which suits large reversible
arithmetic. Native implementations and gate-level bodies are managed
separately: an open oracle with only a native implementation can be used in
the corresponding simulation, but this does not grant gate-level export
capability.

Dynamic operators require a matching PySparQ ABI and a C++17 compiler. For
concrete interface behavior and the reviewed versions see
[Backend compatibility notes](../reference/backend-compatibility.md).

## Reference executor

{obj}`simulate <oracq.infrastructure.execution.simulate>` is a small register
executor that depends only on the standard library. It is convenient for
checking bit order, phases, XOR semantics, and ancilla states. The number of
states, the number of steps, and small-amplitude truncation are all bounded;
do not treat the scale it can run as a real-hardware resource estimate.

## Readout

The RIR core contains no measurement or reset. Final measurement,
post-selection, and statistical processing happen at the host layer; when
dynamic OriginIR is needed, the readout adapter can be used explicitly.
Error-recovery circuits with syndromes ([repetition codes](../zh/manual/algorithms/repetition-codes.html))
keep the error information, and appropriate host processing is required
before these registers are reused.
