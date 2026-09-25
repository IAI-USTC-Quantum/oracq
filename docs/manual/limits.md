# Applicability boundaries and validation status

**English** · [简体中文](../zh/manual/limits.html)

oracq can express open algorithms, assemble modules, and generate executable
small instances. Whether it suits a specific problem also depends on the
input model, the scale, the numerical approximation, and the readout
requirements.

| Layer | What has been validated | What cannot be inferred from it |
|---|---|---|
| [RIR](../reference/rir.md) | Serialization, layout, aliasing, control guards, module calls, and open bindings | Mathematical correctness of arbitrary quantum algorithms |
| [Basic algorithms](../zh/manual/algorithms/index.html) | Small-scale matrix/probability witnesses cross-checked against real-backend complex amplitudes | Large-scale running cost or quantum advantage |
| [Arithmetic and math functions](math-functions.md) | Reversible constructions under finite bit patterns plus some numerical witnesses | Uniform error bounds over all word lengths and approximation intervals |
| [Hamiltonian methods](../zh/manual/algorithms/hamiltonian-simulation.html) | Pauli evolution, Trotter composition, and interface checks | General QSP kernels and automatic error configuration |
| [QLSS/QODE/QPDE](differential-equations.md) | Input adaptation, modular generation, and partial real execution | Full solver accuracy, success channels, and convergence guarantees |
| [QHAM](qham.md) | Algebraic closure at finite HAM truncation and small executions | Automatic convergence certification against the original PDE |

## Scale limits

A single register or fused view is at most 64 bits. Some algorithm
convenience interfaces place the target or signal in one register and
therefore hit this limit first; the RIR itself can use multiple registers.

The [reference executor](backends.md#reference-executor) and the backend
adapters have step-count, state-count, and QRAM materialization budgets.
Downstream OriginIR parsers expand modules. That a program can keep a compact
RIR does not mean it can also run cheaply in the current simulators.

## Mathematical claims

Claims such as `Hermitian`, dissipation, sparsity, element bounds, spectral
bounds, and ancilla re-cleaning are handled by the corresponding algorithms.
The [Python protocols](contracts.md) can check that interfaces exist and that
return shapes match, but cannot prove the mathematical truth of a claim.

`eps` is not a language-core type. Users should keep the generation
configuration and validate, in the concrete application, the errors of spatial
discretization, numerical representation, algorithmic truncation, and
post-selected readout.

## Results and readout

An output-state oracle may carry a success signal. The normalization
direction, the success probability, and the physical norm of the original
vector are different pieces of information. After measurement or
post-selection, the full classical field values cannot be recovered from the
target amplitudes alone.

The interfaces, witness coverage, and documentation build results of the new
algorithms are recorded in the [validation record of this
release](../development/validation.md).
