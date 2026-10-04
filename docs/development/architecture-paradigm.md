# Architecture paradigm and interface conventions

**English** · <a href="../zh/development/architecture-paradigm.html">简体中文</a>

This document is the public, citable statement of oracq's development
paradigm: the layered structure, the stability contract of RIR, the three
boundaries of variation that every public API must respect, and the roadmap
for scientific interface metadata. It records conventions; it does not by
itself change code. The mapping to the external review recommendations is
listed at the end.

## Layered structure

The source tree has three layers, with dependencies pointing downward only.

- `infrastructure`: the register-level intermediate representation (RIR),
  validation, serialization, linking, resource estimation, readout, and
  backend adapters. Nothing here knows what an algorithm is.
- `algorithms`: the algorithm library. `common` holds shared numerical
  kernels and standard transformations (QSP/QSVT, Fourier arithmetic,
  search, estimation); `input_model` holds access models, block-encoding
  algebra, and data-loading constructions; the remaining packages
  (`basics`, `qlss`, `qode`, `qpde`, `qnlss`, `qml`, `optimization`,
  `qec`) hold the solver families.
- `applications`: physical models and application-level composition (QHAM,
  QFVM, the catalog and gallery), which consume algorithm protocols rather
  than concrete implementations.

New algorithms go into the matching category file of `algorithms`, never
into the legacy import-compatibility layer (`_compat.py`). A layer may only
depend on layers below it; a solver family may consume `common` and
`input_model`, and an application may consume solver protocols, but the
reverse directions are forbidden.

## The stable RIR core

RIR's stable core consists of registers, module calls, named resources, and
the structured operations `Repeat`, `Control`, and `Adjoint`. These survive
generation, serialization, and restoration unchanged: register names, bit
widths, and views must survive until backend lowering, and module calls and
`Repeat` structures must not be unconditionally expanded during generation
or text serialization.

Everything that encodes a *plan* rather than a circuit identity lives above
RIR: PDE expressions, QHAM closures, matrix-function approximations, and
precision choices are upper-layer plan objects. They are inputs to
generators, not IR nodes, and they must not be added to the IR.

## The three boundaries of variation

Every public API of the framework must make explicit which of the following
three kinds of variation it supports, and handle it by the corresponding
rule.

1. **Changing an algorithm or approximation parameter triggers
   regeneration.** Truncation orders, error targets, phase sequences, and
   similar choices are baked into the generated program. Because widths are
   concrete, a parameter sweep is a Python-generation task: change the
   parameter, regenerate the program.
2. **Replacing an oracle implementation with another realization of the
   same interface goes through checked binding.** The stored open program
   is not regenerated; the replacement is checked against the declared
   interface and assembly conventions at linking time. Binding may change
   internal circuitry, private workspace, and captured resources; it must
   not silently revise a block-encoding scale, a register format, or an
   approximation parameter already used by the caller.
3. **Changing the backend goes through lowering.** Backend export and
   execution are separated from generation; the core does not depend on
   quantum backends, and native dependencies are imported only at execution
   entry points.

The replaceable phase synthesizer of
{obj}`qsp_phases <oracq.algorithms.common.qsvt.qsp_phases>` is the same rule
applied to a numerical kernel: the synthesizer is a component with a
declared contract (target polynomial plus constraints in, time-ordered
phases in the shared reflection convention out), the default implementation
is bundled, and every synthesized sequence — bundled or replaced — is
checked against the convention by the `qsp_response` round trip before it
may enter a generated program.

## Interface metadata: current state and roadmap

Interface metadata is what a checked binding is allowed to rely on. The
current state and the planned extensions are:

- **Current**: register interfaces (names, bit widths, views), the
  block-encoding scale α, and calling capabilities (controlled and adjoint
  availability, resource captures).
- **Roadmap** (not implemented in this stage): error metrics and the
  norm in which they are stated, valid domains, overflow policies for
  fixed-point layouts, physical layout information, and norm conventions.

Each metadata item must also carry its evidence level, distinguishing
three cases: *caller-declared* (a promise the language does not prove,
e.g. a spectral bound), *numerically verified* (checked on a grid or
against a witness at generation or validation time), and *mathematically
proven* (following from the construction). This stage records the roadmap
only; no code changes are attached to it.

## Relation to the review recommendations

This document answers the framework-level recommendations of the review:
keep the current layering and strengthen the conventions between layers
(the first two sections), make the three boundaries of variation a uniform
public-API rule (the third section), and extend scientific interface
metadata along a staged roadmap (the fourth section). The unified evidence
record and the narrow interoperability adapters recommended by the review
are tracked separately in the future-work plan and are not part of this
stage.
