# Oracles and operator representations

**English** · <a href="../../zh/manual/operators.html">简体中文</a>

When choosing an input representation, first ask how the higher-level algorithm needs to access the data. A BE of a matrix, a sparse position query, and a numeric XOR query provide different capabilities.

| Representation | Readable information | Primary use |
|---|---|---|
| {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` | width, signal_qubits, alpha, capabilities | linear-operator composition, QLSS, and evolution |
| {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` | width, work_width, zero-input and clean-work conventions | initial states, right-hand-side states, and reflections |
| {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` | address_width, data_width | reversible data access |
| {obj}`SparseAccess <oracq.algorithms.input_model.oracles.SparseAccess>` | the two operations position/entry, sparsity, value_width | CKS input and sparse adaptation |
| {obj}`Operation <oracq.infrastructure.builder.Operation>` | public registers, resources, and module dependencies | complete quantum operations |

One object can provide multiple views. For example, a full unitary `U` can prepare `U|0>` and also take part in an LCU as a BE with alpha=1. Algorithms check these capabilities through methods and Python protocols; for the conventions see [the conventions owned by algorithms](contracts.md).

## Normalization constants

If the zero-signal corner block satisfies `⟨0|U_A|0⟩=A/alpha_A`, then alpha is the quantity composition requires. The alpha of a product multiplies; the alpha of an LCU is the sum of `abs(coefficient)*alpha` over the terms.

This constant does not make the executor additionally scale the quantum state. It describes the success corner block of the actual circuit. When modifying alpha, the branch weights and angles depending on it must be regenerated.

## Data words and amplitudes

An XOR database returns bit patterns and does not directly provide amplitude access. The input of {obj}`diagonal_block_encoding <oracq.algorithms.input_model.oracles.diagonal_block_encoding>` is an angle word, and the encoded diagonal entry is `alpha*cos(angle/2)`. General matrix elements need a corresponding numeric-to-amplitude conversion.

Sparse input additionally requires an explicit reversible form of the position oracle. The current CKS position interface permutes index in place; an XOR table that preserves index cannot directly replace it.

## Open oracles

Using {obj}`abstract_block_encoding <oracq.algorithms.input_model.oracles.abstract_block_encoding>`, {obj}`abstract_state_prep <oracq.algorithms.input_model.oracles.abstract_state_prep>`, {obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>`, or {obj}`abstract_sparse_access <oracq.algorithms.input_model.oracles.abstract_sparse_access>`, you can write the algorithm first and provide implementations later. The open state only affects executability or exportability; structural checks remain valid.

For a complete binding example see the [tutorial: replacing oracles](../tutorials/oracle-binding.md); for the API see the [oracle catalog](../api/algorithms/input_model/oracles.rst) and [BE composition](../api/algorithms/input_model/block_encoding.rst).
