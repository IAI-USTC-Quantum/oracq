# Replacing an algorithm's oracle

**English** · [简体中文](../zh/tutorials/oracle-binding.html)

We first write a Bernstein–Vazirani program that does not know its function implementation, then bind a gate implementation to it. This makes the difference between "the algorithm itself is complete" and "the input oracle is still unfinished" visible.

## Declare the input first

```{testcode}
from oracq.algorithms.input_model.oracles import abstract_database
from oracq.algorithms.basics.oracle_algorithms import bernstein_vazirani
from oracq import unresolved

given = abstract_database("BooleanFunction", 3, 1)
opened = bernstein_vazirani(given).program()
print([item.name for item in unresolved(opened)])
assert [item.name for item in unresolved(opened)] == ["BooleanFunction"]
```

```{testoutput}
['BooleanFunction']
```

{obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>` declares an abstract database that has no implementation yet, and {obj}`bernstein_vazirani <oracq.algorithms.basics.oracle_algorithms.bernstein_vazirani>` assembles the algorithm on top of it. The printed list is exactly the names of the still-unbound open slots found by {obj}`unresolved <oracq.infrastructure.linking.unresolved>`: the algorithm body is already complete; the only thing missing is the input implementation named `BooleanFunction`. The declaration provides a three-bit address and a one-bit XOR result. The algorithm assumes the function has the form `f(x)=s·x XOR c`; the declaration itself does not prove that premise.

## Bind a gate implementation

```{testcode}
from oracq.algorithms.basics.oracle_algorithms import affine_boolean_oracle
from oracq import bind, simulate

implementation = affine_boolean_oracle(3, secret=5, bias=1)
closed = bind(opened, {"BooleanFunction": implementation.operation})
assert not unresolved(closed)

state = simulate(closed)
probability = sum(abs(a)**2 for key, a in state.amplitudes.items() if key[0] == 5)
print(state.amplitudes)
print(probability)
assert abs(probability - 1) < 1e-12
```

```{testoutput}
{(5, 0): (-0.7071067811865471+0j), (5, 1): (0.7071067811865471+0j)}
0.9999999999999989
```

The gate implementation is provided by {obj}`affine_boolean_oracle <oracq.algorithms.basics.oracle_algorithms.affine_boolean_oracle>` using ordinary reversible gates; the printed amplitudes are nonzero only on branches where the input reads `5`; the `±1/√2` phase difference on the value bit comes from the affine bias, and the total probability on the second line is approximately `1`. Reading the input yields `5`, the secret bit string interpreted with the least significant bit first. The affine bias changes the phase but not this outcome.

## Swap in a QRAM

The same open slot can be bound to {obj}`qram_database(3,1) <oracq.algorithms.input_model.oracles.qram_database>`. Use {obj}`Binding(..., {"table": "truth"}) <oracq.infrastructure.linking.Binding>` to map its resources to the entry, then supply the `truth` table at run time.

The gate implementation and the QRAM implementation must honor the same XOR semantics. {obj}`bind <oracq.infrastructure.linking.bind>` checks interfaces and capabilities; the mathematical form of the function remains the application's responsibility.

## Related pages

- Manual: [Oracles and operator representations](../manual/operators.md) (the oracle paradigm and views)
- Specification: [Open IR](../reference/open-ir.md) (open declarations, batched binding, and resource capture)
- Algorithm pages: [Bernstein–Vazirani](../zh/manual/algorithms/bernstein-vazirani.html), [XOR database views](../zh/manual/algorithms/xor-database.html)
- API reference: [Oracle declarations and implementations](../api/algorithms/input_model/oracles.rst), [Binding and capability analysis](../api/infrastructure/linking.rst)
- Continue with: [From a paper's access model to implementation comparison](algorithm-research.md)
