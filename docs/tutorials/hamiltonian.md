# Supplying your own Hamiltonian decomposition

**English** · [简体中文](../zh/tutorials/hamiltonian.html)

Suppose you already know the term decomposition of `H=0.3I+0.7X`. You can hand it to the [Trotter implementation](../zh/manual/algorithms/trotter.html) through an ordinary Python method — no dense matrix needs to be built first, and no new language type needs to be registered. The decomposition is written directly as a tuple of {obj}`TrotterTerm <oracq.algorithms.common.hamiltonian.TrotterTerm>`, each factor being a {obj}`PauliOperator <oracq.algorithms.common.hamiltonian.PauliOperator>`; {obj}`hamiltonian_simulation <oracq.algorithms.common.hamiltonian.hamiltonian_simulation>` accepts any object that can produce these terms.

```{testcode}
import cmath
import math
from oracq.algorithms.common.hamiltonian import PauliOperator, TrotterTerm, hamiltonian_simulation
from oracq import simulate

class MyHamiltonian:
    hermitian = True

    def trotter_list(self):
        return (TrotterTerm(0.3, PauliOperator("I")),
                TrotterTerm(0.7, PauliOperator("X")))

evolution = hamiltonian_simulation(MyHamiltonian(), time=0.4, steps=3)
state = simulate(evolution.operation.program())
expected = cmath.exp(-0.12j) * math.cos(0.28)
print(state.amplitudes)
print(evolution.alpha)
assert abs(state.amplitudes[(0, 0)] - expected) < 1e-12
```

```{testoutput}
{(0, 0): (0.9541441386892553-0.11505006784720223j), (1, 0): (-0.03308314468637069-0.274368274461448j)}
1.0
```

The printed `(0, 0)` amplitude is the floating-point form of the analytic value `exp(-0.12j)·cos(0.28)`; the `X` term moves the remaining probability onto `(1, 0)`. The two terms commute here, so the analytic result can be checked against the amplitudes directly. `evolution.alpha=1` because what is generated is the full unitary evolution.

For non-commuting terms, the Trotter step count controls the product-formula approximation. A Hermitian declaration alone is not enough: each term must supply a callable evolution implementation. If the QSP route is chosen, the BE access it needs and an actual QSP kernel must be provided.

This example shows the role of the algorithm contract: the Hamiltonian object provides {obj}`trotter_list() <oracq.algorithms.common.hamiltonian.TrotterizableProtocol.trotter_list>`, and the Trotter algorithm is responsible for checking and calling it; the RIR only receives the modules and gates produced in the end.

## Related pages

- Manual: [Oracles and operator representations](../manual/operators.md) (operator views and protocols for Hamiltonians)
- Algorithm pages: [Hamiltonian evolution](../zh/manual/algorithms/hamiltonian-simulation.html), [Trotter product-formula simulation](../zh/manual/algorithms/trotter.html), [truncated Taylor block encoding](../zh/manual/algorithms/taylor-block-encoding.html)
- API reference: [Hamiltonian evolution](../api/algorithms/common/hamiltonian.rst)
