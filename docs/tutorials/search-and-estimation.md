# Search for an element and estimate the success probability

**English** · <a href="../zh/tutorials/search-and-estimation.html">简体中文</a>

<a href="../zh/manual/algorithms/grover.html">Grover search</a> marks good states with a phase oracle. Below we mark `3` among four basis states and read that state out after one iteration.

```{testcode}
from oracq.algorithms.input_model.oracles import phase_marks
from oracq.algorithms.common.search import grover
from oracq import simulate

result = grover(phase_marks(2, [3]), 2, iterations=1)
state = simulate(result.operation.program())
print(state.amplitudes)
assert abs(state.amplitudes[(3, 0)] - 1) < 1e-12
```

```{testoutput}
{(3, 0): (0.9999999999999996-1.8369701987210287e-16j)}
```

{obj}`grover <oracq.algorithms.common.search.grover>` takes the phase oracle declared by {obj}`phase_marks <oracq.algorithms.input_model.oracles.phase_marks>`. The printed dictionary holds a single basis state: the key `(3, 0)` means the target reads `3` and the signal is zero, and the squared magnitude of the amplitude is approximately `1` (the imaginary part is floating-point noise). The output contains the target and the signal. Here the signal is an empty interface and the result word is zero; a more elaborate input preparation may carry work registers.

## Estimating the probability

Standard <a href="../zh/manual/algorithms/qae.html">amplitude estimation</a> runs phase estimation on the Grover iterate; its entry point is {obj}`amplitude_estimation <oracq.algorithms.common.estimation.amplitude_estimation>`. It takes an initial-state preparation and a set of good states and outputs a phase register. The initial state below is uniform over `0` and `1`, so the good state `1` has probability `1/2`.

```{testcode}
from oracq.algorithms.input_model.oracles import uniform_state
from oracq.algorithms.common.estimation import amplitude_estimation, amplitude_from_phase

operation = amplitude_estimation(uniform_state(1), [1], precision=3)
state = simulate(operation.program())
phases = {key[2] for key in state.amplitudes}
# The decoded values differ from 1/2 only in the last floating-point bit, and trig libraries on
# different platforms may differ in that bit; round to 12 digits before printing so the output
# is verbatim-identical on every platform.
print({value: round(amplitude_from_phase(value, 3), 12) for value in sorted(phases)})
assert phases == {2, 6}
assert all(abs(amplitude_from_phase(value, 3) - 0.5) < 1e-12 for value in phases)
```

```{testoutput}
{2: 0.5, 6: 0.5}
```

The initial state is prepared by {obj}`uniform_state <oracq.algorithms.input_model.oracles.uniform_state>`. The printed dictionary decodes the two phase words back into probabilities via {obj}`amplitude_from_phase <oracq.algorithms.common.estimation.amplitude_from_phase>`: both `2` and `6` give approximately `1/2`; the unrounded decoded values differ from `1/2` only in the last floating-point bit, which comes from the rounding of a finite-precision readout. In general, a finite-precision phase readout yields an approximate estimate; repeated sampling and statistical post-processing happen on the host side.

## Related pages

- Algorithm pages: <a href="../zh/manual/algorithms/grover.html">Grover search</a>, <a href="../zh/manual/algorithms/amplitude-amplification.html">amplitude amplification</a>, <a href="../zh/manual/algorithms/qae.html">amplitude estimation</a>, <a href="../zh/manual/algorithms/qpe.html">phase estimation</a>
- API reference: [Search and amplitude amplification](../api/algorithms/common/search.rst), [Phase, amplitude, and overlap estimation](../api/algorithms/common/estimation.rst), [Oracle declarations and implementations](../api/algorithms/input_model/oracles.rst)
- Continue with: [Running and modifying the algorithm gallery](gallery.md) (further cases covering both algorithms of this example)
