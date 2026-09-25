# Running and modifying the algorithm gallery

**English** · [简体中文](../zh/tutorials/gallery.html)

The gallery places algorithms from different families into one unified run pipeline, making it easy to compare their inputs, registers, and readout. All cases are constructed uniformly by {obj}`algorithm_gallery <oracq.applications.gallery.algorithm_gallery>`.

```bash
uv run python examples/algorithm_gallery.py
```

Each subdirectory contains `closed.rir.yaml`, `modular.originir`, and `toffoli_u3_cz.originir`. Bernstein–Vazirani additionally saves the open description from before binding. `index.json` records the family, the module count, and the readout instructions.

## Choosing which results to check

Do not use the same readout for every algorithm. BV reads the secret string; QPE and amplitude estimation read the phase; the Hadamard test estimates the Z expectation of the probe; repetition-code recovery checks the logical state and keeps the syndrome.

Use the following code to find the cases:

```{testcode}
from oracq.applications.gallery import algorithm_gallery

# Build every gallery case and index them by name for easy lookup.
cases = {case.name: case for case in algorithm_gallery()}
# Print all case names currently covered by the gallery.
print("\n".join(sorted(cases)))
# Confirm the QAOA max-cut case exists (the gallery covers the variational family).
assert "qaoa_maxcut" in cases
# Every case carries a family label: order finding belongs to the number_theory family.
assert cases["order_finding"].family == "number_theory"
# The gallery holds 22 cases; adding or removing one trips this assertion as a reminder to update in sync.
assert len(cases) == 22
```

```{testoutput}
amplitude_amplification
amplitude_estimation
ansatz
bernstein_vazirani
cycle_walk
deutsch_jozsa
fourier_add
grover
hadamard_imag
hadamard_real
hamiltonian_trotter
modular_multiply
order_finding
phase_estimation
qaoa_maxcut
qft
quantum_counting
repetition_bit
repetition_phase
simon
swap_test
vqe_pauli_measurement
```

The 22 printed names are all the cases currently covered by the gallery, spanning search, estimation, arithmetic, walks, error correction, and the variational family; when cases are added or removed, this output must be updated together with the `len(cases) == 22` assertion.

To add application examples, call the corresponding algorithm file from your own script. To extend the public gallery, provide {obj}`GalleryCase <oracq.applications.gallery.GalleryCase>` with an operation and concrete readout instructions, and add an independent mathematical witness.

## Cross-checking against a real backend

```bash
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

The script compares the complex amplitudes of the reference executor, PySparQ, and OriginIR case by case and exits immediately when a tolerance is exceeded. It verifies that the same circuit is semantically consistent across the backends; whether a given parameter suits a real problem still requires application-level analysis.

## Related pages

- Manual: [Algorithm catalog](../zh/manual/algorithms/index.html) (most algorithms covered by the gallery have entry pages), [Export and execution backends](../manual/backends.md)
- API reference: [Algorithm gallery](../api/applications/gallery.rst)
- Continue with: [Search for an element and estimate the success probability](search-and-estimation.md), [Supplying your own Hamiltonian decomposition](hamiltonian.md)
