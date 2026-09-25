"""Invariant assertion library: witness primitives shared across algorithm tests.

All four primitives take a ``unittest.TestCase`` instance as their first
argument and raise ``AssertionError`` on failure. This file name lacks the
``test_`` prefix, so unittest discover does not collect it; the self-tests
live in ``test_witness.py``.
"""

import random

from oracq import ValidationError, bind, simulate, unresolved

# simulate drops amplitudes with |a| < 1e-15; uncomputation checks use a looser zero threshold.
_ZERO_AMPLITUDE = 1e-9


def _register_index(program):
    return {r.name: i for i, r in enumerate(program.main.registers)}


def _basis_indices(program, samples):
    """All basis states (small register space) or a fixed-seed random sample (large space); determinism first."""
    registers = program.main.registers
    width = sum(r.type.width for r in registers)
    size = 1 << width
    if samples is not None:
        return list(samples)
    if size <= 16:
        return list(range(size))
    return sorted(random.Random(0).sample(range(size), 16))


def _initial_for_index(program, index):
    """Split a flat basis-state index into simulate's initial dict in register declaration order (little-endian)."""
    initial = {}
    offset = 0
    for reg in program.main.registers:
        initial[reg.name] = (index >> offset) & ((1 << reg.type.width) - 1)
        offset += reg.type.width
    return initial


def _inner_product(first, second):
    return sum(
        amplitude.conjugate() * second.get(key, 0j) for key, amplitude in first.items()
    )


def assert_unitary(case, program, *, places=9, samples=None):
    """Sampled witness for W†W = I: every basis-state column is normalized (Σ|a|²=1) and pairwise orthogonal.

    samples explicitly provides the list of basis-state indices; by default a
    small register space (≤16 basis states) uses all of them, and a large space
    draws 16 with the fixed seed ``random.Random(0)``.
    """
    indices = _basis_indices(program, samples)
    columns = [
        simulate(program, initial=_initial_for_index(program, index)).amplitudes
        for index in indices
    ]
    for index, amplitudes in zip(indices, columns, strict=True):
        norm = sum(abs(a) ** 2 for a in amplitudes.values())
        case.assertAlmostEqual(
            norm, 1.0, places=places, msg=f"output column of basis state {index} is not normalized: Σ|a|²={norm}"
        )
    for i, first in enumerate(columns):
        for j in range(i + 1, len(columns)):
            inner = _inner_product(first, columns[j])
            case.assertAlmostEqual(
                inner.real,
                0.0,
                places=places,
                msg=f"columns {indices[i]} and {indices[j]} are not orthogonal: inner product {inner}",
            )
            case.assertAlmostEqual(
                inner.imag,
                0.0,
                places=places,
                msg=f"columns {indices[i]} and {indices[j]} are not orthogonal: inner product {inner}",
            )


def assert_uncomputation(case, program, *, initial=None, work_registers=None):
    """Uncomputation witness: local registers must be uncomputed and the named root registers must be 0 in the output.

    simulate forces uncomputation at LocalExit (raising ValidationError when
    not uncomputed), which is converted here into AssertionError; the root
    registers listed in work_registers must evaluate to 0 in every basis state
    with nonzero amplitude. Returns simulate's RegisterState so the caller can
    keep asserting.
    """
    try:
        state = simulate(program, initial=initial)
    except ValidationError as exc:
        case.fail(f"uncomputation failed: {exc}")
    if work_registers:
        index = _register_index(program)
        for name in work_registers:
            position = index[name]
            for key, amplitude in state.amplitudes.items():
                if abs(amplitude) > _ZERO_AMPLITUDE:
                    case.assertEqual(
                        key[position],
                        0,
                        msg=f"work register {name} not uncomputed: basis state {key} amplitude {amplitude}",
                    )
    return state


def assert_bind_invariant(case, abstract_program, bindings, *, tolerance=0.0):
    """Binding invariant: candidate implementations of the same abstract slot must give consistent observable distributions.

    abstract_program must have exactly one unbound slot (resolved via
    unresolved); each (label, Operation) in bindings is treated as one
    candidate implementation of that slot: bind one at a time, simulate, and
    compare the |a|² distributions basis state by basis state. tolerance=0
    means an exact cross-check; bindings with quantization error, such as
    QRAM, pass an error bound (e.g. 0.02, following the existing convention).
    """
    names = [r.name for r in unresolved(abstract_program)]
    case.assertEqual(len(names), 1, msg=f"abstract program must have exactly one unbound slot: {names}")
    slot = names[0]
    reference = None
    reference_label = None
    for label, operation in bindings.items():
        bound = bind(abstract_program, {slot: operation})
        distribution = {
            key: abs(amplitude) ** 2 for key, amplitude in simulate(bound).amplitudes.items()
        }
        if reference is None:
            reference, reference_label = distribution, label
            continue
        for key in set(reference) | set(distribution):
            expected = reference.get(key, 0.0)
            actual = distribution.get(key, 0.0)
            message = (
                f"candidate {label} and {reference_label} disagree on the "
                f"distribution at basis state {key}: {actual} != {expected}"
            )
            if tolerance:
                case.assertAlmostEqual(actual, expected, delta=tolerance, msg=message)
            else:
                case.assertEqual(actual, expected, msg=message)


def block_column(be, column):
    """Column `column` of the BE's (0,0) block (after multiplying by alpha)."""
    state = simulate(be.operation.program(), initial={"target": column})
    size = 1 << be.width
    return [state.amplitudes.get((row, 0), 0) * be.alpha for row in range(size)]


def assert_block_equals(case, be, matrix, *, places=9):
    """Column-by-column cross-check of a block encoding's (0,0) block against a dense matrix (multiplying alpha back in)."""
    for column in range(len(matrix)):
        actual = block_column(be, column)
        for row in range(len(matrix)):
            case.assertAlmostEqual(actual[row], matrix[row][column], places=places)
