"""Independent inverse-system witnesses for the Costa walk success channel."""

import cmath
import math
import unittest

from oracq import Bits, Builder, simulate
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    annotate,
    basis_state,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.algorithms.qlss.qlss import (
    BlockSystem,
    CostaConfig,
    LinearSystem,
    SpectralPromise,
    costa_qlss,
    costa_walk,
    dolph_chebyshev_plan,
    lcu_filter,
    make_costa_qlss,
    schedule,
)


def public_system(*, complex_matrix=False, work_width=0):
    """Use a known singular-value decomposition, independently solve its 2x2 matrix."""
    c, s = math.cos(0.61 / 2), math.sin(0.61 / 2)
    left = ((c, -s), (s, c))
    right = (((cmath.exp(-0.415j), 0), (0, cmath.exp(0.415j))) if complex_matrix
             else ((c, s), (-s, c)))
    singular = (0.4, 0.9)
    matrix = tuple(tuple(sum(left[i][k] * singular[k] * right[k][j] for k in range(2))
                         for j in range(2)) for i in range(2))
    b = Builder("public_svd_block", {"target": Bits(1), "signal": Bits(1)})
    if complex_matrix:
        b.rz(b["target"], 0.83)
    else:
        b.ry(b["target"], -0.61)
    b.z(b["signal"])
    for i, value in enumerate(singular):
        with b.control(b["target"], i):
            b.ry(b["signal"], 2 * math.acos(value))
    b.ry(b["target"], 0.61)
    a = BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))
    rhs = (math.sqrt(0.3), cmath.exp(0.4j) * math.sqrt(0.7))
    prep = gate_state_prep(rhs)
    if work_width:
        b = Builder("rhs_clean_work_extension", {"target": Bits(1), "work": Bits(work_width)},
                    resources_for(("prep", prep.operation)))
        invoke(b, prep.operation, "prep", target=b["target"], work=b["work"][:0])
        # A nontrivial unitary extension agrees only on the clean-work subspace.
        with b.control(b["work"]):
            b.x(b["target"])
        prep = StatePreparation(annotate(b.finish(), "state_prep_isometry",
                                        zero_input=True, clean_work=True))
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    expected = ((matrix[1][1] * rhs[0] - matrix[0][1] * rhs[1]) / determinant,
                (matrix[0][0] * rhs[1] - matrix[1][0] * rhs[0]) / determinant)
    norm = math.sqrt(sum(abs(x) ** 2 for x in expected))
    return a, prep, tuple(x / norm for x in expected), matrix


def conditional_solution(state, expected):
    vector = tuple(state.amplitudes.get((i, 0), 0j) for i in range(2))
    probability = sum(abs(x) ** 2 for x in vector)
    if probability < 1e-15:
        return probability, math.inf
    actual = tuple(x / math.sqrt(probability) for x in vector)
    overlap = sum(x.conjugate() * y for x, y in zip(expected, actual, strict=True))
    phase = overlap.conjugate() / abs(overlap) if overlap else 1
    error = math.sqrt(sum(abs(y * phase - x) ** 2 for x, y in zip(expected, actual, strict=True)))
    return probability, error


class CostaChannelTests(unittest.TestCase):
    def test_declared_zero_signal_channel_solves_even_and_odd_walks(self):
        for complex_matrix, work_width, steps in ((False, 0, 33), (True, 0, 64), (True, 1, 33)):
            with self.subTest(complex_matrix=complex_matrix, work_width=work_width, steps=steps):
                a, prep, expected, matrix = public_system(complex_matrix=complex_matrix,
                                                        work_width=work_width)
                for column in range(2):
                    encoded = simulate(a.operation.program(), initial={"target": column})
                    for row in range(2):
                        self.assertAlmostEqual(encoded.amplitudes.get((row, 0), 0), matrix[row][column], places=12)
                result = costa_qlss(a, prep, CostaConfig(steps=steps, kappa=2.5))
                state = simulate(result.operation.program(), max_steps=10_000_000)
                self.assertAlmostEqual(sum(abs(x) ** 2 for x in state.amplitudes.values()), 1, places=10)
                probability, error = conditional_solution(state, expected)
                self.assertGreater(probability, 0.45)
                self.assertLess(error, 0.01)

    def test_problem_layer_norm_probe_uses_the_correct_solution_channel(self):
        b = Builder("half_identity_block", {"target": Bits(1), "signal": Bits(1)})
        b.z(b["signal"])
        b.ry(b["signal"], 2 * math.acos(0.5))
        encoded = BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))
        problem = LinearSystem(block=BlockSystem(encoded, basis_state(1), SpectralPromise(0.5, 0.5)),
                               rhs_norm=2)
        result = make_costa_qlss(CostaConfig(steps=33))(problem)
        state = simulate(result.operation.program(), max_steps=10_000_000)
        probe = simulate(result.norm_probe.operation.program(), max_steps=10_000_000)
        p_solve = sum(abs(v) ** 2 for key, v in state.amplitudes.items() if key[1] == 0)
        p_joint = sum(abs(v) ** 2 for key, v in probe.amplitudes.items() if key[1] == 0)
        self.assertGreater(p_solve, 0.45)
        self.assertAlmostEqual(result.recover_norm(p_solve, p_joint), 4, places=10)

    def test_coherent_solve_phase_matches_an_independently_prepared_eigenvector(self):
        a, prep, _, _ = public_system(complex_matrix=True)
        config = CostaConfig(steps=3, kappa=2.5)
        kernel = costa_qlss(a, prep, config).operation
        vector = [0j] * 16
        vector[0], vector[12], vector[14] = 2**-0.5, -0.5j, -0.5j
        flags = gate_state_prep(vector)
        flag_state = simulate(flags.operation.program())
        for i, value in enumerate(vector):
            self.assertAlmostEqual(flag_state.amplitudes.get((i, 0), 0), value, places=12)
        walks = [costa_walk(a, prep, schedule((i+1)/config.steps, config.kappa))
                 for i in range(config.steps)]
        filtering = lcu_filter(walks[-1], dolph_chebyshev_plan())
        results = []
        for independent in (False, True):
            registers = {r.name: r.type for r in kernel.module.registers}
            registers["control"] = Bits(1)
            operands = [("kernel", kernel), ("prep", prep.operation), ("flags", flags.operation),
                        ("filter", filtering)] + [(f"walk{i}", w) for i, w in enumerate(walks)]
            b = Builder("coherent_costa_phase", registers, resources_for(*operands))
            b.h(b["control"])
            with b.control(b["control"], 0):
                if independent:
                    invoke(b, prep.operation, "prep", target=b["target"], work=b["signal"][:0])
                    invoke(b, flags.operation, "flags", target=b["signal"][1:5], work=b["signal"][:0])
                    for i, walk in enumerate(walks):
                        invoke(b, walk, f"walk{i}", target=b["target"], signal=b["signal"][:5])
                    invoke(b, filtering, "filter", target=b["target"], signal=b["signal"])
                    b.x(b["signal"][1])
                else:
                    invoke(b, kernel, "kernel", target=b["target"], signal=b["signal"])
            results.append(simulate(b.finish().program(), max_steps=10_000_000))
        for key in results[0].amplitudes.keys() | results[1].amplitudes.keys():
            self.assertAlmostEqual(results[0].amplitudes.get(key, 0),
                                   results[1].amplitudes.get(key, 0), places=11)

    def test_controlled_adjoint_preserves_coherent_zero_and_nonzero_signals(self):
        a, prep, _, _ = public_system(complex_matrix=True, work_width=1)
        op = costa_qlss(a, prep, CostaConfig(steps=3, kappa=2.5)).operation
        registers = {r.name: r.type for r in op.module.registers}
        registers["control"] = Bits(1)
        b = Builder("controlled_costa_inverse", registers, resources_for(("solve", op)))
        b.h(b["target"])
        b.h(b["signal"])
        b.h(b["control"])
        b.rz(b["target"], math.pi / 2)
        b.global_phase(math.pi / 4)
        with b.control(b["control"], 0):
            invoke(b, op, "solve", target=b["target"], signal=b["signal"])
        with b.control(b["control"], 0):
            with b.adjoint():
                invoke(b, op, "solve", target=b["target"], signal=b["signal"])
        state = simulate(b.finish().program(), max_steps=10_000_000)
        amplitude = 1 / math.sqrt(1 << (b["signal"].width + 2))
        for (target, _, _), value in state.amplitudes.items():
            self.assertAlmostEqual(value, amplitude * (1j if target else 1), places=10)


if __name__ == "__main__":
    unittest.main()
