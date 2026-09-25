"""Semantic and resource validation of operator-level optimization (arXiv:2509.08807 appendix D).

Each optimization asserts two things: the encoded object (matrix/state) is unchanged
element by element; the logical resources drop (Toffoli count or controlled-layer
count). merge_similar preserves the encoded matrix alpha*corner; the
sub-normalization alpha itself tightens according to the merging algebra.
"""

import unittest

from oracq import Builder, estimate_resources, simulate
from oracq.algorithms.common.spectral_synthesis import (
    fanout_spectral_diagonal,
    merge_similar,
    uniformly_controlled_prep,
)
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.operators import identity
from oracq.algorithms.input_model.oracles import gate_state_prep
from oracq.algorithms.input_model.spectral import spectral_diagonal


def applied(operation, initial):
    b = Builder("probe", {r.name: r.type for r in operation.module.registers})
    for key, value in initial.items():
        for bit in range(b[key].width):
            if (value >> bit) & 1:
                b.x(b[key][bit])
    b.call(operation, **{r.name: b[r.name] for r in operation.module.registers})
    return simulate(b.finish().program()).amplitudes


def encoded_matrix_corner(be, dim):
    return [
        applied(be.operation, {"target": column})[(column, 0)] * be.alpha
        for column in range(dim)
    ]


class UniformlyControlledPrepTests(unittest.TestCase):
    AMPLITUDES = [0.5, 0.1, -0.2j, 0.4, 0.3, 0.05, 0.2, 0.65]

    def test_unitary_equality(self):
        ucr = applied(uniformly_controlled_prep(self.AMPLITUDES).operation, {})
        plain = applied(gate_state_prep(self.AMPLITUDES).operation, {})
        self.assertEqual(set(ucr), set(plain))
        for key, value in plain.items():
            self.assertAlmostEqual(ucr[key], value, places=12)

    def test_toffoli_reduction(self):
        plain = estimate_resources(gate_state_prep(self.AMPLITUDES).operation.program())
        ucr = estimate_resources(uniformly_controlled_prep(self.AMPLITUDES).operation.program())
        self.assertLess(ucr.atoms.get("toffoli", 0), plain.atoms.get("toffoli", 0))
        self.assertLessEqual(len(ucr.rotations), len(plain.rotations))


class MergeSimilarTests(unittest.TestCase):
    def test_preserves_encoded_matrix(self):
        terms = [(1.0, identity(3)), (-1.0, identity(3)), (0.5, identity(3)), (0.7, identity(3))]
        merged = merge_similar(terms)
        self.assertEqual(len(merged), 1)
        before, after = lcu(terms), lcu(merged)
        self.assertAlmostEqual(before.alpha, 3.2, places=12)
        self.assertAlmostEqual(after.alpha, 1.2, places=12)
        for column in (0, 5, 7):
            raw_before = applied(before.operation, {"target": column})[(column, 0)]
            raw_after = applied(after.operation, {"target": column})[(column, 0)]
            self.assertAlmostEqual(raw_before * before.alpha, raw_after * after.alpha, places=12)

    def test_drops_cancelled_terms(self):
        terms = [(1.0, identity(2)), (-1.0, identity(2)), (0.5, identity(2))]
        self.assertEqual(len(merge_similar(terms)), 1)
        self.assertEqual(merge_similar(terms)[0][0], 0.5)

    def test_reduces_lcu_branches(self):
        terms = [(1.0 + 0.1j * i, identity(4)) for i in range(5)]
        merged = merge_similar(terms)
        before = estimate_resources(lcu(terms).operation.program())
        after = estimate_resources(lcu(merged).operation.program())
        self.assertLess(
            after.atoms.get("toffoli", 0) + len(after.rotations),
            before.atoms.get("toffoli", 0) + len(before.rotations),
        )


class FanoutSpectralDiagonalTests(unittest.TestCase):
    SPECTRUM = {-2: 0.5 + 0.1j, -1: 0.25j, 0: -0.35, 1: 0.4, 2: 0.2}

    def test_corner_equality(self):
        sequential = spectral_diagonal(self.SPECTRUM, 3)
        fanout = fanout_spectral_diagonal(self.SPECTRUM, 3)
        self.assertEqual(sequential.alpha, fanout.alpha)
        for column in range(8):
            plain = applied(sequential.operation, {"target": column})[(column, 0)]
            merged = applied(fanout.operation, {"target": column})[(column, 0)]
            self.assertAlmostEqual(plain, merged, places=11)

    def test_removes_controlled_rotations(self):
        """Structural assertion: the fan-out-form RIR contains no Control nodes wrapping rotation gates."""
        from oracq.infrastructure.ir import Control, Primitive

        def controlled_rotations(body):
            rotations = {"ry", "rx", "rz", "phase"}
            count = 0
            for node in body:
                if isinstance(node, Control):
                    count += sum(
                        1 for sub in node.body if isinstance(sub, Primitive) and sub.op in rotations
                    )
                    count += controlled_rotations(node.body)
            return count

        program = fanout_spectral_diagonal(self.SPECTRUM, 3).operation.program()
        sequential = spectral_diagonal(self.SPECTRUM, 3).operation.program()
        self.assertEqual(controlled_rotations(program.main.body), 0)
        self.assertGreater(controlled_rotations(sequential.main.body), 0)


if __name__ == "__main__":
    unittest.main()
