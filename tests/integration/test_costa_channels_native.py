"""Actual sparse-simulator precision checks against independent NumPy solves."""

import math
import unittest

import numpy as np

from oracq import Bits, Builder, run_pysparq
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import annotate, gate_state_prep, invoke, resources_for
from oracq.algorithms.qlss.qlss import CostaConfig, costa_qlss


class CostaNativeChannelTests(unittest.TestCase):
    def test_hermitian_and_complex_nonsymmetric_solutions(self):
        c, s = math.cos(0.305), math.sin(0.305)
        left = np.array([[c, -s], [s, c]], dtype=complex)
        rhs = np.array([math.sqrt(0.3), np.exp(0.4j) * math.sqrt(0.7)])
        for nonsymmetric in (False, True):
            with self.subTest(nonsymmetric=nonsymmetric):
                right = (np.diag([np.exp(-0.415j), np.exp(0.415j)]) if nonsymmetric
                         else left.T)
                matrix = left @ np.diag([0.4, 0.9]) @ right
                b = Builder("native_costa_svd", {"target": Bits(1), "signal": Bits(1)})
                if nonsymmetric:
                    b.rz(b["target"], 0.83)
                else:
                    b.ry(b["target"], -0.61)
                b.z(b["signal"])
                for i, value in enumerate((0.4, 0.9)):
                    with b.control(b["target"], i):
                        b.ry(b["signal"], 2 * math.acos(value))
                b.ry(b["target"], 0.61)
                encoded = BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))
                for column in range(2):
                    column_builder = Builder("native_block_column",
                                             {"target": Bits(1), "signal": Bits(1)},
                                             resources_for(("a", encoded.operation)))
                    if column:
                        column_builder.x(column_builder["target"])
                    invoke(column_builder, encoded.operation, "a",
                           target=column_builder["target"], signal=column_builder["signal"])
                    state = run_pysparq(column_builder.finish().program())
                    np.testing.assert_allclose([state.amplitudes.get((i, 0), 0) for i in range(2)],
                                               matrix[:, column], atol=1e-12, rtol=0)
                solve = costa_qlss(encoded, gate_state_prep(rhs.tolist()),
                                   CostaConfig(steps=512, kappa=2.5))
                report = {}
                state = run_pysparq(solve.operation.program(), max_steps=50_000_000,
                                    max_states=512, report=report)
                self.assertEqual(report["native_calls"], 0)
                self.assertGreater(report["gate_events"], 20_000)
                self.assertAlmostEqual(sum(abs(x) ** 2 for x in state.amplitudes.values()), 1, places=9)
                actual = np.array([state.amplitudes.get((i, 0), 0j) for i in range(2)])
                probability = float(np.vdot(actual, actual).real)
                self.assertGreater(probability, 0.49)
                self.assertLess(probability, 0.51)
                actual /= math.sqrt(probability)
                expected = np.linalg.solve(matrix, rhs)
                expected /= np.linalg.norm(expected)
                overlap = np.vdot(expected, actual)
                aligned = actual * overlap.conjugate() / abs(overlap)
                self.assertLess(float(np.linalg.norm(aligned - expected)), 1e-3)
                if nonsymmetric:
                    adjoint = np.linalg.solve(matrix.conj().T, rhs)
                    adjoint /= np.linalg.norm(adjoint)
                    self.assertLess(float(abs(np.vdot(adjoint, actual)) ** 2), 0.7)


if __name__ == "__main__":
    unittest.main()
