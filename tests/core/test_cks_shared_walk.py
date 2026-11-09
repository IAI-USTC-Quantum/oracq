"""Compare shared SELECT with the generic CKS LCU on complete coherent inputs."""

import math
import unittest

from oracq import Bits, Builder, ValidationError, estimate_resources, simulate
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import annotate, invoke, qram_database, resources_for
from oracq.algorithms.input_model.sparse import chebyshev_block
from oracq.algorithms.qlss.qlss import CKSConfig, cks_inverse_block


def matrix_block():
    b = Builder("mixed_pauli_involution", {"target": Bits(1), "signal": Bits(1)})
    b.ry(b["signal"], 0.71)
    with b.control(b["signal"], 0):
        b.z(b["target"])
    with b.control(b["signal"]):
        b.x(b["target"])
    b.ry(b["signal"], -0.71)
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0,
                                  self_adjoint_extension=True))


def table_block():
    table = qram_database(1, 1)
    b = Builder("table_involution", {"target": Bits(1), "signal": Bits(1)},
                resources_for(("table", table.operation)))
    invoke(b, table.operation, "table", address=b["target"], data=b["signal"])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0,
                                  self_adjoint_extension=True))


def generic(a, config):
    return lcu([(coefficient, chebyshev_block(a, 2*j+1))
                for j, coefficient in enumerate(config.coefficients())])


def memory(program):
    return {resource.name: [0, 1] for resource in program.main.resources}


class SharedCKSSelectTests(unittest.TestCase):
    def assert_amplitudes(self, expected, actual):
        self.assertEqual(expected.registers, actual.registers)
        for key in expected.amplitudes.keys() | actual.amplitudes.keys():
            self.assertAlmostEqual(expected.amplitudes.get(key, 0), actual.amplitudes.get(key, 0), places=11)

    def test_complete_unitaries_include_padded_and_nonzero_signal_inputs(self):
        a = matrix_block()
        for terms in (1, 2, 3, 4, 5, 7, 8):
            config = CKSConfig(order=8, terms=terms)
            old, new = generic(a, config), cks_inverse_block(a, config)
            self.assertEqual(old.alpha, new.alpha)
            self.assertEqual(old.signal_qubits, new.signal_qubits)
            for target in range(2):
                for signal in range(1 << old.signal_qubits):
                    with self.subTest(terms=terms, target=target, signal=signal):
                        initial = {"target": target, "signal": signal}
                        self.assert_amplitudes(simulate(old.operation.program(), initial=initial),
                                               simulate(new.operation.program(), initial=initial))

    def test_controls_preserve_relative_phase_and_adjoint_cleans_the_flag(self):
        a, config = matrix_block(), CKSConfig(order=5)
        outcomes = []
        for inverse in (generic(a, config), cks_inverse_block(a, config)):
            for undo in (False, True):
                op = inverse.operation
                b = Builder("controlled_cks", {"target": Bits(1), "signal": Bits(inverse.signal_qubits),
                                                "control": Bits(1)})
                for name in ("target", "signal", "control"):
                    b.h(b[name])
                b.global_phase(math.pi/4)
                b.rz(b["target"], math.pi/2)
                with b.control(b["control"], 0):
                    b.call(op, target=b["target"], signal=b["signal"])
                if undo:
                    with b.control(b["control"], 0):
                        with b.adjoint():
                            b.call(op, target=b["target"], signal=b["signal"])
                result = simulate(b.finish().program())
                if undo:
                    amplitude = 1/math.sqrt(1 << (2+inverse.signal_qubits))
                    for (target, _, _), value in result.amplitudes.items():
                        self.assertAlmostEqual(value, amplitude*(1j if target else 1), places=11)
                else:
                    outcomes.append(result)
        self.assert_amplitudes(*outcomes)

    def test_matrix_resource_sharing_keeps_the_same_query_unitary(self):
        a, config = table_block(), CKSConfig(order=3)
        old, new = generic(a, config), cks_inverse_block(a, config)
        self.assertEqual(len(old.operation.module.resources), 3)
        self.assertEqual(len(new.operation.module.resources), 1)
        for target in range(2):
            for signal in range(1 << old.signal_qubits):
                initial = {"target": target, "signal": signal}
                self.assert_amplitudes(
                    simulate(old.operation.program(), memory(old.operation.program()), initial=initial),
                    simulate(new.operation.program(), memory(new.operation.program()), initial=initial))

    def test_query_counts_follow_the_actual_shared_walk_not_a_relabelled_ledger(self):
        a = table_block()
        for order in (2, 4, 8, 16, 32, 64, 128):
            config = CKSConfig(order=order)
            old = estimate_resources(generic(a, config).operation.program())
            new = estimate_resources(cks_inverse_block(a, config).operation.program())
            self.assertEqual(sum(old.qram_queries.values()), order*order)
            self.assertEqual(sum(new.qram_queries.values()), 2*order-1)
            self.assertTrue(old.complete and new.complete)

    def test_unproved_extension_and_invalid_configuration_are_rejected(self):
        b = Builder("not_self_adjoint", {"target": Bits(1), "signal": Bits(1)})
        b.rz(b["target"], math.pi/2)
        a = BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))
        with self.assertRaises(ValidationError):
            cks_inverse_block(a, CKSConfig(order=2))
        with self.assertRaises(ValidationError):
            cks_inverse_block(matrix_block(), CKSConfig(order=129))
        with self.assertRaises(ValidationError):
            cks_inverse_block(matrix_block(), CKSConfig(order=3, terms=4))


if __name__ == "__main__":
    unittest.main()
