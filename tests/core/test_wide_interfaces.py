"""Generic algorithm interfaces accept widths beyond one machine word."""

import unittest
from functools import reduce

from oracq import (
    Builder,
    ValidationError,
    dumps,
    estimate_resources,
    loads,
    simulate,
)
from oracq.algorithms.basics.oracle_algorithms import affine_boolean_oracle, simon_nullspace
from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.common.fourier import fourier_add, qft
from oracq.algorithms.common.hamiltonian import PauliOperator
from oracq.algorithms.input_model.operators import linear_combination, product, zero
from oracq.algorithms.input_model.oracles import gate_database, qram_database
from oracq.algorithms.optimization.gradient import function_phase_oracle, gradient_estimation


class WideInterfaceTests(unittest.TestCase):
    def test_fixed_format_encodes_beyond_one_machine_word(self):
        fmt = FixedFormat(80, 40)
        self.assertEqual(fmt.decode(fmt.encode(-1.5)), -1.5)
        with self.assertRaises(ValidationError):
            FixedFormat(1, 0)

    def test_wide_qft_and_fourier_addition_roundtrip(self):
        for op, qubits in ((qft(65), 65), (fourier_add(65), 130)):
            with self.subTest(qubits=qubits):
                program = loads(dumps(op.program(), format="json"))
                self.assertEqual(estimate_resources(program).qubits, qubits)

    def test_wide_gate_database_loads_sixty_five_bit_words(self):
        word = (1 << 64) + 9
        op = gate_database(2, 65, [0, word, 0, 0]).operation
        b = Builder("main", {"address": op.module.registers[0].type, "data": op.module.registers[1].type})
        b.call(op, address=b["address"], data=b["data"])
        program = b.finish().program()
        self.assertEqual(
            simulate(program, initial={"address": 1, "data": 0}).amplitudes,
            {(1, word): 1 + 0j},
        )

    def test_wide_qram_database_declares_and_queries_sparse_bank(self):
        op = qram_database(65, 70).operation
        b = Builder(
            "main",
            {"address": op.module.registers[0].type, "data": op.module.registers[1].type},
            {"bank": op.module.resources[0].type},
        )
        b.call(op, address=b["address"], data=b["data"], resources={"table": "bank"})
        program = b.finish().program()
        word = (1 << 69) + 5
        self.assertEqual(
            simulate(program, {"bank": {0: word}}, initial={"address": 0, "data": 0}).amplitudes,
            {(0, word): 1 + 0j},
        )

    def test_wide_affine_boolean_oracle_keeps_bit_semantics(self):
        secret = (1 << 64) + 0b101
        op = affine_boolean_oracle(66, secret).operation
        b = Builder("main", {"address": op.module.registers[0].type, "data": op.module.registers[1].type})
        b.call(op, address=b["address"], data=b["data"])
        program = b.finish().program()
        for address, expected in ((0, 0), ((1 << 64) + 0b101, 1), (0b100, 1)):
            with self.subTest(address=address):
                self.assertEqual(
                    simulate(program, initial={"address": address, "data": 0}).amplitudes,
                    {(address, expected): 1 + 0j},
                )

    def test_gradient_packs_coordinates_beyond_one_machine_word(self):
        oracle = function_phase_oracle(
            "def f(x0, x1, x2, x3, x4, x5, x6, x7):\n    return 0.25 * x0\n",
            dimension=8,
            grid_bits=9,
        )
        program = gradient_estimation(oracle, dimension=8, grid_bits=9).program()
        target = next(r for r in program.main.registers if r.name == "target")
        self.assertEqual(target.type.width, 72)
        loads(dumps(program, format="json"))

    def test_block_encoding_signal_merge_beyond_one_machine_word(self):
        # Lightweight stand-in for the large-plate norm probe: merging LCU and
        # product signals past 64 bits must build, validate, and serialize.
        wide = reduce(product, [zero(2)] * 32)
        merged = linear_combination(0.5, wide, 0.5j, reduce(product, [zero(2)] * 32))
        program = merged.operation.program()
        signal = next(r for r in program.main.registers if r.name == "signal")
        self.assertEqual(signal.type.width, 65)
        loads(dumps(program, format="json"))

    def test_wide_classical_helpers_accept_long_words(self):
        self.assertEqual(simon_nullspace([1 << i for i in range(1, 65)], 65), (1,))
        self.assertEqual(len(PauliOperator("X" * 65 + "Z").word), 66)


if __name__ == "__main__":
    unittest.main()
