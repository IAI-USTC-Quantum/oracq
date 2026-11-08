"""Square circuits checked against integer arithmetic and generic multiplication."""
import unittest
from collections import Counter

from oracq import Builder, estimate_resources, export_strict, simulate
from oracq.algorithms.common.arithmetic import BooleanNetwork, FixedFormat, fixed_arithmetic
from oracq.infrastructure.mathfunc import compile_function


class SquareArithmeticTests(unittest.TestCase):
    def test_unsigned_bit_products_and_truncation_exhaustively(self):
        for bits in range(1, 8):
            for width in range(2*bits+3):
                net = BooleanNetwork()
                net.outputs = {"out": net.square(net.input("a", bits), width)}
                for word in range(1 << bits):
                    self.assertEqual(net.evaluate(a=word)["out"], word*word & ((1 << width)-1))

    def test_all_small_fixed_formats_match_mul_and_independent_integer_reference(self):
        for width in range(2, 7):
            for signed in (False, True):
                for fraction in range(width-int(signed)):
                    fmt = FixedFormat(width, fraction, signed)
                    square = fixed_arithmetic("square", fmt)
                    multiply = fixed_arithmetic("mul", fmt)
                    nets = [BooleanNetwork.from_payload(dict(op.module.attributes)["arithmetic_network"])
                            for op in (square, multiply)]
                    for word in range(1 << width):
                        value = word-(1 << width) if signed and word >> (width-1) else word
                        raw = value*value >> fraction
                        expected = {"out": raw & ((1 << width)-1),
                                    "status": 2*bool(raw >> (width-int(signed)))}
                        self.assertEqual(nets[0].evaluate(a=word), expected)
                        self.assertEqual(nets[1].evaluate(a=word, b=word), expected)

    def test_compiled_superposition_preserves_phase_xor_flags_and_clean_dagger(self):
        op = compile_function("def f(x):\n return x*x", fmt=FixedFormat(3, 1)).operation
        b = Builder("square_coherent", {r.name: r.type for r in op.module.registers})
        b.h(b["x"])
        b.x(b["out"])
        b.x(b["status"][0])
        b.call(op, x=b["x"], out=b["out"], status=b["status"])
        state = simulate(b.finish().program())
        self.assertEqual(len(state.amplitudes), 8)
        for (word, output, status), amplitude in state.amplitudes.items():
            value = word if word < 4 else word-8
            raw = value*value >> 1
            self.assertEqual((output, status), ((raw & 7) ^ 7, (2*bool(raw >> 2)) ^ 1))
            self.assertAlmostEqual(amplitude, 8**-0.5)
        inverse = Builder("square_inverse", {r.name: r.type for r in op.module.registers})
        inverse.h(inverse["x"])
        inverse.x(inverse["out"])
        inverse.x(inverse["status"][0])
        args = {r.name: inverse[r.name] for r in op.module.registers}
        inverse.call(op, **args)
        with inverse.adjoint():
            inverse.call(op, **args)
        self.assertEqual({key[1:] for key in simulate(inverse.finish().program()).amplitudes}, {(7, 1)})

    def test_square_constant_specialization_and_actual_export_ledger(self):
        fmt = FixedFormat(4, 1)
        for word in range(16):
            op = fixed_arithmetic("square", fmt, constants=(("a", word),))
            net = BooleanNetwork.from_payload(dict(op.module.attributes)["arithmetic_network"])
            value = word if word < 8 else word-16
            self.assertEqual(net.evaluate()["out"], (value*value >> 1) & 15)
        for width, fraction in ((3, 1), (8, 3), (16, 8)):
            fmt = FixedFormat(width, fraction)
            op = fixed_arithmetic("square", fmt)
            old = estimate_resources(fixed_arithmetic("mul", fmt).program())
            ledger = estimate_resources(op.program())
            artifact = export_strict(op.program())
            atoms = Counter(line.split()[0].lower() for line in artifact.text.splitlines()
                            if line.split()[0] in {"TOFFOLI", "CZ", "H", "X", "Y", "Z", "S", "SDG", "T", "TDG"})
            self.assertEqual(atoms, ledger.atoms)
            self.assertEqual(int(artifact.text.splitlines()[0].split()[1]), ledger.qubits+ledger.mcx_ancilla)
            self.assertLess(ledger.toffoli, old.toffoli)
            self.assertLess(ledger.qubits, old.qubits)
