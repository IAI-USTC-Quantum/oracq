"""Boolean lowering witnesses independent of the network's native fast path."""

import random
import unittest
from collections import Counter

from oracq import Bits, Builder, ValidationError, estimate_resources, export_strict, simulate
from oracq.algorithms.common.arithmetic import BooleanNetwork, FixedFormat, fixed_arithmetic
from oracq.infrastructure.boolean_synthesis import lower_affine
from oracq.infrastructure.mathfunc import compile_function
from oracq.infrastructure.mathfunc.numeric import MathConfig, NumericEmitter


class BooleanSynthesisTests(unittest.TestCase):
    def test_specialization_matches_generic_value_and_status(self):
        for kind in ("add", "sub", "mul", "div", "lt", "eq"):
            original = fixed_arithmetic(kind, FixedFormat(3, 1))
            net = BooleanNetwork.from_payload(dict(original.module.attributes)["arithmetic_network"])
            for constant in range(8):
                for port in ("a", "b"):
                    op = fixed_arithmetic(kind, FixedFormat(3, 1), constants=((port, constant),))
                    self.assertNotIn(port, {r.name for r in op.module.registers})
                    other = "a" if port == "b" else "b"
                    for word in range(8):
                        values = net.evaluate(**{port: constant, other: word})
                        actual = simulate(op.program(), initial={other: word, "out": 1, "status": 3})
                        with self.subTest(kind=kind, port=port, constant=constant, word=word):
                            self.assertEqual(actual.amplitudes, {(word, values["out"] ^ 1, values["status"] ^ 3): 1 + 0j})

    def test_invalid_constant_binding_is_rejected(self):
        net = BooleanNetwork()
        net.outputs = {"out": net.input("a", 2)}
        before = net.payload()
        for constants in ({"a": -1}, {"a": 4}, {"a": True}, {"missing": 0}):
            with self.assertRaises(ValidationError):
                net.specialize(**constants)
            self.assertEqual(net.payload(), before)

    def test_constant_folding_keeps_overflow_and_division_domain_flags(self):
        cases = [
            ("def f(x):\n return x*(3+3)", lambda raw: (-2 * raw) & 15, 2),
        ]
        for source, function, flag in cases:
            op = compile_function(source, fmt=FixedFormat(4, 1)).operation
            for raw in (0, 1, 2):
                result = simulate(op.program(), initial={"x": raw, "out": 5, "status": 3})
                with self.subTest(source=source, raw=raw):
                    self.assertEqual(result.amplitudes, {(raw, function(raw) ^ 5, flag ^ 3): 1 + 0j})
        builder = Builder("constant_domain", {"out": Bits(4), "status": Bits(2)})
        emitter = NumericEmitter(builder, FixedFormat(4, 1), MathConfig())
        value = emitter.arithmetic("div", emitter.constant(1), emitter.constant(0))
        operation = emitter.finish([((builder["out"],), value)], builder["status"])
        self.assertEqual(simulate(operation.program(), initial={"out": 5, "status": 3}).amplitudes, {(5, 2): 1 + 0j})

    def check_basis(self, net, inputs, output):
        operation = lower_affine(net)
        expected = net.evaluate(**inputs)["out"] ^ output
        initial = {**inputs, "out": output}
        result = simulate(operation.program(), initial=initial)
        self.assertEqual(result.amplitudes, {tuple([*inputs.values(), expected]): 1 + 0j})

    def test_linear_function_uses_no_workspace(self):
        net = BooleanNetwork()
        a, b = net.input("a", 3), net.input("b", 3)
        net.outputs = {"out": [net.inv(net.xor(x, y)) for x, y in zip(a, b, strict=True)]}
        self.assertEqual(estimate_resources(lower_affine(net).program()).qubits, 9)
        for x in range(8):
            for y in range(8):
                self.check_basis(net, {"a": x, "b": y}, 5)

    def test_correlated_affine_controls_and_folded_products(self):
        net = BooleanNetwork()
        a, b, c = net.input("a", 3)
        left, right = net.xor(a, b), net.inv(net.xor(b, c))
        product = net.and_(left, right)
        redundant = net.and_(left, net.xor(net.inv(a), b))
        net.outputs = {"out": [product, redundant, net.and_(left, net.xor(a, b))]}
        for word in range(8):
            for output in (0, 3, 7):
                self.check_basis(net, {"a": word}, output)

    def test_dead_nonlinear_nodes_are_not_emitted(self):
        net = BooleanNetwork()
        a, b, c = net.input("a", 3)
        net.and_(net.and_(a, b), c)
        net.outputs = {"out": [net.xor(a, b)]}
        cost = estimate_resources(lower_affine(net).program())
        self.assertEqual(cost.qubits, 4)
        self.assertEqual(cost.toffoli, 0)

    def test_generated_graphs_with_nonzero_outputs(self):
        rng = random.Random(20261108)
        for _ in range(12):
            net = BooleanNetwork()
            nodes = [0, 1, *net.input("a", 3)]
            for _ in range(20):
                left, right = rng.choices(nodes, k=2)
                nodes.append(rng.choice([net.xor, net.and_, net.or_])(left, right))
                if rng.randrange(3) == 0:
                    nodes.append(net.inv(nodes[-1]))
            net.outputs = {"out": rng.choices(nodes, k=3)}
            for word in range(8):
                self.check_basis(net, {"a": word}, rng.randrange(8))

    def test_fixed_point_values_and_status_on_all_small_inputs(self):
        for kind in ("add", "sub", "mul", "div", "reciprocal", "sqrt", "lt", "eq", "select"):
            operation = fixed_arithmetic(kind, FixedFormat(3, 1))
            network = BooleanNetwork.from_payload(dict(operation.module.attributes)["arithmetic_network"])
            rng = random.Random(0)
            for a in range(8):
                for b in range(8) if "b" in network.inputs else (0,):
                    inputs = {"a": a}
                    if "b" in network.inputs:
                        inputs["b"] = b
                    if "select" in network.inputs:
                        inputs["select"] = a % 2
                    initial = {r.name: rng.randrange(1 << r.type.width) for r in operation.module.registers}
                    initial.update(inputs)
                    expected = {**initial}
                    for key, value in network.evaluate(**inputs).items():
                        expected[key] ^= value
                    with self.subTest(kind=kind, inputs=inputs):
                        actual = simulate(operation.program(), initial=initial)
                        self.assertEqual(actual.amplitudes, {tuple(expected[r.name] for r in operation.module.registers): 1 + 0j})

    def test_resource_improvement_against_reference_ssa(self):
        for kind in ("add", "mul", "div", "reciprocal", "sqrt"):
            operation = fixed_arithmetic(kind, FixedFormat(6, 2))
            network = BooleanNetwork.from_payload(dict(operation.module.attributes)["arithmetic_network"])
            before = estimate_resources(network.operation(optimize=False).program())
            after = estimate_resources(operation.program())
            with self.subTest(kind=kind):
                self.assertLess(after.qubits, before.qubits)
                self.assertLess(after.toffoli, before.toffoli)
                self.assertLess(after.clifford, before.clifford)

    def test_superposition_negative_control_adjoint_and_clean_workspace(self):
        net = BooleanNetwork()
        a, b, c = net.input("a", 3)
        net.outputs = {"out": [net.and_(net.inv(net.xor(a, b)), net.xor(b, c))]}
        op = lower_affine(net)
        builder = Builder("coherent", {"a": Bits(3), "out": Bits(1), "control": Bits(1)})
        builder.h(builder["a"])
        builder.h(builder["control"])
        builder.x(builder["out"])
        with builder.control(builder["control"], 0):
            builder.call(op, a=builder["a"], out=builder["out"])
        program = builder.finish().program()
        state = simulate(program)
        for key, amplitude in state.amplitudes.items():
            word, output, control = key
            self.assertEqual(output, 1 ^ (net.evaluate(a=word)["out"] if control == 0 else 0))
            self.assertAlmostEqual(amplitude.real, 0.25)
            self.assertAlmostEqual(amplitude.imag, 0)
        self.assertEqual(len(state.amplitudes), 16)
        builder = Builder("inverse", {"a": Bits(3), "out": Bits(1)})
        builder.h(builder["a"])
        builder.x(builder["out"])
        builder.call(op, a=builder["a"], out=builder["out"])
        with builder.adjoint():
            builder.call(op, a=builder["a"], out=builder["out"])
        self.assertEqual({key[1] for key in simulate(builder.finish().program()).amplitudes}, {1})

    def test_emitted_strict_netlist_matches_counts_and_space(self):
        for kind in ("add", "mul", "div", "sqrt"):
            operation = fixed_arithmetic(kind, FixedFormat(3, 1))
            estimate = estimate_resources(operation.program())
            artifact = export_strict(operation.program())
            atoms = Counter(
                line.split()[0].lower() for line in artifact.text.splitlines()
                if line.split()[0] in {"TOFFOLI", "CZ", "H", "X", "Y", "Z", "S", "SDG", "T", "TDG"}
            )
            with self.subTest(kind=kind):
                self.assertEqual(atoms, estimate.atoms)
                self.assertEqual(int(artifact.text.splitlines()[0].split()[1]), estimate.qubits + estimate.mcx_ancilla)

    def test_alias_pool_repeated_squares_preserve_values_and_flags(self):
        compiled = compile_function(
            "def repeated_square(x):\n y=x\n for i in range(6):\n  y=y*y\n return y",
            fmt=FixedFormat(4, 1),
        )
        builder = Builder("reused_alias", {r.name: r.type for r in compiled.operation.module.registers})
        builder.h(builder["x"])
        builder.x(builder["out"])
        builder.x(builder["status"][0])
        builder.call(compiled.operation, x=builder["x"], out=builder["out"], status=builder["status"])
        state = simulate(builder.finish().program())
        self.assertEqual(len(state.amplitudes), 16)
        for (word, output, flags), amplitude in state.amplitudes.items():
            raw, status = word, 0
            for _ in range(6):
                signed = raw if raw < 8 else raw - 16
                magnitude = (abs(signed) ** 2) >> 1
                status |= 2 * bool(magnitude >> 3)
                raw = magnitude & 15
            self.assertEqual((output, flags), (raw ^ 15, status ^ 1))
            self.assertAlmostEqual(amplitude.real, 0.25)
            self.assertAlmostEqual(amplitude.imag, 0)


if __name__ == "__main__":
    unittest.main()
