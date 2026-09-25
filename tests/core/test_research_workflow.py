"""Boundary regressions for the algorithm research workflow: providers, partial binding, and compositional costs."""

import json
import unittest
from dataclasses import replace

from oracq import (
    AlgorithmContract,
    Binding,
    BindingError,
    Bits,
    Builder,
    ContractError,
    InputRequirement,
    ProtocolContract,
    QLSSProtocol,
    QLSSSolver,
    QODEProtocol,
    QODESolver,
    bind,
    bind_with_report,
    dumps,
    estimate_resources,
    identity,
    loads,
    simulate,
)
from oracq.algorithms.common.state_preparation import apply_be_to_state
from oracq.algorithms.input_model.interfaces import (
    BlockEncodingProtocol,
    as_block_encoding,
)
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import (
    abstract_database,
    abstract_state_prep,
    annotate,
    basis_state,
    declare,
    qram_database,
)


class ProviderBoundaryTests(unittest.TestCase):
    def test_aliases_refer_to_one_implementation(self):
        self.assertIs(ProtocolContract, AlgorithmContract)
        self.assertIs(QLSSProtocol, QLSSSolver)
        self.assertIs(QODEProtocol, QODESolver)

    def test_each_solver_uses_the_checked_view(self):
        class ChangingProvider:
            calls = 0

            def block_encoding(self):
                self.calls += 1
                return identity(self.calls)

        for solver in (
            QODESolver("qode", lambda a, b, t: apply_be_to_state(a, b)),
            QLSSSolver("qlss", "block_encoding", lambda _: None, apply_be_to_state),
        ):
            provider = ChangingProvider()
            args = (provider, basis_state(1))
            result = solver(*args, 0.1) if isinstance(solver, QODESolver) else solver(*args)
            self.assertEqual(provider.calls, 1)
            self.assertEqual(result.width, 1)

    def test_signature_and_noncallable_are_structured(self):
        class WrongSignature:
            def block_encoding(self, required):
                raise AssertionError("the method body must not be entered when the signature is incompatible")

        class NotCallable:
            block_encoding = 42

        for value, code in ((WrongSignature(), "INPUT_SIGNATURE"), (NotCallable(), "INPUT_CALLABLE")):
            with self.subTest(code=code), self.assertRaises(ContractError) as caught:
                as_block_encoding(value)
            self.assertEqual(caught.exception.issues[0].code, code)
            report = QODESolver("solver", lambda *_: None).check(value, basis_state(1))
            self.assertEqual(report.issues[0].code, code)
            self.assertEqual(report.issues[0].path, "solver.generator.block_encoding")

    def test_provider_internal_typeerror_is_preserved(self):
        failure = TypeError("internal error in the provider implementation")

        class BrokenProvider:
            def block_encoding(self):
                raise failure

        with self.assertRaises(TypeError) as caught:
            as_block_encoding(BrokenProvider())
        self.assertIs(caught.exception, failure)

    def test_resolution_checks_before_exposing_inputs(self):
        contract = AlgorithmContract("matrix", (
            InputRequirement("A", (BlockEncodingProtocol,), adapter=as_block_encoding, main_qubit=2),
        ))
        resolved = contract.resolve(A=identity(1))
        with self.assertRaises(ContractError):
            resolved.get("A", BlockEncoding)


def two_bank_program():
    """The same internal resource formal maps to two different entry banks."""
    first = abstract_database("FirstDB", 1, 1)
    second = abstract_database("SecondDB", 1, 1)
    builder = Builder("TwoBanks", {"a": Bits(1), "d": Bits(1)})
    builder.call(first.operation, address=builder["a"], data=builder["d"])
    with builder.repeat(3):
        builder.call(second.operation, address=builder["a"], data=builder["d"])
    implementation = qram_database(1, 1).operation
    return builder.finish().program(), {
        "FirstDB": Binding(implementation, {"table": "bank_a"}),
        "SecondDB": Binding(implementation, {"table": "bank_b"}),
    }


class BindingWorkflowTests(unittest.TestCase):
    def test_custom_role_does_not_require_registration(self):
        slot = declare("Custom", {"q": Bits(1)}, paradigm="research_gradient")
        impl = annotate(Builder("Actual", {"q": Bits(1)}).finish(), "research_gradient")
        self.assertTrue(bind_with_report(slot, {"Custom": impl}).report.ok)

    def test_partial_roundtrip_and_independent_order(self):
        program, bindings = two_bank_program()
        first = bind_with_report(loads(dumps(program)), {"FirstDB": bindings["FirstDB"]})
        self.assertTrue(first.report.ok)
        self.assertEqual([r.name for r in first.report.remaining], ["SecondDB"])
        json.dumps(first.report.to_dict())
        closed = bind(loads(dumps(first.require())), {"SecondDB": bindings["SecondDB"]})
        reverse = bind(bind(program, {"SecondDB": bindings["SecondDB"]}), {"FirstDB": bindings["FirstDB"]})
        self.assertEqual(dumps(closed), dumps(reverse))
        self.assertEqual(dumps(closed), dumps(bind(program, bindings)))

    def test_failed_binding_has_path_and_preserves_source(self):
        slot = abstract_state_prep("State", 1)
        b = Builder("Caller", {"target": Bits(1), "work": Bits(0)})
        b.call(slot.operation, target=b["target"], work=b["work"])
        program = b.finish().program()
        original = dumps(program)
        result = bind_with_report(program, {"State": basis_state(2).operation})
        self.assertFalse(result.report.ok)
        self.assertEqual(result.report.issues[0].code, "BIND_SIGNATURE")
        self.assertEqual(result.report.issues[0].path, ("Caller", "State"))
        self.assertEqual(dumps(program), original)
        with self.assertRaises(BindingError):
            result.require()

    def test_explicit_contradictory_promise_is_rejected(self):
        slot = abstract_state_prep("State", 1)
        wrong = annotate(basis_state(1).operation, "state_prep_isometry", clean_work=False)
        result = bind_with_report(slot.operation, {"State": wrong})
        self.assertEqual(result.report.issues[0].code, "BIND_PROMISE")

    def test_shared_bank_reused_after_partial_roundtrip(self):
        program, bindings = two_bank_program()
        first = bindings["FirstDB"]
        shared = Binding(first.operation, {"table": "shared"})
        partial = loads(dumps(bind(program, {"FirstDB": shared})))
        closed = bind(partial, {"SecondDB": shared})
        self.assertEqual([r.name for r in closed.main.resources], ["shared"])
        self.assertEqual(estimate_resources(closed).qram_queries, {"shared": 4})

    def test_invalid_capture_metadata_fails_roundtrip(self):
        program, bindings = two_bank_program()
        closed = bind(program, bindings)
        broken = replace(closed.main, attributes=(("binding_captures", '{"bank_a": "absent"}'),))
        bad_program = replace(closed, modules=tuple(broken if m.name == closed.entry else m for m in closed.modules))
        with self.assertRaisesRegex(ValueError, "nonexistent resource"):
            dumps(bad_program)


class ResourceWorkflowTests(unittest.TestCase):
    def test_arithmetic_realization_accepts_nonzero_xor_outputs(self):
        from oracq.applications.oracle_study import oracle_study

        _, variants = oracle_study(2)
        for name, (binding, memory) in variants.items():
            for address in range(4):
                for data in range(4):
                    with self.subTest(implementation=name, address=address, data=data):
                        bank = {"table": memory["angle_words"]} if memory else {}
                        actual = simulate(binding.operation.program(), bank, initial={"address": address, "data": data}).amplitudes
                        self.assertEqual(actual, {(address, data ^ ((address + 1) % 4)): 1 + 0j})

    def test_study_implementations_match_analytic_reference(self):
        from oracq.applications.oracle_study import oracle_study, oracle_study_reference

        for width in (1, 2, 3):
            opened, variants = oracle_study(width, 3)
            expected = oracle_study_reference(width, 3)
            for name, (implementation, memory) in variants.items():
                with self.subTest(width=width, implementation=name):
                    closed = bind(loads(dumps(opened)), {"AngleWord": implementation})
                    actual = simulate(closed, memory).amplitudes
                    self.assertLess(max(abs(actual.get(k, 0) - expected.get(k, 0)) for k in actual.keys() | expected.keys()), 1e-12)

    def test_large_rotating_repeat_stays_compact(self):
        b = Builder("HugeRotation", {"q": Bits(1)})
        with b.repeat(2**40):
            b.rz(b["q"], 0.37)
        estimate = estimate_resources(b.finish().program())
        self.assertEqual(estimate.rotations.total, 3 * 2**40)
        self.assertLessEqual(len(estimate.rotations.counts), 3)
        self.assertEqual(estimate.to_dict()["rotations"], 3 * 2**40)
        self.assertEqual(estimate.rotations[-1][0], "rz")

    def test_resource_arguments_are_counted_separately(self):
        program, bindings = two_bank_program()
        estimate = estimate_resources(bind(program, bindings))
        self.assertEqual(estimate.qram_queries, {"bank_a": 1, "bank_b": 3})

    def test_open_call_forms_and_bound_cost(self):
        slot = declare("U", {"q": Bits(1)})

        def build(operation):
            b = Builder("Algorithm", {"q": Bits(1), "c": Bits(2)})
            b.call(operation, q=b["q"])
            with b.adjoint(), b.repeat(3):
                b.call(operation, q=b["q"])
            with b.control(b["c"], 3), b.repeat(5):
                b.call(operation, q=b["q"])
            return b.finish().program()

        program = build(slot)
        estimate = estimate_resources(loads(dumps(program)), require_closed=False)
        self.assertFalse(estimate.complete)
        self.assertIsNone(estimate.qubits)
        self.assertEqual(estimate.qubits_lower_bound, 3)
        forms = {(call.controls, call.adjoint): count for call, count in estimate.oracle_calls.items()}
        self.assertEqual(forms, {(0, False): 1, (0, True): 3, (2, False): 5})
        implementation = Builder("ActualU", {"q": Bits(1)})
        implementation.x(implementation["q"])
        op = implementation.finish()
        actual = estimate_resources(bind(program, {"U": op}))
        direct = estimate_resources(build(op))
        self.assertTrue(actual.complete)
        self.assertEqual(actual.to_dict(), direct.to_dict())
