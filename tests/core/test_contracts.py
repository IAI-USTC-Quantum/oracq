"""Contract acceptance: errors must be located before running the algorithm kernel; open implementations stay composable."""

import json
import unittest
from dataclasses import replace

from oracq import (
    Bits,
    BlockSystem,
    Builder,
    ContractError,
    InputRequirement,
    LinearSystem,
    Operation,
    ProtocolContract,
    QLSSProtocol,
    QODEProblem,
    QODEProtocol,
    SpectralPromise,
    ValidationError,
    bind,
    describe_oracle,
    dumps,
    identity,
    loads,
    scale,
    simulate,
)
from oracq.algorithms.common.state_preparation import apply_be_to_state
from oracq.algorithms.input_model.interfaces import (
    BlockEncodingProtocol,
    StatePreparationProtocol,
)
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    abstract_block_encoding,
    abstract_database,
    abstract_sparse_access,
    abstract_state_prep,
    annotate,
    basis_state,
    declare,
)
from oracq.algorithms.qlss.qlss import CostaConfig, make_costa_qlss
from oracq.algorithms.qnlss.carleman import PolynomialODE
from oracq.algorithms.qode.cbmd import ContourPlan
from oracq.algorithms.qode.lchs import QuadraturePlan
from oracq.algorithms.qode.ode import linear_qode
from oracq.algorithms.qode.schrodingerization import SchrodingerPlan


class OracleContractTests(unittest.TestCase):
    def test_descriptors_come_from_ir_and_survive_roundtrip(self):
        a = abstract_block_encoding("GivenA", 3, 2, 7.0)
        self.assertEqual(
            (a.type, a.main_qubit, a.anc_qubit, a.alpha), ("block_encoding", 3, 2, 7.0)
        )
        self.assertEqual(a.spec.implementation, "open")
        restored = type(a)(Operation.from_program(loads(dumps(a.operation.program()))))
        self.assertEqual(a.spec, restored.spec)
        self.assertEqual(json.loads(json.dumps(a.spec.to_dict()))["alpha"], 7)
        self.assertEqual(identity(1).anc_qubit, 0)

    def test_sparse_is_a_bundle_not_a_fictitious_unitary(self):
        a = abstract_sparse_access("Sparse", 3, 8, 4, work_width=2)
        self.assertEqual(a.type, "cks_sparse")
        self.assertIsNone(a.anc_qubit)
        self.assertEqual(dict(a.spec.parameters)["sparsity"], 4)
        self.assertEqual(dict(a.spec.components)["position"].anc_qubit, 2)
        self.assertEqual(dict(a.spec.components)["entry"].type, "sparse_entry_xor")

    def test_aggregate_report_and_kernel_not_executed(self):
        calls = []

        def kernel(*inputs):
            calls.append(inputs)
            return apply_be_to_state(inputs[0], inputs[1])

        solver = QODEProtocol("custom", kernel)
        a = abstract_database("WrongA", 2, 3)
        b = abstract_state_prep("LimitedB", 1, reversible=False)
        report = solver.check(a, b)
        self.assertFalse(report.ok)
        self.assertTrue({"INPUT_PROTOCOL", "INPUT_CAPABILITY"} <= {i.code for i in report.issues})
        with self.assertRaises(ContractError):
            solver(a, b, 0.1)
        self.assertEqual(calls, [])

    def test_explicit_unitary_state_prep_adapter(self):
        b = Builder("HadamardInput", {"q": Bits(1)})
        b.h(b["q"])
        unitary = b.finish()
        prep = StatePreparation.from_unitary(unitary, target="q")
        self.assertTrue(prep.capabilities.adjoint)
        self.assertTrue(dict(prep.spec.parameters)["clean_work"])
        state = simulate(prep.operation.program())
        self.assertAlmostEqual(state.amplitudes[(0, 0)], 2**-0.5)
        self.assertAlmostEqual(state.amplitudes[(1, 0)], 2**-0.5)
        req = ProtocolContract(
            "state_consumer", (InputRequirement("b", (StatePreparationProtocol,)),)
        )
        self.assertTrue(req.check(b=unitary).ok)
        self.assertTrue(req.check(b=prep).ok)

    def test_unitary_adapter_requires_work_promise(self):
        op = Builder("WithWork", {"q": Bits(1), "tmp": Bits(2)}).finish()
        with self.assertRaisesRegex(ContractError, "clean_work"):
            StatePreparation.from_unitary(op, target="q", work="tmp")
        prep = StatePreparation.from_unitary(op, target="q", work="tmp", clean_work=True)
        self.assertEqual(prep.anc_qubit, 2)

    def test_capability_restrictions_survive_annotation_and_wrapping(self):
        limited = declare(
            "ForwardOnly", {"q": Bits(1)}, supports_adjoint=False, supports_controlled=False
        )
        tagged = annotate(limited, "unitary", description="forward only")
        self.assertFalse(describe_oracle(tagged).capabilities.adjoint)
        with self.assertRaises(ContractError):
            annotate(limited, "unitary", supports_adjoint=True)
        prep = StatePreparation.from_unitary(limited, target="q")
        self.assertFalse(prep.capabilities.adjoint)
        self.assertFalse(prep.capabilities.controlled)

    def test_dirty_work_is_rejected_by_algorithm(self):
        bad = StatePreparation(
            annotate(
                basis_state(1, work_width=1).operation, "state_prep_isometry", clean_work=False
            )
        )
        report = linear_qode("lchs").check(identity(1), bad)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.path.endswith("clean_work") for i in report.issues))

    def test_changed_alpha_requires_regeneration(self):
        slot = abstract_block_encoding("A", 1, 0, 2.0)
        with self.assertRaisesRegex(ValidationError, "be_alpha"):
            bind(slot.operation, {"A": identity(1).operation})
        p = bind(slot.operation, {"A": scale(2, identity(1)).operation})
        self.assertEqual(describe_oracle(Operation.from_program(p)).implementation, "closed")

    def test_wrong_wrappers_fail_with_validation_error(self):
        for cls in (StatePreparation, StateOracle):
            with self.subTest(cls=cls.__name__), self.assertRaises(ValidationError):
                cls(object())
        with self.assertRaises(ValidationError):
            BlockSystem(object(), basis_state(1), SpectralPromise(1, 1))

    def test_contract_reports_missing_extra_and_exact_layout(self):
        contract = ProtocolContract(
            "example",
            (InputRequirement("A", (BlockEncodingProtocol,), main_qubit=2, anc_qubit=1, alpha=4),),
        )
        report = contract.check(A=identity(1), unused=basis_state(1))
        self.assertEqual(len(report.issues), 4)
        self.assertEqual(contract.check().issues[0].code, "INPUT_MISSING")

    def test_qlss_has_inspectable_requirements(self):
        a, b = identity(1), abstract_state_prep("OnlyForward", 1, reversible=False)
        problem = LinearSystem(block=BlockSystem(a, b, SpectralPromise(1, 1)))
        solver = make_costa_qlss(CostaConfig(steps=1))
        self.assertIn(BlockEncodingProtocol, solver.contract.inputs[0].protocols)
        self.assertFalse(solver.check(problem).ok)
        self.assertTrue(any(i.path.endswith("b.adjoint") for i in solver.check(problem).issues))

    def test_problem_level_qode_enforces_math_declaration(self):
        calls = []

        def kernel(a, b, t):
            calls.append(t)
            return apply_be_to_state(a, b)

        solver = QODEProtocol("lchs_contract_witness", kernel, requires_dissipative=True)
        p = QODEProblem(scale(-1, identity(1)), basis_state(1))
        with self.assertRaisesRegex(ContractError, "dissipative"):
            solver.solve(p, 0.1)
        self.assertFalse(calls)
        out = solver.solve(replace(p, dissipative=True, initial_norm=2.0), 0.1)
        self.assertEqual(calls, [0.1])
        self.assertEqual(dict(out.operation.module.attributes)["qode_initial_norm"], 2)
        self.assertTrue(dict(out.operation.module.attributes)["qode_dissipative_promise"])

    def test_output_contract_is_checked(self):
        solver = QODEProtocol("bad_result", lambda *_: basis_state(1))
        with self.assertRaisesRegex(ContractError, "output"):
            solver(identity(1), basis_state(1), 0.1)
        solver = QODEProtocol(
            "bad_width", lambda *_: apply_be_to_state(identity(2), basis_state(2))
        )
        with self.assertRaisesRegex(ContractError, "OUTPUT_LAYOUT"):
            solver(identity(1), basis_state(1), 0.1)
        invalid_qlss = QLSSProtocol("bad_qlss", "block_encoding", lambda _: basis_state(1))
        with self.assertRaisesRegex(ContractError, "output"):
            invalid_qlss.solve(
                LinearSystem(block=BlockSystem(identity(1), basis_state(1), SpectralPromise(1, 1)))
            )

    def test_configs_reject_nan_bool_fractional_counts(self):
        invalid = [
            lambda: CostaConfig(steps=True),
            lambda: CostaConfig(kappa=float("nan")),
            lambda: QuadraturePlan((0,), (complex(float("nan"), 0),)),
            lambda: QuadraturePlan(("bad",), (1,)),
            lambda: QuadraturePlan((0,), (0,)),
            lambda: QuadraturePlan.cauchy(cutoff=1.2),
            lambda: SchrodingerPlan(auxiliary_width=1.5),
            lambda: SchrodingerPlan(period=float("nan")),
            lambda: SchrodingerPlan(selected_index=True),
            lambda: ContourPlan(cutoff=True),
            lambda: ContourPlan(poles=(complex(0, float("inf")),)),
            lambda: SpectralPromise(True, 0.1),
            lambda: PolynomialODE(1, ((1, identity(1)), (1, identity(1))), basis_state(1)),
            lambda: PolynomialODE(1, ((1, identity(1)),), basis_state(1), float("nan")),
            lambda: linear_qode("lchs", bogus=1),
            lambda: linear_qode("lchs", plan=SchrodingerPlan()),
        ]
        for case in invalid:
            with self.subTest(case=invalid.index(case)), self.assertRaises(ValidationError):
                case()

    def test_legacy_state_prep_clean_work_in_bind(self):
        prep = basis_state(1)
        op = replace(
            prep.operation,
            module=replace(
                prep.operation.module,
                attributes=tuple(
                    (k, v) for k, v in prep.operation.module.attributes if k != "clean_work"
                ),
            ),
        )
        bound = bind(abstract_state_prep("B", 1).operation, {"B": op})
        self.assertEqual(describe_oracle(Operation.from_program(bound)).implementation, "closed")
