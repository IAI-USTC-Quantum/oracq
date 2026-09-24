"""算法私有结构协议、多角色操作与实际 Hamiltonian 组合。"""

import cmath
import math
import unittest
from typing import Protocol, runtime_checkable

from oracq import (
    Bits,
    BlockSystem,
    Builder,
    SpectralPromise,
    ValidationError,
    identity,
    requires,
    simulate,
)
from oracq.algorithms.common.hamiltonian import (
    EncodedOperator,
    PauliHamiltonian,
    PauliOperator,
    TrotterizableProtocol,
    TrotterTerm,
    hamiltonian_simulation,
)
from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.interfaces import (
    BlockEncodingProtocol,
    StatePreparationProtocol,
    UnitaryProtocol,
)
from oracq.algorithms.input_model.oracles import basis_state
from oracq.algorithms.qode.ode import linear_qode


class AlgorithmProtocolTests(unittest.TestCase):
    def test_one_gate_satisfies_multiple_protocols(self):
        b = Builder("OneGate", {"q": Bits(1)})
        b.x(b["q"])
        gate = b.finish()
        for interface in (UnitaryProtocol, StatePreparationProtocol, BlockEncodingProtocol):
            self.assertIs(requires(gate, interface), gate)
        problem = BlockSystem(identity(1), gate, SpectralPromise(1, 1))
        self.assertEqual(simulate(problem.rhs.operation.program()).amplitudes, {(1, 0): 1 + 0j})
        encoded = lcu(((1, gate), (1, identity(1))))
        state = simulate(encoded.operation.program())
        for output in (0, 1):
            self.assertAlmostEqual(state.amplitudes.get((output, 0), 0), 0.5)

    def test_full_unitary_register_layout_has_no_hidden_clean_work(self):
        b = Builder("TwoRegisterUnitary", {"a": Bits(1), "scratch": Bits(1)})
        b.x(b["scratch"])
        state = b.finish().state_preparation()
        self.assertEqual(state.width, 2)
        self.assertEqual(simulate(state.operation.program()).amplitudes, {(2, 0): 1 + 0j})

    def test_user_can_define_protocol_and_adapter_without_registration(self):
        @runtime_checkable
        class HasDiagonal(Protocol):
            def diagonal_values(self): ...

        class UserOperator:
            def diagonal_values(self):
                return (1, 1)

            def block_encoding(self):
                return identity(1)

        a = UserOperator()
        self.assertEqual(requires(a, HasDiagonal).diagonal_values(), (1, 1))
        self.assertTrue(linear_qode("schrodingerization").check(a, basis_state(1)).ok)

    def test_method_presence_is_not_a_valid_return_contract(self):
        class BrokenOperator:
            def block_encoding(self):
                return "not an encoding"

        report = linear_qode("lchs").check(BrokenOperator(), basis_state(1))
        self.assertFalse(report.ok)
        self.assertEqual(report.issues[0].code, "INPUT_ADAPTER")

    def test_trotter_protocol_keeps_phase_and_repeat(self):
        h = PauliHamiltonian(((0.3, "I"), (0.7, "X")))
        self.assertIsInstance(h, TrotterizableProtocol)
        result = hamiltonian_simulation(h, 0.4, steps=3)
        state = simulate(result.operation.program())
        phase = cmath.exp(-0.12j)
        self.assertAlmostEqual(state.amplitudes[(0, 0)], phase * math.cos(0.28), places=11)
        self.assertAlmostEqual(state.amplitudes[(1, 0)], -1j * phase * math.sin(0.28), places=11)
        self.assertEqual(result.alpha, 1)
        self.assertEqual(result.operation.module.body[0].count, 3)

    def test_trotter_only_input_does_not_need_block_encoding(self):
        class UserHamiltonian:
            hermitian = True

            def trotter_list(self):
                return (TrotterTerm(0.5, PauliOperator("Z")),)

        h = UserHamiltonian()
        self.assertNotIsInstance(h, BlockEncodingProtocol)
        result = hamiltonian_simulation(h, 0.2)
        self.assertAlmostEqual(
            simulate(result.operation.program()).amplitudes[(0, 0)], cmath.exp(-0.1j)
        )

    def test_nonhermitian_and_missing_qsp_are_distinct_errors(self):
        with self.assertRaisesRegex(ValidationError, "Hermitian"):
            hamiltonian_simulation(EncodedOperator(identity(1), hermitian=False), 0.1)
        with self.assertRaisesRegex(ValidationError, "实际 qsp"):
            hamiltonian_simulation(EncodedOperator(identity(1), hermitian=True), 0.1)

    def test_trotter_validates_each_term(self):
        with self.assertRaises(ValidationError):
            TrotterTerm(1, object())
        with self.assertRaises(ValidationError):
            PauliHamiltonian(((1, "X"), (1, "ZZ")))
