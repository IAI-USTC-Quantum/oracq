"""Numerical witnesses for the PREPARE–SELECT standard decomposition and alias-sampling PREPARE."""

import cmath
import math
import unittest

from witness import assert_block_equals

from oracq import ValidationError, bind, simulate, unresolved
from oracq.algorithms.common.hamiltonian import PauliHamiltonian
from oracq.algorithms.common.prepare_select import (
    abstract_prepare,
    alias_prepare,
    alias_table,
    gate_prepare,
    lcu_prepare_select,
    qram_prepare,
    select_pauli,
)
from oracq.algorithms.common.transforms import qubitization_walk
from oracq.algorithms.input_model.oracles import gate_database
from oracq.infrastructure.linking import Binding

TERMS = ((0.6, "XZ"), (-0.8, "ZI"), (0.3j, "YY"), (-0.5, "IX"))
COEFFICIENTS = tuple(c for c, _ in TERMS)
ALPHA = sum(abs(c) for c in COEFFICIENTS)


def dense_hamiltonian(terms):
    """Expand Σ c_i P_i following the repository bit convention (bit i acts on target[i], little-endian)."""
    size = 1 << len(terms[0][1])
    matrix = [[0j] * size for _ in range(size)]
    for coefficient, word in terms:
        for column in range(size):
            row, phase = column, 1 + 0j
            for bit, letter in enumerate(word):
                value = (column >> bit) & 1
                if letter in "XY":
                    row ^= 1 << bit
                if letter == "Y":
                    phase *= -1j if value else 1j
                elif letter == "Z":
                    phase *= -1 if value else 1
            matrix[row][column] += coefficient * phase
    return matrix


def marginal(state, index):
    result = {}
    for key, amplitude in state.amplitudes.items():
        result[key[index]] = result.get(key[index], 0) + abs(amplitude) ** 2
    return result


class PrepareSelectTests(unittest.TestCase):
    def test_gate_and_qram_prepare_bindings_agree(self):
        expected = [math.sqrt(abs(c) / ALPHA) for c in COEFFICIENTS]
        gate = gate_prepare(COEFFICIENTS)
        state = simulate(gate.operation.program())
        for index, amplitude in enumerate(expected):
            self.assertAlmostEqual(state.amplitudes.get((index, 0), 0), amplitude, places=11)
        slot = abstract_prepare(COEFFICIENTS)
        name = slot.operation.module.name
        bound = bind(slot.operation, {name: gate.operation})
        self.assertFalse(unresolved(bound))
        self.assertEqual(
            simulate(bound).amplitudes, simulate(gate.operation.program()).amplitudes
        )
        qram = qram_prepare(COEFFICIENTS, angle_width=10)
        slot = abstract_prepare(COEFFICIENTS, work_width=qram.preparation.work_width)
        bound = bind(
            slot.operation,
            {slot.operation.module.name: Binding(qram.preparation.operation, {"angles": "angles"})},
        )
        state = simulate(bound, qram.memory)
        self.assertTrue(all(key[1] == 0 for key in state.amplitudes))
        for index, amplitude in enumerate(expected):
            self.assertAlmostEqual(
                abs(state.amplitudes.get((index, 0), 0)), amplitude, delta=0.02
            )

    def test_select_applies_indexed_pauli_word_with_phase(self):
        select = select_pauli(TERMS)
        for index, (coefficient, word) in enumerate(TERMS):
            dense = dense_hamiltonian([(1, word)])
            phase = cmath.exp(1j * cmath.phase(coefficient))
            for column in range(4):
                state = simulate(
                    select.program(), initial={"selector": index, "target": column}
                )
                for row in range(4):
                    self.assertAlmostEqual(
                        state.amplitudes.get((index, row), 0),
                        dense[row][column] * phase,
                        places=11,
                    )

    def test_block_encoding_recovers_hamiltonian_over_alpha(self):
        be = lcu_prepare_select(TERMS)
        attributes = dict(be.operation.module.attributes)
        self.assertAlmostEqual(be.alpha, ALPHA)
        self.assertEqual(attributes["be_form"], "prepare_select")
        self.assertEqual(attributes["lcu_terms"], 4)
        self.assertEqual(attributes["selector_width"], 2)
        self.assertEqual(be.signal_qubits, 2)
        matrix = dense_hamiltonian(TERMS)
        assert_block_equals(self, be, matrix, places=10)

    def test_abstract_prepare_binds_inside_block_encoding(self):
        slot = abstract_prepare(COEFFICIENTS)
        be = lcu_prepare_select(TERMS, prepare=slot)
        self.assertEqual(
            {r.name for r in unresolved(be.operation.program())},
            {slot.operation.module.name},
        )
        bound = bind(
            be.operation.program(),
            {slot.operation.module.name: gate_prepare(COEFFICIENTS).operation},
        )
        matrix = dense_hamiltonian(TERMS)
        for column in range(4):
            state = simulate(bound, initial={"target": column})
            for row in range(4):
                self.assertAlmostEqual(
                    state.amplitudes.get((row, 0), 0) * be.alpha,
                    matrix[row][column],
                    places=10,
                )

    def test_alias_table_reproduces_normalized_distribution(self):
        table = alias_table(COEFFICIENTS, precision=8)
        size = 1 << 2
        for index, probability in enumerate(table.probabilities):
            exact = (
                table.keep[index]
                + sum(1 - table.keep[j] for j in range(size) if table.alt[j] == index)
            ) / size
            self.assertAlmostEqual(exact, probability, places=12)
            self.assertAlmostEqual(table.distribution()[index], probability, delta=0.02)
        for word in table.table.values():
            self.assertLess(word, 1 << 10)

    def test_alias_preparation_samples_target_distribution(self):
        table = alias_table(COEFFICIENTS, precision=8)
        database = gate_database(2, 10, dict(table.table))
        alias = alias_prepare(COEFFICIENTS, precision=8, database=database)
        self.assertEqual(alias.table, table)
        attributes = dict(alias.preparation.operation.module.attributes)
        self.assertEqual(attributes["implementation"], "alias_sampling")
        self.assertFalse(attributes["clean_work"])
        state = simulate(alias.preparation.operation.program())
        for index, probability in enumerate(table.distribution()):
            self.assertAlmostEqual(marginal(state, 0).get(index, 0), probability, places=10)
        qram_alias = alias_prepare(COEFFICIENTS, precision=8)
        state = simulate(qram_alias.preparation.operation.program(), qram_alias.memory)
        for index, probability in enumerate(table.distribution()):
            self.assertAlmostEqual(marginal(state, 0).get(index, 0), probability, places=10)

    def test_alias_block_encoding_approaches_hamiltonian(self):
        table = alias_table(COEFFICIENTS, precision=8)
        alias = alias_prepare(
            COEFFICIENTS, precision=8, database=gate_database(2, 10, dict(table.table))
        )
        be = lcu_prepare_select(TERMS, prepare=alias)
        self.assertEqual(be.signal_qubits, 2 + alias.preparation.work_width)
        matrix = dense_hamiltonian(TERMS)
        for column in range(4):
            state = simulate(be.operation.program(), initial={"target": column})
            for row in range(4):
                self.assertAlmostEqual(
                    state.amplitudes.get((row, 0), 0) * be.alpha,
                    matrix[row][column],
                    delta=0.02,
                )

    def test_qram_prepare_block_encoding_approaches_hamiltonian(self):
        qram = qram_prepare(COEFFICIENTS, angle_width=10)
        be = lcu_prepare_select(TERMS, prepare=qram)
        matrix = dense_hamiltonian(TERMS)
        memory = {"prep__angles": qram.memory["angles"]}
        for column in range(4):
            state = simulate(be.operation.program(), memory, initial={"target": column})
            for row in range(4):
                self.assertAlmostEqual(
                    state.amplitudes.get((row, 0), 0) * be.alpha,
                    matrix[row][column],
                    delta=0.02,
                )

    def test_qubitization_walk_composition(self):
        be = lcu_prepare_select(TERMS)
        walk = qubitization_walk(be)
        program = walk.program()
        self.assertFalse(unresolved(program))
        self.assertEqual(dict(walk.module.attributes)["algorithm"], "qubitization_walk")
        matrix = dense_hamiltonian(TERMS)
        # walk = (2|0><0| - I)·BE: the amplitudes on the signal=0 block are exactly H/alpha.
        for column in range(4):
            state = simulate(program, initial={"target": column})
            for row in range(4):
                self.assertAlmostEqual(
                    state.amplitudes.get((row, 0), 0) * be.alpha,
                    matrix[row][column],
                    places=10,
                )

    def test_pauli_hamiltonian_input_and_single_term_shortcut(self):
        hamiltonian = PauliHamiltonian(((1.0, "X"), (2.0, "Z")))
        be = lcu_prepare_select(hamiltonian)
        self.assertAlmostEqual(be.alpha, 3.0)
        matrix = dense_hamiltonian(hamiltonian.terms)
        assert_block_equals(self, be, matrix, places=10)
        single = lcu_prepare_select(((2.0, "Z"),))
        self.assertAlmostEqual(single.alpha, 2.0)
        state = simulate(single.operation.program(), initial={"target": 1})
        self.assertAlmostEqual(state.amplitudes[(1, 0)], -1)

    def test_bad_inputs_fail_at_generation(self):
        cases = [
            lambda: gate_prepare([1.0]),
            lambda: gate_prepare([0.0, 0.0]),
            lambda: gate_prepare([math.inf, 1.0]),
            lambda: select_pauli([(1.0, "X")]),
            lambda: select_pauli([(1.0, "X"), (1.0, "XZ")]),
            lambda: select_pauli([(1.0, "XA"), (1.0, "ZZ")]),
            lambda: lcu_prepare_select([]),
            lambda: lcu_prepare_select(TERMS, prepare=gate_prepare([1.0, 1.0])),
            lambda: alias_prepare(COEFFICIENTS, database=gate_database(2, 4, [0] * 4)),
        ]
        for case in cases:
            with self.assertRaises(ValidationError):
                case()


if __name__ == "__main__":
    unittest.main()
