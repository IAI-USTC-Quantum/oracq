"""Numerical witnesses for the graph-oracle input model and Szegedy/MNRS walk search."""

import math
import unittest

from oracq import Bits, Builder, ValidationError, bind, declare, simulate, unresolved
from oracq.algorithms.input_model.graph_walks import (
    abstract_adjacency,
    as_adjacency,
    gate_adjacency,
    hitting_times,
    qram_adjacency,
    quantum_walk_search,
    suggest_steps,
    szegedy_setup,
    szegedy_walk,
    transition_matrix,
)
from oracq.algorithms.input_model.oracles import (
    abstract_state_prep,
    invoke,
    phase_marks,
    resources_for,
)


def cycle_table(n):
    """Alternating edge-coloring neighbor table of an even cycle, guaranteeing the involution N(N(v,j),j)=v."""
    return [
        [(v + 1) % n if v % 2 == 0 else (v - 1) % n, (v - 1) % n if v % 2 == 0 else (v + 1) % n]
        for v in range(n)
    ]


HYPERCUBE_Q3 = [[v ^ 1, v ^ 2, v ^ 4, v] for v in range(8)]
COMPLETE_K4 = [[(v + 1 + j) % 4 for j in range(3)] + [v] for v in range(4)]


def vertex_probability(state, vertex_bits, marked):
    mask = (1 << vertex_bits) - 1
    return sum(
        abs(a) ** 2 for key, a in state.amplitudes.items() if (key[0] & mask) in set(marked)
    )


class AdjacencyOracleTests(unittest.TestCase):
    def test_gate_adjacency_implements_neighbor_table(self):
        oracle = gate_adjacency(HYPERCUBE_Q3)
        self.assertEqual(oracle.vertex_bits, 3)
        self.assertEqual(oracle.degree_bits, 2)
        for v in range(8):
            for j in range(4):
                state = simulate(oracle.operation.program(), initial={"vertex": v, "index": j})
                self.assertAlmostEqual(
                    state.amplitudes.get((v, j, HYPERCUBE_Q3[v][j]), 0), 1, places=10
                )

    def test_qram_adjacency_matches_gate_implementation(self):
        v_bits, g_bits = 3, 2
        memory = {v | (j << v_bits): HYPERCUBE_Q3[v][j] for v in range(8) for j in range(4)}
        oracle = qram_adjacency(v_bits, g_bits)
        for v in (0, 3, 7):
            for j in range(4):
                state = simulate(
                    oracle.operation.program(),
                    memory={"table": memory},
                    initial={"vertex": v, "index": j},
                )
                self.assertAlmostEqual(
                    state.amplitudes.get((v, j, HYPERCUBE_Q3[v][j]), 0), 1, places=10
                )

    def test_xor_database_adapter(self):
        database = gate_adjacency(HYPERCUBE_Q3).xor_database()
        self.assertEqual(database.address_width, 5)
        self.assertEqual(database.data_width, 3)
        state = simulate(database.operation.program(), initial={"address": 3 | (1 << 3)})
        self.assertAlmostEqual(state.amplitudes.get((3 | (1 << 3), HYPERCUBE_Q3[3][1]), 0), 1)

    def test_bad_tables_fail_at_generation(self):
        with self.assertRaises(ValidationError):
            gate_adjacency([[1, 0], [1]])  # not rectangular
        with self.assertRaises(ValidationError):
            gate_adjacency([[1], [3]])  # vertex out of range
        with self.assertRaises(ValidationError):
            hitting_times(transition_matrix([[1], [0], [2]]), {0})  # vertex 2 cannot reach the marked set

    def test_as_adjacency_accepts_bare_operation(self):
        adjacency = gate_adjacency(HYPERCUBE_Q3)
        recovered = as_adjacency(adjacency.operation)
        self.assertEqual(recovered, adjacency)
        szegedy_setup(recovered)


class SzegedyWalkTests(unittest.TestCase):
    def test_walk_step_is_unitary(self):
        walk = szegedy_walk(gate_adjacency(HYPERCUBE_Q3))
        b = Builder(
            "WalkUnitaryWitness",
            {"current": Bits(3), "peer": Bits(3), "index": Bits(2)},
            resources_for(("walk", walk)),
        )
        invoke(b, walk, "walk", current=b["current"], peer=b["peer"], index=b["index"])
        with b.adjoint():
            invoke(b, walk, "walk", current=b["current"], peer=b["peer"], index=b["index"])
        operation = b.finish()
        for initial in ((0, 0, 0), (3, 5, 1), (7, 7, 3)):
            state = simulate(
                operation.program(),
                initial={"current": initial[0], "peer": initial[1], "index": initial[2]},
            )
            self.assertAlmostEqual(state.amplitudes.get(initial, 0), 1, places=10)

    def test_stationary_state_is_fixed_point(self):
        adjacency = gate_adjacency(cycle_table(8))
        setup, walk = szegedy_setup(adjacency), szegedy_walk(adjacency)
        b = Builder(
            "WalkFixedPoint",
            {"target": Bits(7), "work": Bits(0)},
            resources_for(("setup", setup.operation), ("walk", walk)),
        )
        invoke(b, setup.operation, "setup", target=b["target"], work=b["work"])
        invoke(
            b,
            walk,
            "walk",
            current=b["target"][:3],
            peer=b["target"][3:6],
            index=b["target"][6:],
        )
        walked = simulate(b.finish().program())
        initial = simulate(setup.operation.program())
        self.assertEqual(set(walked.amplitudes), set(initial.amplitudes))
        for key, amplitude in initial.amplitudes.items():
            self.assertAlmostEqual(walked.amplitudes[key], amplitude, places=10)


class QuantumWalkSearchTests(unittest.TestCase):
    def test_search_amplifies_marked_vertex_on_hypercube(self):
        adjacency = gate_adjacency(HYPERCUBE_Q3)
        steps = suggest_steps(transition_matrix(HYPERCUBE_Q3), {0})
        self.assertEqual(steps, 3)
        search = quantum_walk_search(
            szegedy_setup(adjacency), szegedy_walk(adjacency), phase_marks(3, (0,)), steps
        )
        self.assertEqual(dict(search.module.attributes)["algorithm"], "quantum_walk_search")
        baseline = quantum_walk_search(
            szegedy_setup(adjacency), szegedy_walk(adjacency), phase_marks(3, (0,)), 0
        )
        p0 = vertex_probability(simulate(baseline.program()), 3, (0,))
        p1 = vertex_probability(simulate(search.program()), 3, (0,))
        self.assertAlmostEqual(p0, 0.125, places=10)
        self.assertAlmostEqual(p1, 0.78125, places=9)

    def test_search_recovers_grover_limit_on_complete_graph(self):
        adjacency = gate_adjacency(COMPLETE_K4)
        steps = suggest_steps(transition_matrix(COMPLETE_K4), {0})
        self.assertEqual(steps, 2)
        search = quantum_walk_search(
            szegedy_setup(adjacency), szegedy_walk(adjacency), phase_marks(2, (0,)), steps
        )
        probability = vertex_probability(simulate(search.program()), 2, (0,))
        self.assertAlmostEqual(probability, 1.0, places=10)

    def test_abstract_handles_bind_to_gate_implementations(self):
        adjacency = abstract_adjacency(3, 2, name="G")
        setup, walk = szegedy_setup(adjacency), szegedy_walk(adjacency)
        marked = declare("M", {"target": Bits(3)}, paradigm="phase_oracle")
        search = quantum_walk_search(setup, walk, marked, 3)
        program = search.program()
        self.assertEqual({r.name for r in unresolved(program)}, {"G", "M"})
        bound = bind(
            program,
            {"G": gate_adjacency(HYPERCUBE_Q3).operation, "M": phase_marks(3, (0,))},
        )
        self.assertFalse(unresolved(bound))
        probability = vertex_probability(simulate(bound), 3, (0,))
        self.assertAlmostEqual(probability, 0.78125, places=9)

    def test_walk_handle_may_stay_abstract(self):
        walk = declare(
            "W", {"current": Bits(3), "peer": Bits(3), "index": Bits(2)}, paradigm="unitary"
        )
        search = quantum_walk_search(
            abstract_state_prep("S", 8), walk, phase_marks(3, (0,)), 2
        )
        self.assertEqual({r.name for r in unresolved(search.program())}, {"S", "W"})
        bound = bind(
            search.program(),
            {
                "S": szegedy_setup(gate_adjacency(HYPERCUBE_Q3)).operation,
                "W": szegedy_walk(gate_adjacency(HYPERCUBE_Q3)),
            },
        )
        self.assertFalse(unresolved(bound))

    def test_search_rejects_mismatched_handles(self):
        adjacency = gate_adjacency(HYPERCUBE_Q3)
        with self.assertRaises(ValidationError):
            quantum_walk_search(
                szegedy_setup(adjacency), szegedy_walk(adjacency), phase_marks(3, (0,)), -1
            )
        with self.assertRaises(ValidationError):
            quantum_walk_search(
                abstract_state_prep("Narrow", 6),
                szegedy_walk(adjacency),
                phase_marks(3, (0,)),
                1,
            )


class HittingTimeTests(unittest.TestCase):
    def test_cycle_hitting_times_match_closed_form(self):
        hits = hitting_times(transition_matrix(cycle_table(8)), {0})
        for k in range(8):
            distance = min(k, 8 - k)
            self.assertAlmostEqual(hits[k], distance * (8 - distance), places=9)

    def test_suggest_steps_matches_grover_scaling(self):
        transition = transition_matrix(COMPLETE_K4)
        for v in (1, 2, 3):
            self.assertAlmostEqual(hitting_times(transition, {0})[v], 4.0, places=9)
        self.assertEqual(suggest_steps(transition, {0}), 2)
        self.assertTrue(math.isfinite(suggest_steps(transition_matrix(cycle_table(8)), {0})))

    def test_all_marked_needs_no_steps(self):
        self.assertEqual(suggest_steps(transition_matrix(COMPLETE_K4), {0, 1, 2, 3}), 0)


if __name__ == "__main__":
    unittest.main()
