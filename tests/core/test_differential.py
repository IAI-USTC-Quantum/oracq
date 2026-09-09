"""方法之间的结构与开放输入差异；不是算法数值正确性测试。"""

import unittest
from functools import partial

from pyqecclang import dumps, loads
from pyqecclang.algorithms.arithmetic import FixedFormat
from pyqecclang.algorithms.carleman import PolynomialODE, carleman_qode
from pyqecclang.algorithms.cbmd import ContourPlan
from pyqecclang.algorithms.hamiltonian import taylor_hamiltonian
from pyqecclang.algorithms.ode import linear_qode
from pyqecclang.algorithms.oracles import abstract_block_encoding, abstract_state_prep
from pyqecclang.applications.flow_data import RoeFlowData
from pyqecclang.applications.qfvm import (
    bind_qfvm,
    geometry_cells,
    roe_qfvm_block_encoding,
    roe_qfvm_inputs,
)
from pyqecclang.infrastructure.linking import unresolved


class DifferentialStructureTests(unittest.TestCase):
    def test_four_methods_keep_input_oracles(self):
        a = abstract_block_encoding("input_A", 1, 1, 1.0)
        initial = abstract_state_prep("input_b", 1)
        function = partial(taylor_hamiltonian, degree=1)
        for method in ("lchs", "cbmd", "schrodingerization"):
            with self.subTest(method=method):
                state = linear_qode(method, hamiltonian_function=function)(a, initial, 0.1)
                program = state.operation.program()
                self.assertEqual({r.name for r in unresolved(program)}, {"input_A", "input_b"})
                self.assertIn(method, dict(program.main.attributes)["algorithm"])
                self.assertEqual(loads(dumps(program)), program)
        f2 = abstract_block_encoding("input_F2", 2, 1, 1.0)
        problem = PolynomialODE(1, ((1, a), (2, f2)), initial)
        state = carleman_qode(
            problem,
            0.1,
            linear_qode("cbmd", hamiltonian_function=function, plan=ContourPlan(cutoff=0)),
        )
        self.assertEqual(
            {r.name for r in unresolved(state.operation.program())},
            {"input_A", "input_b", "input_F2"},
        )
        placements = [
            dict(m.attributes)
            for m in state.operation.program().modules
            if "carleman_column_level" in dict(m.attributes)
        ]
        self.assertTrue(
            any(
                m["carleman_column_level"] == 2 and m["carleman_row_level"] == 1 for m in placements
            )
        )

    def test_cbmd_plan_records_omitted_terms(self):
        import json

        plan = ContourPlan(cutoff=2)
        metadata = json.loads(plan.metadata())
        self.assertEqual(len(plan.weights), 5)
        self.assertEqual(len(plan.auxiliary_coefficients), len(plan.poles))
        self.assertIn("auxiliary_nonhermitian_evolutions", metadata["omitted"])

    def test_qfvm_uses_raw_flow_and_modular_arithmetic(self):
        inputs = roe_qfvm_inputs(fmt=FixedFormat(4, 1))
        be = roe_qfvm_block_encoding(inputs, amax=4)
        program = be.operation.program()
        self.assertEqual(be.alpha, 36)
        self.assertEqual(len(unresolved(program)), 4)
        closed = bind_qfvm(program, inputs)
        self.assertFalse(unresolved(closed))
        self.assertEqual(
            {r.name for r in closed.main.resources},
            {"rho", "momentum", "energy", "geometry"},
        )
        kinds = {dict(m.attributes).get("arithmetic_kind") for m in closed.modules}
        self.assertTrue({"sqrt", "div", "mul", "select"} <= kinds)
        self.assertEqual(len(geometry_cells(inputs)), 1 << (inputs.width + 4))
        isometry = next(
            m
            for m in program.modules
            if dict(m.attributes).get("algorithm") == "cks_real_symmetric_isometry"
        )
        self.assertEqual(dict(isometry.attributes)["sparsity"], 9)
        self.assertTrue(dict(be.operation.module.attributes)["self_adjoint_extension"])

    def test_local_riemann_updates_only_neighbor_faces_and_tree(self):
        flow = RoeFlowData([(1, 0.125, 2)] * 16, fmt=FixedFormat(10, 5))
        before = flow.store.snapshot()
        patch = flow.update({5: (1.25, 0.25, 2.5)})
        self.assertEqual(patch.recomputed_faces, (4, 5))
        self.assertEqual(patch.recomputed_cells, (4, 5, 6))
        for bank in ("rho", "momentum", "energy"):
            self.assertEqual(set(patch.changes[bank]), {5})
            for cell in range(16):
                if cell != 5:
                    self.assertEqual(flow.store.banks[bank].get(cell), before[bank].get(cell))
        self.assertFalse(any("matrix" in bank for bank in flow.store.banks))
        self.assertLess(len(patch.changes["rhs_angles"]), 4 * 16 - 1)
