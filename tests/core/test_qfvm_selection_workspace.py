"""Probe the emitted selection prefix against raw QRAM geometry, before Roe arithmetic."""

import unittest
from dataclasses import replace

from oracq import Bits, Builder, FixedFormat, simulate
from oracq.applications.qfvm import bind_qfvm, geometry_cells, qfvm_sparse_access, roe_qfvm_inputs
from oracq.infrastructure.builder import Operation
from oracq.infrastructure.ir import Adjoint, Call, Primitive


class GeometryWorkspaceTests(unittest.TestCase):
    def test_selection_and_uncompute_match_classical_geometry(self):
        inputs = roe_qfvm_inputs(fmt=FixedFormat(4, 1))
        entry = qfvm_sparse_access(inputs).entry
        body = entry.module.body
        end = next(i for i, node in enumerate(body)
                   if isinstance(node, Call) and node.module.startswith("roe_entry_"))
        prefix = body[:end]
        width = inputs.geometry_width
        probe = Builder("selection_workspace_probe", {
            **{r.name: r.type for r in entry.module.registers}, "selected_out": Bits(width),
        })
        for register in entry.module.locals:
            probe.local(register.name, register.type)
        output = probe["selected_out"]
        selected = probe["selected_geometry"]
        module = replace(entry.module, name=probe.name, registers=probe.registers,
                         body=(*prefix, Primitive("xor", (selected, output)), Adjoint(prefix)))
        program = bind_qfvm(Operation(module, entry.dependencies).program(), inputs)
        geometry = geometry_cells(inputs)
        memory = {resource.name: {} for resource in program.main.resources}
        memory["geometry"] = geometry
        n = inputs.width
        for column in (0, 1, 3, 4, 7, 15, 31):
            rows = {0, 3, 31}
            rows.update(geometry[column + (rank << n)] & ((1 << n) - 1) for rank in range(9))
            for row in sorted(rows):
                expected = 7
                for rank in range(9):
                    word = geometry[column + (rank << n)]
                    if word >> (width - 1) and word & ((1 << n) - 1) == row:
                        expected ^= word
                result = simulate(program, memory, initial={"row": row, "column": column, "selected_out": 7})
                self.assertEqual(result.amplitudes, {(row, column, 0, expected): 1 + 0j})


if __name__ == "__main__":
    unittest.main()
