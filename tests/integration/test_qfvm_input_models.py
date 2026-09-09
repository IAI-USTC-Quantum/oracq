"""QFVM 位置置换/填充元素及替换描述的真实后端见证。"""

import unittest

from pyqecclang import Bits, Builder, FixedFormat, arithmetic_native_registry, run_pysparq
from pyqecclang.applications.qfvm import (
    bind_qfvm,
    geometry_cells,
    qfvm_sparse_access,
    roe_qfvm_inputs,
)


class QfvmInputNativeTests(unittest.TestCase):
    def test_location_full_permutation_and_inverse(self):
        inputs = roe_qfvm_inputs(fmt=FixedFormat(4, 1))
        location = qfvm_sparse_access(inputs).location
        b = Builder(
            "location_witness",
            {"column": Bits(5), "index": Bits(5), "record": Bits(5), "work": Bits(0)},
        )
        b.h(b["column"][:2])
        b.h(b["index"])
        b.xor(b["index"], b["record"])
        b.call(location, column=b["column"], index=b["index"], work=b["work"])
        program = bind_qfvm(b.finish().program(), inputs)
        geometry = geometry_cells(inputs)
        result = run_pysparq(
            program, {"geometry": geometry}, native_registry=arithmetic_native_registry(program)
        )
        rows = {}
        for (column, index, record, work), amplitude in result.amplitudes.items():
            self.assertEqual(work, 0)
            self.assertGreater(abs(amplitude), 0)
            rows.setdefault(column, {})[record] = index
        for column, mapping in rows.items():
            self.assertEqual(set(mapping.values()), set(range(32)))
            for slot in range(9):
                expected = (
                    (column + slot) % 32 if column % 4 == 3 else geometry[column + (slot << 5)] & 31
                )
                self.assertEqual(mapping[slot], expected)
        inverse = Builder(
            "location_inverse", {"column": Bits(5), "index": Bits(5), "work": Bits(0)}
        )
        inverse.h(inverse["column"][:2])
        inverse.h(inverse["index"])
        inverse.call(
            location, column=inverse["column"], index=inverse["index"], work=inverse["work"]
        )
        with inverse.adjoint():
            inverse.call(
                location, column=inverse["column"], index=inverse["index"], work=inverse["work"]
            )
        p = bind_qfvm(inverse.finish().program(), inputs)
        restored = run_pysparq(
            p, {"geometry": geometry}, native_registry=arithmetic_native_registry(p)
        )
        self.assertEqual(len(restored.amplitudes), 128)
        self.assertTrue(all(abs(v - 128**-0.5) < 1e-11 for v in restored.amplitudes.values()))

    def test_padding_entry_xor_with_qram_raw_fields(self):
        inputs = roe_qfvm_inputs(fmt=FixedFormat(4, 1))
        entry = qfvm_sparse_access(inputs).entry
        b = Builder("padding_entry", {r.name: r.type for r in entry.module.registers})
        b.x(b["row"][:2])
        b.x(b["column"][:2])
        b.x(b["data"][:2])
        b.call(entry, row=b["row"], column=b["column"], data=b["data"])
        p = bind_qfvm(b.finish().program(), inputs)
        memory = {
            "geometry": geometry_cells(inputs),
            "rho": {i: 2 for i in range(4)},
            "momentum": {},
            "energy": {i: 4 for i in range(4)},
        }
        report = {}
        result = run_pysparq(
            p, memory, native_registry=arithmetic_native_registry(p), report=report
        )
        self.assertEqual(result.amplitudes, {(3, 3, 1): 1 + 0j})
        self.assertGreater(report["native_calls"], 100)

    def test_two_solver_descriptions_parse_in_actual_originir(self):
        from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser

        from pyqecclang import export_toffoli_u3_cz
        from pyqecclang.algorithms.oracles import (
            SparseAccess,
            basis_state,
            gate_database,
            sparse_entry,
            sparse_location_gate,
        )
        from pyqecclang.algorithms.qlss import (
            CostaConfig,
            LinearSystem,
            SparseSystem,
            SpectralPromise,
            make_cks_qlss,
            make_costa_qlss,
        )

        fmt = FixedFormat(4, 2)
        access = SparseAccess(
            sparse_location_gate(1, [[0, 1], [0, 1]], work_width=0),
            sparse_entry(gate_database(2, 4, {0: 2, 1: 15, 2: 15, 3: 2}), 1),
            1,
            4,
            2,
        )
        problem = LinearSystem(
            sparse=SparseSystem(
                access,
                fmt,
                1,
                basis_state(1),
                SpectralPromise(0.75, 0.25),
                diagonal_nonnegative=True,
                hermitian=True,
            ),
            rhs_norm=1,
        )
        for solver in (make_cks_qlss(), make_costa_qlss(CostaConfig(steps=1))):
            result = solver(problem)
            parser = OriginIR_BaseParser()
            parser.parse(export_toffoli_u3_cz(result.operation.program()).text)
            self.assertGreater(parser.n_qubit, 1)
