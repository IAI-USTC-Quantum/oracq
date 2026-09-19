"""QRAM 指针式读写（qmem）的语义测试：寻址、多维视图、随机写与台账。"""

import json
import unittest
from pathlib import Path

from jsonschema import Draft202012Validator

from pyqecclang import (
    QRAM,
    Builder,
    QMem,
    UInt,
    ValidationError,
    dumps,
    estimate_resources,
    export_originir,
    export_strict,
    loads,
    quantikz,
    simulate,
)

ROOT = Path(__file__).resolve().parents[2]
DATA = list(range(16))


class PointerLoadTests(unittest.TestCase):
    def test_constant_offset_load(self):
        b = Builder("const_off", {"out": UInt(4)}, {"rom": QRAM(4, 4)})
        mem = QMem(b, "rom")
        (mem.ptr() + 5).load(b["out"])
        state = simulate(b.finish().program(), {"rom": DATA})
        self.assertEqual(state.amplitudes, {(5,): 1 + 0j})

    def test_negative_offset_wraps_mod_address_space(self):
        b = Builder("wrap", {"out": UInt(4)}, {"rom": QRAM(4, 4)})
        mem = QMem(b, "rom")
        (mem.ptr(2) - 5).load(b["out"])
        state = simulate(b.finish().program(), {"rom": DATA})
        self.assertEqual(state.amplitudes, {(13,): 1 + 0j})

    def test_quantum_pointer_register(self):
        b = Builder("qptr", {"idx": UInt(2), "out": UInt(4)}, {"rom": QRAM(2, 4)})
        mem = QMem(b, "rom")
        (mem.ptr(b["idx"]) + 1).load(b["out"])
        state = simulate(b.finish().program(), {"rom": [3, 5, 7, 9]}, initial={"idx": 2})
        self.assertEqual(state.amplitudes, {(2, 9): 1 + 0j})

    def test_pointer_chain_uses_ripple_adder(self):
        b = Builder("chain", {"i": UInt(2), "j": UInt(2), "out": UInt(4)}, {"rom": QRAM(2, 4)})
        mem = QMem(b, "rom")
        (mem.ptr() + b["i"] + b["j"]).load(b["out"])
        operation = b.finish()
        state = simulate(operation.program(), {"rom": [3, 5, 7, 9]}, initial={"i": 1, "j": 2})
        self.assertEqual(state.amplitudes, {(1, 2, 9): 1 + 0j})
        estimate = operation.estimate()
        self.assertGreater(estimate.toffoli, 0)
        self.assertEqual(dict(estimate.qram_queries), {"rom": 1})

    def test_ripple_adder_correct_in_coherent_superposition(self):
        table = [3, 5, 7, 9]
        b = Builder("sweep", {"i": UInt(2), "j": UInt(2), "out": UInt(4)}, {"rom": QRAM(2, 4)})
        mem = QMem(b, "rom")
        b.h(b["i"])
        b.h(b["j"])
        (mem.ptr() + b["i"] + b["j"]).load(b["out"])
        state = simulate(b.finish().program(), {"rom": table})
        self.assertEqual(len(state.amplitudes), 16)
        for (i, j, out), amplitude in state.amplitudes.items():
            self.assertEqual(out, table[(i + j) % 4])
            self.assertAlmostEqual(abs(amplitude), 0.25)

    def test_pointer_chain_all_basis_combinations(self):
        table = [3, 5, 7, 9]
        for i in range(4):
            for j in range(4):
                b = Builder("det", {"i": UInt(2), "j": UInt(2), "o": UInt(4)}, {"rom": QRAM(2, 4)})
                mem = QMem(b, "rom")
                (mem.ptr() + b["i"] + b["j"] + 3).load(b["o"])
                state = simulate(
                    b.finish().program(), {"rom": table}, initial={"i": i, "j": j}
                )
                self.assertEqual(list(state.amplitudes), [(i, j, table[(i + j + 3) % 4])])

    def test_superposed_address_loads_every_cell(self):
        b = Builder("sup", {"idx": UInt(2), "out": UInt(4)}, {"rom": QRAM(2, 4)})
        mem = QMem(b, "rom")
        b.h(b["idx"])
        mem[b["idx"]].load(b["out"])
        state = simulate(b.finish().program(), {"rom": [1, 2, 4, 8]})
        expected = {(0, 1): 0.25, (1, 2): 0.25, (2, 4): 0.25, (3, 8): 0.25}
        for key, probability in expected.items():
            self.assertAlmostEqual(abs(state.amplitudes[key]) ** 2, probability)

    def test_direct_address_emits_single_load(self):
        b = Builder("direct", {"idx": UInt(2), "out": UInt(4)}, {"rom": QRAM(2, 4)})
        QMem(b, "rom")[b["idx"]].load(b["out"])
        body = b.finish().program().main.body
        self.assertEqual(len(body), 1)
        self.assertEqual(body[0].resource, "rom")
        self.assertEqual(body[0].address.parts[0].register, "idx")


class MultiDimensionalTests(unittest.TestCase):
    def test_two_dimensional_quantum_row(self):
        b = Builder("grid", {"row": UInt(2), "out": UInt(4)}, {"rom": QRAM(4, 4)})
        grid = QMem(b, "rom", shape=(4, 4))
        grid[b["row"], 3].load(b["out"])
        state = simulate(b.finish().program(), {"rom": DATA}, initial={"row": 2})
        self.assertEqual(state.amplitudes, {(2, 11): 1 + 0j})

    def test_two_dimensional_quantum_row_and_col(self):
        b = Builder("grid2", {"row": UInt(2), "col": UInt(2), "out": UInt(4)}, {"rom": QRAM(4, 4)})
        grid = QMem(b, "rom", shape=(4, 4))
        b.h(b["row"])
        b.h(b["col"])
        grid[b["row"], b["col"]].load(b["out"])
        state = simulate(b.finish().program(), {"rom": DATA})
        self.assertEqual(len(state.amplitudes), 16)
        for (row, col, out), amplitude in state.amplitudes.items():
            self.assertEqual(out, DATA[row * 4 + col])
            self.assertAlmostEqual(abs(amplitude), 0.25)

    def test_disjoint_shifts_avoid_ripple_adder(self):
        b = Builder("noadder", {"row": UInt(2), "col": UInt(2), "out": UInt(4)}, {"rom": QRAM(4, 4)})
        grid = QMem(b, "rom", shape=(4, 4))
        grid[b["row"], b["col"]].load(b["out"])
        estimate = b.finish().estimate()
        self.assertEqual(estimate.toffoli, 0)
        self.assertEqual(estimate.qram_total, 1)

    def test_slice_view_rebases_addresses(self):
        b = Builder("slice", {"col": UInt(2), "out": UInt(4)}, {"rom": QRAM(4, 4)})
        grid = QMem(b, "rom", shape=(4, 4))
        view = grid[2:]
        self.assertEqual(view.shape, (2, 4))
        view[1, b["col"]].load(b["out"])
        state = simulate(b.finish().program(), {"rom": DATA}, initial={"col": 1})
        self.assertEqual(state.amplitudes, {(1, 13): 1 + 0j})

    def test_shape_validation(self):
        b = Builder("shape", {"out": UInt(4)}, {"rom": QRAM(4, 4)})
        with self.assertRaisesRegex(ValidationError, "地址空间"):
            QMem(b, "rom", shape=(5, 5))
        mem = QMem(b, "rom", shape=(4, 4))
        with self.assertRaisesRegex(ValidationError, "越界"):
            mem[4]
        with self.assertRaisesRegex(ValidationError, "维量子下标宽度"):
            mem[b["out"]]
        with self.assertRaisesRegex(ValidationError, "没有声明"):
            QMem(b, "missing")


class StoreTests(unittest.TestCase):
    def test_store_then_load_roundtrip(self):
        b = Builder("wr", {"addr": UInt(2), "val": UInt(4), "out": UInt(4)}, {"ram": QRAM(2, 4)})
        mem = QMem(b, "ram")
        mem[b["addr"]].store(b["val"])
        mem[b["addr"]].load(b["out"])
        state = simulate(b.finish().program(), {"ram": [0, 0, 0, 0]}, initial={"addr": 2, "val": 13})
        self.assertEqual(state.amplitudes, {(2, 13, 13): 1 + 0j})

    def test_store_overwrites_and_other_cells_untouched(self):
        b = Builder("ow", {"addr": UInt(2), "val": UInt(4), "out0": UInt(4), "out1": UInt(4)}, {"ram": QRAM(2, 4)})
        mem = QMem(b, "ram")
        mem[b["addr"]].store(b["val"])
        mem.ptr(0).load(b["out0"])
        mem.ptr(1).load(b["out1"])
        program = b.finish().program()
        state = simulate(program, {"ram": [7, 7, 7, 7]}, initial={"addr": 1, "val": 0})
        self.assertEqual(state.amplitudes, {(1, 0, 7, 0): 1 + 0j})

    def test_superposed_store_is_rejected_at_runtime(self):
        b = Builder("supstore", {"a": UInt(2), "v": UInt(4)}, {"ram": QRAM(2, 4)})
        mem = QMem(b, "ram")
        b.h(b["a"])
        mem[b["a"]].store(b["v"])
        with self.assertRaisesRegex(ValidationError, "确定基矢"):
            simulate(b.finish().program(), {"ram": [0, 0, 0, 0]})

    def test_store_rejected_inside_control_and_adjoint(self):
        for block in ("control", "adjoint"):
            b = Builder(f"bad_{block}", {"a": UInt(2), "v": UInt(4), "c": UInt(1)}, {"ram": QRAM(2, 4)})
            mem = QMem(b, "ram")
            context = b.control(b["c"]) if block == "control" else b.adjoint()
            with context:
                mem[b["a"]].store(b["v"])
            with self.assertRaisesRegex(ValidationError, "非酉副作用"):
                b.finish()

    def test_store_capability_blocks_controlled_call(self):
        inner = Builder("inner", {"a": UInt(2), "v": UInt(4)}, {"ram": QRAM(2, 4)})
        QMem(inner, "ram")[inner["a"]].store(inner["v"])
        target = inner.finish()
        b = Builder("outer", {"a": UInt(2), "v": UInt(4), "c": UInt(1)}, {"ram": QRAM(2, 4)})
        with b.control(b["c"]):
            b.call(target, a=b["a"], v=b["v"], resources={"ram": "ram"})
        with self.assertRaisesRegex(ValidationError, "supports_controlled"):
            b.finish()

    def test_store_allows_repeat(self):
        b = Builder("rep", {"a": UInt(2), "v": UInt(4)}, {"ram": QRAM(2, 4)})
        mem = QMem(b, "ram")
        with b.repeat(2):
            mem[b["a"]].store(b["v"])
        state = simulate(b.finish().program(), {"ram": [0, 0, 0, 0]}, initial={"a": 3, "v": 5})
        self.assertEqual(state.amplitudes, {(3, 5): 1 + 0j})


class LedgerAndBackendTests(unittest.TestCase):
    def test_store_is_counted_without_gate_cost(self):
        b = Builder("led", {"a": UInt(2), "v": UInt(4), "o": UInt(4)}, {"ram": QRAM(2, 4)})
        mem = QMem(b, "ram")
        mem[b["a"]].store(b["v"])
        mem[b["a"]].load(b["o"])
        estimate = b.finish().estimate()
        self.assertEqual(dict(estimate.qram_writes), {"ram": 1})
        self.assertEqual(dict(estimate.qram_queries), {"ram": 1})
        self.assertEqual(estimate.gate_total, 0)

    def test_serialization_and_schema_roundtrip(self):
        b = Builder("ser", {"a": UInt(2), "v": UInt(4)}, {"ram": QRAM(2, 4)})
        QMem(b, "ram")[b["a"]].store(b["v"])
        program = b.finish().program()
        text = dumps(program)
        self.assertIn('"tag": "Store"', text)
        self.assertEqual(loads(text), program)
        schema = json.loads((ROOT / "docs/reference/schemas/rir.schema.json").read_text())
        Draft202012Validator.check_schema(schema)
        Draft202012Validator(schema).validate(json.loads(text))

    def test_originir_emits_qramwrite_and_strict_counts_it(self):
        b = Builder("exp", {"a": UInt(2), "v": UInt(4), "o": UInt(4)}, {"ram": QRAM(2, 4)})
        mem = QMem(b, "ram")
        mem[b["a"]].store(b["v"])
        mem[b["a"]].load(b["o"])
        program = b.finish().program()
        artifact = export_originir(program)
        write_lines = [line for line in artifact.text.splitlines() if line.startswith("QRAMWRITE")]
        self.assertEqual(len(write_lines), 1)
        self.assertTrue(write_lines[0].startswith("QRAMWRITE ram_ram "))
        strict = export_strict(program)
        self.assertEqual(dict(strict.qram_writes), {"ram_ram": 1})
        self.assertEqual(dict(strict.qram_queries), {"ram_ram": 1})
        self.assertIn("QRAM", quantikz(program))

    def test_estimate_resources_module_function(self):
        b = Builder("fn", {"a": UInt(2), "v": UInt(4)}, {"ram": QRAM(2, 4)})
        QMem(b, "ram")[b["a"]].store(b["v"])
        report = estimate_resources(b.finish().program())
        self.assertEqual(report.to_dict()["qram_writes"], {"ram": 1})
        self.assertEqual(report.to_dict()["qram_write_total"], 1)


if __name__ == "__main__":
    unittest.main()
