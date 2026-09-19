"""QFVM-on-QMem 与既有 QFVM 路径的等价性对拍（simulate 逐振幅）。

数据访问层（几何表 fuse 编址 vs 二维指针、状态表 vs 三个独立库、周期邻居指针）
用叠加探针全覆盖；位置 oracle 与残差态做完整对拍；含 Roe 算术的物理层
（roe_face 编译电路，展开约 1.6e7 步）在单点上做深度对拍。
"""

import unittest

from pyqecclang import QRAM, Builder, simulate
from pyqecclang.algorithms.arithmetic import FixedFormat
from pyqecclang.algorithms.oracles import invoke
from pyqecclang.applications.flow_data import RoeFlowData
from pyqecclang.applications.qfvm import (
    bind_qfvm,
    geometry_cells,
    qfvm_sparse_access,
    rhs_qram_preparation,
    roe_entry,
    roe_qfvm_inputs,
)
from pyqecclang.applications.qfvm_qmem import (
    RoeQmemData,
    qfvm_qmem_location,
    qfvm_qmem_physical,
    qfvm_qmem_rhs,
)
from pyqecclang.infrastructure.ir import Bits, fuse
from pyqecclang.infrastructure.qmem import QMem

FMT = FixedFormat(5, 2)
DELTA = 0.5
STATES = [(1.0, 0.0, 1.0), (1.0, 0.25, 1.0), (1.25, 0.0, 1.0), (1.25, -0.25, 1.0)]
AW = 6


def rounded(state):
    return {k: round(a.real, 12) for k, a in state.amplitudes.items()}


class DataAccessProbeTests(unittest.TestCase):
    """叠加探针：QMem 二维寻址与 legacy fuse/invoke 编址在真实数据上逐字一致。"""

    def setUp(self):
        self.inputs = roe_qfvm_inputs(fmt=FMT, angle_width=AW)
        self.geometry = geometry_cells(self.inputs)

    def test_geometry_probe_matches_fuse_addressing(self):
        width = self.inputs.geometry_width

        def legacy_probe():
            b = Builder("legacy_geom_probe", {"column": Bits(self.inputs.width), "acc": Bits(width)})
            slot = b.local("slot", Bits(4))
            word = b.local("word", Bits(width))
            for rank in range(9):
                for bit in range(4):
                    if (rank >> bit) & 1:
                        b.x(slot[bit])
                invoke(b, self.inputs.geometry.operation, "geometry",
                       address=fuse(b["column"], slot), data=word)
                b.xor(word, b["acc"])
                invoke(b, self.inputs.geometry.operation, "geometry",
                       address=fuse(b["column"], slot), data=word)
                for bit in range(4):
                    if (rank >> bit) & 1:
                        b.x(slot[bit])
            b.h(b["column"])
            return bind_qfvm(b.finish().program(), self.inputs)

        b = Builder(
            "qmem_geom_probe",
            {"column": Bits(self.inputs.width), "acc": Bits(width)},
            {"geometry": QRAM(self.inputs.width + 4, width)},
        )
        geometry = QMem(b, "geometry", shape=(16, 1 << self.inputs.width))
        slot = b.local("slot", Bits(4))
        word = b.local("word", Bits(width))
        for rank in range(9):
            for bit in range(4):
                if (rank >> bit) & 1:
                    b.x(slot[bit])
            cell = geometry[slot, b["column"]]
            cell.load(word)
            b.xor(word, b["acc"])
            cell.load(word)
            for bit in range(4):
                if (rank >> bit) & 1:
                    b.x(slot[bit])
        b.h(b["column"])
        legacy = simulate(legacy_probe(), {"geometry": self.geometry})
        modern = simulate(b.finish().program(), {"geometry": self.geometry})
        self.assertEqual(rounded(legacy), rounded(modern))

    def test_state_probe_matches_field_databases(self):
        data = RoeQmemData(STATES, fmt=FMT, angle_width=AW, entropy_delta=DELTA)
        cw = data.cell_width

        def legacy_probe():
            b = Builder("legacy_state_probe", {"cell": Bits(cw), "word0": Bits(FMT.width),
                                               "word1": Bits(FMT.width), "word2": Bits(FMT.width)})
            b.h(b["cell"])
            for name, out in zip(
                ("rho", "momentum", "energy"), ("word0", "word1", "word2"), strict=True
            ):
                invoke(b, getattr(self.inputs, name).operation, name,
                       address=b["cell"], data=b[out])
            return bind_qfvm(b.finish().program(), self.inputs)

        b = Builder(
            "qmem_state_probe",
            {"cell": Bits(cw), "word0": Bits(FMT.width),
             "word1": Bits(FMT.width), "word2": Bits(FMT.width)},
            {"state": QRAM(cw + 2, FMT.width)},
        )
        state = QMem(b, "state", shape=(3, 1 << cw))
        b.h(b["cell"])
        for field, out in enumerate(("word0", "word1", "word2")):
            state[field, b["cell"]].load(b[out])
        legacy = simulate(
            legacy_probe(),
            {key: data.flow.store.snapshot()[key] for key in ("rho", "momentum", "energy")},
        )
        modern = simulate(b.finish().program(), {"state": data.state_bank})
        self.assertEqual(rounded(legacy), rounded(modern))

    def test_periodic_neighbor_pointer_matches_legacy_scratch(self):
        data = RoeQmemData(STATES, fmt=FMT, angle_width=AW, entropy_delta=DELTA)
        cw = data.cell_width

        def legacy_probe():
            b = Builder("legacy_shift_probe", {"cell": Bits(cw), "word": Bits(FMT.width)})
            b.h(b["cell"])
            addr = b.local("addr", Bits(cw))
            b.xor(b["cell"], addr)
            b.add_const(addr.reinterpret("uint"), 1)
            invoke(b, self.inputs.rho.operation, "rho", address=addr, data=b["word"])
            b.add_const(addr.reinterpret("uint"), (-1) % (1 << cw))
            b.xor(b["cell"], addr)
            return bind_qfvm(b.finish().program(), self.inputs)

        b = Builder(
            "qmem_shift_probe",
            {"cell": Bits(cw), "word": Bits(FMT.width)},
            {"state": QRAM(cw + 2, FMT.width)},
        )
        state = QMem(b, "state", shape=(3, 1 << cw))
        b.h(b["cell"])
        addr = b.local("addr", Bits(cw))
        b.xor(b["cell"], addr)
        b.add_const(addr.reinterpret("uint"), 1)
        state[0, addr].load(b["word"])
        b.add_const(addr.reinterpret("uint"), (-1) % (1 << cw))
        b.xor(b["cell"], addr)
        legacy = simulate(legacy_probe(), {"rho": data.flow.store.snapshot()["rho"]})
        modern = simulate(b.finish().program(), {"state": data.state_bank})
        self.assertEqual(rounded(legacy), rounded(modern))

    def test_state_bank_layout(self):
        data = RoeQmemData(STATES, fmt=FMT, angle_width=AW, entropy_delta=DELTA)
        banks = data.flow.store.snapshot()
        cw = data.cell_width
        for field, key in enumerate(("rho", "momentum", "energy")):
            for address, word in banks[key].items():
                self.assertEqual(data.state_bank[(field << cw) | address], word)


class LocationEquivalenceTests(unittest.TestCase):
    def test_location_permutation_matches_legacy(self):
        inputs = roe_qfvm_inputs(fmt=FMT, angle_width=AW)
        geometry = geometry_cells(inputs)
        legacy_location = qfvm_sparse_access(inputs, entropy_delta=DELTA).location
        legacy = Builder(
            "legacy_location_driver",
            {r.name: r.type for r in legacy_location.module.registers},
        )
        legacy.h(legacy["column"])
        legacy.call(
            legacy_location, column=legacy["column"], index=legacy["index"], work=legacy["work"]
        )
        legacy_state = simulate(bind_qfvm(legacy.finish().program(), inputs), {"geometry": geometry})

        modern_location = qfvm_qmem_location(inputs)
        modern = Builder(
            "qmem_location_driver",
            {r.name: r.type for r in modern_location.module.registers},
            {r.name: r.type for r in modern_location.module.resources},
        )
        modern.h(modern["column"])
        modern.call(
            modern_location,
            column=modern["column"],
            index=modern["index"],
            work=modern["work"],
            resources={r.name: r.name for r in modern_location.module.resources},
        )
        modern_state = simulate(modern.finish().program(), {"geometry": geometry})
        self.assertEqual(rounded(legacy_state), rounded(modern_state))


class ResidualEquivalenceTests(unittest.TestCase):
    def test_rhs_amplitudes_match_legacy(self):
        inputs = roe_qfvm_inputs(fmt=FMT, angle_width=AW)
        flow = RoeFlowData(STATES, fmt=FMT, angle_width=AW, entropy_delta=DELTA)
        prep = rhs_qram_preparation(inputs)
        legacy = Builder(
            "legacy_rhs_driver",
            {r.name: r.type for r in prep.operation.module.registers},
            {
                "rhs_angles": QRAM(inputs.cell_width + 2, AW),
                "rhs_sign": QRAM(inputs.cell_width + 2, 1),
            },
        )
        legacy.call(
            prep.operation,
            target=legacy["target"],
            work=legacy["work"],
            resources={"prep__angles": "rhs_angles", "sign__table": "rhs_sign"},
        )
        banks = flow.store.snapshot()
        legacy_state = simulate(
            legacy.finish().program(),
            {"rhs_angles": banks["rhs_angles"], "rhs_sign": banks["rhs_sign"]},
        )

        data = RoeQmemData(STATES, fmt=FMT, angle_width=AW, entropy_delta=DELTA)
        modern_prep = qfvm_qmem_rhs(data)
        modern_state = simulate(modern_prep.operation.program(), data.vector.snapshot(signed=True))

        def marginal(state, width):
            return {key[0] & ((1 << width) - 1): amplitude for key, amplitude in state.amplitudes.items()}

        left = marginal(legacy_state, inputs.cell_width + 2)
        right = marginal(modern_state, data.cell_width + 2)
        self.assertEqual(set(left), set(right))
        for key in left:
            self.assertAlmostEqual(left[key], right[key], delta=1e-9)


class PhysicalEquivalenceTests(unittest.TestCase):
    def test_single_point_matches_legacy(self):
        """含编译 Roe 算术的物理层深对拍：source=1、row=col=0、band=1（中心带）。"""
        inputs = roe_qfvm_inputs(fmt=FMT, angle_width=AW)
        data = RoeQmemData(STATES, fmt=FMT, angle_width=AW, entropy_delta=DELTA)

        legacy_physical = roe_entry(inputs, entropy_delta=DELTA)
        legacy = Builder(
            "legacy_physical_driver",
            {r.name: r.type for r in legacy_physical.module.registers},
        )
        legacy.x(legacy["source"][0])
        legacy.call(
            legacy_physical,
            source=legacy["source"],
            row=legacy["row"],
            col=legacy["col"],
            band=legacy["band"],
            value=legacy["value"],
            status=legacy["status"],
        )
        banks = data.flow.store.snapshot()
        legacy_state = simulate(
            bind_qfvm(legacy.finish().program(), inputs),
            {key: banks[key] for key in ("rho", "momentum", "energy")},
            max_steps=60_000_000,
        )

        modern_physical = qfvm_qmem_physical(inputs, entropy_delta=DELTA)
        modern = Builder(
            "qmem_physical_driver",
            {r.name: r.type for r in modern_physical.module.registers},
            {r.name: r.type for r in modern_physical.module.resources},
        )
        modern.x(modern["source"][0])
        modern.call(
            modern_physical,
            source=modern["source"],
            row=modern["row"],
            col=modern["col"],
            band=modern["band"],
            value=modern["value"],
            status=modern["status"],
            resources={"state": "state"},
        )
        modern_state = simulate(
            modern.finish().program(), {"state": data.state_bank}, max_steps=60_000_000
        )
        self.assertEqual(rounded(legacy_state), rounded(modern_state))


if __name__ == "__main__":
    unittest.main()
