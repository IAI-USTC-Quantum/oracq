"""qdata 量子数据结构读路径的真实后端对拍；用具有 uniqc 和 pysparq 的解释器运行。"""

import unittest

from pyqecclang import (
    Builder,
    run_originir,
    run_pysparq,
    run_pysparq_rir,
    simulate,
)
from pyqecclang.algorithms.common.arithmetic import FixedFormat
from pyqecclang.algorithms.input_model.qdata import QMatrix, QVector

FMT = FixedFormat(8, 4)
AW = 10


class QDataBackendTests(unittest.TestCase):
    def compare_sparse(self, program, memory):
        """高工作区案例保留参考执行器与两个真实 PySparQ 路径的对拍。"""
        expected = simulate(program, memory).amplitudes
        for native in (
            run_pysparq(program, memory).amplitudes,
            run_pysparq_rir(program, memory).amplitudes,
        ):
            self.assertEqual(set(expected), set(native))
            for key in expected:
                self.assertAlmostEqual(expected[key], native[key], places=9)

    def compare(self, program, memory):
        """预算内案例额外覆盖真实 OriginIR 稠密后端。"""
        self.compare_sparse(program, memory)
        vector = simulate(program, memory).statevector()
        origin = list(run_originir(program, memory))
        for index, (a, b) in enumerate(zip(vector, origin, strict=True)):
            with self.subTest(index=index):
                self.assertAlmostEqual(a, complex(b), places=9)

    def test_vector_preparation_unsigned(self):
        vector = QVector((1.0, 2.0, 0.0, 3.0), fmt=FMT, angle_width=AW)
        self.compare(vector.preparation().operation.program(), vector.snapshot())

    def test_vector_preparation_signed(self):
        vector = QVector((1.0, -2.0, 0.5, 0.0), fmt=FMT, angle_width=AW)
        self.compare(
            vector.preparation(signed=True).operation.program(), vector.snapshot(signed=True)
        )

    def test_matrix_query_superposed(self):
        matrix = QMatrix(
            [
                [0.5, 0.5, 0.25, 0.0],
                [0.5, 0.0, 0.25, 0.5],
                [0.25, 0.5, 0.0, 0.5],
                [0.0, 0.0, 0.5, 0.5],
            ],
            fmt=FMT,
            angle_width=AW,
        )
        b = Builder(
            "query_driver",
            {r.name: r.type for r in matrix.query().operation.module.registers},
            {r.name: r.type for r in matrix.query().operation.module.resources},
        )
        b.h(b["address"])
        query = matrix.query().operation
        b.call(query, address=b["address"], data=b["data"], resources={"entries": "entries"})
        self.compare(b.finish().program(), {"entries": matrix.snapshot()["entries"]})

    def test_row_and_amplitude_preparation(self):
        # 地址计算私有工作区固定占 16 位；4 位角字的行制备共 24 位。
        self.check_row_and_amplitude_preparation(4, self.compare)

    def test_wide_row_and_amplitude_preparation_sparse_backends(self):
        # 保留原 10 位角字及全部输入；行制备共 30 位，采用真实稀疏后端。
        self.check_row_and_amplitude_preparation(AW, self.compare_sparse)

    def check_row_and_amplitude_preparation(self, angle_width, compare):
        matrix = QMatrix(
            [
                [0.5, 0.5, 0.25, 0.0],
                [0.5, 0.0, 0.25, 0.5],
                [0.25, 0.5, 0.0, 0.5],
                [0.0, 0.0, 0.5, 0.5],
            ],
            fmt=FMT,
            angle_width=angle_width,
        )
        snapshot = matrix.snapshot()
        row = matrix.row_preparation()
        b = Builder(
            "row_driver",
            {r.name: r.type for r in row.module.registers},
            {r.name: r.type for r in row.module.resources},
        )
        b.h(b["row"])
        b.call(row, row=b["row"], item=b["item"], work=b["work"], resources={"row_angles": "row_angles"})
        compare(b.finish().program(), {"row_angles": snapshot["row_angles"]})

        amp = matrix.amplitude_preparation()
        c = Builder(
            "amp_driver",
            {r.name: r.type for r in amp.module.registers},
            {r.name: r.type for r in amp.module.resources},
        )
        c.h(c["item"])
        c.call(amp, row=c["row"], item=c["item"], work=c["work"], resources={"root_angles": "root_angles"})
        compare(c.finish().program(), {"root_angles": snapshot["root_angles"]})


if __name__ == "__main__":
    unittest.main()
