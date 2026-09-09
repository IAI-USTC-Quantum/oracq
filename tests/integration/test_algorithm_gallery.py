"""算法分类目录通过真实 PySparQ / OriginIR 的复幅度对拍。"""

import unittest

from pyqecclang import run_originir, run_pysparq, simulate
from pyqecclang.applications.gallery import algorithm_gallery


class AlgorithmGalleryNativeTests(unittest.TestCase):
    def test_every_gallery_case_on_both_real_backends(self):
        for case in algorithm_gallery():
            with self.subTest(algorithm=case.name):
                program = case.operation.program()
                expected = simulate(program).amplitudes
                native = run_pysparq(program).amplitudes
                for key in expected.keys() | native.keys():
                    self.assertAlmostEqual(expected.get(key, 0), native.get(key, 0), places=10)
                for index, amplitude in enumerate(run_originir(program)):
                    key, cursor = [], 0
                    for reg in program.main.registers:
                        key.append((index >> cursor) & ((1 << reg.type.width) - 1))
                        cursor += reg.type.width
                    self.assertAlmostEqual(
                        expected.get(tuple(key), 0), complex(amplitude), places=10
                    )
