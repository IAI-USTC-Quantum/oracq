"""新应用只做真实 OriginIR-ext 语法验收，不进行算法数值认证。"""

import unittest

from oracq.applications.catalog import build_case


class CatalogSyntaxTests(unittest.TestCase):
    def test_application_descriptions_parse_in_real_backend(self):
        from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser

        for name in (
            "dj_qram",
            "grover_qram",
            "costa_qram",
            "qfvm_qram",
            "qham_qode",
            "qham_qpde",
            "sparse_qram",
            "measure_reset",
            "banked_qram",
        ):
            with self.subTest(case=name):
                case = build_case(name)
                artifact = case.artifact()
                parser = OriginIR_BaseParser()
                parser.parse(artifact.text)
                if case.readout:
                    from uniqc.circuit_builder.classical_program import parse_originir_ext_dynamic

                    from oracq.infrastructure.readout import export_with_readout

                    execution = export_with_readout(case.closed(), case.readout)
                    self.assertIsNotNone(parse_originir_ext_dynamic(execution.text))
                self.assertEqual(
                    parser.n_qubit, sum(r.type.width for r in case.closed().main.registers)
                )
