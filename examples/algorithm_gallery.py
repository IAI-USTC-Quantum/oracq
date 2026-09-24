"""生成算法展示目录；--native 使用真实后端比较完整复幅度。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from oracq import (
    dumps,
    export_originir,
    export_toffoli_u3_cz,
    run_originir,
    run_pysparq,
    simulate,
)
from oracq.applications.gallery import algorithm_gallery


def main() -> None:
    """遍历算法展示目录，写出各案例的导出产物并汇总索引。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path, default=Path("out/algorithm-gallery"))
    parser.add_argument("--native", action="store_true")
    args = parser.parse_args()
    records = []
    for case in algorithm_gallery():
        directory = args.output / case.name
        directory.mkdir(parents=True, exist_ok=True)
        program = case.operation.program()
        (directory / "closed.rir.yaml").write_text(dumps(program), encoding="utf-8")
        if case.opened is not None:
            (directory / "open.rir.yaml").write_text(dumps(case.opened.program()), encoding="utf-8")
        (directory / "modular.originir").write_text(export_originir(program).text, encoding="utf-8")
        (directory / "toffoli_u3_cz.originir").write_text(
            export_toffoli_u3_cz(program).text, encoding="utf-8"
        )
        expected = simulate(program).amplitudes
        difference = None
        if args.native:
            native = run_pysparq(program).amplitudes
            difference = max(
                abs(expected.get(k, 0) - native.get(k, 0)) for k in expected.keys() | native.keys()
            )
            vector = run_originir(program)
            for index, amplitude in enumerate(vector):
                cursor, key = 0, []
                for reg in program.main.registers:
                    key.append((index >> cursor) & ((1 << reg.type.width) - 1))
                    cursor += reg.type.width
                difference = max(difference, abs(complex(amplitude) - expected.get(tuple(key), 0)))
            if difference > 1e-10:
                raise AssertionError((case.name, difference))
        records.append(
            {
                "name": case.name,
                "family": case.family,
                "readout": case.readout,
                "modules": len(program.modules),
                "native_max_difference": difference,
            }
        )
        print(case.name, "passed", flush=True)
    (args.output / "index.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
