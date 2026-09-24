"""求解器栈 14 个输入模型案例的 Toffoli+Clifford+T+QRAM 资源台账。

复用 examples/input_models.py 生成的闭合 RIR 产物（out/input-models/），
逐案例做组合式资源估计：门级成本与 QRAM 查询数在同一问题、不同数据
路径（结构化门端口 / 谱嵌入 / QRAM）下的对照。若产物缺失先运行
`PYTHONPATH=src python examples/input_models.py`。

运行：PYTHONPATH=src python tools/build_qram_queries.py [--cases out/input-models]
"""

import argparse
import json
from pathlib import Path

from oracq import estimate_resources, loads

CASES = Path("out/input-models")
OUT = Path("out/resource-estimates")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=Path, default=CASES, help="闭合 RIR 案例目录")
    parser.add_argument(
        "--out",
        type=Path,
        default=OUT / "solver_qram_queries.json",
        help="台账产物路径",
    )
    args = parser.parse_args()
    if not args.cases.is_dir():
        raise SystemExit(f"缺少 {args.cases} 产物，请先运行 examples/input_models.py")
    records = []
    for folder in sorted(args.cases.iterdir()):
        closed = folder / "closed.rir.yaml"
        if not closed.is_file():
            continue
        program = loads(closed.read_text(encoding="utf-8"))
        report = json.loads((folder / "report.json").read_text(encoding="utf-8"))
        estimate = estimate_resources(program)
        record = {
            "case": folder.name,
            "family": folder.name.split("_", 1)[0],
            "resources": report["resources"],
            "modules": report["modules"],
            **estimate.to_dict(),
        }
        records.append(record)
        print(
            f"{folder.name:28s} qubits={estimate.qubits:4d} toffoli={estimate.toffoli:7d} "
            f"t_exact={estimate.t_exact:4d} rot={len(estimate.rotations):5d} "
            f"qram={estimate.qram_total:3d} {dict(estimate.qram_queries)}",
            flush=True,
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "group": "solver_qram_queries",
        "note": "out/input-models 14 个闭合案例的资源台账；同一问题不同数据路径的对照",
        "records": records,
    }
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"[artifact] {args.out}", flush=True)


if __name__ == "__main__":
    main()
