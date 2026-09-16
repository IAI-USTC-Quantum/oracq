"""求解器栈 14 个输入模型案例的 Toffoli+Clifford+T+QRAM 资源台账。

复用 examples/input_models.py 生成的闭合 RIR 产物（out/input-models/），
逐案例做组合式资源估计：门级成本与 QRAM 查询数在同一问题、不同数据
路径（结构化门端口 / 谱嵌入 / QRAM）下的对照。若产物缺失先运行
`PYTHONPATH=src python examples/input_models.py`。

运行：PYTHONPATH=src python tools/build_qram_queries.py
"""

import json
from pathlib import Path

from pyqecclang import estimate_resources, loads

CASES = Path("out/input-models")
OUT = Path("out/resource-estimates")


def main():
    if not CASES.is_dir():
        raise SystemExit("缺少 out/input-models 产物，请先运行 examples/input_models.py")
    records = []
    for folder in sorted(CASES.iterdir()):
        closed = folder / "closed.rir.json"
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
    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "group": "solver_qram_queries",
        "note": "out/input-models 14 个闭合案例的资源台账；同一问题不同数据路径的对照",
        "records": records,
    }
    (OUT / "solver_qram_queries.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"[artifact] {OUT / 'solver_qram_queries.json'}", flush=True)


if __name__ == "__main__":
    main()
