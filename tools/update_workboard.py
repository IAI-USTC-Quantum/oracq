"""根据实现证据与案例报告刷新工作面板。"""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    path = ROOT / "docs/archive/workboard.json"
    board = json.loads(path.read_text())
    board["updated"] = "2026-09-09"
    evidence = {
        "P0": ["docs/archive/coverage.json", "docs/archive/coverage.md"],
        "P1": ["src/oracq/linking.py", "docs/reference/open-ir.md", "tests/core/test_open_ir.py"],
        "P2": ["src/oracq/oracles.py", "docs/archive/oracle-paradigms.md"],
        "P3": ["src/oracq/access.py", "src/oracq/combinators.py"],
        "P4": ["src/oracq/algorithms/basics/oracle_algorithms.py"],
        "P5": ["src/oracq/algorithms/qlss/qlss.py"],
        "P6": ["src/oracq/applications.py", "out/catalog/qfvm_qram/program.originir"],
        "P7": ["src/oracq/algorithms/qode/legacy.py", "out/catalog/qham_qode/program.originir"],
        "P8": [
            "tools/build_catalog.py",
            "tests/core/test_workloads.py",
            "tests/integration/test_catalog.py",
        ],
    }
    for task in board["tasks"]:
        task["status"] = "done" if task["id"] != "P8" else "in_progress"
        task["evidence"] = evidence[task["id"]]
    index_path = ROOT / "out/catalog/index.json"
    cases = json.loads(index_path.read_text())["cases"] if index_path.exists() else []
    verification_path = ROOT / "docs/archive/phase-validation.json"
    verification = json.loads(verification_path.read_text()) if verification_path.exists() else {}
    if (
        cases
        and all(row["native_parsed"] for row in cases)
        and verification.get("status") == "passed"
    ):
        board["tasks"][-1]["status"] = "done"
        board["tasks"][-1]["evidence"].append("docs/archive/phase-validation.json")
    board["case_count"] = len(cases)
    board["native_parsed_count"] = sum(row["native_parsed"] for row in cases)
    path.write_text(json.dumps(board, ensure_ascii=False, indent=2) + "\n")
    render(board, cases)


def render(board, cases):
    lines = [
        "# 范式实现工作面板",
        "",
        f"更新时间：{board['updated']}。本轮验收表达、开放 IR、绑定和后端描述可达性；数学正确性留待下一阶段。",
        "",
        "| 编号 | 工作项 | 状态 | 证据 |",
        "|---|---|---|---|",
    ]
    for task in board["tasks"]:
        evidence = "、".join(f"[{Path(p).name}](../{p})" for p in task.get("evidence", []))
        lines.append(f"| {task['id']} | {task['title']} | {task['status']} | {evidence} |")
    lines += [
        "",
        "## 案例产物",
        "",
        "每个目录提供 open.rir.yaml、partial.rir.yaml、closed.rir.yaml、bindings.json、memory.json 和 program.originir。"
        "它们是可重新生成的描述产物，保存在被 Git 忽略的 out/catalog。",
        "",
        "| 案例 | 开放槽 | 绑定后模块 | OriginIR DEF | 原生解析 | 产物 |",
        "|---|---:|---:|---:|---|---|",
    ]
    for row in cases:
        name = row["case"]
        links = " / ".join(
            f"[{label}](../out/catalog/{name}/{file})"
            for label, file in [
                ("开放", "open.rir.yaml"),
                ("部分绑定", "partial.rir.yaml"),
                ("闭合", "closed.rir.yaml"),
                ("OriginIR", "program.originir"),
                ("报告", "report.json"),
            ]
        )
        lines.append(
            f"| {name} | {len(row['unresolved'])} | {row['closed_modules']} | "
            f"{row['originir_definitions']} | {'通过' if row['native_parsed'] else '待验收'} | {links} |"
        )
    lines += [
        "",
        "## 覆盖与限制",
        "",
        "- [旧案例逐项映射](coverage.md)包含 61 个正例和 16 个负例；这里不是旧源码/CLIR 的等价证明。",
        "- QFVM 具体样例使用玩具物理查表核；完整 Roe 可以继续保留为开放声明。",
        "- QHAM 当前具体组装是 m=1 提升；更高阶数学构造留在下一阶段。",
        "- 当前绑定针对固定宽度和 alpha 的接口；更换这些常量需要重新运行 Python 生成器。",
        "- 所有应用均标记 correctness=not_assessed；没有因为尚未完成精度或收敛证明而阻断组装。",
        "",
        "## 下一阶段验证顺序",
        "",
        "先检查普通 oracle 的矩阵和位语义，再检查 BE 与 state prep，随后是 Costa 初态/反射/filtering，最后验证 QODE/PDE、Roe 物理核与 QHAM 的条件输出和外层行为。",
        "",
    ]
    (ROOT / "docs/archive/workboard.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
