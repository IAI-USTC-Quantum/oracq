"""Refresh the workboard from implementation evidence and case reports."""

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
        "# Paradigm Implementation Workboard",
        "",
        f"Updated: {board['updated']}. This round accepts expressiveness, open IR, "
        "binding, and backend-description reachability; mathematical correctness "
        "is deferred to the next stage.",
        "",
        "| ID | Work item | Status | Evidence |",
        "|---|---|---|---|",
    ]
    for task in board["tasks"]:
        evidence = ", ".join(f"[{Path(p).name}](../{p})" for p in task.get("evidence", []))
        lines.append(f"| {task['id']} | {task['title']} | {task['status']} | {evidence} |")
    lines += [
        "",
        "## Case artifacts",
        "",
        "Each directory provides open.rir.yaml, partial.rir.yaml, closed.rir.yaml, "
        "bindings.json, memory.json, and program.originir. They are regenerable "
        "description artifacts stored in the Git-ignored out/catalog.",
        "",
        "| Case | Open slots | Bound modules | OriginIR DEF | Native parse | Artifacts |",
        "|---|---:|---:|---:|---|---|",
    ]
    for row in cases:
        name = row["case"]
        links = " / ".join(
            f"[{label}](../out/catalog/{name}/{file})"
            for label, file in [
                ("open", "open.rir.yaml"),
                ("partial", "partial.rir.yaml"),
                ("closed", "closed.rir.yaml"),
                ("OriginIR", "program.originir"),
                ("report", "report.json"),
            ]
        )
        lines.append(
            f"| {name} | {len(row['unresolved'])} | {row['closed_modules']} | "
            f"{row['originir_definitions']} | {'passed' if row['native_parsed'] else 'pending'} | {links} |"
        )
    lines += [
        "",
        "## Coverage and limitations",
        "",
        "- The [legacy case mapping](coverage.md) covers 61 positive and 16 negative cases; it is not a proof of equivalence with the legacy source/CLIR.",
        "- The concrete QFVM samples use toy physics lookup kernels; the full Roe operator can remain an open declaration.",
        "- The current concrete QHAM assembly is the m=1 lifting; higher-order mathematical constructions are deferred to the next stage.",
        "- Current bindings target interfaces with fixed widths and alpha; changing those constants requires rerunning the Python generators.",
        "- All applications are marked correctness=not_assessed; assembly is not blocked by pending accuracy or convergence proofs.",
        "",
        "## Next-stage validation order",
        "",
        "Check matrix and bit semantics of plain oracles first, then block encodings "
        "and state preparation, followed by the Costa initial state/reflection/filtering, "
        "and finally validate QODE/PDE, the Roe physics kernels, and QHAM's conditional "
        "outputs and outer-loop behavior.",
        "",
    ]
    (ROOT / "docs/archive/workboard.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
