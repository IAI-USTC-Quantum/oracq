"""集中运行工程检查；真实后端必须由显式解释器提供，失败即退出。"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def backend_interpreter(path):
    """保留虚拟环境解释器的路径；resolve 会把它解引用为基础 Python。"""
    return os.path.abspath(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend-python", type=Path)
    parser.add_argument("--docs", action="store_true", help="构建 Sphinx HTML 并执行教程")
    parser.add_argument("--output", type=Path, default=ROOT / "out/checks/latest")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    commands = [
        ("lint", [sys.executable, "-m", "ruff", "check", "src", "tests", "examples", "tools"]),
        ("core_schema", [sys.executable, "-m", "pytest", "tests/core", "tests/schema", "-q"]),
        (
            "contracts_example",
            [sys.executable, "examples/algorithm_contracts.py", "-o", str(output / "contracts")],
        ),
        (
            "ode_examples",
            [sys.executable, "examples/ode_input_models.py", "-o", str(output / "ode")],
        ),
        (
            "input_models_example",
            [sys.executable, "examples/input_models.py", "-o", str(output / "input-models")],
        ),
        (
            "algorithm_gallery",
            [sys.executable, "examples/algorithm_gallery.py", "-o", str(output / "gallery")],
        ),
        (
            "resource_estimates",
            [sys.executable, "tools/build_resource_estimates.py"],
        ),
        (
            "qram_queries",
            [sys.executable, "tools/build_qram_queries.py"],
        ),
    ]
    if args.docs:
        commands.extend(
            [
                (
                    "docs_search",
                    [sys.executable, "-m", "unittest", "discover", "-s", "tests/docs", "-v"],
                ),
                (
                    "sphinx_html",
                    [
                        sys.executable,
                        "-m",
                        "sphinx",
                        "-W",
                        "--keep-going",
                        "-b",
                        "html",
                        "docs",
                        str(output / "docs/html"),
                    ],
                ),
                (
                    "sphinx_doctest",
                    [
                        sys.executable,
                        "-m",
                        "sphinx",
                        "-W",
                        "--keep-going",
                        "-b",
                        "doctest",
                        "docs",
                        str(output / "docs/doctest"),
                    ],
                ),
            ]
        )
    if args.backend_python is not None:
        commands.append(
            (
                "native",
                [
                    backend_interpreter(args.backend_python),
                    "-B",
                    "-m",
                    "unittest",
                    "discover",
                    "-s",
                    "tests/integration",
                    "-v",
                ],
            )
        )
        commands.append(
            (
                "verification",
                [
                    backend_interpreter(args.backend_python),
                    "-B",
                    "tools/run_verification.py",
                ],
            )
        )
    results = {
        "native_requested": args.backend_python is not None,
        "docs_requested": args.docs,
        "checks": [],
    }
    for label, command in commands:
        print("Checking", label, flush=True)
        result = subprocess.run(command, cwd=ROOT, env=env, check=False)
        results["checks"].append({"name": label, "exit_code": result.returncode})
        (output / "report.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
        if result.returncode:
            raise SystemExit(result.returncode)
    print("Checks passed; native requested:", results["native_requested"])


if __name__ == "__main__":
    main()
