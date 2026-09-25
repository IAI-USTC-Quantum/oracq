"""Run all paper-grade numerical validation groups; artifacts are written to out/verification/.

Requires an interpreter with both pysparq and uniqc installed (a real
backend, no mock substitutes):

    PYTHONPATH=src <python with pysparq+uniqc> tools/run_verification.py [--group arithmetic]...
"""

import argparse
import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERIFY_DIR = ROOT / "tests" / "verification"


def discover():
    return sorted(path.stem.removeprefix("verify_") for path in VERIFY_DIR.glob("verify_*.py"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", action="append", choices=discover(), help="run only the specified groups")
    args = parser.parse_args()
    groups = args.group or discover()
    sys.path.insert(0, str(ROOT / "src"))
    sys.path.insert(0, str(VERIFY_DIR))
    reports = []
    for group in groups:
        print(f"== verify_{group}", flush=True)
        module = importlib.import_module(f"verify_{group}")
        reports.append(module.run())
    total = sum(len(report.cases) for report in reports)
    print(f"verification complete: {len(reports)} groups, {total} cases", flush=True)


if __name__ == "__main__":
    main()
