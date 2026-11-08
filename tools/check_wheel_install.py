"""Install one real wheel in a clean environment and test its compiler studio.

Usage: python tools/check_wheel_install.py dist/oracq-*.whl
Requires the same uv executable as the project's build workflow. No quantum
backend is installed or replaced by this artifact test.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

PACKAGE_CHECK = """
import importlib.metadata as metadata
import importlib.resources as resources
import pathlib, sys
import oracq
root = pathlib.Path(sys.argv[1])
assert metadata.version('oracq') == sys.argv[2]
assert not pathlib.Path(oracq.__file__).resolve().is_relative_to(root)
expected = (root/'src/oracq/applications/studio/index.html').read_bytes()
actual = resources.files('oracq.applications.studio').joinpath('index.html').read_bytes()
assert actual == expected, 'installed compiler UI differs from source'
print('PASS installed OracQ', metadata.version('oracq'), 'from', oracq.__file__)
"""


def main() -> None:
    """Use an isolated interpreter and cwd so the source checkout cannot shadow it."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    expected = tomllib.loads((root/"pyproject.toml").read_text())["project"]["version"]
    wheel = args.wheel.resolve(strict=True)
    if wheel.suffix != ".whl":
        parser.error("expected one built .whl file")
    env = os.environ.copy()
    for name in ("PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV"):
        env.pop(name, None)
    with tempfile.TemporaryDirectory(prefix="oracq-wheel-") as directory:
        work = Path(directory)
        target = work/"env"
        subprocess.run(["uv", "venv", str(target), "--python", sys.executable], check=True, env=env)
        binary = target/("Scripts" if os.name == "nt" else "bin")
        python = binary/("python.exe" if os.name == "nt" else "python")
        subprocess.run(["uv", "pip", "install", "--python", str(python), str(wheel)],
                       check=True, env=env)
        subprocess.run([str(python), "-c", PACKAGE_CHECK, str(root), expected],
                       cwd=work, env=env, check=True)
        subprocess.run([str(binary/("oracq.exe" if os.name == "nt" else "oracq")), "--help"],
                       cwd=work, env=env, check=True, stdout=subprocess.DEVNULL)
        subprocess.run([str(python), "-m", "unittest", "discover", "-s",
                        str(root/"tests/core"), "-p", "test_studio*.py", "-v"],
                       cwd=work, env=env, check=True)


if __name__ == "__main__":
    main()
