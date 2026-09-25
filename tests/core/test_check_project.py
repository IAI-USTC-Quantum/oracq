"""The acceptance tool must keep virtual-environment paths instead of dereferencing to the base interpreter."""

import importlib.util
import subprocess
import sys
import tempfile
import unittest
import venv
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "project_checker", Path(__file__).resolve().parents[2] / "tools/check_project.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
backend_interpreter = module.backend_interpreter


class ProjectCheckTests(unittest.TestCase):
    def test_backend_interpreter_keeps_virtual_environment(self):
        with tempfile.TemporaryDirectory() as directory:
            env = Path(directory) / "backend"
            venv.EnvBuilder(with_pip=False, symlinks=True).create(env)
            python = env / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
            command = backend_interpreter(python)
            result = subprocess.run(
                [command, "-c", "import sys; print(sys.prefix)"],
                check=True,
                text=True,
                capture_output=True,
            )
            self.assertEqual(Path(result.stdout.strip()), env)
