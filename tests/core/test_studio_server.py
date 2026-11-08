"""Exercise the local packaged page and isolated real compiler through HTTP."""
import json
import re
import subprocess
import sys
import threading
import unittest
import urllib.error
import urllib.request

from oracq import loads
from oracq.applications.studio.server import StudioServer


class StudioServerTests(unittest.TestCase):
    def test_local_session_and_real_compiler_worker(self):
        with StudioServer() as server:
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            base = f"http://127.0.0.1:{server.server_port}"
            try:
                with urllib.request.urlopen(base) as response:
                    html = response.read().decode()
                    self.assertIn("frame-ancestors 'none'", response.headers["Content-Security-Policy"])
                token = re.search(r"const session='([^']+)'", html).group(1)
                data = json.dumps({"source": "def f(x):\n return x+1", "width": 3,
                                   "fraction": 0, "degree": 2, "checkpoint_size": 2}).encode()

                def post(headers):
                    request = urllib.request.Request(base+"/api/compile", data=data,
                        headers={"Content-Type": "application/json", **headers})
                    return urllib.request.urlopen(request)

                with self.assertRaises(urllib.error.HTTPError) as missing:
                    post({})
                self.assertEqual(missing.exception.code, 403)
                with self.assertRaises(urllib.error.HTTPError) as foreign:
                    post({"X-Studio-Session": token, "Origin": "https://example.invalid"})
                self.assertEqual(foreign.exception.code, 403)
                with post({"X-Studio-Session": token}) as response:
                    report = json.load(response)
                self.assertEqual(report["schema"], "oracq.studio.v1")
                self.assertTrue(report["variants"]["selected"]["complete"])
                self.assertTrue(loads(report["artifacts"]["program.rir.json"]).main.registers)
                self.assertFalse(report["scope"]["hardware_run"])
            finally:
                server.shutdown()
                thread.join(timeout=3)

    def test_cli_requires_paths_and_rejects_invalid_gui_port(self):
        for args in (["estimate"], ["gui", "--port", "-1"], ["gui", "unexpected"]):
            with self.subTest(args=args):
                result = subprocess.run([sys.executable, "-m", "oracq", *args],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 2)
                self.assertNotIn("Traceback", result.stderr)
