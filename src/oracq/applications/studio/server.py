"""Local compiler studio transport; importing this module starts no service."""
from __future__ import annotations

import hmac
import json
import secrets
import subprocess
import sys
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib.metadata import version
from importlib.resources import files
from typing import cast

from oracq.applications.studio.compiler import CompilationRequest


class Handler(BaseHTTPRequestHandler):
    """Serve the packaged page and compile a validated local-session request."""

    @property
    def studio(self) -> StudioServer:
        """Return the server's typed session state."""
        return cast("StudioServer", self.server)

    def log_message(self, format: str, *args: object) -> None:
        """Keep routine browser requests out of compiler output."""

    def handle(self) -> None:
        """A browser may disconnect while its compilation is finishing."""
        try:
            super().handle()
        except ConnectionError:
            pass

    def local(self) -> bool:
        """Require the exact loopback origin, including its selected port."""
        host = f"127.0.0.1:{self.studio.server_port}"
        return (self.headers.get("Host") == host
                and self.headers.get("Origin", f"http://{host}") == f"http://{host}")

    def reply(self, status: int, payload: object, mime: str = "application/json") -> None:
        """Write a bounded response with a session-specific script policy."""
        data = payload if isinstance(payload, bytes) else json.dumps(payload).encode()
        self.send_response(status)
        self.send_header("Content-Type", mime+"; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy", "default-src 'none'; connect-src 'self'; "
                         f"script-src 'nonce-{self.studio.token}'; style-src 'unsafe-inline'; "
                         "frame-ancestors 'none'; base-uri 'none'")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self) -> None:
        """Read the packaged UI or current software version."""
        if not self.local():
            self.reply(403, {"error": "local origin required"})
        elif self.path == "/":
            html = files("oracq.applications.studio").joinpath("index.html").read_text(encoding="utf-8")
            self.reply(200, html.replace("__SESSION__", self.studio.token).encode(), "text/html")
        elif self.path == "/api/info":
            self.reply(200, {"software": "OracQ", "version": version("oracq")})
        else:
            self.reply(404, {"error": "not found"})

    def do_POST(self) -> None:
        """Compile in one isolated worker; never execute user source as Python."""
        if (not self.local() or not hmac.compare_digest(
                self.headers.get("X-Studio-Session", ""), self.studio.token)):
            self.reply(403, {"error": "local studio session required"})
            return
        if self.path != "/api/compile":
            self.reply(404, {"error": "not found"})
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= 16384 or self.headers.get("Content-Type") != "application/json":
                raise ValueError("a JSON request of at most 16 KiB is required")
            data = self.rfile.read(length)
            CompilationRequest(**json.loads(data)).validate()
        except (ValueError, TypeError) as exc:
            self.reply(400, {"error": str(exc)})
            return
        if not self.studio.busy.acquire(blocking=False):
            self.reply(409, {"error": "a compilation is already running"})
            return
        try:
            result = subprocess.run(
                [sys.executable, "-m", "oracq.applications.studio.worker"], input=data,
                capture_output=True, timeout=60,
            )
            if result.returncode:
                self.reply(422, {"error": result.stderr.decode(errors="replace")[-2048:]})
            else:
                self.reply(200, json.loads(result.stdout))
        except ConnectionError:
            return
        except subprocess.TimeoutExpired:
            self.reply(504, {"error": "compilation exceeded 60 seconds; use compile-function"})
        except (OSError, ValueError) as exc:
            self.reply(500, {"error": str(exc)})
        finally:
            self.studio.busy.release()


class StudioServer(ThreadingHTTPServer):
    """A loopback-only session with one compilation slot."""

    def __init__(self, port: int = 0):
        if type(port) is not int or not 0 <= port <= 65535:
            raise ValueError("port must be an integer in [0,65535]")
        self.token = secrets.token_urlsafe(32)
        self.busy = threading.Lock()
        super().__init__(("127.0.0.1", port), Handler)


def serve(*, port: int = 0, open_browser: bool = True) -> None:
    """Start the independent local compiler UI until interrupted."""
    with StudioServer(port) as server:
        url = f"http://127.0.0.1:{server.server_port}/"
        print(url, flush=True)
        if open_browser:
            webbrowser.open(url)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
