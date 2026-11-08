"""Isolated bounded compilation worker used by the local studio server."""
from __future__ import annotations

import json
import sys

from oracq.applications.studio.compiler import CompilationRequest, compile_comparison


def main() -> None:
    """Read one JSON request from stdin and emit one compiler report."""
    try:
        import resource
    except ImportError:  # Windows still has the server's wall-clock deadline.
        pass
    else:
        resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    try:
        data = sys.stdin.buffer.read(16385)
        if len(data) > 16384:
            raise ValueError("studio request exceeds 16 KiB")
        request = CompilationRequest(**json.loads(data))
        print(json.dumps(compile_comparison(request), allow_nan=False))
    except (ValueError, TypeError, OSError, MemoryError) as exc:
        print(f"studio compilation failed: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
