"""描述导出与可选执行后端的公共入口。"""

from pyqecclang.infrastructure.backends.originir import (
    OriginIRArtifact,
    export_originir,
    run_originir,
)

__all__ = ["OriginIRArtifact", "export_originir", "run_originir"]
