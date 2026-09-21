"模块私有工作区的顺序复用布局。"

from __future__ import annotations

from pyqecclang.infrastructure.ir import Program
from pyqecclang.infrastructure.linking import calls


def workspace_table(program: Program) -> dict[str, int]:
    """计算每个入口可达模块所需的私有工作区位数。

    位数为模块自身 ``locals`` 宽度之和加上其被调模块工作区位数的最大值；
    顺序调用复用同一段物理工作区。

    Args:
        program: 调用图无环的 ``Program``。

    Returns:
        dict: 模块名到位数的映射；入口不可达的模块不出现在表中。
    """
    result: dict[str, int]
    modules, result = program.module_map, {}

    def visit(key: str) -> int:
        """递归求模块自身 ``locals`` 宽度与被调链工作区最大值之和。"""
        if key in result:
            return result[key]
        module = modules[key]
        own = sum(r.type.width for r in module.locals)
        child = max((visit(call.module) for call in calls(module.body)), default=0)
        result[key] = own + child
        return result[key]

    visit(program.entry)
    return result
