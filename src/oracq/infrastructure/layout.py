"Sequential-reuse layout of module-private workspaces."

from __future__ import annotations

from oracq.infrastructure.ir import Program
from oracq.infrastructure.linking import calls


def workspace_table(program: Program) -> dict[str, int]:
    """Compute the private workspace bit count needed by each entry-reachable module.

    The count is the sum of the module's own ``locals`` widths plus the maximum
    workspace bit count of its callees; sequential calls reuse the same stretch
    of physical workspace.

    Args:
        program: ``Program`` whose call graph is acyclic.

    Returns:
        dict: mapping from module names to bit counts; modules unreachable from the entry do not appear.
    """
    result: dict[str, int]
    modules, result = program.module_map, {}

    def visit(key: str) -> int:
        """Recursively compute the module's own ``locals`` width plus the maximum workspace of its callee chain."""
        if key in result:
            return result[key]
        module = modules[key]
        own = sum(r.type.width for r in module.locals)
        child = max((visit(call.module) for call in calls(module.body)), default=0)
        result[key] = own + child
        return result[key]

    visit(program.entry)
    return result
