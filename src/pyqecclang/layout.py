"""模块私有工作区的顺序复用布局。"""

from .linking import calls


def workspace_table(program):
    modules, result = program.module_map, {}

    def visit(key):
        if key in result:
            return result[key]
        module = modules[key]
        own = sum(r.type.width for r in module.locals)
        child = max((visit(call.module) for call in calls(module.body)), default=0)
        result[key] = own + child
        return result[key]

    visit(program.entry)
    return result
