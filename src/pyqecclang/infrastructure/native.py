"PySparQ 模块级实现注册；运行能力不冒充门级闭合。"

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from pyqecclang.infrastructure.builder import Operation
from pyqecclang.infrastructure.ir import Module, Ref, ValidationError
from pyqecclang.infrastructure.linking import calls
from pyqecclang.infrastructure.validation import validate


@dataclass(frozen=True)
class NativeSite:
    module: Module
    arguments: tuple[Ref, ...]
    resources: tuple[str, ...]


@dataclass(frozen=True)
class NativeEntry:
    module: Module
    factory: Callable
    label: str
    kind: str


class NativeRegistry:
    def __init__(self):
        self.entries = {}

    def register(
        self, operation: Operation | Module, factory, *, label=None, kind="pysparq_custom"
    ):
        module = operation.module if isinstance(operation, Operation) else operation
        entry = NativeEntry(module, factory, label or module.name, kind)
        if module.name in self.entries and self.entries[module.name] != entry:
            raise ValidationError("原生实现名字重复")
        self.entries[module.name] = entry
        return self

    def matching(self, program):
        result = set()
        for module in program.modules:
            if module.name in self.entries:
                if self.entries[module.name].module != module:
                    raise ValidationError(f"原生实现与当前模块描述不匹配：{module.name}")
                result.add(module.name)
        return frozenset(result)

    def missing(self, program):
        validate(program)
        supported = self.matching(program)
        result, visited = [], set()

        def visit(key):
            if key in visited or key in supported:
                return
            visited.add(key)
            module = program.module_map[key]
            if module.body is None:
                result.append(key)
            else:
                for call in calls(module.body):
                    visit(call.module)

        visit(program.entry)
        return tuple(sorted(result))


@dataclass
class NativeContext:
    ps: object
    site: NativeSite
    names: dict
    qrams: dict
    memories: dict

    def register_id(self, parameter):
        index = next(i for i, r in enumerate(self.site.module.registers) if r.name == parameter)
        ref = self.site.arguments[index]
        if len(ref.parts) != 1:
            raise ValidationError("当前原生工厂需要完整寄存器参数")
        span = ref.parts[0]
        name = self.names[span.register]
        if span.start != 0 or span.width != self.ps.System.size_of(name):
            raise ValidationError("当前原生工厂不接收切片；请使用完整寄存器或门级实现")
        return self.ps.System.get_id(name)

    def qram(self, parameter):
        index = next(i for i, r in enumerate(self.site.module.resources) if r.name == parameter)
        return self.qrams[self.site.resources[index]]


class DynamicCppFactory:
    """按需调用真实 compile_operator，C++ 代码不写入 RIR。"""

    def __init__(self, name, source, parameters, *, cache_dir, base_class="SelfAdjointOperator"):
        self.name, self.source = name, source
        self.parameters = tuple(parameters)
        self.cache_dir, self.base_class = str(cache_dir), base_class
        self._class = None

    def __call__(self, context):
        if self._class is None:
            from pysparq.dynamic_operator import compile_operator

            self._class = compile_operator(
                name=self.name,
                cpp_code=self.source,
                base_class=self.base_class,
                constructor_args=[("size_t", name) for name in self.parameters],
                cache_dir=self.cache_dir,
            )
        return self._class(**{name: context.register_id(name) for name in self.parameters})
