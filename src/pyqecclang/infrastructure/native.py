"PySparQ 模块级实现注册；运行能力不冒充门级闭合。"

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from types import ModuleType
from typing import cast

from pyqecclang.infrastructure.builder import Operation
from pyqecclang.infrastructure.ir import Module, Program, Ref, ValidationError
from pyqecclang.infrastructure.linking import calls
from pyqecclang.infrastructure.validation import validate


@dataclass(frozen=True)
class NativeSite:
    """执行流中由原生实现接管的一个模块调用点。

    Attributes:
        module: 被调用的模块声明。
        arguments: 与模块寄存器一一对应的实参视图。
        resources: 与模块资源一一对应的资源绑定名。
    """

    module: Module
    arguments: tuple[Ref, ...]
    resources: tuple[str, ...]


@dataclass(frozen=True)
class NativeEntry:
    """一个模块名下的原生实现登记项。

    Attributes:
        module: 登记时声明的模块描述，用于与程序中的同名模块比对。
        factory: 执行期才以 ``NativeContext`` 为参调用的实现工厂。
        label: 执行报告中显示的实现名。
        kind: 实现类别标记，默认 ``pysparq_custom``。
    """

    module: Module
    factory: Callable
    label: str
    kind: str


class NativeRegistry:
    """按模块名登记 PySparQ 模块级原生实现，并对照程序核对覆盖情况。"""

    def __init__(self) -> None:
        """初始化空的模块名到原生实现登记表。"""
        self.entries: dict[str, NativeEntry] = {}

    def register(
        self,
        operation: Operation | Module,
        factory: Callable[[NativeContext], object],
        *,
        label: str | None = None,
        kind: str = "pysparq_custom",
    ) -> NativeRegistry:
        """登记一个模块的原生实现。

        Args:
            operation: 提供模块声明的 ``Operation`` 或 ``Module``。
            factory: 执行期以 ``NativeContext`` 为参调用的工厂；C++ 编译等
                重活在工厂内部按需延迟进行，注册阶段不做。
            label: 执行报告中的实现显示名；省略时用模块名。
            kind: 实现类别标记。

        Returns:
            NativeRegistry: 返回自身，支持链式登记。

        Raises:
            ValidationError: 模块名已登记且登记内容不同。
        """
        module = operation.module if isinstance(operation, Operation) else operation
        entry = NativeEntry(module, factory, label or module.name, kind)
        if module.name in self.entries and self.entries[module.name] != entry:
            raise ValidationError("原生实现名字重复")
        self.entries[module.name] = entry
        return self

    def matching(self, program: Program) -> frozenset[str]:
        """求程序中被原生实现接管的模块名集合。

        Args:
            program: 待核对的程序。

        Returns:
            frozenset[str]: 与登记条目同名且模块声明一致的模块名。

        Raises:
            ValidationError: 同名模块与登记时的声明不匹配。
        """
        result = set()
        for module in program.modules:
            if module.name in self.entries:
                if self.entries[module.name].module != module:
                    raise ValidationError(f"原生实现与当前模块描述不匹配：{module.name}")
                result.add(module.name)
        return frozenset(result)

    def missing(self, program: Program) -> tuple[str, ...]:
        """列出从入口可达、尚未被原生实现覆盖的开放声明模块名。

        Args:
            program: 待检查的程序。

        Returns:
            tuple[str, ...]: 按名字排序的开放声明模块名；无缺口时为空元组。

        Raises:
            ValidationError: 程序未通过 RIR 结构校验。
        """
        validate(program)
        supported = self.matching(program)
        result, visited = [], set()

        def visit(key: str) -> None:
            """跳过已接管的模块，递归收集入口可达的开放声明名。"""
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
    """传给原生工厂的执行环境，衔接 RIR 形参与 PySparQ 运行时对象。

    Attributes:
        ps: 当前已导入的 ``pysparq`` 模块。
        site: 正在执行的原生调用点。
        names: RIR 寄存器名到 PySparQ 寄存器名的映射。
        qrams: 资源名到 PySparQ QRAM 对象的映射。
        memories: 资源名到稀疏字数据字典的映射。
    """

    ps: object
    site: NativeSite
    names: dict
    qrams: dict
    memories: dict

    def register_id(self, parameter: str) -> int:
        """把模块的量子形参解析为 PySparQ 的整数寄存器 id。

        Args:
            parameter: 模块寄存器形参名。

        Returns:
            PySparQ ``System.get_id`` 对应寄存器给出的 id。

        Raises:
            ValidationError: 实参不是恰好覆盖整个寄存器的单一视图（切片或拼接）。
        """
        index = next(i for i, r in enumerate(self.site.module.registers) if r.name == parameter)
        ref = self.site.arguments[index]
        if len(ref.parts) != 1:
            raise ValidationError("当前原生工厂需要完整寄存器参数")
        span = ref.parts[0]
        name = self.names[span.register]
        # ps 为已导入的 pysparq 模块对象，System 属性由动态模块提供。
        if span.start != 0 or span.width != cast(ModuleType, self.ps).System.size_of(name):
            raise ValidationError("当前原生工厂不接收切片；请使用完整寄存器或门级实现")
        return cast(ModuleType, self.ps).System.get_id(name)

    def qram(self, parameter: str) -> object:
        """把模块的资源形参解析为该调用点绑定的 QRAM 对象。

        Args:
            parameter: 模块资源形参名。

        Returns:
            绑定到该形参的 PySparQ QRAM 对象。
        """
        index = next(i for i, r in enumerate(self.site.module.resources) if r.name == parameter)
        return self.qrams[self.site.resources[index]]


class DynamicCppFactory:
    """按需调用真实 compile_operator，C++ 代码不写入 RIR。"""

    def __init__(
        self,
        name: str,
        source: str,
        parameters: Iterable[str],
        *,
        cache_dir: str,
        base_class: str = "SelfAdjointOperator",
    ) -> None:
        """记录算子名、C++ 源与参数，延迟到首次调用再编译。"""
        self.name: str = name
        self.source: str = source
        self.parameters: tuple[str, ...] = tuple(parameters)
        self.cache_dir: str = str(cache_dir)
        self.base_class: str = base_class
        self._class: type | None = None

    def __call__(self, context: NativeContext) -> object:
        """首次调用时延迟编译 C++ 算子源码，再按形参寄存器 id 构造实例。"""
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
