"Python 生成阶段构造器。生成结果不包含 Python 可调用对象。"

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass

from pyqecclang.infrastructure.ir import (
    QRAM,
    Adjoint,
    Call,
    Control,
    Load,
    Module,
    Primitive,
    Program,
    Ref,
    Register,
    RegType,
    Repeat,
    Resource,
    Span,
    ValidationError,
)
from pyqecclang.infrastructure.validation import validate


@dataclass(frozen=True)
class Operation:
    module: Module
    dependencies: tuple[Module, ...] = ()

    def unitary(self):
        """RIR 操作的全寄存器酉作用；不表示某个目标矩阵本身是酉矩阵。"""
        return self

    def state_preparation(self):
        """把完整公开寄存器空间作为目标，制备 U|0>。"""
        from pyqecclang.algorithms.oracles import StatePreparation

        return StatePreparation.from_unitary(self)

    def block_encoding(self):
        """完整酉矩阵自身的 alpha=1、零信号 BE；保留原模块调用。"""
        from pyqecclang.algorithms.operators import _name, block_encoding
        from pyqecclang.algorithms.oracles import invoke, resources_for
        from pyqecclang.infrastructure.ir import Bits

        prep = self.state_preparation()
        b = Builder(
            _name("unitary_as_be", self),
            {"target": Bits(prep.width), "signal": Bits(0)},
            resources_for(("unitary", prep.operation)),
        )
        invoke(b, prep.operation, "unitary", target=b["target"], work=b["signal"])
        return block_encoding(b.finish(), alpha=1.0)

    @classmethod
    def from_program(cls, program: Program):
        """由序列化程序恢复入口操作及模块依赖，保留开放声明。"""
        validate(program)
        return cls(program.main, tuple(m for m in program.modules if m.name != program.entry))

    def program(self) -> Program:
        modules = {}
        for item in (*self.dependencies, self.module):
            if item.name in modules and modules[item.name] != item:
                raise ValidationError(f"同名模块定义冲突：{item.name}")
            modules[item.name] = item
        return validate(Program(self.module.name, tuple(modules[k] for k in sorted(modules))))

    def quantikz(self, **kwargs):
        """把入口模块绘制成 quantikz 代码；结构保持 RIR，不展开调用。"""
        from pyqecclang.infrastructure.backends.quantikz import quantikz

        return quantikz(self, **kwargs)

    def estimate(self, **kwargs):
        """Toffoli+Clifford+T+QRAM 级别的组合式资源估计（Repeat 符号相乘）。"""
        from pyqecclang.infrastructure.estimate import estimate_resources

        return estimate_resources(self.program(), **kwargs)


class Builder:
    def __init__(
        self,
        name: str,
        registers: dict[str, RegType],
        resources: dict[str, QRAM] | None = None,
        *,
        attributes: dict | None = None,
    ):
        self.name = name
        self.attributes = tuple(sorted((attributes or {}).items()))
        self.registers = tuple(Register(k, v) for k, v in registers.items())
        self.resources = tuple(Resource(k, v) for k, v in (resources or {}).items())
        self._refs = {r.name: Ref((Span(r.name, 0, r.type.width),), r.type) for r in self.registers}
        self.locals = []
        self._frames = [[]]
        self._dependencies = {}
        self._closed = False

    def local(self, name: str, type: RegType) -> Ref:
        if self._closed or name in self._refs or name in {r.name for r in self.resources}:
            raise ValidationError("局部寄存器名字冲突或构造器已结束")
        self.locals.append(Register(name, type))
        ref = Ref((Span(name, 0, type.width),), type)
        self._refs[name] = ref
        return ref

    def __getitem__(self, name: str) -> Ref:
        return self._refs[name]

    def emit(self, instruction):
        if self._closed:
            raise ValidationError("构造器已经结束")
        self._frames[-1].append(instruction)

    def gate(self, name: str, target: Ref, angle: float | None = None):
        self.emit(Primitive(name, (target,), angle=angle))

    def h(self, target):
        self.gate("h", target)

    def x(self, target):
        self.gate("x", target)

    def z(self, target):
        self.gate("z", target)

    def ry(self, target, angle):
        self.gate("ry", target, angle)

    def rz(self, target, angle):
        self.gate("rz", target, angle)

    def xor(self, source: Ref, target: Ref):
        self.emit(Primitive("xor", (source, target)))

    def swap(self, first: Ref, second: Ref):
        self.emit(Primitive("swap", (first, second)))

    def add_const(self, target: Ref, value: int):
        self.emit(Primitive("add_const", (target,), value=value))

    def global_phase(self, angle: float):
        self.emit(Primitive("gphase", (), angle=angle))

    def qram(self, resource: str, address: Ref, data: Ref):
        self.emit(Load(resource, address, data))

    def call(
        self, operation: Operation, *, resources: dict[str, str] | None = None, **arguments: Ref
    ):
        target = operation.module
        if set(arguments) != {r.name for r in target.registers}:
            raise ValidationError("模块调用的量子参数名字不匹配")
        resources = resources or {}
        if set(resources) != {r.name for r in target.resources}:
            raise ValidationError("模块调用的资源参数名字不匹配")
        for module in (*operation.dependencies, target):
            existing = self._dependencies.get(module.name)
            if existing is not None and existing != module:
                raise ValidationError(f"同名模块定义冲突：{module.name}")
            self._dependencies[module.name] = module
        self.emit(
            Call(
                target.name,
                tuple(arguments[r.name] for r in target.registers),
                tuple(resources[r.name] for r in target.resources),
            )
        )

    @contextmanager
    def _block(self, factory):
        if self._closed:
            raise ValidationError("构造器已经结束")
        self._frames.append([])
        try:
            yield self
        except BaseException:
            self._frames.pop()
            raise
        else:
            nodes = tuple(self._frames.pop())
            self.emit(factory(nodes))

    def repeat(self, count: int):
        return self._block(lambda body: Repeat(count, body))

    def control(self, register: Ref, value: int | None = None):
        if value is None:
            value = (1 << register.width) - 1
        return self._block(lambda body: Control(register, value, body))

    def adjoint(self):
        return self._block(Adjoint)

    def finish(self) -> Operation:
        if self._closed or len(self._frames) != 1:
            raise ValidationError("构造器已结束或仍有未关闭的控制块")
        operation = Operation(
            Module(
                self.name,
                self.registers,
                self.resources,
                tuple(self._frames[0]),
                self.attributes,
                tuple(self.locals),
            ),
            tuple(self._dependencies[k] for k in sorted(self._dependencies)),
        )
        operation.program()
        self._closed = True
        return operation
