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
    Store,
    ValidationError,
)
from pyqecclang.infrastructure.validation import validate


@dataclass(frozen=True)
class Operation:
    """入口模块连同其依赖模块构成的完整 RIR 操作。

    Attributes:
        module: 入口模块定义。
        dependencies: 入口之外的全部被依赖模块定义。
    """

    module: Module
    dependencies: tuple[Module, ...] = ()

    def unitary(self):
        """RIR 操作的全寄存器酉作用；不表示某个目标矩阵本身是酉矩阵。"""
        return self

    def state_preparation(self):
        """把完整公开寄存器空间作为目标，制备 U|0>。"""
        from pyqecclang.algorithms.input_model.oracles import StatePreparation

        return StatePreparation.from_unitary(self)

    def block_encoding(self):
        """完整酉矩阵自身的 alpha=1、零信号 BE；保留原模块调用。"""
        from pyqecclang.algorithms.input_model.operators import _name, block_encoding
        from pyqecclang.algorithms.input_model.oracles import invoke, resources_for
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
        """把入口模块与依赖合并成模块按名字排序并经校验的 ``Program``。

        Returns:
            Program: 以 ``module.name`` 为入口的程序。

        Raises:
            ValidationError: 同名依赖模块存在冲突定义，或合并结果违反 RIR 规则。
        """
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
    """逐指令构造 RIR 模块的生成阶段构造器；生成结果不含 Python 可调用对象。

    Args:
        name: 生成的模块名。
        registers: 公开寄存器名到位宽类型的映射。
        resources: 公开 QRAM 资源名到规格的映射；省略时模块无资源。
        attributes: 追加到模块的属性键值对；省略时为空。

    公开寄存器经 ``builder[name]`` 取 ``Ref`` 视图，``local`` 声明模块私有
    工作寄存器，``finish`` 关闭构造器并返回经校验的 ``Operation``。
    """

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
        """声明一个模块私有工作寄存器并返回覆盖它的完整视图。

        Args:
            name: 寄存器名，不得与公开寄存器、资源或其他局部寄存器重名。
            type: 寄存器的类型与位宽。

        Returns:
            Ref: 覆盖整个新寄存器的引用。

        Raises:
            ValidationError: 名字冲突，或构造器已经结束。
        """
        if self._closed or name in self._refs or name in {r.name for r in self.resources}:
            raise ValidationError("局部寄存器名字冲突或构造器已结束")
        self.locals.append(Register(name, type))
        ref = Ref((Span(name, 0, type.width),), type)
        self._refs[name] = ref
        return ref

    def __getitem__(self, name: str) -> Ref:
        return self._refs[name]

    def emit(self, instruction):
        """把一条 RIR 指令追加到当前最内层作用域。

        Args:
            instruction: RIR ``Instruction`` 节点。

        Raises:
            ValidationError: 构造器已经结束。
        """
        if self._closed:
            raise ValidationError("构造器已经结束")
        self._frames[-1].append(instruction)

    def gate(self, name: str, target: Ref, angle: float | None = None):
        """发射单目标 ``Primitive`` 门。

        Args:
            name: 门名。
            target: 目标视图。
            angle: 旋转类门的弧度角；非旋转门省略。
        """
        self.emit(Primitive(name, (target,), angle=angle))

    def h(self, target):
        """发射 Hadamard 门。"""
        self.gate("h", target)

    def x(self, target):
        """发射 Pauli-X 门。"""
        self.gate("x", target)

    def z(self, target):
        """发射 Pauli-Z 门。"""
        self.gate("z", target)

    def ry(self, target, angle):
        """发射绕 Y 轴旋转 ``angle`` 弧度的门。"""
        self.gate("ry", target, angle)

    def rz(self, target, angle):
        """发射绕 Z 轴旋转 ``angle`` 弧度的门。"""
        self.gate("rz", target, angle)

    def xor(self, source: Ref, target: Ref):
        """发射按位异或门，把 ``source`` 的值异或进等宽的 ``target``。"""
        self.emit(Primitive("xor", (source, target)))

    def swap(self, first: Ref, second: Ref):
        """发射 ``first`` 与 ``second`` 两个等宽视图的交换门。"""
        self.emit(Primitive("swap", (first, second)))

    def add_const(self, target: Ref, value: int):
        """发射常量加法门，把 ``target`` 加上 ``value`` 并按位宽回绕。"""
        self.emit(Primitive("add_const", (target,), value=value))

    def global_phase(self, angle: float):
        """发射全局相位门，整体幅值乘以 ``exp(i*angle)``。"""
        self.emit(Primitive("gphase", (), angle=angle))

    def qram(self, resource: str, address: Ref, data: Ref):
        """发射 QRAM 读取：按 ``address`` 查询 ``resource``，把命中的字异或进 ``data``。"""
        self.emit(Load(resource, address, data))

    def store(self, resource: str, address: Ref, data: Ref):
        """发射 QRAM 写入：把 ``data`` 存入 ``resource`` 的 ``address`` 单元。"""
        self.emit(Store(resource, address, data))

    def call(
        self, operation: Operation, *, resources: dict[str, str] | None = None, **arguments: Ref
    ):
        """原样调用另一个 ``Operation``；被调模块保留为 ``Call`` 节点，不在生成期展开。

        Args:
            operation: 被调用的操作；其模块与依赖一并登记为本模块的依赖。
            resources: 被调模块资源形参名到本模块资源名的映射。
            **arguments: 以被调模块寄存器名为关键字的实参视图。

        Raises:
            ValidationError: 量子或资源参数名不匹配，或同名依赖模块定义冲突。
        """
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
        """打开重复 ``count`` 次的作用域；体内指令聚合为 ``Repeat`` 节点，计数保持符号，不在生成期展开。

        Returns:
            ``with`` 语句使用的作用域上下文；块须按 ``with`` 层级全部关闭后才能 ``finish``。
        """
        return self._block(lambda body: Repeat(count, body))

    def control(self, register: Ref, value: int | None = None):
        """打开受控作用域；体内指令聚合为 ``Control`` 节点。

        Args:
            register: 控制视图。
            value: 生效的整数掩码；省略时取寄存器全一掩码。

        Returns:
            ``with`` 语句使用的作用域上下文；块须按 ``with`` 层级全部关闭后才能 ``finish``。
        """
        if value is None:
            value = (1 << register.width) - 1
        return self._block(lambda body: Control(register, value, body))

    def adjoint(self):
        """打开逆作用域；体内指令聚合为 ``Adjoint`` 块，执行时整体逆序取逆。

        Returns:
            ``with`` 语句使用的作用域上下文；块须按 ``with`` 层级全部关闭后才能 ``finish``。
        """
        return self._block(Adjoint)

    def finish(self) -> Operation:
        """结束构造：组装模块与全部登记依赖，经校验后返回 ``Operation`` 并关闭构造器。

        Returns:
            Operation: 入口为本次构造的模块，依赖为登记过的被调模块。

        Raises:
            ValidationError: 构造器已结束、仍有未关闭的作用域，或组装结果违反 RIR 规则。
        """
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
