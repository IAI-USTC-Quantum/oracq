"MIR → 模块化 RIR 的自动 compute/XOR/uncompute。"

from __future__ import annotations

import json
from dataclasses import dataclass, replace

from pyqecclang.algorithms.common.arithmetic import FixedFormat
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError
from pyqecclang.infrastructure.mathfunc.numeric import Numeric, NumericEmitter


def ports(name, kind, fmt, index_width=0):
    """按数学类型展开一个端口的物理寄存器名与位宽。

    complex 拆为实部/虚部两个定点寄存器，real 占一个定点字，
    bool 占一位，index 占 ``index_width`` 位的无符号寄存器。

    Args:
        name: 端口基名；复数端口在其上追加 ``_real``/``_imag`` 后缀。
        kind: 端口类型，为 real、complex、bool 或 index。
        fmt: ``FixedFormat`` 定点格式，决定数值端口的寄存器位宽。
        index_width: index 类型的无符号位宽。

    Returns:
        list: ``(寄存器名, Bits)`` 对组成的列表。
    """
    if kind == "complex":
        return [(name + "_real", Bits(fmt.width)), (name + "_imag", Bits(fmt.width))]
    return [(name, Bits(1 if kind == "bool" else index_width if kind == "index" else fmt.width))]


@dataclass(frozen=True)
class CompiledFunction:
    """一次数学函数编译的完整结果。

    Attributes:
        operation: 降低得到的可逆 ``Operation``；输入保持不变，输出与 status 以 XOR 写回。
        math_ir: 生成该模块的 ``MathProgram`` 原对象。
        fmt: 生成使用的 ``FixedFormat`` 定点格式。
        output_names: 入口各返回值对应的输出端口名。
        input_layout: 各逻辑输入到其展开物理寄存器名的对应表，每项为 (参数名, (寄存器名, ...))。
        output_layout: 各逻辑输出到其展开物理寄存器名的对应表，结构与 ``input_layout`` 相同。
    """

    operation: object
    math_ir: object
    fmt: FixedFormat
    output_names: tuple[str, ...]
    input_layout: tuple[tuple[str, tuple[str, ...]], ...]
    output_layout: tuple[tuple[str, tuple[str, ...]], ...]

    def program(self):
        """返回编译结果的完整 RIR 程序。

        Returns:
            Program: 入口模块与依赖模块合并后的程序；语义同 ``Operation.program``。
        """
        return self.operation.program()


class Lowerer:
    """把 MIR 数学函数图降低为模块化 RIR 的可逆量子模块。

    每个函数（含 helper）降低为独立 RIR 模块，调用保留为 Call，
    不在降低阶段展开。输入端口保持不变，结果与状态位以 XOR 写回，
    私有中间量最终反算清零。

    Args:
        program: 待降低的 ``MathProgram``；构造时立即做完整校验。
        fmt: ``FixedFormat`` 定点格式。
        config: ``MathConfig`` 数学核阶数与近似区间配置。
        output_names: 入口各返回值的输出寄存器名；省略时单返回值为 out，多返回值依次为 out_0、out_1 等。
    """

    def __init__(self, program, fmt, config, output_names=None):
        self.program = program.validate()
        self.fmt, self.config = fmt, config
        self.outputs = output_names
        self.cache = {}

    def lower(self, key):
        """降低指定函数并缓存结果，同一函数只编译一次。

        Args:
            key: ``program.function_map`` 中的函数符号名，通常为入口。

        Returns:
            CompiledFunction: 该函数的模块、输入/输出布局与生成配置。

        Raises:
            ValidationError: 输出端口名与返回值个数不符、展开后的寄存器重名，或定点格式容不下 index 参数的完整范围。
        """
        if key in self.cache:
            return self.cache[key]
        graph = self.program.function_map[key]
        names = (
            tuple(self.outputs)
            if key == self.program.entry and self.outputs is not None
            else ("out",)
            if len(graph.returns) == 1
            else tuple("out_" + str(i) for i in range(len(graph.returns)))
        )
        if len(names) != len(graph.returns):
            raise ValidationError("output_names 与返回值个数不符")
        in_layout = [(p.name, ports(p.name, p.kind, self.fmt, p.width)) for p in graph.parameters]
        out_layout = [
            (name, ports(name, graph.nodes[index].kind, self.fmt))
            for name, index in zip(names, graph.returns, strict=True)
        ]
        physical = [item for _, items in in_layout + out_layout for item in items] + [
            ("status", Bits(2))
        ]
        if len({name for name, _ in physical}) != len(physical):
            raise ValidationError("展开后的输入/输出/status 寄存器重名，请指定 output_names")
        b = Builder(
            _name("function", key, self.fmt, self.config, names),
            dict(physical),
            attributes={
                "oracle_paradigm": "reversible_function",
                "math_function": graph.label,
                "math_ir_version": "0.1",
                "fixed_width": self.fmt.width,
                "fixed_fraction": self.fmt.fraction,
                "math_config": json.dumps(
                    {"degree": self.config.degree, "intervals": self.config.intervals}
                ),
                "update_semantics": "inputs preserved; outputs XOR; status XOR; locals zero restored",
                "correctness": "pending",
            },
        )
        e = NumericEmitter(b, self.fmt, self.config)
        parameters = {}
        for p, (_, layout) in zip(graph.parameters, in_layout, strict=True):
            refs = tuple(b[name] for name, _ in layout)
            if p.kind == "index":
                if p.width + self.fmt.fraction >= self.fmt.width:
                    raise ValidationError("定点格式容不下 Index 的完整无符号范围")
                word = e.local()
                b.xor(refs[0], word[self.fmt.fraction : self.fmt.fraction + p.width])
                refs = (word,)
            parameters[p.name] = Numeric(
                "real" if p.kind == "index" else p.kind, refs, e.zero_status
            )
        values, called = {}, {}

        def value(index):
            if index in values:
                return values[index]
            node = graph.nodes[index]
            args = [value(i) for i in node.args]
            if node.op == "input":
                result = parameters[node.data[0]]
            elif node.op == "const":
                result = e.constant(
                    complex(*node.data) if node.kind == "complex" else node.data[0], node.kind
                )
            elif node.op in {"neg", "not", "real", "imag", "conj", "abs"}:
                result = e.unary(node.op, args[0])
            elif node.op == "select":
                result = e.choose(*args)
            elif node.op == "complex":
                result = e.complex(*args)
            elif node.op == "intrinsic":
                result = e.intrinsic(node.data[0], args, node.kind)
            elif node.op == "pow":
                exponent = graph.nodes[node.args[1]]
                if (
                    exponent.op == "const"
                    and exponent.kind == "real"
                    and int(exponent.data[0]) == exponent.data[0]
                ):
                    result = e.integer_power(args[0], int(exponent.data[0]))
                else:
                    logarithm = e.intrinsic("log", args[:1], node.kind)
                    result = e.intrinsic("exp", [e.binary("mul", args[1], logarithm)], node.kind)
            elif node.op == "call":
                callee, output = node.data
                callkey = (callee, node.args)
                if callkey not in called:
                    child = self.lower(callee)
                    child_graph = self.program.function_map[callee]
                    actual = {}
                    for (_, names_), arg in zip(child.input_layout, args, strict=True):
                        if len(names_) != len(arg.parts):
                            raise ValidationError("helper 输入类型展开不匹配")
                        actual.update(zip(names_, arg.parts, strict=True))
                    returned = []
                    flag = e.local(2)
                    for (_, names_), ret in zip(
                        child.output_layout, child_graph.returns, strict=True
                    ):
                        kind = child_graph.nodes[ret].kind
                        refs = tuple(e.local(1 if kind == "bool" else None) for _ in names_)
                        actual.update(zip(names_, refs, strict=True))
                        returned.append((kind, refs))
                    e.invoke(child.operation, {**actual, "status": flag})
                    combined = e.flags(*(a.status for a in args), flag)
                    called[callkey] = [Numeric(kind, refs, combined) for kind, refs in returned]
                result = called[callkey][output]
            else:
                result = e.binary(node.op, *args)
            values[index] = result
            return result

        outputs = [
            (tuple(b[name] for name, _ in layout), value(ret))
            for (_, layout), ret in zip(out_layout, graph.returns, strict=True)
        ]
        operation = e.finish(outputs, b["status"])
        if key == self.program.entry:
            attrs = dict(operation.module.attributes)
            attrs["math_ir"] = self.program.dumps()
            operation = replace(
                operation, module=replace(operation.module, attributes=tuple(sorted(attrs.items())))
            )
        result = CompiledFunction(
            operation,
            self.program,
            self.fmt,
            names,
            tuple((name, tuple(key for key, _ in items)) for name, items in in_layout),
            tuple((name, tuple(key for key, _ in items)) for name, items in out_layout),
        )
        self.cache[key] = result
        return result
