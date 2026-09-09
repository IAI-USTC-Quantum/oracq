"""MIR → 模块化 RIR 的自动 compute/XOR/uncompute。"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace

from ..arithmetic import FixedFormat
from ..builder import Builder
from ..ir import Bits, ValidationError
from ..library import _name
from .numeric import Numeric, NumericEmitter


def ports(name, kind, fmt, index_width=0):
    if kind == "complex":
        return [(name + "_real", Bits(fmt.width)), (name + "_imag", Bits(fmt.width))]
    return [(name, Bits(1 if kind == "bool" else index_width if kind == "index" else fmt.width))]


@dataclass(frozen=True)
class CompiledFunction:
    operation: object
    math_ir: object
    fmt: FixedFormat
    output_names: tuple[str, ...]
    input_layout: tuple[tuple[str, tuple[str, ...]], ...]
    output_layout: tuple[tuple[str, tuple[str, ...]], ...]

    def program(self):
        return self.operation.program()


class Lowerer:
    def __init__(self, program, fmt, config, output_names=None):
        self.program = program.validate()
        self.fmt, self.config = fmt, config
        self.outputs = output_names
        self.cache = {}

    def lower(self, key):
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
