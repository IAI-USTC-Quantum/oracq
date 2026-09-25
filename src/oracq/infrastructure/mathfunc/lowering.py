"Automatic compute/XOR/uncompute lowering from MIR to modular RIR."

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass, replace

from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.input_model.operators import _name
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Program, Ref, RegType, ValidationError
from oracq.infrastructure.mathfunc.graph import MathProgram
from oracq.infrastructure.mathfunc.numeric import MathConfig, Numeric, NumericEmitter


def ports(
    name: str, kind: str, fmt: FixedFormat, index_width: int = 0
) -> list[tuple[str, RegType]]:
    """Expand one port's physical register names and bit widths by math kind.

    complex splits into real/imaginary fixed-point registers, real occupies
    one fixed-point word, bool one bit, and index an unsigned register of
    ``index_width`` bits.

    Args:
        name: The port base name; complex ports append the ``_real``/``_imag`` suffixes to it.
        kind: The port kind: real, complex, bool or index.
        fmt: The ``FixedFormat`` fixed-point format, deciding numeric ports' register widths.
        index_width: The unsigned bit width of the index kind.

    Returns:
        list: A list of ``(register name, Bits)`` pairs.
    """
    if kind == "complex":
        return [(name + "_real", Bits(fmt.width)), (name + "_imag", Bits(fmt.width))]
    return [(name, Bits(1 if kind == "bool" else index_width if kind == "index" else fmt.width))]


@dataclass(frozen=True)
class CompiledFunction:
    """The complete result of one math function compilation.

    Attributes:
        operation: The lowered reversible ``Operation``; inputs stay unchanged, outputs and status are XOR-written back.
        math_ir: The original ``MathProgram`` object that produced the module.
        fmt: The ``FixedFormat`` fixed-point format used for generation.
        output_names: The output port names of the entry's return values.
        input_layout: Mapping from each logical input to its expanded physical register names; each entry is (parameter name, (register name, ...)).
        output_layout: Mapping from each logical output to its expanded physical register names, structured like ``input_layout``.
    """

    operation: Operation
    math_ir: MathProgram
    fmt: FixedFormat
    output_names: tuple[str, ...]
    input_layout: tuple[tuple[str, tuple[str, ...]], ...]
    output_layout: tuple[tuple[str, tuple[str, ...]], ...]

    def program(self) -> Program:
        """Return the full RIR program of the compilation result.

        Returns:
            Program: The program merging the entry module and its dependency modules; semantics match ``Operation.program``.
        """
        return self.operation.program()


class Lowerer:
    """Lowers an MIR math function graph into reversible quantum modules of modular RIR.

    Every function, helpers included, lowers to an independent RIR module;
    calls are kept as Call and are not expanded during lowering. Input
    ports stay unchanged, results and status bits are XOR-written back,
    and private intermediates are finally uncomputed to zero.

    Args:
        program: The ``MathProgram`` to lower; fully validated immediately at construction.
        fmt: The ``FixedFormat`` fixed-point format.
        config: The ``MathConfig`` math kernel degree and approximation interval configuration.
        output_names: Output register names of the entry's return values; when omitted a single return is out and multiple returns are out_0, out_1 and so on.
    """

    def __init__(
        self,
        program: MathProgram,
        fmt: FixedFormat,
        config: MathConfig,
        output_names: Sequence[str] | None = None,
    ) -> None:
        """Validate and bind the MIR program and lowering configuration, and initialize the function lowering cache."""
        self.program: MathProgram = program.validate()
        self.fmt: FixedFormat = fmt
        self.config: MathConfig = config
        self.outputs: Sequence[str] | None = output_names
        self.cache: dict[str, CompiledFunction] = {}

    def lower(self, key: str) -> CompiledFunction:
        """Lower the named function and cache the result; each function is compiled once.

        Args:
            key: Function symbol name in ``program.function_map``, usually the entry.

        Returns:
            CompiledFunction: The function's module, input/output layouts and generation configuration.

        Raises:
            ValidationError: The output port names do not match the number of return values, expanded registers collide, or the fixed-point format cannot hold an index parameter's full range.
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
            raise ValidationError("output_names does not match the number of return values")
        in_layout = [(p.name, ports(p.name, p.kind, self.fmt, p.width)) for p in graph.parameters]
        out_layout = [
            (name, ports(name, graph.nodes[index].kind, self.fmt))
            for name, index in zip(names, graph.returns, strict=True)
        ]
        physical = [item for _, items in in_layout + out_layout for item in items] + [
            ("status", Bits(2))
        ]
        if len({name for name, _ in physical}) != len(physical):
            raise ValidationError("expanded input/output/status registers collide; specify output_names")
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
        parameters: dict[str, Numeric] = {}
        for p, (_, layout) in zip(graph.parameters, in_layout, strict=True):
            refs = tuple(b[name] for name, _ in layout)
            if p.kind == "index":
                if p.width + self.fmt.fraction >= self.fmt.width:
                    raise ValidationError("the fixed-point format cannot hold the full unsigned range of Index")
                word = e.local()
                b.xor(refs[0], word[self.fmt.fraction : self.fmt.fraction + p.width])
                refs = (word,)
            parameters[p.name] = Numeric(
                "real" if p.kind == "index" else p.kind, refs, e.zero_status
            )
        values: dict[int, Numeric] = {}
        called: dict[tuple[str, tuple[int, ...]], list[Numeric]] = {}

        def value(index: int) -> Numeric:
            """Recursively evaluate the named node on demand and cache its circuit representation."""
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
                    actual: dict[str, Ref] = {}
                    for (_, names_), arg in zip(child.input_layout, args, strict=True):
                        if len(names_) != len(arg.parts):
                            raise ValidationError("helper input type expansion mismatch")
                        actual.update(zip(names_, arg.parts, strict=True))
                    returned: list[tuple[str, tuple[Ref, ...]]] = []
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
