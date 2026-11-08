"Automatic compute/XOR/uncompute lowering from MIR to modular RIR."

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, replace

from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.input_model.operators import _name
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Adjoint, Bits, Program, Ref, RegType, ValidationError
from oracq.infrastructure.mathfunc.graph import MathNode, MathProgram
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
                    {"degree": self.config.degree, "intervals": self.config.intervals,
                     "checkpoint_size": self.config.checkpoint_size}
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

        def emit_node(
            node: MathNode, args: list[Numeric], emitter: NumericEmitter,
            parameter_values: dict[str, Numeric],
            call_cache: dict[tuple[str, tuple[int, ...]], list[Numeric]],
        ) -> Numeric:
            """Emit one MIR node with unchanged numeric and path-status semantics."""
            if node.op == "input":
                result = parameter_values[node.data[0]]
            elif node.op == "const":
                result = emitter.constant(
                    complex(*node.data) if node.kind == "complex" else node.data[0], node.kind
                )
            elif node.op in {"neg", "not", "real", "imag", "conj", "abs"}:
                result = emitter.unary(node.op, args[0])
            elif node.op == "select":
                result = emitter.choose(*args)
            elif node.op == "complex":
                result = emitter.complex(*args)
            elif node.op == "intrinsic":
                result = emitter.intrinsic(node.data[0], args, node.kind)
            elif node.op == "pow":
                exponent = graph.nodes[node.args[1]]
                if (
                    exponent.op == "const"
                    and exponent.kind == "real"
                    and int(exponent.data[0]) == exponent.data[0]
                ):
                    result = emitter.integer_power(args[0], int(exponent.data[0]))
                else:
                    logarithm = emitter.intrinsic("log", args[:1], node.kind)
                    result = emitter.intrinsic("exp", [emitter.binary("mul", args[1], logarithm)], node.kind)
            elif node.op == "call":
                callee, output = node.data
                callkey = (callee, node.args)
                if callkey not in call_cache:
                    child = self.lower(callee)
                    child_graph = self.program.function_map[callee]
                    actual: dict[str, Ref] = {}
                    for (_, names_), arg in zip(child.input_layout, args, strict=True):
                        if len(names_) != len(arg.parts):
                            raise ValidationError("helper input type expansion mismatch")
                        actual.update(zip(names_, arg.parts, strict=True))
                    returned: list[tuple[str, tuple[Ref, ...]]] = []
                    flag = emitter.local(2)
                    for (_, names_), ret in zip(
                        child.output_layout, child_graph.returns, strict=True
                    ):
                        kind = child_graph.nodes[ret].kind
                        refs = tuple(emitter.local(1 if kind == "bool" else None) for _ in names_)
                        actual.update(zip(names_, refs, strict=True))
                        returned.append((kind, refs))
                    emitter.invoke(child.operation, {**actual, "status": flag})
                    combined = emitter.flags(*(a.status for a in args), flag)
                    call_cache[callkey] = [Numeric(kind, refs, combined) for kind, refs in returned]
                result = call_cache[callkey][output]
            else:
                result = emitter.binary(node.op, *args)
            return result

        reachable: set[int] = set()
        pending = list(graph.returns)
        while pending:
            current = pending.pop()
            if current not in reachable:
                reachable.add(current)
                pending.extend(graph.nodes[current].args)
        computational = [i for i in sorted(reachable) if graph.nodes[i].op not in {"input", "const"}]
        size = self.config.checkpoint_size
        checkpointed = size > 0 and len(computational) > size
        consumers = Counter(operand for current in reachable for operand in graph.nodes[current].args)
        consumers.update(graph.returns)
        cuts = {current for current in computational if consumers[current] > 1}
        cuts.update(current for current in graph.returns if current in computational)
        costs: dict[int, int] = {}
        for current in computational:
            children = {operand for operand in graph.nodes[current].args if operand in costs and operand not in cuts}
            cost = 1 + sum(costs[operand] for operand in children)
            # Preserve shared values, but inline single-use expressions until
            # the region budget would be exceeded. This avoids artificial
            # checkpoints at every fixed SSA index boundary.
            for operand in sorted(children, key=lambda item: costs[item], reverse=True):
                if cost <= max(1, size):
                    break
                cuts.add(operand)
                cost -= costs[operand]
            costs[current] = cost

        def region(index: int) -> Numeric:
            """Compute one region, copying only its result and path status out.

            Boundary status words are explicit read-only inputs. In particular,
            they are not indiscriminately ORed at the call site: a select must
            still mask errors in its unselected branch.
            """
            boundary: set[int] = set()
            visited: set[int] = set()

            def visit(current: int) -> None:
                if current in visited:
                    return
                visited.add(current)
                node = graph.nodes[current]
                if node.op == "input" or (current != index and current in cuts):
                    boundary.add(current)
                    return
                for operand in node.args:
                    visit(operand)

            outputs_in_region = [index]
            for output in outputs_in_region:
                visit(output)
            arguments = {current: value(current) for current in sorted(boundary)}
            register_types: dict[str, RegType] = {}
            for current, argument in arguments.items():
                for part, ref in enumerate(argument.parts):
                    register_types[f"arg_{current}_{part}"] = ref.type
                register_types[f"flag_{current}"] = Bits(2)
            out_layouts = {output: ports(f"out_{output}", graph.nodes[output].kind, self.fmt)
                           for output in outputs_in_region}
            for output, out_ports in out_layouts.items():
                register_types.update(out_ports)
                register_types[f"status_{output}"] = Bits(2)
            rb = Builder(
                _name("math_region", key, tuple(outputs_in_region), tuple(arguments), self.fmt, self.config),
                register_types,
                attributes={"workspace_contract": "zero_in_zero_out",
                            "math_checkpoint": index, "math_checkpoint_size": size},
            )
            re = NumericEmitter(rb, self.fmt, self.config)
            region_values = {
                current: Numeric(argument.kind,
                    tuple(rb[f"arg_{current}_{part}"] for part in range(len(argument.parts))),
                    rb[f"flag_{current}"])
                for current, argument in arguments.items()
            }
            region_calls: dict[tuple[str, tuple[int, ...]], list[Numeric]] = {}

            def local_value(current: int) -> Numeric:
                if current not in region_values:
                    node = graph.nodes[current]
                    region_values[current] = emit_node(
                        node, [local_value(operand) for operand in node.args], re, {}, region_calls,
                    )
                return region_values[current]

            results = {output: local_value(output) for output in outputs_in_region}
            forward = tuple(rb._frames[0])
            for output, result in results.items():
                for (name, _), ref in zip(out_layouts[output], result.parts, strict=True):
                    rb.xor(ref, rb[name])
                rb.xor(result.status, rb[f"status_{output}"])
            rb.emit(Adjoint(forward))
            operation = rb.finish()
            actual = {f"arg_{current}_{part}": ref
                      for current, argument in arguments.items()
                      for part, ref in enumerate(argument.parts)}
            actual.update({f"flag_{current}": argument.status for current, argument in arguments.items()})
            saved: dict[int, Numeric] = {}
            for output, out_ports in out_layouts.items():
                outputs = tuple(e.local(type_.width) for _, type_ in out_ports)
                flag = e.local(2)
                actual.update({name: ref for (name, _), ref in zip(out_ports, outputs, strict=True)})
                actual[f"status_{output}"] = flag
                saved[output] = Numeric(graph.nodes[output].kind, outputs, flag)
            e.invoke(operation, actual)
            values.update(saved)
            return values[index]

        def value(index: int) -> Numeric:
            """Cache region boundaries, or retain the reference SSA computation."""
            if index in values:
                return values[index]
            node = graph.nodes[index]
            if checkpointed and node.op not in {"input", "const"}:
                result = region(index)
            else:
                result = emit_node(node, [value(i) for i in node.args], e, parameters, called)
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
