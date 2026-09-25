"Pure-function compilation frontend that reads the Python AST; it never executes the function being compiled."

from __future__ import annotations

import ast
import builtins
import cmath
import hashlib
import inspect
import math
import textwrap
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import NoReturn, cast

from oracq.infrastructure.ir import ValidationError
from oracq.infrastructure.mathfunc.graph import (
    Index,
    MathFunction,
    MathNode,
    MathProgram,
    Parameter,
)

ELEMENTARY = {
    "sqrt",
    "exp",
    "log",
    "log10",
    "sin",
    "cos",
    "tan",
    "asin",
    "acos",
    "atan",
    "sinh",
    "cosh",
    "tanh",
    "asinh",
    "acosh",
    "atanh",
    "phase",
    "polar",
    "rect",
}
"""Set of elementary math function names compilable into intrinsic nodes; phase, polar and rect are cmath only."""
REAL_EXTRA = {"atan2", "hypot"}
"""Set of additional binary math function names supported on the real domain only."""
BUILTINS = {"abs", "complex", "min", "max"}
"""Set of Python builtin function names allowed to be called directly by name."""
KINDS = {"float": "real", "complex": "complex", "bool": "bool", "int": "real"}
"""Mapping from Python annotation or type names to MIR value kinds; unannotated parameters are treated as float."""


class FunctionCompileError(ValidationError):
    """Exception raised when the math function frontend fails to parse or compile restricted Python source."""


@dataclass(frozen=True)
class Value:
    """Compile-time SSA value: a pair of node index and value kind.

    Attributes:
        id: Position of the node backing this value within the current function's node sequence.
        kind: Value kind: real, complex or bool.
    """

    id: int
    kind: str


@dataclass
class Source:
    """Location record of a function to compile: its AST definition and evaluation namespace.

    Attributes:
        tree: The function's ``ast.FunctionDef`` definition node.
        namespace: Mapping of global and closure variables used to resolve name references.
    """

    tree: ast.FunctionDef
    namespace: dict[str, object]


class Frontend:
    """Compiles restricted pure Python function source into an MIR math function graph.

    Reads and interprets the AST of an ordinary Python function without
    executing the function being compiled. In source-string form, every def
    other than the entry is kept as a pure helper in its own MIR function;
    calls are represented as call nodes and are not expanded at compile time.

    Args:
        source: An ordinary Python function object, or a source string containing several defs and math/cmath imports; the module docstring and ``__future__`` imports are ignored.
        inputs: Mapping from parameter names to real/complex/bool or ``Index(width)``; inferred from annotations when omitted.
        constants: Mapping from parameter names to finite numeric values, bound as generation-time constants.
        helpers: Helper functions appended to the entry namespace.
        max_unroll: Maximum number of iterations a static range loop may unroll.
        entry: Entry function name for the source-string form; the last def is used when omitted.

    Raises:
        FunctionCompileError: The source has a syntax error, contains unsupported statements, or has no usable function definition.
    """

    def __init__(
        self,
        source: str | Callable[..., object] | Source,
        *,
        inputs: Mapping[str, Index | type | str] | None = None,
        constants: Mapping[str, bool | int | float | complex] | None = None,
        helpers: Mapping[str, Callable[..., object]] | None = None,
        max_unroll: int = 128,
        entry: str | None = None,
    ) -> None:
        """Parse the entry source and initialize the function table, compile cache and active stack."""
        self.functions: dict[str, MathFunction]
        self.cache: dict[tuple, str]
        self.active: set[tuple]
        self.functions, self.cache, self.active = {}, {}, set()
        self.max_unroll: int = max_unroll
        self.sources: dict[str, Source] = {}
        if isinstance(source, str):
            try:
                tree = ast.parse(textwrap.dedent(source))
            except SyntaxError as exc:
                raise FunctionCompileError(
                    f"Python syntax error at line {exc.lineno}: {exc.msg}"
                ) from exc
            ns: dict[str, object] = {"math": math, "cmath": cmath}
            for stmt in tree.body:
                if (
                    isinstance(stmt, ast.Expr)
                    and isinstance(stmt.value, ast.Constant)
                    and isinstance(stmt.value.value, str)
                ):
                    continue
                if isinstance(stmt, ast.ImportFrom) and stmt.module == "__future__":
                    # __future__ only affects real Python compilation and has no semantics while the frontend interprets the AST.
                    continue
                if isinstance(stmt, (ast.Import, ast.ImportFrom)):
                    self.imports(stmt, ns)
                elif isinstance(stmt, ast.FunctionDef):
                    self.sources[stmt.name] = Source(stmt, ns)
                else:
                    raise FunctionCompileError("the source may only contain math/cmath imports and pure function definitions")
            if not self.sources:
                raise FunctionCompileError("the source contains no functions")
            for key, value in self.sources.items():
                ns[key] = value
            if entry is not None and entry not in self.sources:
                raise FunctionCompileError("the source contains no entry function: " + entry)
            entry = self.sources[entry] if entry is not None else list(self.sources.values())[-1]  # type: ignore[assignment]
        else:
            entry = self.source(source)  # type: ignore[assignment]
        if helpers:
            cast(Source, entry).namespace.update(helpers)
        self.entry: str = self.compile(cast(Source, entry), inputs, constants or {})

    def source(self, function: object) -> Source:
        """Normalize a function object or Source into a location record.

        Reads the source of ``function``, locates the def of the same name
        uniquely, and builds the evaluation namespace from its globals and
        closure variables.

        Args:
            function: An ordinary Python function object, or an already constructed Source.

        Returns:
            Source: The location record, ready to be handed to compile.

        Raises:
            FunctionCompileError: The input is not an ordinary function with readable source, the source is unavailable, or the definition cannot be located uniquely.
        """
        if isinstance(function, Source):
            return function
        if not inspect.isfunction(function):
            raise FunctionCompileError(
                "only plain Python functions with readable source are supported; a def source string may also be passed"
            )
        try:
            tree = ast.parse(textwrap.dedent(inspect.getsource(function)))
        except (OSError, TypeError, IndentationError) as exc:
            raise FunctionCompileError("cannot read the function source; pass a def source string instead") from exc
        definitions = [
            n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == function.__name__
        ]
        if len(definitions) != 1:
            raise FunctionCompileError("cannot uniquely locate the function definition; lambda is not supported")
        closure = inspect.getclosurevars(function)
        return Source(definitions[0], {**function.__globals__, **closure.nonlocals})

    @staticmethod
    def imports(stmt: ast.Import | ast.ImportFrom, namespace: dict[str, object]) -> None:
        """Handle one math/cmath import statement and write it into the namespace.

        Args:
            stmt: An ``ast.Import`` or ``ast.ImportFrom`` node.
            namespace: The mapping that receives the imported names.

        Raises:
            FunctionCompileError: The import target is not math/cmath, or an unsupported math name was imported.
        """
        if isinstance(stmt, ast.Import):
            for alias in stmt.names:
                if alias.name not in {"math", "cmath"}:
                    raise FunctionCompileError("only math/cmath may be imported")
                namespace[alias.asname or alias.name] = {"math": math, "cmath": cmath}[alias.name]
        else:
            if stmt.module not in {"math", "cmath"} or stmt.level:
                raise FunctionCompileError("only imports from math/cmath are allowed")
            module = {"math": math, "cmath": cmath}[stmt.module]
            for alias in stmt.names:
                if alias.name not in ELEMENTARY | REAL_EXTRA | {"pi", "e", "tau"}:
                    raise FunctionCompileError("unsupported math import: " + alias.name)
                namespace[alias.asname or alias.name] = getattr(module, alias.name)

    def fail(self, node: ast.AST, message: str) -> NoReturn:
        """Raise an error carrying the source location.

        Args:
            node: The AST node marking the error location; a question mark replaces the line number when no location info exists.
            message: Description of the error.

        Returns:
            NoReturn: Never returns; always ends by raising ``FunctionCompileError``.

        Raises:
            FunctionCompileError: Always raised, with the message prefixed by the current function label and line number.
        """
        raise FunctionCompileError(f"{self.label}:{getattr(node, 'lineno', '?')}: {message}")

    def compile(
        self,
        source: Callable[..., object] | Source,
        inputs: Mapping[str, Index | type | str] | None = None,
        constants: Mapping[str, bool | int | float | complex] | None = None,
    ) -> str:
        """Compile one function definition to MIR and return its symbol name.

        Infers dynamic input kinds, interprets the function body as an SSA
        node sequence, and caches repeated compilations by AST, inputs,
        constants and captured values; the symbol name is derived from a
        hash of the function graph content.

        Args:
            source: A Python function object, Source or source string.
            inputs: Mapping from parameter names to kinds; when None, inferred from annotations and defaulting to real when unannotated.
            constants: Mapping from parameter names to finite numeric values, bound as generation-time constants.

        Returns:
            str: The function symbol name registered in ``self.functions``.

        Raises:
            FunctionCompileError: A parameter declaration, statement or expression leaves the restricted subset, or a recursive helper is detected.
        """
        source = self.source(source)
        constants = constants or {}
        tree = source.tree
        if tree.args.vararg or tree.args.kwarg or tree.args.kwonlyargs or tree.args.posonlyargs:
            raise FunctionCompileError("functions require plain named parameters; varargs, kwargs and keyword-only parameters are not supported")
        parameters = tree.args.args
        defaults = (
            dict(
                zip(
                    [p.arg for p in parameters][-len(tree.args.defaults) :],
                    tree.args.defaults,
                    strict=True,
                )
            )
            if tree.args.defaults
            else {}
        )
        if inputs is None:
            inputs = {
                p.arg: KINDS.get(getattr(p.annotation, "id", "float"), "real")
                for p in parameters
                if p.arg not in constants and p.arg not in defaults
            }
        if set(inputs) & set(constants) & {p.arg for p in parameters}:
            raise FunctionCompileError("a parameter cannot be both a dynamic input and a generation-time constant")
        if set(inputs) - {p.arg for p in parameters}:
            raise FunctionCompileError("inputs contains unknown parameters")
        dynamic: list[Parameter] = []
        for p in parameters:
            if p.arg in inputs:
                spec: Index | type | str | Parameter = inputs[p.arg]
                if isinstance(spec, Index):
                    dynamic.append(Parameter(p.arg, "index", spec.width))
                else:
                    kind = KINDS.get(getattr(spec, "__name__", spec), spec)  # type: ignore[arg-type]
                    if kind not in {"real", "complex", "bool"}:
                        raise FunctionCompileError("input types must be real/complex/bool or an Index instance")
                    dynamic.append(Parameter(p.arg, kind))
        captured = tuple(
            sorted(
                (k, v)
                for k, v in source.namespace.items()
                if type(v) in (bool, int, float, complex)
            )
        )
        key = (
            ast.dump(tree, include_attributes=False),
            tuple(dynamic),
            repr(sorted(constants.items())),
            captured,
        )
        if key in self.cache:
            return self.cache[key]
        marker = (tree.name, ast.dump(tree, include_attributes=False))
        if marker in self.active:
            raise FunctionCompileError("recursive helpers are not supported: " + tree.name)
        self.active.add(marker)
        saved = {k: getattr(self, k, None) for k in ("nodes", "intern", "namespace", "label")}
        self.nodes: list[MathNode]
        self.intern: dict[MathNode, int]
        self.namespace: dict[str, object]
        self.label: str
        self.nodes, self.intern = [], {}
        self.namespace = {**source.namespace, **constants}
        self.label = tree.name
        env: dict[str, Value | tuple[Value, ...]] = {}
        try:
            for p in parameters:
                if p.arg in inputs:
                    spec = next(x for x in dynamic if x.name == p.arg)
                    env[p.arg] = self.node(
                        "input", "real" if spec.kind == "index" else spec.kind, data=(p.arg,)
                    )
                elif p.arg in constants:
                    env[p.arg] = self.constant(constants[p.arg], p)
                elif p.arg in defaults:
                    env[p.arg] = self.expr(defaults[p.arg], env)
                else:
                    self.fail(p, "parameter is not declared as a dynamic input or generation-time constant: " + p.arg)
            _, result = self.statements(tree.body, env)
            if result is None:
                self.fail(tree, "all control paths must return a value")
            values = result if isinstance(result, tuple) else (result,)
            if any(not isinstance(v, Value) for v in values):
                self.fail(tree, "return values must be scalars or a flat tuple")
            digest = hashlib.sha256(repr((key, self.nodes, values)).encode()).hexdigest()[:20]
            symbol = "math_" + tree.name + "_" + digest
            function = MathFunction(
                symbol, tree.name, tuple(dynamic), tuple(self.nodes), tuple(v.id for v in values)
            )
            self.functions[symbol] = function
            self.cache[key] = symbol
        finally:
            for k, v in saved.items():
                setattr(self, k, v)
            self.active.remove(marker)
        return symbol

    def node(
        self,
        op: str,
        kind: str,
        args: Sequence[Value] = (),
        data: Sequence[bool | int | float | str] = (),
    ) -> Value:
        """Construct or reuse one SSA node and return its value.

        Only one node with identical op, kind, inputs and data is kept per
        current function.

        Args:
            op: Node operation name, such as add, select or intrinsic.
            kind: Result value kind: real/complex/bool.
            args: List of operand Values.
            data: Extra data attached to the node, such as a constant value or an intrinsic name.

        Returns:
            Value: The compile-time value pointing to the node.
        """
        values = tuple(v.id for v in args)
        key = MathNode(op, kind, values, tuple(data))
        if key not in self.intern:
            self.intern[key] = len(self.nodes)
            self.nodes.append(key)
        return Value(self.intern[key], kind)

    def constant(self, value: bool | int | float | complex, node: ast.AST) -> Value:
        """Freeze one finite numeric value into a const node.

        Args:
            value: A bool, int, float or complex constant.
            node: The AST node used for error location.

        Returns:
            Value: The constant value with a kind matching the constant.

        Raises:
            FunctionCompileError: The numeric type is unsupported, or the value is non-finite such as NaN/Inf.
        """
        if type(value) not in (bool, int, float, complex):
            self.fail(node, "only finite numeric constants are allowed")
        if not (math.isfinite(value.real) and math.isfinite(value.imag)):
            self.fail(node, "fixed-point cannot represent NaN/Inf")
        kind = "bool" if type(value) is bool else "complex" if type(value) is complex else "real"
        return self.node(
            "const",
            kind,
            data=(value.real, value.imag)
            if kind == "complex"
            else (cast("bool | int | float", value),),
        )

    def literal(self, value: Value | tuple[Value, ...], node: ast.AST) -> bool | int | float | complex:
        """Extract the Python numeric value carried by a compile-time constant node.

        Args:
            value: A Value expected to be a const node.
            node: The AST node used for error location.

        Returns:
            The raw Python numeric value of the constant; complex values are represented as complex.

        Raises:
            FunctionCompileError: The value is not a constant node.
        """
        if isinstance(value, Value) and self.nodes[value.id].op == "const":
            data = self.nodes[value.id].data
            return complex(*data) if value.kind == "complex" else data[0]
        self.fail(node, "a generation-time numeric constant is required here")

    def binary(
        self, op: str, a: Value | tuple[Value, ...], b: Value | tuple[Value, ...], node: ast.AST
    ) -> Value:
        """Construct a binary operation node and infer the result kind.

        Arithmetic between real/complex is typed by the promotion rules,
        comparisons and boolean operations produce bool; when both sides are
        constants the operation is evaluated at generation time and folded
        into a single constant node.

        Args:
            op: Operation name, such as add, sub, mul, div, pow, lt, eq, and, or.
            a: Left operand Value.
            b: Right operand Value.
            node: The AST node used for error location.

        Returns:
            Value: The operation result value.

        Raises:
            FunctionCompileError: An operand is not a scalar, the type combination is illegal, or constant folding hits an undefined operation.
        """
        if not isinstance(a, Value) or not isinstance(b, Value):
            self.fail(node, "arithmetic only accepts scalars")
        if op in {"and", "or"}:
            if a.kind != "bool" or b.kind != "bool":
                self.fail(node, "invalid boolean operand types")
            kind = "bool"
        elif op in {"lt", "eq"}:
            if op == "lt" and "complex" in (a.kind, b.kind):
                self.fail(node, "ordering comparison is not supported for complex values")
            kind = "bool"
        else:
            if "bool" in (a.kind, b.kind):
                self.fail(node, "boolean values require an explicit select to become numeric")
            kind = "complex" if "complex" in (a.kind, b.kind) else "real"
        if self.nodes[a.id].op == "const" and self.nodes[b.id].op == "const":
            av, bv = self.literal(a, node), self.literal(b, node)
            functions: dict[str, Callable[[], bool | int | float | complex]] = {
                "add": lambda: av + bv,
                "sub": lambda: av - bv,
                "mul": lambda: av * bv,
                "div": lambda: av / bv,
                "pow": lambda: av**bv,
                "lt": lambda: cast("bool | int | float", av) < cast("bool | int | float", bv),
                "eq": lambda: av == bv,
                "and": lambda: av and bv,
                "or": lambda: av or bv,
            }
            try:
                return self.constant(functions[op](), node)
            except (ArithmeticError, ValueError) as exc:
                self.fail(node, "constant expression is undefined: " + str(exc))
        return self.node(op, kind, (a, b))

    def choose(
        self,
        test: Value | tuple[Value, ...],
        yes: Value | tuple[Value, ...],
        no: Value | tuple[Value, ...],
        node: ast.AST,
    ) -> Value | tuple[Value, ...]:
        """Construct a select node choosing branch values by a boolean condition.

        Tuple branches are selected element-wise; real and complex branches
        are unified by the promotion rules.

        Args:
            test: The condition value of bool kind.
            yes: The value, or tuple of values, chosen when the condition is true.
            no: The value, or tuple of values, chosen when the condition is false.
            node: The AST node used for error location.

        Returns:
            Value or tuple: The selection result, shaped like the branches.

        Raises:
            FunctionCompileError: The condition is not boolean, or the two branches disagree in structure or type.
        """
        if not isinstance(test, Value) or test.kind != "bool":
            self.fail(node, "the quantum condition must be a boolean expression")
        if isinstance(yes, tuple) and isinstance(no, tuple) and len(yes) == len(no):
            return tuple(cast(Value, self.choose(test, a, b, node)) for a, b in zip(yes, no, strict=True))
        if not isinstance(yes, Value) or not isinstance(no, Value):
            self.fail(node, "branch return or assignment structure is inconsistent")
        if {yes.kind, no.kind} <= {"real", "complex"}:
            kind = "complex" if "complex" in (yes.kind, no.kind) else "real"
        elif yes.kind == no.kind:
            kind = yes.kind
        else:
            self.fail(node, "branch types are inconsistent")
        return self.node("select", kind, (test, yes, no))

    def expr(
        self, node: ast.expr, env: dict[str, Value | tuple[Value, ...]]
    ) -> Value | tuple[Value, ...]:
        """Interpret one expression AST as a compile-time value.

        Supports constants, names, tuple/list, generation-time tuple
        subscripts, math/cmath constant attributes, arithmetic, comparisons,
        boolean short-circuiting, conditional expressions and function calls.

        Args:
            node: The expression AST node.
            env: Mapping from local names to Values or tuples of them.

        Returns:
            Value or tuple: The value of the expression; tuple/list expressions produce a tuple of Values.

        Raises:
            FunctionCompileError: The expression leaves the restricted subset.
        """
        if isinstance(node, ast.Constant):
            return self.constant(node.value, node)  # type: ignore[arg-type]
        if isinstance(node, ast.Name):
            if node.id in env:
                return env[node.id]
            value = self.namespace.get(node.id)
            return self.constant(value, node)  # type: ignore[arg-type]
        if isinstance(node, (ast.Tuple, ast.List)):
            return tuple(self.expr(x, env) for x in node.elts)  # type: ignore[misc]
        if isinstance(node, ast.Subscript):
            value = self.expr(node.value, env)
            index = self.literal(self.expr(node.slice, env), node)
            if not isinstance(value, tuple) or type(index) is not int:
                self.fail(node, "only generation-time tuple subscripts are supported")
            try:
                return value[index]
            except IndexError:
                self.fail(node, "tuple index out of range")
        if isinstance(node, ast.Attribute):
            if (
                isinstance(node.value, ast.Name)
                and (
                    self.namespace.get(node.value.id) is math
                    or self.namespace.get(node.value.id) is cmath
                )
                and node.attr in {"pi", "e", "tau"}
            ):
                return self.constant(getattr(math, node.attr), node)
            value = self.expr(node.value, env)
            if not isinstance(value, Value) or node.attr not in {"real", "imag"}:
                self.fail(node, "only the real and imag numeric attributes are supported")
            return self.node(node.attr, "real", (value,))
        if isinstance(node, ast.BinOp):
            operations = {
                ast.Add: "add",
                ast.Sub: "sub",
                ast.Mult: "mul",
                ast.Div: "div",
                ast.Pow: "pow",
            }
            if type(node.op) not in operations:
                self.fail(node, "unsupported binary operation")
            return self.binary(
                operations[type(node.op)],
                self.expr(node.left, env),
                self.expr(node.right, env),
                node,
            )
        if isinstance(node, ast.UnaryOp):
            value = self.expr(node.operand, env)
            if isinstance(node.op, ast.UAdd):
                return value
            if isinstance(node.op, ast.USub):
                if self.nodes[cast(Value, value).id].op == "const":
                    return self.constant(-self.literal(value, node), node)
                return self.node("neg", cast(Value, value).kind, (cast(Value, value),))
            if isinstance(node.op, ast.Not) and cast(Value, value).kind == "bool":
                return self.node("not", "bool", (cast(Value, value),))
            self.fail(node, "unsupported unary operation")
        if isinstance(node, ast.Compare):
            pairs: list[Value] = []
            left = self.expr(node.left, env)
            for operator, right_node in zip(node.ops, node.comparators, strict=True):
                right = self.expr(right_node, env)
                if isinstance(operator, (ast.Eq, ast.NotEq)):
                    test = self.binary("eq", left, right, node)
                elif isinstance(operator, (ast.Lt, ast.GtE)):
                    test = self.binary("lt", left, right, node)
                elif isinstance(operator, (ast.Gt, ast.LtE)):
                    test = self.binary("lt", right, left, node)
                else:
                    self.fail(node, "unsupported comparison")
                if isinstance(operator, (ast.NotEq, ast.GtE, ast.LtE)):
                    test = self.node("not", "bool", (test,))
                pairs.append(test)
                left = right
            result: Value | tuple[Value, ...] = pairs[0]
            for test in pairs[1:]:
                result = self.node("and", "bool", (cast(Value, result), test))
            return result
        if isinstance(node, ast.BoolOp):
            values = [self.expr(v, env) for v in node.values]
            result = values[0]
            for value in values[1:]:
                # Python short-circuit semantics also apply to the status path.
                result = (
                    self.choose(result, value, self.constant(False, node), node)
                    if isinstance(node.op, ast.And)
                    else self.choose(result, self.constant(True, node), value, node)
                )
            return result
        if isinstance(node, ast.IfExp):
            return self.choose(
                self.expr(node.test, env),
                self.expr(node.body, env),
                self.expr(node.orelse, env),
                node,
            )
        if isinstance(node, ast.Call):
            return self.call(node, env)
        self.fail(node, "unsupported expression: " + type(node).__name__)

    def call(
        self, node: ast.Call, env: dict[str, Value | tuple[Value, ...]]
    ) -> Value | tuple[Value, ...]:
        """Interpret a function call expression.

        Handles the ``conjugate`` method, named pure helper calls, the
        whitelisted builtins complex/abs/min/max, and math/cmath intrinsics;
        helper calls are kept as call nodes and are not expanded at compile
        time.

        Args:
            node: The ``ast.Call`` node.
            env: Mapping from local names to values.

        Returns:
            Value or tuple: The call result; multi-return helpers give a tuple of Values.

        Raises:
            FunctionCompileError: The call target, argument count/types or keyword usage is unsupported.
        """
        if (
            isinstance(node.func, ast.Attribute)
            and node.func.attr == "conjugate"
            and not node.args
            and not node.keywords
        ):
            value = self.expr(node.func.value, env)
            return self.node("conj", cast(Value, value).kind, (cast(Value, value),))
        module: object | None = None
        if isinstance(node.func, ast.Name):
            symbol = node.func.id
            target = self.namespace.get(symbol, getattr(builtins, symbol, None))
            for candidate in (math, cmath):
                for key in ELEMENTARY | REAL_EXTRA:
                    if target is getattr(candidate, key, None) and target is not None:
                        module, symbol = candidate, key
            if module is None and symbol in BUILTINS and target is getattr(builtins, symbol):
                target = None
        elif isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name):
            module = self.namespace.get(node.func.value.id)
            symbol = node.func.attr
            if module is not math and module is not cmath:
                self.fail(node, "arbitrary object methods are not pure math functions")
            target = None
        else:
            self.fail(node, "only named pure functions and math/cmath function calls are allowed")
        args = [self.expr(x, env) for x in node.args]
        if target is not None and module is None:
            source = self.source(target)
            names = [p.arg for p in source.tree.args.args]
            kwargs = {kw.arg: self.expr(kw.value, env) for kw in node.keywords}
            if None in kwargs or len(args) > len(names):
                self.fail(node, "invalid helper arguments")
            actual = dict(zip(names, args, strict=False))
            if set(actual) & set(kwargs) or set(kwargs) - set(names):
                self.fail(node, "duplicate or unknown helper arguments")
            actual.update(kwargs)  # type: ignore[arg-type]
            if any(not isinstance(x, Value) for x in actual.values()):
                self.fail(node, "helper arguments must be scalars")
            key = self.compile(source, {k: cast(Value, v).kind for k, v in actual.items()})
            child = self.functions[key]
            ordered = tuple(cast(Value, actual[p.name]) for p in child.parameters)
            values = tuple(
                self.node("call", child.nodes[ret].kind, ordered, (key, i))
                for i, ret in enumerate(child.returns)
            )
            return values[0] if len(values) == 1 else values
        if node.keywords:
            self.fail(node, "math intrinsics currently use positional arguments")
        if symbol == "complex" and module is None and len(args) in (1, 2):
            return self.node(
                "complex",
                "complex",
                (
                    cast(Value, args[0]),
                    cast(Value, args[1]) if len(args) > 1 else self.constant(0, node),
                ),
            )
        if symbol == "abs" and module is None and len(args) == 1:
            return self.node("abs", "real", cast("list[Value]", args))
        if symbol in {"min", "max"} and module is None and len(args) >= 2:
            result = args[0]
            for other in args[1:]:
                condition = self.binary("lt", result, other, node)
                result = (
                    self.choose(condition, result, other, node)
                    if symbol == "min"
                    else self.choose(condition, other, result, node)
                )
            return result
        if (module is not math and module is not cmath) or symbol not in ELEMENTARY | REAL_EXTRA:
            self.fail(node, "unsupported pure math call: " + symbol)
        if symbol in {"phase", "polar", "rect"} and module is not cmath:
            self.fail(node, "this function belongs to cmath")
        allowed = (
            (1, 2) if symbol == "log" else (2,) if symbol in {"atan2", "hypot", "rect"} else (1,)
        )
        if len(args) not in allowed or any(
            not isinstance(x, Value) or x.kind == "bool" for x in args
        ):
            self.fail(node, "invalid number or type of math function arguments")
        if module is math and any(cast(Value, x).kind == "complex" for x in args):
            self.fail(node, "math functions do not accept complex values; use cmath instead")
        if symbol == "polar":
            return (
                self.node("abs", "real", cast("list[Value]", args)),
                self.node("intrinsic", "real", cast("list[Value]", args), ("phase",)),
            )
        if symbol in {"atan2", "hypot", "rect"} and any(cast(Value, x).kind != "real" for x in args):
            self.fail(node, "this math function requires real arguments")
        kind = "real" if symbol in {"phase", "atan2", "hypot"} or module is math else "complex"
        return self.node("intrinsic", kind, cast("list[Value]", args), (symbol,))

    def assign(
        self,
        target: ast.expr,
        value: Value | tuple[Value, ...],
        env: dict[str, Value | tuple[Value, ...]],
    ) -> None:
        """Bind a value to an assignment target.

        Supports a single local name and tuple/list unpacking; mutation
        targets such as attributes are forbidden.

        Args:
            target: The assignment target AST node.
            value: The Value, or tuple of Values, to bind.
            env: The local name mapping to update.

        Raises:
            FunctionCompileError: The target is not a local name, or the unpacking structure does not match.
        """
        if isinstance(target, ast.Name):
            env[target.id] = value
        elif (
            isinstance(target, (ast.Tuple, ast.List))
            and isinstance(value, tuple)
            and len(target.elts) == len(value)
        ):
            for dest, item in zip(target.elts, value, strict=True):
                self.assign(dest, item, env)
        else:
            self.fail(target, "only local name assignment or tuple unpacking is allowed; object mutation is forbidden")

    def statements(
        self, body: list[ast.stmt], env: dict[str, Value | tuple[Value, ...]]
    ) -> tuple[dict[str, Value | tuple[Value, ...]], Value | tuple[Value, ...] | None]:
        """Interpret a statement block sequentially, returning the resulting environment and return value.

        Supports return, assignment/unpacking, augmented assignment, imports,
        structured if (the environments and return values of both branches
        are merged via select) and full unrolling of static bounded range
        loops.

        Args:
            body: List of statement AST nodes.
            env: The local name mapping when entering the block.

        Returns:
            tuple: ``(environment, return value)``; the return value is None when there is no return.

        Raises:
            FunctionCompileError: A statement leaves the restricted subset, an if return path is incomplete, or a loop exceeds the unroll limit.
        """
        env = dict(env)
        for i, node in enumerate(body):
            if (
                isinstance(node, ast.Expr)
                and isinstance(node.value, ast.Constant)
                and isinstance(node.value.value, str)
            ):
                continue
            if isinstance(node, ast.Return):
                return env, self.expr(node.value, env)  # type: ignore[arg-type]
            if isinstance(node, (ast.Assign, ast.AnnAssign)):
                value = self.expr(node.value, env)  # type: ignore[arg-type]
                for target in node.targets if isinstance(node, ast.Assign) else (node.target,):
                    self.assign(target, value, env)
            elif isinstance(node, ast.AugAssign) and isinstance(node.target, ast.Name):
                fake = ast.BinOp(left=node.target, op=node.op, right=node.value)
                ast.copy_location(fake, node)
                self.assign(node.target, self.expr(fake, env), env)
            elif isinstance(node, (ast.Import, ast.ImportFrom)):
                self.imports(node, self.namespace)
            elif isinstance(node, ast.If):
                test = self.expr(node.test, env)
                yes, yret = self.statements(node.body, env)
                no, nret = self.statements(node.orelse, env)
                if yret is not None or nret is not None:
                    if yret is None:
                        yes, yret = self.statements(body[i + 1 :], yes)
                    if nret is None:
                        no, nret = self.statements(body[i + 1 :], no)
                    if yret is None or nret is None:
                        self.fail(node, "incomplete return path in if statement")
                    return env, self.choose(test, yret, nret, node)
                common = yes.keys() & no.keys()
                env = {key: self.choose(test, yes[key], no[key], node) for key in sorted(common)}
            elif isinstance(node, ast.For):
                if (
                    not isinstance(node.iter, ast.Call)
                    or not isinstance(node.iter.func, ast.Name)
                    or node.iter.func.id != "range"
                    or node.iter.keywords
                    or node.orelse
                ):
                    self.fail(node, "only static range loops are supported")
                values = [self.literal(self.expr(arg, env), node) for arg in node.iter.args]
                if any(type(v) is not int for v in values):
                    self.fail(node, "range bounds must be generation-time integer arguments")
                indices = range(*cast("list[int]", values))
                if len(indices) > self.max_unroll:
                    self.fail(node, "static loop exceeds the unroll limit")
                for j in indices:
                    self.assign(node.target, self.constant(j, node), env)
                    env, result = self.statements(node.body, env)
                    if result is not None:
                        self.fail(node, "return is not supported inside a static loop")
            else:
                self.fail(node, "side effects and unsupported statements are forbidden: " + type(node).__name__)
        return env, None

    def program(self) -> MathProgram:
        """Collect the compilation results and return the validated MIR program.

        Returns:
            ``MathProgram``: The entry points at the compiled symbol name, and the function table is sorted by symbol name.

        Raises:
            ValidationError: The produced function graph failed MIR structural validation.
        """
        return MathProgram(
            self.entry, tuple(self.functions[k] for k in sorted(self.functions))
        ).validate()
