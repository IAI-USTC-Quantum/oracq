"读取 Python AST 的纯函数编译前端；不执行待编译函数。"

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

from pyqecclang.infrastructure.ir import ValidationError
from pyqecclang.infrastructure.mathfunc.graph import (
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
"""可编译为 intrinsic 节点的初等数学函数名集合；phase、polar、rect 仅限 cmath。"""
REAL_EXTRA = {"atan2", "hypot"}
"""仅实数域支持的补充二元数学函数名集合。"""
BUILTINS = {"abs", "complex", "min", "max"}
"""允许按名字直接调用的 Python 内建函数名集合。"""
KINDS = {"float": "real", "complex": "complex", "bool": "bool", "int": "real"}
"""Python 注解/类型名到 MIR 值类型的映射；未注解参数按 float 处理。"""


class FunctionCompileError(ValidationError):
    """数学函数前端在解析或编译受限 Python 源码失败时抛出的异常。"""


@dataclass(frozen=True)
class Value:
    """编译期 SSA 值：节点编号与其值类型的配对。

    Attributes:
        id: 该值对应节点在当前函数节点序列中的位置。
        kind: 值类型，为 real、complex 或 bool。
    """

    id: int
    kind: str


@dataclass
class Source:
    """待编译函数的定位记录：AST 定义及其求值命名空间。

    Attributes:
        tree: 函数的 ``ast.FunctionDef`` 定义节点。
        namespace: 解析名字引用时使用的全局与闭包变量映射。
    """

    tree: ast.FunctionDef
    namespace: dict[str, object]


class Frontend:
    """把受限 Python 纯函数源码编译为 MIR 数学函数图。

    读取并解释普通 Python 函数的 AST，不执行待编译函数。源码字符串中除
    入口外的其余 def 作为纯 helper 保留为独立 MIR 函数，调用以 call 节点
    表示，不在编译期展开。

    Args:
        source: 普通 Python 函数对象，或含若干 def 与 math/cmath 导入的源码字符串。
        inputs: 参数名到 real/complex/bool 或 ``Index(width)`` 的映射；省略时按注解推导。
        constants: 参数名到有限数值的映射，作为生成期常量绑定。
        helpers: 追加到入口命名空间的辅助函数。
        max_unroll: 静态 range 循环允许展开的最大迭代数。
        entry: 源码字符串形式下的入口函数名；省略时使用最后一个 def。

    Raises:
        FunctionCompileError: 源码有语法错误、包含不支持的语句或没有可用函数定义。
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
        """解析入口源并初始化函数表、编译缓存与活动栈。"""
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
                    f"Python 语法错误，第 {exc.lineno} 行：{exc.msg}"
                ) from exc
            ns: dict[str, object] = {"math": math, "cmath": cmath}
            for stmt in tree.body:
                if (
                    isinstance(stmt, ast.Expr)
                    and isinstance(stmt.value, ast.Constant)
                    and isinstance(stmt.value.value, str)
                ):
                    continue
                if isinstance(stmt, (ast.Import, ast.ImportFrom)):
                    self.imports(stmt, ns)
                elif isinstance(stmt, ast.FunctionDef):
                    self.sources[stmt.name] = Source(stmt, ns)
                else:
                    raise FunctionCompileError("源码只允许 math/cmath 导入和纯函数定义")
            if not self.sources:
                raise FunctionCompileError("源码中没有函数")
            for key, value in self.sources.items():
                ns[key] = value
            if entry is not None and entry not in self.sources:
                raise FunctionCompileError("源码中没有入口函数：" + entry)
            entry = self.sources[entry] if entry is not None else list(self.sources.values())[-1]  # type: ignore[assignment]
        else:
            entry = self.source(source)  # type: ignore[assignment]
        if helpers:
            cast(Source, entry).namespace.update(helpers)
        self.entry: str = self.compile(cast(Source, entry), inputs, constants or {})

    def source(self, function: object) -> Source:
        """把函数对象或 Source 归一化为定位记录。

        读取 ``function`` 的源码并唯一定位同名 def，结合其全局变量与
        闭包变量构成求值命名空间。

        Args:
            function: 普通 Python 函数对象，或已构造的 Source。

        Returns:
            Source: 可交给 compile 的定位记录。

        Raises:
            FunctionCompileError: 输入不是可读源码的普通函数、源码不可用或定义无法唯一定位。
        """
        if isinstance(function, Source):
            return function
        if not inspect.isfunction(function):
            raise FunctionCompileError(
                "只支持可读取源码的普通 Python 函数；也可传入 def 源码字符串"
            )
        try:
            tree = ast.parse(textwrap.dedent(inspect.getsource(function)))
        except (OSError, TypeError, IndentationError) as exc:
            raise FunctionCompileError("无法读取函数源码；请传入 def 源码字符串") from exc
        definitions = [
            n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == function.__name__
        ]
        if len(definitions) != 1:
            raise FunctionCompileError("无法唯一定位函数定义；不支持 lambda")
        closure = inspect.getclosurevars(function)
        return Source(definitions[0], {**function.__globals__, **closure.nonlocals})

    @staticmethod
    def imports(stmt: ast.Import | ast.ImportFrom, namespace: dict[str, object]) -> None:
        """处理一条 math/cmath 导入语句并写入命名空间。

        Args:
            stmt: ``ast.Import`` 或 ``ast.ImportFrom`` 节点。
            namespace: 接收导入名字的映射。

        Raises:
            FunctionCompileError: 导入目标不是 math/cmath，或导入了不支持的数学名字。
        """
        if isinstance(stmt, ast.Import):
            for alias in stmt.names:
                if alias.name not in {"math", "cmath"}:
                    raise FunctionCompileError("只允许导入 math/cmath")
                namespace[alias.asname or alias.name] = {"math": math, "cmath": cmath}[alias.name]
        else:
            if stmt.module not in {"math", "cmath"} or stmt.level:
                raise FunctionCompileError("只允许从 math/cmath 导入")
            module = {"math": math, "cmath": cmath}[stmt.module]
            for alias in stmt.names:
                if alias.name not in ELEMENTARY | REAL_EXTRA | {"pi", "e", "tau"}:
                    raise FunctionCompileError("未支持的数学导入：" + alias.name)
                namespace[alias.asname or alias.name] = getattr(module, alias.name)

    def fail(self, node: ast.AST, message: str) -> NoReturn:
        """抛出带源码位置的错误。

        Args:
            node: 出错位置对应的 AST 节点；无位置信息时以问号代替行号。
            message: 错误说明。

        Returns:
            NoReturn: 不返回；总是以抛出 ``FunctionCompileError`` 结束。

        Raises:
            FunctionCompileError: 总是抛出，消息带当前函数标签与行号前缀。
        """
        raise FunctionCompileError(f"{self.label}:{getattr(node, 'lineno', '?')}：{message}")

    def compile(
        self,
        source: Callable[..., object] | Source,
        inputs: Mapping[str, Index | type | str] | None = None,
        constants: Mapping[str, bool | int | float | complex] | None = None,
    ) -> str:
        """编译一个函数定义为 MIR 并返回其符号名。

        推导动态输入类型，把函数体解释为 SSA 节点序列，并按 AST、输入、
        常量与捕获数值缓存重复编译；符号名由函数图内容的哈希得到。

        Args:
            source: Python 函数对象、Source 或源码字符串。
            inputs: 参数名到类型的映射；为 None 时按注解推导，未注解默认 real。
            constants: 参数名到有限数值的映射，绑定为生成期常量。

        Returns:
            str: 已登记到 ``self.functions`` 的函数符号名。

        Raises:
            FunctionCompileError: 参数声明、语句或表达式超出受限子集，或检测到递归 helper。
        """
        source = self.source(source)
        constants = constants or {}
        tree = source.tree
        if tree.args.vararg or tree.args.kwarg or tree.args.kwonlyargs or tree.args.posonlyargs:
            raise FunctionCompileError("函数需要普通具名参数；不支持 *args/**kwargs/keyword-only")
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
            raise FunctionCompileError("参数不能同时为动态输入和生成期常量")
        if set(inputs) - {p.arg for p in parameters}:
            raise FunctionCompileError("inputs 包含未知参数")
        dynamic: list[Parameter] = []
        for p in parameters:
            if p.arg in inputs:
                spec: Index | type | str | Parameter = inputs[p.arg]
                if isinstance(spec, Index):
                    dynamic.append(Parameter(p.arg, "index", spec.width))
                else:
                    kind = KINDS.get(getattr(spec, "__name__", spec), spec)  # type: ignore[arg-type]
                    if kind not in {"real", "complex", "bool"}:
                        raise FunctionCompileError("输入类型应为 real/complex/bool 或 Index(width)")
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
            raise FunctionCompileError("不支持递归 helper：" + tree.name)
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
                    self.fail(p, "参数未声明为动态输入或生成期常量：" + p.arg)
            _, result = self.statements(tree.body, env)
            if result is None:
                self.fail(tree, "所有控制路径必须返回数值")
            values = result if isinstance(result, tuple) else (result,)
            if any(not isinstance(v, Value) for v in values):
                self.fail(tree, "返回值必须是标量或一层 tuple")
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
        """构造（或复用）一个 SSA 节点并返回其值。

        操作、类型、输入与 data 都相同的节点在当前函数内只保留一份。

        Args:
            op: 节点操作名，如 add、select 或 intrinsic。
            kind: 结果值类型 real/complex/bool。
            args: 操作数 Value 列表。
            data: 附着于节点的额外数据，如常量数值或 intrinsic 名。

        Returns:
            Value: 指向该节点的编译期值。
        """
        values = tuple(v.id for v in args)
        key = MathNode(op, kind, values, tuple(data))
        if key not in self.intern:
            self.intern[key] = len(self.nodes)
            self.nodes.append(key)
        return Value(self.intern[key], kind)

    def constant(self, value: bool | int | float | complex, node: ast.AST) -> Value:
        """把一个有限数值固化为 const 节点。

        Args:
            value: bool、int、float 或 complex 常量。
            node: 报错时定位用的 AST 节点。

        Returns:
            Value: 与常量类型一致的常量值。

        Raises:
            FunctionCompileError: 数值类型不受支持，或为 NaN/Inf 等非有限值。
        """
        if type(value) not in (bool, int, float, complex):
            self.fail(node, "只允许有限数值常量")
        if not (math.isfinite(value.real) and math.isfinite(value.imag)):
            self.fail(node, "固定点不表示 NaN/Inf")
        kind = "bool" if type(value) is bool else "complex" if type(value) is complex else "real"
        return self.node(
            "const",
            kind,
            data=(value.real, value.imag)
            if kind == "complex"
            else (cast("bool | int | float", value),),
        )

    def literal(self, value: Value | tuple[Value, ...], node: ast.AST) -> bool | int | float | complex:
        """取出编译期常量节点承载的 Python 数值。

        Args:
            value: 应为 const 节点的 Value。
            node: 报错时定位用的 AST 节点。

        Returns:
            该常量的原始 Python 数值；复数以 complex 表示。

        Raises:
            FunctionCompileError: 该值不是常量节点。
        """
        if isinstance(value, Value) and self.nodes[value.id].op == "const":
            data = self.nodes[value.id].data
            return complex(*data) if value.kind == "complex" else data[0]
        self.fail(node, "这里需要生成期数值常量")

    def binary(
        self, op: str, a: Value | tuple[Value, ...], b: Value | tuple[Value, ...], node: ast.AST
    ) -> Value:
        """构造二元运算节点并推导结果类型。

        算术在 real/complex 间按提升规则定型，比较与布尔运算产出 bool；
        两侧均为常量时在生成期直接求值并折叠为一个常量节点。

        Args:
            op: 运算名，如 add、sub、mul、div、pow、lt、eq、and、or。
            a: 左操作数 Value。
            b: 右操作数 Value。
            node: 报错时定位用的 AST 节点。

        Returns:
            Value: 运算结果值。

        Raises:
            FunctionCompileError: 操作数不是标量、类型组合非法或常量折叠遇到无定义运算。
        """
        if not isinstance(a, Value) or not isinstance(b, Value):
            self.fail(node, "算术只接收标量")
        if op in {"and", "or"}:
            if a.kind != "bool" or b.kind != "bool":
                self.fail(node, "布尔操作数类型无效")
            kind = "bool"
        elif op in {"lt", "eq"}:
            if op == "lt" and "complex" in (a.kind, b.kind):
                self.fail(node, "复数不支持顺序比较")
            kind = "bool"
        else:
            if "bool" in (a.kind, b.kind):
                self.fail(node, "布尔值需要显式选择为数值")
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
                self.fail(node, "常量表达式无定义：" + str(exc))
        return self.node(op, kind, (a, b))

    def choose(
        self,
        test: Value | tuple[Value, ...],
        yes: Value | tuple[Value, ...],
        no: Value | tuple[Value, ...],
        node: ast.AST,
    ) -> Value | tuple[Value, ...]:
        """构造按布尔条件选择分支值的 select 节点。

        分支为 tuple 时逐元素选择；real 与 complex 分支按提升规则统一类型。

        Args:
            test: bool 类型的条件值。
            yes: 条件为真时的值或值的 tuple。
            no: 条件为假时的值或值的 tuple。
            node: 报错时定位用的 AST 节点。

        Returns:
            Value 或 tuple: 选择结果，形状与分支一致。

        Raises:
            FunctionCompileError: 条件不是布尔，或两分支结构与类型不一致。
        """
        if not isinstance(test, Value) or test.kind != "bool":
            self.fail(node, "量子条件必须是布尔表达式")
        if isinstance(yes, tuple) and isinstance(no, tuple) and len(yes) == len(no):
            return tuple(cast(Value, self.choose(test, a, b, node)) for a, b in zip(yes, no, strict=True))
        if not isinstance(yes, Value) or not isinstance(no, Value):
            self.fail(node, "分支返回/赋值结构不一致")
        if {yes.kind, no.kind} <= {"real", "complex"}:
            kind = "complex" if "complex" in (yes.kind, no.kind) else "real"
        elif yes.kind == no.kind:
            kind = yes.kind
        else:
            self.fail(node, "分支类型不一致")
        return self.node("select", kind, (test, yes, no))

    def expr(
        self, node: ast.expr, env: dict[str, Value | tuple[Value, ...]]
    ) -> Value | tuple[Value, ...]:
        """把一个表达式 AST 解释为编译期值。

        支持常量、名字、tuple/list、生成期 tuple 下标、math/cmath 常数
        属性、算术、比较、布尔短路、条件表达式与函数调用。

        Args:
            node: 表达式 AST 节点。
            env: 局部名字到 Value 或其 tuple 的映射。

        Returns:
            Value 或 tuple: 表达式的值；tuple/list 表达式产出 Value 的 tuple。

        Raises:
            FunctionCompileError: 表达式超出受限子集。
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
                self.fail(node, "只支持生成期 tuple 下标")
            try:
                return value[index]
            except IndexError:
                self.fail(node, "tuple 下标越界")
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
                self.fail(node, "只支持 .real/.imag 数值属性")
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
                self.fail(node, "未支持的二元运算")
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
            self.fail(node, "未支持的单目运算")
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
                    self.fail(node, "未支持的比较")
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
                # Python 短路语义也应用于状态路径。
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
        self.fail(node, "未支持的表达式：" + type(node).__name__)

    def call(
        self, node: ast.Call, env: dict[str, Value | tuple[Value, ...]]
    ) -> Value | tuple[Value, ...]:
        """解释函数调用表达式。

        处理 ``conjugate`` 方法、具名纯 helper 调用、白名单内建
        complex/abs/min/max，以及 math/cmath 数学 intrinsic；helper 调用
        保留为 call 节点，不在编译期展开。

        Args:
            node: ``ast.Call`` 节点。
            env: 局部名字到值的映射。

        Returns:
            Value 或 tuple: 调用结果；多返回 helper 为 Value 的 tuple。

        Raises:
            FunctionCompileError: 调用目标、参数数目/类型或关键字用法不受支持。
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
                self.fail(node, "任意对象方法不属于纯数学函数")
            target = None
        else:
            self.fail(node, "只允许具名纯函数和 math/cmath 函数调用")
        args = [self.expr(x, env) for x in node.args]
        if target is not None and module is None:
            source = self.source(target)
            names = [p.arg for p in source.tree.args.args]
            kwargs = {kw.arg: self.expr(kw.value, env) for kw in node.keywords}
            if None in kwargs or len(args) > len(names):
                self.fail(node, "helper 参数无效")
            actual = dict(zip(names, args, strict=False))
            if set(actual) & set(kwargs) or set(kwargs) - set(names):
                self.fail(node, "helper 参数重复或未知")
            actual.update(kwargs)  # type: ignore[arg-type]
            if any(not isinstance(x, Value) for x in actual.values()):
                self.fail(node, "helper 参数必须为标量")
            key = self.compile(source, {k: cast(Value, v).kind for k, v in actual.items()})
            child = self.functions[key]
            ordered = tuple(cast(Value, actual[p.name]) for p in child.parameters)
            values = tuple(
                self.node("call", child.nodes[ret].kind, ordered, (key, i))
                for i, ret in enumerate(child.returns)
            )
            return values[0] if len(values) == 1 else values
        if node.keywords:
            self.fail(node, "当前数学 intrinsic 使用位置参数")
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
            self.fail(node, "未支持的纯数学调用：" + symbol)
        if symbol in {"phase", "polar", "rect"} and module is not cmath:
            self.fail(node, "该函数属于 cmath")
        allowed = (
            (1, 2) if symbol == "log" else (2,) if symbol in {"atan2", "hypot", "rect"} else (1,)
        )
        if len(args) not in allowed or any(
            not isinstance(x, Value) or x.kind == "bool" for x in args
        ):
            self.fail(node, "数学函数参数数目/类型无效")
        if module is math and any(cast(Value, x).kind == "complex" for x in args):
            self.fail(node, "math 接口不接收复数；请使用 cmath")
        if symbol == "polar":
            return (
                self.node("abs", "real", cast("list[Value]", args)),
                self.node("intrinsic", "real", cast("list[Value]", args), ("phase",)),
            )
        if symbol in {"atan2", "hypot", "rect"} and any(cast(Value, x).kind != "real" for x in args):
            self.fail(node, "此数学函数要求实数参数")
        kind = "real" if symbol in {"phase", "atan2", "hypot"} or module is math else "complex"
        return self.node("intrinsic", kind, cast("list[Value]", args), (symbol,))

    def assign(
        self,
        target: ast.expr,
        value: Value | tuple[Value, ...],
        env: dict[str, Value | tuple[Value, ...]],
    ) -> None:
        """把一个值绑定到赋值目标。

        支持单个局部名字与 tuple/list 解包，禁止属性等突变目标。

        Args:
            target: 赋值目标 AST 节点。
            value: 待绑定的 Value 或其 tuple。
            env: 待更新的局部名字映射。

        Raises:
            FunctionCompileError: 目标不是局部名字，或解包结构不匹配。
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
            self.fail(target, "只允许局部名字赋值或 tuple 解包，禁止对象突变")

    def statements(
        self, body: list[ast.stmt], env: dict[str, Value | tuple[Value, ...]]
    ) -> tuple[dict[str, Value | tuple[Value, ...]], Value | tuple[Value, ...] | None]:
        """顺序解释一个语句块，返回执行后的环境与返回值。

        支持 return、赋值/解包、增量赋值、导入、结构化 if（两分支的环境
        与返回值经 select 合并）以及静态有界 range 循环的完全展开。

        Args:
            body: 语句 AST 节点列表。
            env: 进入该块时的局部名字映射。

        Returns:
            tuple: ``(环境, 返回值)``；没有 return 时返回值为 None。

        Raises:
            FunctionCompileError: 语句超出受限子集、if 返回路径不完整或循环超出展开上限。
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
                        self.fail(node, "if 返回路径不完整")
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
                    self.fail(node, "只支持静态 range 循环")
                values = [self.literal(self.expr(arg, env), node) for arg in node.iter.args]
                if any(type(v) is not int for v in values):
                    self.fail(node, "range 边界必须为整数生成参数")
                indices = range(*cast("list[int]", values))
                if len(indices) > self.max_unroll:
                    self.fail(node, "静态循环超过展开上限")
                for j in indices:
                    self.assign(node.target, self.constant(j, node), env)
                    env, result = self.statements(node.body, env)
                    if result is not None:
                        self.fail(node, "静态循环内部不支持 return")
            else:
                self.fail(node, "禁止副作用或未支持的语句：" + type(node).__name__)
        return env, None

    def program(self) -> MathProgram:
        """汇总编译结果并返回经校验的 MIR 程序。

        Returns:
            ``MathProgram``：入口指向编译得到的符号名，函数表按符号名排序。

        Raises:
            ValidationError: 生成的函数图未通过 MIR 结构校验。
        """
        return MathProgram(
            self.entry, tuple(self.functions[k] for k in sorted(self.functions))
        ).validate()
