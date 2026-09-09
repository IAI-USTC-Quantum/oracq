"""读取 Python AST 的纯函数编译前端；不执行待编译函数。"""

from __future__ import annotations

import ast
import builtins
import cmath
import hashlib
import inspect
import math
import textwrap
from dataclasses import dataclass

from ..ir import ValidationError
from .graph import Index, MathFunction, MathNode, MathProgram, Parameter

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
REAL_EXTRA = {"atan2", "hypot"}
BUILTINS = {"abs", "complex", "min", "max"}
KINDS = {"float": "real", "complex": "complex", "bool": "bool", "int": "real"}


class FunctionCompileError(ValidationError):
    pass


@dataclass(frozen=True)
class Value:
    id: int
    kind: str


@dataclass
class Source:
    tree: object
    namespace: dict


class Frontend:
    def __init__(
        self, source, *, inputs=None, constants=None, helpers=None, max_unroll=128, entry=None
    ):
        self.functions, self.cache, self.active = {}, {}, set()
        self.max_unroll = max_unroll
        self.sources = {}
        if isinstance(source, str):
            try:
                tree = ast.parse(textwrap.dedent(source))
            except SyntaxError as exc:
                raise FunctionCompileError(
                    f"Python 语法错误，第 {exc.lineno} 行：{exc.msg}"
                ) from exc
            ns = {"math": math, "cmath": cmath}
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
            entry = self.sources[entry] if entry is not None else list(self.sources.values())[-1]
        else:
            entry = self.source(source)
        if helpers:
            entry.namespace.update(helpers)
        self.entry = self.compile(entry, inputs, constants or {})

    def source(self, function):
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
    def imports(stmt, namespace):
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

    def fail(self, node, message):
        raise FunctionCompileError(f"{self.label}:{getattr(node, 'lineno', '?')}：{message}")

    def compile(self, source, inputs=None, constants=None):
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
        dynamic = []
        for p in parameters:
            if p.arg in inputs:
                spec = inputs[p.arg]
                if isinstance(spec, Index):
                    dynamic.append(Parameter(p.arg, "index", spec.width))
                else:
                    kind = KINDS.get(getattr(spec, "__name__", spec), spec)
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
        self.nodes, self.intern = [], {}
        self.namespace = {**source.namespace, **constants}
        self.label = tree.name
        env = {}
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

    def node(self, op, kind, args=(), data=()):
        values = tuple(v.id for v in args)
        key = MathNode(op, kind, values, tuple(data))
        if key not in self.intern:
            self.intern[key] = len(self.nodes)
            self.nodes.append(key)
        return Value(self.intern[key], kind)

    def constant(self, value, node):
        if type(value) not in (bool, int, float, complex):
            self.fail(node, "只允许有限数值常量")
        if not (math.isfinite(value.real) and math.isfinite(value.imag)):
            self.fail(node, "固定点不表示 NaN/Inf")
        kind = "bool" if type(value) is bool else "complex" if type(value) is complex else "real"
        return self.node(
            "const", kind, data=(value.real, value.imag) if kind == "complex" else (value,)
        )

    def literal(self, value, node):
        if isinstance(value, Value) and self.nodes[value.id].op == "const":
            data = self.nodes[value.id].data
            return complex(*data) if value.kind == "complex" else data[0]
        self.fail(node, "这里需要生成期数值常量")

    def binary(self, op, a, b, node):
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
            functions = {
                "add": lambda: av + bv,
                "sub": lambda: av - bv,
                "mul": lambda: av * bv,
                "div": lambda: av / bv,
                "pow": lambda: av**bv,
                "lt": lambda: av < bv,
                "eq": lambda: av == bv,
                "and": lambda: av and bv,
                "or": lambda: av or bv,
            }
            try:
                return self.constant(functions[op](), node)
            except (ArithmeticError, ValueError) as exc:
                self.fail(node, "常量表达式无定义：" + str(exc))
        return self.node(op, kind, (a, b))

    def choose(self, test, yes, no, node):
        if not isinstance(test, Value) or test.kind != "bool":
            self.fail(node, "量子条件必须是布尔表达式")
        if isinstance(yes, tuple) and isinstance(no, tuple) and len(yes) == len(no):
            return tuple(self.choose(test, a, b, node) for a, b in zip(yes, no, strict=True))
        if not isinstance(yes, Value) or not isinstance(no, Value):
            self.fail(node, "分支返回/赋值结构不一致")
        if {yes.kind, no.kind} <= {"real", "complex"}:
            kind = "complex" if "complex" in (yes.kind, no.kind) else "real"
        elif yes.kind == no.kind:
            kind = yes.kind
        else:
            self.fail(node, "分支类型不一致")
        return self.node("select", kind, (test, yes, no))

    def expr(self, node, env):
        if isinstance(node, ast.Constant):
            return self.constant(node.value, node)
        if isinstance(node, ast.Name):
            if node.id in env:
                return env[node.id]
            value = self.namespace.get(node.id)
            return self.constant(value, node)
        if isinstance(node, (ast.Tuple, ast.List)):
            return tuple(self.expr(x, env) for x in node.elts)
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
                if self.nodes[value.id].op == "const":
                    return self.constant(-self.literal(value, node), node)
                return self.node("neg", value.kind, (value,))
            if isinstance(node.op, ast.Not) and value.kind == "bool":
                return self.node("not", "bool", (value,))
            self.fail(node, "未支持的单目运算")
        if isinstance(node, ast.Compare):
            pairs = []
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
            result = pairs[0]
            for test in pairs[1:]:
                result = self.node("and", "bool", (result, test))
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

    def call(self, node, env):
        if (
            isinstance(node.func, ast.Attribute)
            and node.func.attr == "conjugate"
            and not node.args
            and not node.keywords
        ):
            value = self.expr(node.func.value, env)
            return self.node("conj", value.kind, (value,))
        module = None
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
            actual.update(kwargs)
            if any(not isinstance(x, Value) for x in actual.values()):
                self.fail(node, "helper 参数必须为标量")
            key = self.compile(source, {k: v.kind for k, v in actual.items()})
            child = self.functions[key]
            ordered = tuple(actual[p.name] for p in child.parameters)
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
                (args[0], args[1] if len(args) > 1 else self.constant(0, node)),
            )
        if symbol == "abs" and module is None and len(args) == 1:
            return self.node("abs", "real", args)
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
        if module is math and any(x.kind == "complex" for x in args):
            self.fail(node, "math 接口不接收复数；请使用 cmath")
        if symbol == "polar":
            return (
                self.node("abs", "real", args),
                self.node("intrinsic", "real", args, ("phase",)),
            )
        if symbol in {"atan2", "hypot", "rect"} and any(x.kind != "real" for x in args):
            self.fail(node, "此数学函数要求实数参数")
        kind = "real" if symbol in {"phase", "atan2", "hypot"} or module is math else "complex"
        return self.node("intrinsic", kind, args, (symbol,))

    def assign(self, target, value, env):
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

    def statements(self, body, env):
        env = dict(env)
        for i, node in enumerate(body):
            if (
                isinstance(node, ast.Expr)
                and isinstance(node.value, ast.Constant)
                and isinstance(node.value.value, str)
            ):
                continue
            if isinstance(node, ast.Return):
                return env, self.expr(node.value, env)
            if isinstance(node, (ast.Assign, ast.AnnAssign)):
                value = self.expr(node.value, env)
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
                indices = range(*values)
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

    def program(self):
        return MathProgram(
            self.entry, tuple(self.functions[k] for k in sorted(self.functions))
        ).validate()
