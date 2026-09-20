"MIR 0.1：纯数学函数的有类型 SSA 图；不含 Python 可调用对象。"

from __future__ import annotations

import json
from dataclasses import asdict, dataclass

from pyqecclang.infrastructure.ir import ValidationError
from pyqecclang.infrastructure.validation import name


@dataclass(frozen=True)
class Index:
    """无符号整数索引参数的位宽说明，用于 ``inputs`` 映射。

    为行列索引等无符号整数提供较短的公开寄存器；进入数学计算时
    转换为当前定点表示，格式必须容纳其完整范围。

    Attributes:
        width: 无符号位宽，范围为 1..64。
    """

    width: int

    def __post_init__(self):
        if type(self.width) is not int or not 1 <= self.width <= 64:
            raise ValidationError("Index 位宽必须为 1..64")


@dataclass(frozen=True)
class Parameter:
    """MIR 函数形参的名字与数学类型。

    Attributes:
        name: 参数名，在所在函数内唯一。
        kind: 参数类型，为 real、complex、bool 或 index。
        width: 仅 index 参数使用的无符号位宽；其他类型必须为 0。
    """

    name: str
    kind: str = "real"
    width: int = 0


@dataclass(frozen=True)
class MathNode:
    """MIR 顺序 SSA 图中的一个值节点。

    节点编号为其在函数 ``nodes`` 序列中的位置；``args`` 只能引用
    编号更小的节点。

    Attributes:
        op: 操作名，如 ``input``、``const``、四则运算或 ``call``。
        kind: 结果值类型，为 real、complex 或 bool。
        args: 操作数节点编号组成的 tuple。
        data: 附加数据：input 存参数名，const 存常量编码，intrinsic 存数学函数名，call 存目标函数符号与返回索引。
    """

    op: str
    kind: str
    args: tuple[int, ...] = ()
    data: tuple = ()


@dataclass(frozen=True)
class MathFunction:
    """MIR 中的一个纯数学函数：形参表与顺序 SSA 节点序列。

    Attributes:
        name: 函数符号名，在程序内唯一。
        label: 源码中的函数名，用于展示与模块属性。
        parameters: 按形参顺序排列的参数表。
        nodes: 顺序 SSA 节点序列。
        returns: 非空的返回节点编号 tuple；多个编号表示多个独立结果。
    """

    name: str
    label: str
    parameters: tuple[Parameter, ...]
    nodes: tuple[MathNode, ...]
    returns: tuple[int, ...]


@dataclass(frozen=True)
class MathProgram:
    """MIR 0.1 顶层对象：入口函数符号与全部函数表。

    Attributes:
        entry: 入口函数的符号名，必须出现在 ``functions`` 中。
        functions: 全部 ``MathFunction`` 组成的不可变 tuple。
        version: MIR 版本号，当前为 0.1。
    """

    entry: str
    functions: tuple[MathFunction, ...]
    version: str = "0.1"

    @property
    def function_map(self):
        """函数符号名到 ``MathFunction`` 的映射。"""
        return {f.name: f for f in self.functions}

    def validate(self):
        """核对程序的结构与类型约束并返回自身。

        覆盖版本号、函数表不可变性与入口有效性、参数唯一性与类型、
        顺序 SSA 引用、操作元数与类型提升、常量有限性、调用接口
        一致性，以及函数图无递归。

        Returns:
            MathProgram: 校验通过的自身。

        Raises:
            ValidationError: 任一结构或类型约束不满足。
        """
        if self.version != "0.1":
            raise ValidationError("未知 MIR 版本")
        if type(self.functions) is not tuple:
            raise ValidationError("MIR 函数表必须不可变")
        functions = self.function_map
        if len(functions) != len(self.functions) or self.entry not in functions:
            raise ValidationError("MIR 入口或函数表无效")
        for f in self.functions:
            name(f.name)
            if not isinstance(f.label, str) or any(
                type(v) is not tuple for v in (f.parameters, f.nodes, f.returns)
            ):
                raise ValidationError("MIR 函数记录必须具有字符串标签和不可变数组")
            params = {p.name: p for p in f.parameters}
            if len(params) != len(f.parameters):
                raise ValidationError("MIR 参数重名")
            for p in f.parameters:
                name(p.name)
                if type(p.width) is not int or (p.kind != "index" and p.width != 0):
                    raise ValidationError("MIR 非 index 参数 width 必须为 0")
                if p.kind not in {"real", "complex", "bool", "index"} or (
                    p.kind == "index" and not 1 <= p.width <= 64
                ):
                    raise ValidationError("MIR 参数类型无效")
            if not f.returns or any(
                type(i) is not int or not 0 <= i < len(f.nodes) for i in f.returns
            ):
                raise ValidationError("MIR 返回值无效")
            for i, node in enumerate(f.nodes):
                if type(node.args) is not tuple or type(node.data) is not tuple:
                    raise ValidationError("MIR 节点必须不可变")
                if node.kind not in {"real", "complex", "bool"}:
                    raise ValidationError("MIR 值类型无效")
                if any(type(a) is not int or not 0 <= a < i for a in node.args):
                    raise ValidationError("MIR 不是顺序 SSA 图")
                arities = {
                    "input": 0,
                    "const": 0,
                    "neg": 1,
                    "not": 1,
                    "real": 1,
                    "imag": 1,
                    "conj": 1,
                    "abs": 1,
                    "add": 2,
                    "sub": 2,
                    "mul": 2,
                    "div": 2,
                    "pow": 2,
                    "lt": 2,
                    "eq": 2,
                    "and": 2,
                    "or": 2,
                    "complex": 2,
                    "select": 3,
                }
                if node.op in arities and len(node.args) != arities[node.op]:
                    raise ValidationError("MIR 操作元数无效")
                if node.op not in {*arities, "intrinsic", "call"}:
                    raise ValidationError("未知 MIR 操作")
                if node.op == "input" and (len(node.data) != 1 or node.data[0] not in params):
                    raise ValidationError("MIR 输入引用无效")
                if node.op == "const":
                    if len(node.data) != (2 if node.kind == "complex" else 1):
                        raise ValidationError("MIR 常量编码无效")
                    import math

                    if any(
                        type(x) not in (bool, int, float) or not math.isfinite(x) for x in node.data
                    ):
                        raise ValidationError("MIR 只接收有限数值常量")
                if node.op == "call":
                    if len(node.data) != 2 or node.data[0] not in functions:
                        raise ValidationError("MIR 调用目标无效")
                    child = functions[node.data[0]]
                    if len(node.args) != len(child.parameters) or not 0 <= node.data[1] < len(
                        child.returns
                    ):
                        raise ValidationError("MIR 调用布局无效")
                    if node.kind != child.nodes[child.returns[node.data[1]]].kind:
                        raise ValidationError("MIR 调用返回类型无效")
        intrinsic_names = {
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
            "rect",
            "atan2",
            "hypot",
        }
        for f in self.functions:
            params = {p.name: p for p in f.parameters}
            for node in f.nodes:
                kinds = [f.nodes[i].kind for i in node.args]
                numeric = {"real", "complex"}
                valid = True
                if node.op in {"add", "sub", "mul", "div", "pow"}:
                    valid = set(kinds) <= numeric and node.kind == (
                        "complex" if "complex" in kinds else "real"
                    )
                elif node.op == "lt":
                    valid = kinds == ["real", "real"] and node.kind == "bool"
                elif node.op == "eq":
                    valid = (
                        set(kinds) <= numeric or kinds == ["bool", "bool"]
                    ) and node.kind == "bool"
                elif node.op in {"and", "or", "not"}:
                    valid = set(kinds) == {"bool"} and node.kind == "bool"
                elif node.op in {"neg", "conj"}:
                    valid = kinds[0] in numeric and node.kind == kinds[0]
                elif node.op in {"real", "imag", "abs"}:
                    valid = kinds[0] in numeric and node.kind == "real"
                elif node.op == "complex":
                    valid = kinds == ["real", "real"] and node.kind == "complex"
                elif node.op == "select":
                    promoted = "complex" if "complex" in kinds[1:] else kinds[1]
                    valid = (
                        kinds[0] == "bool"
                        and (kinds[1] == kinds[2] or set(kinds[1:]) <= numeric)
                        and node.kind == promoted
                    )
                elif node.op == "input":
                    parameter = params[node.data[0]]
                    valid = node.kind == ("real" if parameter.kind == "index" else parameter.kind)
                elif node.op == "const" and node.kind == "bool":
                    valid = type(node.data[0]) is bool
                elif node.op == "intrinsic":
                    if len(node.data) != 1 or node.data[0] not in intrinsic_names:
                        raise ValidationError("未知 MIR 数学 intrinsic")
                    function = node.data[0]
                    counts = (
                        {1, 2}
                        if function == "log"
                        else {2}
                        if function in {"rect", "atan2", "hypot"}
                        else {1}
                    )
                    valid = len(kinds) in counts and set(kinds) <= numeric
                    if function in {"rect", "atan2", "hypot"}:
                        valid &= set(kinds) == {"real"}
                    if function in {"phase", "atan2", "hypot"}:
                        valid &= node.kind == "real"
                    elif function == "rect":
                        valid &= node.kind == "complex"
                    else:
                        valid &= node.kind in numeric and (
                            node.kind == "complex" or set(kinds) == {"real"}
                        )
                elif node.op == "call":
                    child = functions[node.data[0]]
                    valid = kinds == [
                        ("real" if p.kind == "index" else p.kind) for p in child.parameters
                    ]
                if node.op not in {"input", "const", "call", "intrinsic"} and node.data:
                    valid = False
                if not valid:
                    raise ValidationError("MIR 操作类型/数据约束不满足：" + node.op)
        active, visited = set(), set()

        def visit(key):
            if key in active:
                raise ValidationError("MIR 不支持递归")
            if key in visited:
                return
            active.add(key)
            for node in functions[key].nodes:
                if node.op == "call":
                    visit(node.data[0])
            active.remove(key)
            visited.add(key)

        for key in functions:
            visit(key)
        return self

    def dumps(self):
        """先校验再把程序序列化为规范化 JSON 文本。

        Returns:
            str: 按键排序、两格缩进、以换行结尾的 JSON 文本；非 ASCII 字符原样保留。

        Raises:
            ValidationError: 程序未通过校验。
        """
        self.validate()
        return (
            json.dumps(asdict(self), ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
            + "\n"
        )

    @classmethod
    def loads(cls, text):
        """从 JSON 文本反序列化并校验 MIR 程序。

        Args:
            text: ``MathProgram.dumps`` 生成的 JSON 文本。

        Returns:
            MathProgram: 校验通过的程序对象。

        Raises:
            ValidationError: 文本不是合法的 MIR JSON，或未通过校验。
        """
        try:
            raw = json.loads(text)
            if set(raw) != {"entry", "functions", "version"}:
                raise ValidationError("MIR 顶层字段无效")
            functions = []
            for f in raw["functions"]:
                if set(f) != {"name", "label", "parameters", "nodes", "returns"} or any(
                    set(n) != {"op", "kind", "args", "data"} for n in f["nodes"]
                ):
                    raise ValidationError("MIR 函数或节点包含未知字段")
                functions.append(
                    MathFunction(
                        f["name"],
                        f["label"],
                        tuple(Parameter(**p) for p in f["parameters"]),
                        tuple(
                            MathNode(n["op"], n["kind"], tuple(n["args"]), tuple(n["data"]))
                            for n in f["nodes"]
                        ),
                        tuple(f["returns"]),
                    )
                )
            return cls(raw["entry"], tuple(functions), raw["version"]).validate()
        except (KeyError, TypeError, ValueError, IndexError, OverflowError, RecursionError) as exc:
            raise ValidationError(f"非法 MIR JSON：{exc}") from exc
