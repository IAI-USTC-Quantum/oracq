"MIR 0.1: a typed SSA graph of pure math functions; it holds no Python callables."

from __future__ import annotations

import json
from dataclasses import asdict, dataclass

from oracq.infrastructure.ir import ValidationError
from oracq.infrastructure.validation import name


@dataclass(frozen=True)
class Index:
    """Bit width declaration of an unsigned integer index parameter, for the ``inputs`` mapping.

    Provides a shorter public register for unsigned integers such as row
    and column indices; it is converted to the current fixed-point
    representation when entering math computation, and the format must
    hold its full range.

    Attributes:
        width: Unsigned bit width, in the range 1..64.
    """

    width: int

    def __post_init__(self) -> None:
        """Validate that the bit width is in 1..64."""
        if type(self.width) is not int or not 1 <= self.width <= 64:
            raise ValidationError("Index bit width must be between 1 and 64")


@dataclass(frozen=True)
class Parameter:
    """Name and math kind of an MIR function formal parameter.

    Attributes:
        name: Parameter name, unique within its function.
        kind: Parameter kind: real, complex, bool or index.
        width: Unsigned bit width used only by index parameters; must be 0 for other kinds.
    """

    name: str
    kind: str = "real"
    width: int = 0


@dataclass(frozen=True)
class MathNode:
    """One value node of the MIR sequential SSA graph.

    The node index is its position in the function's ``nodes`` sequence;
    ``args`` may only reference nodes with smaller indices.

    Attributes:
        op: Operation name, such as ``input``, ``const``, an arithmetic operation or ``call``.
        kind: Result value kind: real, complex or bool.
        args: Tuple of operand node indices.
        data: Attached data: input stores the parameter name, const stores the constant encoding, intrinsic stores the math function name, call stores the callee symbol and return index.
    """

    op: str
    kind: str
    args: tuple[int, ...] = ()
    data: tuple = ()


@dataclass(frozen=True)
class MathFunction:
    """One pure math function in MIR: a formal parameter table and a sequential SSA node sequence.

    Attributes:
        name: Function symbol name, unique within the program.
        label: The function name in the source, used for display and module attributes.
        parameters: The parameter table in formal parameter order.
        nodes: The sequential SSA node sequence.
        returns: Non-empty tuple of return node indices; multiple indices mean multiple independent results.
    """

    name: str
    label: str
    parameters: tuple[Parameter, ...]
    nodes: tuple[MathNode, ...]
    returns: tuple[int, ...]


@dataclass(frozen=True)
class MathProgram:
    """MIR 0.1 top-level object: the entry function symbol and the full function table.

    Attributes:
        entry: The entry function's symbol name, which must appear in ``functions``.
        functions: Immutable tuple of all ``MathFunction`` objects.
        version: The MIR version number, currently 0.1.
    """

    entry: str
    functions: tuple[MathFunction, ...]
    version: str = "0.1"

    @property
    def function_map(self) -> dict[str, MathFunction]:
        """Mapping from function symbol names to ``MathFunction``."""
        return {f.name: f for f in self.functions}

    def validate(self) -> MathProgram:
        """Check the program's structural and type constraints and return self.

        Covers the version number, function table immutability and entry
        validity, parameter uniqueness and kinds, sequential SSA references,
        operation arities and type promotion, constant finiteness, call
        interface consistency, and acyclicity of the function graph.

        Returns:
            MathProgram: Self, once validation has passed.

        Raises:
            ValidationError: Any structural or type constraint is violated.
        """
        if self.version != "0.1":
            raise ValidationError("unknown MIR version")
        if type(self.functions) is not tuple:
            raise ValidationError("the MIR function table must be immutable")
        functions = self.function_map
        if len(functions) != len(self.functions) or self.entry not in functions:
            raise ValidationError("invalid MIR entry or function table")
        for f in self.functions:
            name(f.name)
            if not isinstance(f.label, str) or any(
                type(v) is not tuple for v in (f.parameters, f.nodes, f.returns)
            ):
                raise ValidationError("MIR function records must have a string label and immutable arrays")
            params = {p.name: p for p in f.parameters}
            if len(params) != len(f.parameters):
                raise ValidationError("duplicate MIR parameter names")
            for p in f.parameters:
                name(p.name)
                if type(p.width) is not int or (p.kind != "index" and p.width != 0):
                    raise ValidationError("width must be 0 for non-index MIR parameters")
                if p.kind not in {"real", "complex", "bool", "index"} or (
                    p.kind == "index" and not 1 <= p.width <= 64
                ):
                    raise ValidationError("invalid MIR parameter kind")
            if not f.returns or any(
                type(i) is not int or not 0 <= i < len(f.nodes) for i in f.returns
            ):
                raise ValidationError("invalid MIR return values")
            for i, node in enumerate(f.nodes):
                if type(node.args) is not tuple or type(node.data) is not tuple:
                    raise ValidationError("MIR nodes must be immutable")
                if node.kind not in {"real", "complex", "bool"}:
                    raise ValidationError("invalid MIR value kind")
                if any(type(a) is not int or not 0 <= a < i for a in node.args):
                    raise ValidationError("MIR is not a sequential SSA graph")
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
                    raise ValidationError("invalid MIR operation arity")
                if node.op not in {*arities, "intrinsic", "call"}:
                    raise ValidationError("unknown MIR operation")
                if node.op == "input" and (len(node.data) != 1 or node.data[0] not in params):
                    raise ValidationError("invalid MIR input reference")
                if node.op == "const":
                    if len(node.data) != (2 if node.kind == "complex" else 1):
                        raise ValidationError("invalid MIR constant encoding")
                    import math

                    if any(
                        type(x) not in (bool, int, float) or not math.isfinite(x) for x in node.data
                    ):
                        raise ValidationError("MIR only accepts finite numeric constants")
                if node.op == "call":
                    if len(node.data) != 2 or node.data[0] not in functions:
                        raise ValidationError("invalid MIR call target")
                    child = functions[node.data[0]]
                    if len(node.args) != len(child.parameters) or not 0 <= node.data[1] < len(
                        child.returns
                    ):
                        raise ValidationError("invalid MIR call layout")
                    if node.kind != child.nodes[child.returns[node.data[1]]].kind:
                        raise ValidationError("invalid MIR call return kind")
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
                        raise ValidationError("unknown MIR math intrinsic")
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
                    raise ValidationError("MIR operation type or data constraint violated: " + node.op)
        active: set[str] = set()
        visited: set[str] = set()

        def visit(key: str) -> None:
            """Traverse the call graph depth-first; a back edge means recursion."""
            if key in active:
                raise ValidationError("recursion is not supported in MIR")
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

    def dumps(self) -> str:
        """Validate first, then serialize the program into canonical JSON text.

        Returns:
            str: JSON text with sorted keys, two-space indentation and a trailing newline; non-ASCII characters are kept verbatim.

        Raises:
            ValidationError: The program failed validation.
        """
        self.validate()
        return (
            json.dumps(asdict(self), ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
            + "\n"
        )

    @classmethod
    def loads(cls, text: str) -> MathProgram:
        """Deserialize and validate an MIR program from JSON text.

        Args:
            text: JSON text produced by ``MathProgram.dumps``.

        Returns:
            MathProgram: The validated program object.

        Raises:
            ValidationError: The text is not valid MIR JSON, or failed validation.
        """
        try:
            raw = json.loads(text)
            if set(raw) != {"entry", "functions", "version"}:
                raise ValidationError("invalid MIR top-level fields")
            functions: list[MathFunction] = []
            for f in raw["functions"]:
                if set(f) != {"name", "label", "parameters", "nodes", "returns"} or any(
                    set(n) != {"op", "kind", "args", "data"} for n in f["nodes"]
                ):
                    raise ValidationError("MIR function or node contains unknown fields")
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
            raise ValidationError(f"invalid MIR JSON: {exc}") from exc
