"生成阶段可检查的输入契约；结构/能力与数学真实性分别处理。"

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

from pyqecclang.infrastructure.ir import Register, Resource, ValidationError


@dataclass(frozen=True)
class OracleCapabilities:
    adjoint: bool
    controlled: bool


@dataclass(frozen=True)
class OracleSpec:
    """从实际操作读取的快照；CKS 是操作集合，没有统一 anc 总宽度。"""

    type: str
    main_qubit: int | None
    anc_qubit: int | None
    capabilities: OracleCapabilities
    implementation: str
    registers: tuple[Register, ...] = ()
    resources: tuple[Resource, ...] = ()
    alpha: float | None = None
    parameters: tuple[tuple[str, str | int | float | bool], ...] = ()
    components: tuple[tuple[str, OracleSpec], ...] = ()

    def to_dict(self):
        return asdict(self)


class OracleView:
    """包装类共享只读查询；宽度始终以 RIR 寄存器签名为准。"""

    @property
    def spec(self) -> OracleSpec:
        return self.describe()

    def describe(self):
        return _operation_spec(self.operation, self.oracle_kind)

    @property
    def type(self) -> str:
        return self.spec.type

    @property
    def main_qubit(self) -> int | None:
        return self.spec.main_qubit

    @property
    def anc_qubit(self) -> int | None:
        return self.spec.anc_qubit

    @property
    def capabilities(self) -> OracleCapabilities:
        return self.spec.capabilities


@dataclass(frozen=True)
class ContractIssue:
    code: str
    path: str
    expected: object
    actual: object
    message: str


class ContractError(ValidationError):
    def __init__(self, issues):
        self.issues = tuple(issues)
        super().__init__("; ".join(f"{i.code} {i.path}: {i.message}" for i in self.issues))


def fail(code, path, expected, actual, message):
    raise ContractError((ContractIssue(code, path, expected, actual, message),))


def require_instance(value, cls, path):
    if not isinstance(value, cls):
        fail("INPUT_TYPE", path, cls.__name__, type(value).__name__, f"需要 {cls.__name__}")


def finite_real(value, path, *, minimum=None, strict=False):
    if type(value) not in (int, float) or not math.isfinite(value):
        fail("CONFIG_VALUE", path, "finite real", repr(value), "需要有限实数，不能是 bool")
    if minimum is not None and (value <= minimum if strict else value < minimum):
        fail("CONFIG_VALUE", path, f"{'>' if strict else '>='}{minimum}", value, "超出允许范围")


def positive_integer(value, path, *, minimum=1, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        fail("CONFIG_VALUE", path, f"integer {minimum}..{maximum}", repr(value), "整数配置超出范围")


def validate_signature(operation, fields, path):
    from pyqecclang.infrastructure.builder import Operation

    require_instance(operation, Operation, path)
    operation.program()
    actual = {r.name: r.type for r in operation.module.registers}
    if set(actual) != set(fields) or any(t.kind != "bits" for t in actual.values()):
        fail("ORACLE_SIGNATURE", path, list(fields), list(actual), "需要约定的 bits 寄存器签名")


def _operation_spec(operation, kind=None):
    from pyqecclang.infrastructure.linking import capabilities, unresolved

    program = operation.program()
    module = operation.module
    attrs = dict(module.attributes)
    caps = capabilities(program)
    regs = {r.name: r.type.width for r in module.registers}
    kind = kind or attrs.get("oracle_paradigm", "unitary")
    main = regs.get("target", regs.get("address", regs.get("column")))
    anc = regs.get("signal", regs.get("work"))
    return OracleSpec(
        str(kind),
        main,
        anc,
        OracleCapabilities(caps["supports_adjoint"], caps["supports_controlled"]),
        "open" if unresolved(program) else "closed",
        module.registers,
        module.resources,
        attrs.get("be_alpha"),
        module.attributes,
    )


def describe_oracle(value) -> OracleSpec:
    from pyqecclang.infrastructure.builder import Operation

    if isinstance(value, Operation):
        return _operation_spec(value)
    describe = getattr(value, "describe", None)
    if callable(describe):
        result = describe()
        require_instance(result, OracleSpec, "describe.output")
        return result
    # 自定义算法对象不需要向语言注册；没有描述方法时只给出宿主类型。
    return OracleSpec(
        type(value).__name__,
        getattr(value, "main_qubit", None),
        getattr(value, "anc_qubit", None),
        OracleCapabilities(False, False),
        "host",
    )


def requires(value, protocol, *, path="input"):
    """算法侧的结构检查；协议由调用者定义，不查询任何全局类型表。"""
    if not isinstance(value, protocol):
        fail(
            "INPUT_PROTOCOL",
            path,
            protocol.__name__,
            type(value).__name__,
            "需要满足 " + protocol.__name__,
        )
    return value


@dataclass(frozen=True)
class InputRequirement:
    name: str
    protocols: tuple[type, ...]
    adapter: object = None
    main_qubit: int | None = None
    anc_qubit: int | None = None
    alpha: float | None = None
    adjoint: bool = False
    controlled: bool = False
    zero_input: bool = False
    clean_work: bool = False

    def __post_init__(self):
        if (
            not self.name
            or type(self.protocols) is not tuple
            or not self.protocols
            or any(not isinstance(t, type) for t in self.protocols)
        ):
            raise ValidationError("InputRequirement 需要名字与非空的协议类元组")
        if self.adapter is not None and not callable(self.adapter):
            raise ValidationError("InputRequirement.adapter 必须可调用")
        for key in ("main_qubit", "anc_qubit"):
            value = getattr(self, key)
            if value is not None:
                positive_integer(value, key, minimum=0)
        if self.alpha is not None:
            finite_real(self.alpha, "alpha", minimum=0, strict=True)
        for key in ("adjoint", "controlled", "zero_input", "clean_work"):
            if type(getattr(self, key)) is not bool:
                raise ValidationError(f"InputRequirement.{key} 必须是 bool")

    def to_dict(self):
        return {
            "name": self.name,
            "protocols": [p.__name__ for p in self.protocols],
            "adapter": getattr(self.adapter, "__name__", None),
            **{
                key: getattr(self, key)
                for key in (
                    "main_qubit",
                    "anc_qubit",
                    "alpha",
                    "adjoint",
                    "controlled",
                    "zero_input",
                    "clean_work",
                )
            },
        }

    def inspect(self, value, *, prefix=""):
        path = f"{prefix}.{self.name}" if prefix else self.name
        expected_names = tuple(p.__name__ for p in self.protocols)
        if not any(isinstance(value, p) for p in self.protocols):
            return None, (
                ContractIssue(
                    "INPUT_PROTOCOL",
                    path,
                    expected_names,
                    type(value).__name__,
                    f"需要满足至少一个协议：{expected_names}",
                ),
            )
        try:
            value = self.adapter(value) if self.adapter is not None else value
            spec = describe_oracle(value)
        except ValidationError as exc:
            return None, (
                ContractIssue(
                    "INPUT_ADAPTER", path, expected_names, type(value).__name__, str(exc)
                ),
            )
        issues = []

        def issue(code, field, expected, actual, message):
            issues.append(ContractIssue(code, path + "." + field, expected, actual, message))

        for field in ("main_qubit", "anc_qubit", "alpha"):
            expected, actual = getattr(self, field), getattr(spec, field)
            if expected is not None and expected != actual:
                issue(
                    "INPUT_LAYOUT" if field != "alpha" else "INPUT_ALPHA",
                    field,
                    expected,
                    actual,
                    f"需要 {expected}，收到 {actual}",
                )
        for field in ("adjoint", "controlled"):
            if getattr(self, field) and not getattr(spec.capabilities, field):
                issue("INPUT_CAPABILITY", field, True, False, f"缺少 supports_{field}")
        attrs = dict(spec.parameters)
        for field in ("zero_input", "clean_work"):
            if getattr(self, field) and attrs.get(field) is not True:
                issue("INPUT_PROMISE", field, True, attrs.get(field), f"需要显式的 {field} 契约")
        return spec, tuple(issues)


@dataclass(frozen=True)
class ContractReport:
    protocol: str
    inputs: tuple[tuple[str, OracleSpec], ...]
    issues: tuple[ContractIssue, ...] = ()
    assumptions: tuple[str, ...] = ()
    adapters: tuple[str, ...] = ()

    @property
    def ok(self):
        return not self.issues

    def require(self):
        if self.issues:
            raise ContractError(self.issues)
        return self

    def to_dict(self):
        return {"ok": self.ok, **asdict(self)}


@dataclass(frozen=True)
class ProtocolContract:
    name: str
    inputs: tuple[InputRequirement, ...]
    same_width: tuple[tuple[str, str], ...] = ()
    assumptions: tuple[str, ...] = ()

    def __post_init__(self):
        if type(self.inputs) is not tuple or any(
            not isinstance(r, InputRequirement) for r in self.inputs
        ):
            raise ValidationError("ProtocolContract.inputs 需要不可变的 InputRequirement 元组")
        names = [r.name for r in self.inputs]
        if not self.name or len(names) != len(set(names)):
            raise ValidationError("protocol 名称无效或输入名重复")
        if any(len(pair) != 2 or any(n not in names for n in pair) for pair in self.same_width):
            raise ValidationError("同宽约束引用了不存在的输入")

    def check(self, **values) -> ContractReport:
        specs, issues = {}, []
        expected = {r.name for r in self.inputs}
        for key in sorted(values.keys() - expected):
            issues.append(
                ContractIssue("INPUT_UNEXPECTED", f"{self.name}.{key}", None, key, "多余输入")
            )
        for req in self.inputs:
            if req.name not in values:
                issues.append(
                    ContractIssue(
                        "INPUT_MISSING",
                        f"{self.name}.{req.name}",
                        tuple(p.__name__ for p in req.protocols),
                        None,
                        "缺少输入",
                    )
                )
                continue
            spec, found = req.inspect(values[req.name], prefix=self.name)
            if spec is not None:
                specs[req.name] = spec
            issues.extend(found)
        for left, right in self.same_width:
            if (
                left in specs
                and right in specs
                and specs[left].main_qubit != specs[right].main_qubit
            ):
                issues.append(
                    ContractIssue(
                        "INPUT_WIDTH",
                        f"{self.name}.{right}.main_qubit",
                        specs[left].main_qubit,
                        specs[right].main_qubit,
                        f"需要与 {left} 的目标宽度一致",
                    )
                )
        return ContractReport(self.name, tuple(specs.items()), tuple(issues), self.assumptions)

    def to_dict(self):
        return {
            "name": self.name,
            "inputs": [r.to_dict() for r in self.inputs],
            "same_width": self.same_width,
            "assumptions": self.assumptions,
        }
