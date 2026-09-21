"生成阶段可检查的输入契约；结构/能力与数学真实性分别处理。"

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

from pyqecclang.infrastructure.ir import Register, Resource, ValidationError


@dataclass(frozen=True)
class OracleCapabilities:
    """oracle 支持的逆与受控调用能力，取值由依赖图保守合取推导。

    Attributes:
        adjoint: 是否支持逆调用。
        controlled: 是否支持受控调用。
    """

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
        """把快照转成可存 JSON 的普通字典。"""
        return asdict(self)


class OracleView:
    """包装类共享只读查询；宽度始终以 RIR 寄存器签名为准。"""

    @property
    def spec(self) -> OracleSpec:
        """该视图当前的 ``OracleSpec`` 快照；每次访问重新推导。"""
        return self.describe()

    def describe(self):
        """从 RIR 模块签名、属性与依赖图推导该视图的 ``OracleSpec``。"""
        return _operation_spec(self.operation, self.oracle_kind)

    @property
    def type(self) -> str:
        """描述性范式角色名，如 ``block_encoding``；不作为唯一类型或分派依据。"""
        return self.spec.type

    @property
    def main_qubit(self) -> int | None:
        """主寄存器（``target``/``address``/``column``）的位数；无法识别时为 ``None``。"""
        return self.spec.main_qubit

    @property
    def anc_qubit(self) -> int | None:
        """公开辅助寄存器（``signal``/``work``）的位数；操作集合没有统一宽度时为 ``None``。"""
        return self.spec.anc_qubit

    @property
    def capabilities(self) -> OracleCapabilities:
        """从整个依赖图推导的有效调用能力。"""
        return self.spec.capabilities


@dataclass(frozen=True)
class ContractIssue:
    """一次契约检查发现的具体问题。

    Attributes:
        code: 问题类别码，如 ``INPUT_PROTOCOL``、``INPUT_WIDTH``。
        path: 问题所在的输入路径，如 ``protocol.input.field``。
        expected: 契约期望的值。
        actual: 实际收到的值。
        message: 面向用户的问题描述。
    """

    code: str
    path: str
    expected: object
    actual: object
    message: str


class ContractError(ValidationError):
    """一次抛出全部已发现契约问题的异常，仍是 ``ValidationError`` 的子类。

    Args:
        issues: ``ContractIssue`` 可迭代对象，按 ``code path: message``
            逐条拼进异常消息。

    Attributes:
        issues: 结构化问题的不可变元组。
    """

    def __init__(self, issues):
        self.issues = tuple(issues)
        super().__init__("; ".join(f"{i.code} {i.path}: {i.message}" for i in self.issues))


def fail(code, path, expected, actual, message):
    """立即以单个问题抛出 ``ContractError``。

    Args:
        code: 问题类别码。
        path: 问题所在的输入路径。
        expected: 契约期望的值。
        actual: 实际收到的值。
        message: 面向用户的问题描述。

    Raises:
        ContractError: 总是抛出，仅携带这一个问题。
    """
    raise ContractError((ContractIssue(code, path, expected, actual, message),))


def require_instance(value, cls, path):
    """要求 ``value`` 是 ``cls`` 的实例，否则抛出 ``ContractError``。

    Args:
        value: 待检查的值。
        cls: 期望的类型。
        path: 问题上报时使用的输入路径。

    Raises:
        ContractError: 实例检查失败时以 ``INPUT_TYPE`` 上报。
    """
    if not isinstance(value, cls):
        fail("INPUT_TYPE", path, cls.__name__, type(value).__name__, f"需要 {cls.__name__}")


def finite_real(value, path, *, minimum=None, strict=False):
    """要求 ``value`` 是有限实数并可选下界；``bool`` 不被接受。

    Args:
        value: 待检查的数值。
        path: 问题上报时使用的输入路径。
        minimum: 允许的最小值；省略时不检查下界。
        strict: 为 ``True`` 时要求严格大于 ``minimum``，否则允许相等。

    Raises:
        ContractError: 值不是有限实数或低于下界时以 ``CONFIG_VALUE`` 上报。
    """
    if type(value) not in (int, float) or not math.isfinite(value):
        fail("CONFIG_VALUE", path, "finite real", repr(value), "需要有限实数，不能是 bool")
    if minimum is not None and (value <= minimum if strict else value < minimum):
        fail("CONFIG_VALUE", path, f"{'>' if strict else '>='}{minimum}", value, "超出允许范围")


def positive_integer(value, path, *, minimum=1, maximum=None):
    """要求 ``value`` 是范围内的整数；``bool`` 不被接受。

    Args:
        value: 待检查的整数。
        path: 问题上报时使用的输入路径。
        minimum: 允许的最小值，默认为 1。
        maximum: 允许的最大值；省略时不检查上界。

    Raises:
        ContractError: 值不是整数或越界时以 ``CONFIG_VALUE`` 上报。
    """
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        fail("CONFIG_VALUE", path, f"integer {minimum}..{maximum}", repr(value), "整数配置超出范围")


def validate_signature(operation, fields, path):
    """校验操作的公开寄存器恰好是约定的 ``bits`` 寄存器集合。

    Args:
        operation: 待检查的 ``Operation``。
        fields: 约定的寄存器名集合。
        path: 问题上报时使用的输入路径。

    Raises:
        ContractError: 输入不是 ``Operation``，或寄存器名与 ``fields``
            不一致、存在非 ``bits`` 寄存器时抛出。
        ValidationError: ``operation`` 的程序合并校验失败时透传。
    """
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
    """取得任意输入对象的 ``OracleSpec`` 描述快照。

    ``Operation`` 从模块属性与依赖图直接推导；带可调用 ``describe`` 方法
    的对象调用该方法并要求返回 ``OracleSpec``；其余宿主对象只报告类型名
    与可选的 ``main_qubit``/``anc_qubit`` 属性，不查询任何全局类型表。

    Args:
        value: ``Operation``、实现 ``describe`` 的对象或任意宿主对象。

    Returns:
        OracleSpec: 输入的只读描述快照。

    Raises:
        ContractError: ``describe`` 返回值不是 ``OracleSpec`` 时以
            ``INPUT_TYPE`` 上报。
        ValidationError: ``Operation`` 分支的程序合并校验失败时透传。
    """
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
    """协议对单个命名输入的需求声明。

    ``protocols`` 列出可接受的协议类（满足其一即可），可选的 ``adapter``
    把原始值转换成具体视图后再检查布局、能力与承诺。

    Attributes:
        name: 输入名，在所属 ``ProtocolContract`` 内唯一。
        protocols: 可接受的协议类元组。
        adapter: 可选适配函数，如 ``as_block_encoding``；省略时直接检查原值。
        main_qubit: 期望的主寄存器位数；``None`` 表示不检查。
        anc_qubit: 期望的辅助寄存器位数；``None`` 表示不检查。
        alpha: 期望的 BE 归一化常数；``None`` 表示不检查。
        adjoint: 要求输入支持逆调用。
        controlled: 要求输入支持受控调用。
        zero_input: 要求输入显式承诺零输入制备。
        clean_work: 要求输入显式承诺工作区复净。

    Raises:
        ValidationError: 名字为空、协议元组无效、adapter 不可调用或字段
            取值非法时在构造时抛出。
    """

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
        """转成可存 JSON 的字典；协议记录类名，adapter 记录函数名。"""
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
        """检查候选值是否满足该需求，返回描述与问题而不直接抛出。

        Args:
            value: 待检查的输入值。
            prefix: 报告路径前缀，通常为所属协议的名字。

        Returns:
            tuple: ``(spec, issues)``。协议不匹配或取得描述失败时
            ``spec`` 为 ``None``；布局、能力与承诺问题逐项追加进
            ``issues``，不阻断后续检查。

        Raises:
            ContractError: 仅当 adapter 抛出 ``ValidationError`` 以外的
                异常时透传。
        """
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
    """一次契约检查的结果报告，只覆盖结构与声明，不证明数学性质。

    Attributes:
        protocol: 被检查协议的名字。
        inputs: ``(输入名, OracleSpec)`` 元组，记录成功取得的输入描述。
        issues: 检查发现的全部问题；为空表示通过。
        assumptions: 协议声明的假设，原样保留。
        adapters: 检查涉及的适配器名；默认为空。
    """

    protocol: str
    inputs: tuple[tuple[str, OracleSpec], ...]
    issues: tuple[ContractIssue, ...] = ()
    assumptions: tuple[str, ...] = ()
    adapters: tuple[str, ...] = ()

    @property
    def ok(self):
        """没有发现任何问题时为 ``True``。"""
        return not self.issues

    def require(self):
        """通过时返回自身，否则一次抛出全部问题。

        Returns:
            ContractReport: 自身，便于链式调用。

        Raises:
            ContractError: 存在任何未解决的问题时抛出。
        """
        if self.issues:
            raise ContractError(self.issues)
        return self

    def to_dict(self):
        """转成可存 JSON 的字典，附加顶层 ``ok`` 字段。"""
        return {"ok": self.ok, **asdict(self)}


@dataclass(frozen=True)
class ProtocolContract:
    """算法对全部输入的契约声明，可整体检查并生成聚合报告。

    Attributes:
        name: 协议名。
        inputs: 按名字区分的 ``InputRequirement`` 元组。
        same_width: 要求 ``main_qubit`` 相同的输入名二元组。
        assumptions: 原样透传到报告中的假设声明。

    Raises:
        ValidationError: 名字为空、输入名重复、inputs 元素非法或同宽约束
            引用不存在的输入时在构造时抛出。
    """

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
        """按名字逐项检查输入并生成聚合报告。

        只调用输入对象的方法取得描述，不运行算法内核或量子模拟。

        Args:
            **values: 以需求名为键的输入值；多余或缺失都记为问题。

        Returns:
            ContractReport: 聚合所有输入的检查结果。
        """
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
        """转成可存 JSON 的字典；各输入递归使用自身 ``to_dict``。"""
        return {
            "name": self.name,
            "inputs": [r.to_dict() for r in self.inputs],
            "same_width": self.same_width,
            "assumptions": self.assumptions,
        }
