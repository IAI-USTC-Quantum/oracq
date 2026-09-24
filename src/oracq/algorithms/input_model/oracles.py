"Oracle 范式目录及普通 gate/QRAM 实现。数学正确性认证不属于此模块。"

from __future__ import annotations

import cmath
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import cast

from oracq.algorithms.input_model.contracts import (
    OracleSpec,
    OracleView,
    fail,
    positive_integer,
    require_instance,
    validate_signature,
)
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import (
    QRAM,
    Bits,
    Module,
    Ref,
    Register,
    RegType,
    Resource,
    ValidationError,
    fuse,
)

PARADIGMS = {
    "unitary": "完整寄存器空间上的量子操作接口",
    "block_encoding": "零信号投影角块与 be_alpha 的组合约定",
    "database_xor": "|address,data> -> |address,data XOR memory[address]>",
    "state_prep_isometry": "从零态子空间制备目标态；逆和控制要求可逆扩张",
    "sparse_location_inplace": "CKS P_A(1)：|column,index> -> |column,nu(column,index)>",
    "sparse_entry_xor": "CKS P_A(2)：|row,column,data> -> |row,column,data XOR A[row,column]>",
    "phase_oracle": "按已定义谓词给目标基态施加相位",
    "reversible_function": "输入保留、输出可逆更新的领域计算接口",
    "algorithm_stage": "具有明确寄存器接口的未完成算法阶段",
}
"""内置 oracle 角色约定；应用可使用自定义标识符，无需修改此目录。"""


def declare(
    name: str,
    registers: Mapping[str, RegType],
    *,
    paradigm: str = "unitary",
    resources: Mapping[str, QRAM] | None = None,
    attributes: Mapping[str, str | int | float | bool] | None = None,
    supports_adjoint: bool = True,
    supports_controlled: bool = True,
) -> Operation:
    """声明一个开放的 oracle 操作。

    生成体为空、``implementation_status`` 标记为 ``unresolved`` 的 ``Operation``；
    结构在生成期校验，实现可经 ``annotate`` 或绑定闭合。

    Args:
        name: 模块名。
        registers: 寄存器名到 ``Bits`` 等类型的有序映射。
        paradigm: 内置或应用定义的角色标识符；角色名本身不证明数学性质。
        resources: 资源名到 ``QRAM`` 等资源类型的映射。
        attributes: 并入模块的额外属性。
        supports_adjoint: 是否声明支持逆操作。
        supports_controlled: 是否声明支持受控操作。

    Returns:
        Operation: 未解析的开放声明。

    Raises:
        ValidationError: 角色名非法，或生成的模块未通过结构校验。
    """
    from oracq.infrastructure.validation import name as validate_name

    validate_name(paradigm)
    attrs = dict(attributes or {})
    attrs.update(
        oracle_paradigm=paradigm,
        supports_adjoint=supports_adjoint,
        supports_controlled=supports_controlled,
        implementation_status="unresolved",
    )
    operation = Operation(
        Module(
            name,
            tuple(Register(k, v) for k, v in registers.items()),
            tuple(Resource(k, v) for k, v in (resources or {}).items()),
            None,
            tuple(sorted(attrs.items())),
        )
    )
    operation.program()
    return operation


def annotate(
    operation: Operation,
    paradigm: str,
    **attributes: str | int | float | bool,
) -> Operation:
    """为已生成的实现操作登记范式与属性。

    把 ``implementation_status`` 置为 ``constructed`` 并写入 ``oracle_paradigm``；
    已被显式限制为 ``False`` 的 adjoint/controlled 能力不能经标注重新授予。范式为
    ``state_prep_isometry`` 时缺省补 ``zero_input`` 与 ``clean_work`` 承诺。

    Args:
        operation: 已生成的 ``Operation``。
        paradigm: 内置或应用定义的角色标识符。
        **attributes: 并入模块属性的其他键值对。

    Returns:
        Operation: 替换属性后的新 ``Operation``。

    Raises:
        ValidationError: 输入类型或范式名无效，或试图恢复已被限制的能力。
    """
    require_instance(operation, Operation, "annotate.operation")
    from oracq.infrastructure.validation import name as validate_name

    validate_name(paradigm)
    attrs = dict(operation.module.attributes)
    for cap in ("supports_adjoint", "supports_controlled"):
        if attrs.get(cap) is False and attributes.get(cap) is True:
            fail("INPUT_CAPABILITY", "annotate." + cap, False, True, "标注不能授予已被限制的能力")
    attrs.update(
        oracle_paradigm=paradigm,
        implementation_status="constructed",
        **attributes,
    )
    # 标注描述不会授予已被限制的变换能力。
    attrs.setdefault("supports_adjoint", True)
    attrs.setdefault("supports_controlled", True)
    if paradigm == "state_prep_isometry":
        attrs.setdefault("zero_input", True)
        attrs.setdefault("clean_work", True)
    return Operation(
        replace(operation.module, attributes=tuple(sorted(attrs.items()))), operation.dependencies
    )


def resources_for(*items: tuple[str, Operation]) -> dict[str, QRAM]:
    """汇总多个操作的资源并以 ``前缀__资源名`` 重命名。

    Args:
        *items: 形如 ``(prefix, operation)`` 的二元组，``operation`` 为带资源
            声明的 ``Operation``。

    Returns:
        dict: 键为 ``prefix__资源名``、值为资源类型的映射，供嵌套调用的
        ``Builder`` 声明资源。
    """
    return {f"{prefix}__{r.name}": r.type for prefix, op in items for r in op.module.resources}


def invoke(builder: Builder, operation: Operation, prefix: str = "", **arguments: Ref) -> None:
    """在构建器中以带前缀的资源映射调用一个操作。

    Args:
        builder: 当前 ``Builder``，其资源须已按同一前缀声明。
        operation: 被调用的 ``Operation``。
        prefix: 资源名前缀；空串表示资源名原样直通。
        **arguments: 被调模块寄存器名到实参视图的绑定。

    Raises:
        ValidationError: 寄存器或资源参数名与被调模块不匹配，或依赖定义冲突。
    """
    mapping = {
        r.name: f"{prefix}__{r.name}" if prefix else r.name for r in operation.module.resources
    }
    builder.call(operation, resources=mapping, **arguments)


@dataclass(frozen=True)
class XorDatabase(OracleView):
    """XOR 数据库视图：``|address,data> -> |address,data XOR memory[address]>``。

    经典查询表的可逆量子访问接口；包装的 ``Operation`` 恰含 ``address`` 与
    ``data`` 两个 bits 寄存器，XOR 语义在任意初值下自逆。

    Attributes:
        operation: 被包装的 ``Operation``。
    """

    oracle_kind = "database_xor"
    operation: Operation

    def xor_database(self) -> XorDatabase:
        """返回自身；实现 ``XorDatabaseProtocol`` 的视图适配方法。

        Returns:
            XorDatabase: 该视图自身。
        """
        return self

    def __post_init__(self) -> None:
        """校验包装操作恰有 ``address`` 与 ``data`` 两个 bits 寄存器。"""
        validate_signature(self.operation, ("address", "data"), "XorDatabase")
        regs = {r.name: r.type for r in self.operation.module.registers}
        if set(regs) != {"address", "data"} or any(r.kind != "bits" for r in regs.values()):
            raise ValidationError("XOR database 需要 address/data bits 接口")

    @property
    def address_width(self) -> int:
        """``address`` 寄存器的位宽。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "address")

    @property
    def data_width(self) -> int:
        """``data`` 寄存器的位宽。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "data")


@dataclass(frozen=True)
class StatePreparation(OracleView):
    """零输入态制备视图：从零态子空间出发的等距。

    约定 ``V|0,0> = |psi,0>``，即 ``work`` 寄存器复净；逆与受控操作要求存在
    可逆扩张。包装的 ``Operation`` 恰含 ``target`` 与 ``work`` 两个 bits 寄存器。

    Attributes:
        operation: 被包装的 ``Operation``。
    """

    oracle_kind = "state_prep_isometry"
    operation: Operation

    def state_preparation(self) -> StatePreparation:
        """返回自身；实现 ``StatePreparationProtocol`` 的视图适配方法。

        Returns:
            StatePreparation: 该视图自身。
        """
        return self

    def __post_init__(self) -> None:
        """校验包装操作恰有 ``target`` 与 ``work`` 两个 bits 寄存器。"""
        validate_signature(self.operation, ("target", "work"), "StatePreparation")
        regs = {r.name: r.type for r in self.operation.module.registers}
        if set(regs) != {"target", "work"} or any(r.kind != "bits" for r in regs.values()):
            raise ValidationError("StatePreparation 需要 target/work bits 接口")

    @classmethod
    def from_unitary(
        cls,
        operation: Operation,
        *,
        target: str | None = None,
        work: str | None = None,
        clean_work: bool = False,
    ) -> StatePreparation:
        """显式赋予 U|0> 初态角色；有工作寄存器时由调用方承诺复净。

        Args:
            operation: 完整酉 ``Operation``。
            target: 承担初态制备的寄存器名；缺省时把全部寄存器拼接为 target。
            work: 工作寄存器名；仅在显式指定 target 时可指定。
            clean_work: work 非零宽时必须为 True，承诺制备后 work 复净。

        Returns:
            StatePreparation: 以零输入角色包装该操作的制备视图。
        """
        require_instance(operation, Operation, "StatePreparation.from_unitary")
        if target is None:
            if work is not None:
                raise ValidationError("指定 work 时必须同时指定 target")
            operation.program()
            width = sum(r.type.width for r in operation.module.registers)
            b = Builder(
                _name("unitary_zero_input", operation),
                {"target": Bits(width), "work": Bits(0)},
                resources_for(("unitary", operation)),
            )
            cursor = 0
            arguments: dict[str, Ref] = {}
            for reg in operation.module.registers:
                arguments[reg.name] = b["target"][cursor : cursor + reg.type.width].reinterpret(
                    reg.type.kind
                )
                cursor += reg.type.width
            invoke(b, operation, "unitary", **arguments)
            return cls(
                annotate(b.finish(), "state_prep_isometry", zero_input=True, clean_work=True)
            )
        expected = (target,) if work is None else (target, work)
        if len(set(expected)) != len(expected):
            raise ValidationError("target/work 不能指向同一个寄存器")
        validate_signature(operation, expected, "StatePreparation.from_unitary")
        widths = {r.name: r.type.width for r in operation.module.registers}
        ww = 0 if work is None else widths[work]
        if ww and clean_work is not True:
            fail(
                "INPUT_PROMISE",
                "StatePreparation.from_unitary.clean_work",
                True,
                clean_work,
                "需要承诺 U|0,0> 中 work 复净",
            )
        b = Builder(
            _name("as_state_preparation", operation, target, work),
            {"target": Bits(widths[target]), "work": Bits(ww)},
            resources_for(("unitary", operation)),
        )
        arguments = {target: b["target"]}
        if work is not None:
            arguments[work] = b["work"]
        invoke(b, operation, "unitary", **arguments)
        return cls(annotate(b.finish(), "state_prep_isometry", zero_input=True, clean_work=True))

    @property
    def width(self) -> int:
        """``target`` 寄存器的位宽。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def work_width(self) -> int:
        """``work`` 寄存器的位宽。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "work")


@dataclass(frozen=True)
class StateOracle(OracleView):
    """态输出视图：成功子空间约定为 ``signal`` 全零。

    求解内核等算法阶段的输出载体；包装的 ``Operation`` 恰含 ``target`` 与
    ``signal`` 两个 bits 寄存器，读取 ``signal == 0`` 分支得到目标态。

    Attributes:
        operation: 被包装的 ``Operation``。
    """

    oracle_kind = "state_oracle"
    operation: Operation

    def state_oracle(self) -> StateOracle:
        """返回自身；实现 ``StateOracleProtocol`` 的视图适配方法。

        Returns:
            StateOracle: 该视图自身。
        """
        return self

    def __post_init__(self) -> None:
        """校验包装操作恰有 ``target`` 与 ``signal`` 寄存器。"""
        validate_signature(self.operation, ("target", "signal"), "StateOracle")

    @property
    def width(self) -> int:
        """``target`` 寄存器的位宽。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def signal_qubits(self) -> int:
        """``signal`` 寄存器的位宽。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "signal")


@dataclass(frozen=True)
class SparseAccess(OracleView):
    """CKS 稀疏访问的束视图：位置与元素两个查询操作的组合，不伪装成总酉。

    Attributes:
        location: 位置操作 ``P_A(1)``，原地置换 ``index``，签名为 ``column``、
            ``index``、``work``。
        entry: 元素操作 ``P_A(2)``，任意行列可查，签名为 ``row``、``column``、
            ``data``。
        width: 矩阵维度位数，范围为 1..64。
        value_width: 元素值字宽，范围为 1..64。
        sparsity: 每列非零元上界，范围为 1..2**width。
    """

    location: Operation
    entry: Operation
    width: int
    value_width: int
    sparsity: int

    def sparse_access(self) -> SparseAccess:
        """返回自身；实现 ``CKSSparseProtocol`` 的视图适配方法。

        Returns:
            SparseAccess: 该视图自身。
        """
        return self

    def describe(self) -> OracleSpec:
        """汇总两个组件得到束整体的只读描述。

        类型为 ``cks_sparse`` 且 ``anc_qubit`` 恒为 ``None``；adjoint 与
        controlled 能力取两组件的与，全部组件闭合才标记为 ``closed``，并附
        ``sparsity`` 与 ``value_width`` 参数。

        Returns:
            OracleSpec: 束整体的描述，组件以 ``position``、``entry`` 命名。
        """
        from oracq.algorithms.input_model.contracts import (
            OracleCapabilities,
            OracleSpec,
            describe_oracle,
        )

        components = (
            ("position", describe_oracle(self.location)),
            ("entry", describe_oracle(self.entry)),
        )
        return OracleSpec(
            "cks_sparse",
            self.width,
            None,
            OracleCapabilities(
                all(s.capabilities.adjoint for _, s in components),
                all(s.capabilities.controlled for _, s in components),
            ),
            "closed" if all(s.implementation == "closed" for _, s in components) else "open",
            parameters=(("sparsity", self.sparsity), ("value_width", self.value_width)),
            components=components,
        )

    def __post_init__(self) -> None:
        """校验位宽、稀疏度取值与两个组件操作的寄存器签名。"""
        positive_integer(self.width, "SparseAccess.width", maximum=64)
        positive_integer(self.value_width, "SparseAccess.value_width", maximum=64)
        positive_integer(self.sparsity, "SparseAccess.sparsity", maximum=1 << self.width)
        validate_signature(self.location, ("column", "index", "work"), "SparseAccess.position")
        validate_signature(self.entry, ("row", "column", "data"), "SparseAccess.entry")
        if (
            not 1 <= self.width <= 64
            or not 1 <= self.value_width <= 64
            or not 1 <= self.sparsity <= 1 << self.width
        ):
            raise ValidationError("SparseAccess 位宽/稀疏度无效")
        location = {r.name: r.type for r in self.location.module.registers}
        entry = {r.name: r.type for r in self.entry.module.registers}
        if (
            set(location) != {"column", "index", "work"}
            or location["column"] != Bits(self.width)
            or location["index"] != Bits(self.width)
        ):
            raise ValidationError("SparseAccess 需要原地 column/index/work 位置接口")
        if entry != {
            "row": Bits(self.width),
            "column": Bits(self.width),
            "data": Bits(self.value_width),
        }:
            raise ValidationError("SparseAccess 需要任意 row/column 的元素 XOR 接口")


def abstract_database(name: str, address_width: int, data_width: int) -> XorDatabase:
    """返回未绑定实现的开放 XOR 数据库声明。

    Args:
        name: 声明名。
        address_width: ``address`` 寄存器位宽。
        data_width: ``data`` 寄存器位宽。

    Returns:
        XorDatabase: 体为空的开放声明，经实现绑定后闭合。
    """
    return XorDatabase(
        declare(
            name,
            {"address": Bits(address_width), "data": Bits(data_width)},
            paradigm="database_xor",
        )
    )


def abstract_state_prep(
    name: str,
    width: int,
    work_width: int = 0,
    *,
    reversible: bool = True,
) -> StatePreparation:
    """返回未绑定实现的开放态制备声明。

    Args:
        name: 声明名。
        width: ``target`` 寄存器位宽。
        work_width: ``work`` 寄存器位宽，缺省为零。
        reversible: 为 ``False`` 时声明不支持 adjoint 与 controlled。

    Returns:
        StatePreparation: 带 ``zero_input`` 与 ``clean_work`` 承诺的开放声明。
    """
    return StatePreparation(
        declare(
            name,
            {"target": Bits(width), "work": Bits(work_width)},
            paradigm="state_prep_isometry",
            attributes={"zero_input": True, "clean_work": True},
            supports_adjoint=reversible,
            supports_controlled=reversible,
        )
    )


def abstract_block_encoding(
    name: str,
    width: int,
    signal_width: int,
    alpha: float,
) -> BlockEncoding:
    """返回未绑定实现的开放块编码声明。

    Args:
        name: 声明名。
        width: ``target`` 寄存器位宽。
        signal_width: ``signal`` 寄存器位宽。
        alpha: 归一化常数，写入 ``be_alpha`` 属性。

    Returns:
        BlockEncoding: 体为空的开放声明。
    """
    return BlockEncoding(
        declare(
            name,
            {"target": Bits(width), "signal": Bits(signal_width)},
            paradigm="block_encoding",
            attributes={"be_alpha": alpha},
        )
    )


def abstract_sparse_access(
    name: str,
    width: int,
    value_width: int,
    sparsity: int,
    work_width: int | None = None,
) -> SparseAccess:
    """返回未绑定实现的开放稀疏访问束。

    生成 ``name_position`` 与 ``name_entry`` 两个开放操作；位置操作带
    ``sparsity`` 与完整置换扩张属性。

    Args:
        name: 声明名前缀。
        width: 矩阵维度位数。
        value_width: 元素值字宽。
        sparsity: 每列非零元上界，须不超过 ``2**width``。
        work_width: 位置操作的 ``work`` 位宽，缺省取 ``width``。

    Returns:
        SparseAccess: 位置与元素两个开放声明组成的束。

    Raises:
        ValidationError: ``sparsity`` 超出维度。
    """
    if not 1 <= sparsity <= 1 << width:
        raise ValidationError("sparsity 超出维度")
    work_width = width if work_width is None else work_width
    location = declare(
        name + "_position",
        {"column": Bits(width), "index": Bits(width), "work": Bits(work_width)},
        paradigm="sparse_location_inplace",
        attributes={"sparsity": sparsity, "full_permutation_extension": True},
    )
    entry = declare(
        name + "_entry",
        {"row": Bits(width), "column": Bits(width), "data": Bits(value_width)},
        paradigm="sparse_entry_xor",
    )
    return SparseAccess(location, entry, width, value_width, sparsity)


def gate_database(
    address_width: int,
    data_width: int,
    table: Mapping[int, int] | Sequence[int],
    *,
    name: str | None = None,
) -> XorDatabase:
    """以门级真值表实现 XOR 数据库。

    非零表字在地址控制下逐位施加 X 门；门数随表规模增长，适合小实例见证。

    Args:
        address_width: ``address`` 寄存器位宽。
        data_width: ``data`` 寄存器位宽。
        table: 地址到字的映射字典，或按地址枚举的字序列。
        name: 模块名；缺省按内容确定性生成。

    Returns:
        XorDatabase: 标注 ``implementation="gate_truth_table"`` 的门级实现。

    Raises:
        ValidationError: 表项不是整数，或地址、字超出位宽范围。
    """
    items = tuple(sorted(table.items())) if isinstance(table, dict) else tuple(enumerate(table))
    b = Builder(
        name or _name("xor_table", address_width, data_width, items),
        {"address": Bits(address_width), "data": Bits(data_width)},
    )
    for address, value in items:
        if (
            type(address) is not int
            or type(value) is not int
            or not (0 <= address < 1 << address_width and 0 <= value < 1 << data_width)
        ):
            raise ValidationError("查询表地址或字越界")
        if value:
            with b.control(b["address"], address):
                for bit in range(data_width):
                    if (value >> bit) & 1:
                        b.x(b["data"][bit])
    return XorDatabase(annotate(b.finish(), "database_xor", implementation="gate_truth_table"))


def qram_database(address_width: int, data_width: int, *, name: str | None = None) -> XorDatabase:
    """以 QRAM 资源实现 XOR 数据库。

    声明名为 ``table`` 的 ``QRAM(address_width, data_width)`` 资源并发出一条
    查询原语；数据表在执行期以 memory 形式提供，不进入 IR。

    Args:
        address_width: ``address`` 寄存器位宽。
        data_width: ``data`` 寄存器位宽。
        name: 模块名；缺省按位宽生成。

    Returns:
        XorDatabase: 标注 ``implementation="qram"`` 的资源版实现。
    """
    b = Builder(
        name or f"qram_xor_{address_width}_{data_width}",
        {"address": Bits(address_width), "data": Bits(data_width)},
        {"table": QRAM(address_width, data_width)},
    )
    b.qram("table", b["address"], b["data"])
    return XorDatabase(annotate(b.finish(), "database_xor", implementation="qram"))


def basis_state(width: int, value: int = 0, *, work_width: int = 0) -> StatePreparation:
    """以逐位 X 门制备计算基态 ``|value>``。

    Args:
        width: ``target`` 寄存器位宽。
        value: 目标基态的整数编号，范围为 ``0 .. 2**width - 1``。
        work_width: ``work`` 寄存器位宽，缺省为零。

    Returns:
        StatePreparation: 带 ``zero_input`` 承诺的门级制备。

    Raises:
        ValidationError: ``value`` 不是整数或越界。
    """
    if type(value) is not int or not 0 <= value < 1 << width:
        raise ValidationError("基态值越界")
    b = Builder(
        _name("basis", width, value, work_width), {"target": Bits(width), "work": Bits(work_width)}
    )
    for bit in range(width):
        if (value >> bit) & 1:
            b.x(b["target"][bit])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def uniform_state(width: int, *, work_width: int = 0) -> StatePreparation:
    """对整个 ``target`` 施加 Hadamard，制备均匀叠加态。

    Args:
        width: ``target`` 寄存器位宽。
        work_width: ``work`` 寄存器位宽，缺省为零。

    Returns:
        StatePreparation: 带 ``zero_input`` 承诺的门级制备。
    """
    b = Builder(
        _name("uniform", width, work_width), {"target": Bits(width), "work": Bits(work_width)}
    )
    b.h(b["target"])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def _state_angles(
    amplitudes: Iterable[complex],
) -> tuple[tuple[complex, ...], int, list[tuple[int, int, int, float]]]:
    """校验幅度向量并推导 Ry 旋转树各节点的旋转角。"""
    values = tuple(complex(v) for v in amplitudes)
    if not values or len(values) & (len(values) - 1):
        raise ValidationError("态向量长度必须是二的幂")
    if not all(math.isfinite(v.real) and math.isfinite(v.imag) for v in values):
        raise ValidationError("态幅度必须有限")
    if sum(abs(v) ** 2 for v in values) == 0:
        raise ValidationError("零向量不能制备为归一化态")
    n = (len(values) - 1).bit_length()
    nodes: list[tuple[int, int, int, float]] = []
    for depth in range(n):
        bit = n - depth - 1
        for prefix in range(1 << depth):
            start, size = prefix << (bit + 1), 1 << bit
            left = sum(abs(v) ** 2 for v in values[start : start + size])
            right = sum(abs(v) ** 2 for v in values[start + size : start + 2 * size])
            nodes.append((depth, prefix, bit, 2 * math.atan2(math.sqrt(right), math.sqrt(left))))
    return values, n, nodes


def gate_state_prep(
    amplitudes: Iterable[complex],
    *,
    work_width: int = 0,
    name: str | None = None,
) -> StatePreparation:
    """以多路复用 Ry 旋转树制备任意复幅度态。

    幅度按二叉树分解为受控 ``Ry`` 角度，非零相位分量经受控 ``global_phase``
    逐点补偿，因此幅度可为复数。

    Args:
        amplitudes: 长度为二的幂、范数非零的有限复幅度序列。
        work_width: ``work`` 寄存器位宽，缺省为零。
        name: 模块名；缺省按幅度内容确定性生成。

    Returns:
        StatePreparation: 标注 ``implementation="multiplexed_rotations"`` 的门级
        制备。

    Raises:
        ValidationError: 幅度向量非法，或仅有一个幅度而无目标位。
    """
    values, n, nodes = _state_angles(amplitudes)
    if n == 0:
        raise ValidationError("当前态制备需要至少一个目标位")
    b = Builder(
        name or _name("state", values, work_width), {"target": Bits(n), "work": Bits(work_width)}
    )
    for depth, prefix, bit, angle in nodes:
        if not depth:
            b.ry(b["target"][bit], angle)
        else:
            with b.control(b["target"][bit + 1 :], prefix):
                b.ry(b["target"][bit], angle)
    for index, value in enumerate(values):
        if value and cmath.phase(value):
            with b.control(b["target"], index):
                b.global_phase(cmath.phase(value))
    return StatePreparation(
        annotate(
            b.finish(),
            "state_prep_isometry",
            zero_input=True,
            implementation="multiplexed_rotations",
        )
    )


def qram_state_prep(width: int, angle_width: int = 8) -> StatePreparation:
    """以 QRAM 角度表实现旋转树态制备。

    逐层把 ``target`` 高位前缀换入地址并查询角表，按角度字的各个位受控施加
    Ry 旋转后复净地址，每层查询两次；角度分辨率 ``2*pi/2**angle_width`` 是量化
    误差的来源。

    Args:
        width: ``target`` 寄存器位宽。
        angle_width: 角度字的位宽。

    Returns:
        StatePreparation: 标注 ``implementation="qram_rotation_tree"`` 的资源版
        制备；``work`` 为 ``max(1, width) + angle_width`` 位，属性附 ``qram_queries``。
    """
    address_width = max(1, width)
    b = Builder(
        f"qram_state_{width}_{angle_width}",
        {"target": Bits(width), "work": Bits(address_width + angle_width)},
        {"angles": QRAM(address_width, angle_width)},
    )
    addr, angle = b["work"][:address_width], b["work"][address_width:]
    for depth in range(width):
        bit = width - depth - 1
        if depth:
            b.xor(b["target"][bit + 1 :], addr[:depth])
        offset = (1 << depth) - 1
        b.add_const(addr.reinterpret("uint"), offset)
        b.qram("angles", addr, angle)
        for k in range(angle_width):
            with b.control(angle[k]):
                b.ry(b["target"][bit], 2 * math.pi * (1 << k) / (1 << angle_width))
        b.qram("angles", addr, angle)
        b.add_const(addr.reinterpret("uint"), (-offset) % (1 << address_width))
        if depth:
            b.xor(b["target"][bit + 1 :], addr[:depth])
    return StatePreparation(
        annotate(
            b.finish(),
            "state_prep_isometry",
            zero_input=True,
            implementation="qram_rotation_tree",
            clean_work=True,
            qram_queries=2 * width,
        )
    )


def qram_state_angles(amplitudes: Iterable[complex], angle_width: int = 8) -> dict[int, int]:
    """由幅度向量生成 QRAM 旋转树的角表。

    键为 ``(1 << depth) - 1 + prefix`` 形式的树节点下标，值为量化到
    ``2*pi/2**angle_width`` 网格的 Ry 角度字，执行期作为角表 memory 提供。

    Args:
        amplitudes: 长度为二的幂、范数非零的非负实幅度序列。
        angle_width: 角度字的位宽。

    Returns:
        dict: 树节点下标到角度字的映射。

    Raises:
        ValidationError: 幅度向量非法，或含负实部、非零虚部的幅度。
    """
    values, _, nodes = _state_angles(amplitudes)
    if any(v.imag or v.real < 0 for v in values):
        raise ValidationError("当前 QRAM 角表实现接收非负实幅度")
    scale = (1 << angle_width) / (2 * math.pi)
    return {
        ((1 << depth) - 1 + prefix): round(angle * scale) % (1 << angle_width)
        for depth, prefix, _, angle in nodes
    }


def phase_marks(width: int, marked: Iterable[int]) -> Operation:
    """生成对指定计算基态施加相位 ``pi`` 的相位 oracle。

    Args:
        width: ``target`` 寄存器位宽。
        marked: 需要标记的基态整数编号的可迭代对象，重复项自动去重。

    Returns:
        Operation: 标注为 ``phase_oracle`` 范式的操作。
    """
    marked = tuple(sorted(set(marked)))
    b = Builder(_name("phase_marks", width, marked), {"target": Bits(width)})
    for value in marked:
        with b.control(b["target"], value):
            b.global_phase(math.pi)
    return annotate(b.finish(), "phase_oracle")


def _transposition(b: Builder, ref: Ref, first: int, second: int) -> None:
    """把 ``first`` 与 ``second`` 的对换编译为受控 X 的往返门序列。"""
    path = [first]
    for bit in range(ref.width):
        if ((first ^ second) >> bit) & 1:
            path.append(path[-1] ^ (1 << bit))
    edges = list(zip(path, path[1:], strict=False))
    for left, right in edges + list(reversed(edges[:-1])):
        bit = (left ^ right).bit_length() - 1
        controls = fuse(ref[:bit], ref[bit + 1 :])
        value = (left & ((1 << bit) - 1)) | ((left >> (bit + 1)) << bit)
        if controls.width:
            with b.control(controls, value):
                b.x(ref[bit])
        else:
            b.x(ref[bit])


def sparse_location_gate(
    width: int,
    permutations: Iterable[Iterable[int]],
    *,
    work_width: int | None = None,
    name: str | None = None,
) -> Operation:
    """以门级置换网络实现 CKS 原地位置查询。

    每列在列控制下把置换按循环分解为对换，对换再化为受控 X 的往返路径；
    门数随 ``2**width`` 增长，仅适合小实例见证。

    Args:
        width: ``column`` 与 ``index`` 寄存器位宽。
        permutations: 按列给出的 ``2**width`` 个置换，每个须为
            ``0 .. 2**width - 1`` 的完整排列。
        work_width: ``work`` 寄存器位宽，缺省取 ``width``。
        name: 模块名；缺省按内容确定性生成。

    Returns:
        Operation: 标注 ``sparse_location_inplace`` 范式与完整置换扩张属性的操作。

    Raises:
        ValidationError: 列数不是 ``2**width``，或某列不是完整置换。
    """
    work_width = width if work_width is None else work_width
    permutations = tuple(tuple(row) for row in permutations)
    dimension = 1 << width
    if len(permutations) != dimension or any(
        sorted(row) != list(range(dimension)) for row in permutations
    ):
        raise ValidationError("CKS 原地位置实现需要每列的完整置换扩张")
    b = Builder(
        name or _name("sparse_position", width, permutations, work_width),
        {"column": Bits(width), "index": Bits(width), "work": Bits(work_width)},
    )
    for column, permutation in enumerate(permutations):
        visited: set[int] = set()
        with b.control(b["column"], column):
            for first in range(dimension):
                if first in visited:
                    continue
                cycle: list[int] = []
                current = first
                while current not in visited:
                    visited.add(current)
                    cycle.append(current)
                    current = cast("tuple[int, ...]", permutation)[current]
                for second in cycle[1:]:
                    _transposition(b, b["index"], first, second)
    return annotate(b.finish(), "sparse_location_inplace", full_permutation_extension=True)


def sparse_location_qram(width: int) -> Operation:
    """以正反两张 QRAM 表实现 CKS 原地位置查询。

    三步完成：``work ^= forward[column, index]``、``swap(index, work)``、
    ``work ^= inverse[column, index]``；末步把 ``work`` 复净为零，``index``
    持有新值。

    Args:
        width: ``column``、``index`` 与 ``work`` 寄存器位宽。

    Returns:
        Operation: 标注 ``sparse_location_inplace`` 范式的操作，声明 ``forward``
        与 ``inverse`` 两张 ``QRAM(2*width, width)`` 资源。
    """
    b = Builder(
        f"sparse_position_qram_{width}",
        {"column": Bits(width), "index": Bits(width), "work": Bits(width)},
        {"forward": QRAM(2 * width, width), "inverse": QRAM(2 * width, width)},
    )
    b.qram("forward", fuse(b["column"], b["index"]), b["work"])
    b.swap(b["index"], b["work"])
    b.qram("inverse", fuse(b["column"], b["index"]), b["work"])
    return annotate(b.finish(), "sparse_location_inplace", full_permutation_extension=True)


def sparse_entry(database: XorDatabase, width: int) -> Operation:
    """把 XOR 数据库适配为任意行列的矩阵元素查询。

    ``row`` 与 ``column`` 拼接成 ``2*width`` 位地址驱动数据库，``data`` 保存
    XOR 结果。

    Args:
        database: ``address_width`` 恰为 ``2*width`` 的 XOR 数据库。
        width: 行、列各自的位宽。

    Returns:
        Operation: 标注 ``sparse_entry_xor`` 范式的元素查询操作。

    Raises:
        ValidationError: 数据库地址宽度不是 ``2*width``。
    """
    if database.address_width != 2 * width:
        raise ValidationError("矩阵条目查询需要 row 和 column 两组地址位")
    b = Builder(
        _name("matrix_entry", database.operation),
        {"row": Bits(width), "column": Bits(width), "data": Bits(database.data_width)},
        resources_for(("db", database.operation)),
    )
    invoke(b, database.operation, "db", address=fuse(b["row"], b["column"]), data=b["data"])
    return annotate(b.finish(), "sparse_entry_xor")


def diagonal_block_encoding(
    database: XorDatabase,
    *,
    alpha: float = 1.0,
    angle_scale: float | None = None,
) -> BlockEncoding:
    """按查询字控制信号 Ry，编码由该角表决定的对角矩阵。

    Args:
        database: 存放旋转角字的 XOR 数据库，地址驱动 target。
        alpha: 块编码归一化常数。
        angle_scale: 数据字第 k 位的角步长；缺省为 ``2π/2^数据位宽``。

    Returns:
        BlockEncoding: 对角元为 ``alpha*cos(encoded_angle/2)`` 的块编码。
    """
    angle_scale = 2 * math.pi / (1 << database.data_width) if angle_scale is None else angle_scale
    n, w = database.address_width, database.data_width
    b = Builder(
        _name("diagonal_be", database.operation, alpha, angle_scale),
        {"target": Bits(n), "signal": Bits(w + 1)},
        resources_for(("db", database.operation)),
    )
    word, flag = b["signal"][:w], b["signal"][w:]
    invoke(b, database.operation, "db", address=b["target"], data=word)
    for k in range(w):
        with b.control(word[k]):
            b.ry(flag, angle_scale * (1 << k))
    with b.adjoint():
        invoke(b, database.operation, "db", address=b["target"], data=word)
    return BlockEncoding(
        annotate(
            b.finish(),
            "block_encoding",
            be_alpha=alpha,
            matrix_interpretation="alpha*cos(encoded_angle/2) on diagonal",
        )
    )


def banked_database(
    address_width: int,
    data_width: int,
    *,
    name: str = "BankedData",
    abstract: bool = False,
) -> Operation:
    """将一个宽数据字拆成不超过 64 位的寄存器与 QRAM bank。

    Args:
        address_width: 地址字的位宽。
        data_width: 逻辑数据字的总位宽，须为正整数。
        name: 生成操作的名称，缺省为 ``BankedData``。
        abstract: 为 True 时生成体为空的开放声明而非 QRAM 实现。

    Returns:
        Operation: 按小端序拆分为多个不超过 64 位 bank 的数据库操作。
    """
    if type(data_width) is not int or data_width < 1:
        raise ValidationError("banked 数据宽度必须为正整数")
    chunks = tuple(min(64, data_width - offset) for offset in range(0, data_width, 64))
    registers = {
        "address": Bits(address_width),
        **{f"data{i}": Bits(w) for i, w in enumerate(chunks)},
    }
    if abstract:
        return declare(
            name,
            registers,
            paradigm="database_xor",
            attributes={"logical_data_width": data_width, "word_order": "little_endian_banks"},
        )
    b = Builder(
        _name("banked_qram", address_width, data_width),
        registers,
        {f"bank{i}": QRAM(address_width, w) for i, w in enumerate(chunks)},
    )
    for i, _ in enumerate(chunks):
        b.qram(f"bank{i}", b["address"], b[f"data{i}"])
    return annotate(
        b.finish(), "database_xor", logical_data_width=data_width, word_order="little_endian_banks"
    )
