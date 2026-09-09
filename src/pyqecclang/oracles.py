"""Oracle 范式目录及普通 gate/QRAM 实现。数学正确性认证不属于此模块。"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass, replace

from .builder import Builder, Operation
from .ir import QRAM, Bits, Module, Register, Resource, ValidationError, fuse
from .library import BlockEncoding, _name

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


def declare(
    name,
    registers,
    *,
    paradigm="unitary",
    resources=None,
    attributes=None,
    supports_adjoint=True,
    supports_controlled=True,
) -> Operation:
    if paradigm not in PARADIGMS:
        raise ValidationError(f"未知 oracle paradigm：{paradigm}")
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


def annotate(operation, paradigm, **attributes):
    attrs = dict(operation.module.attributes)
    attrs.update(
        oracle_paradigm=paradigm,
        implementation_status="constructed",
        supports_adjoint=True,
        supports_controlled=True,
        **attributes,
    )
    return Operation(
        replace(operation.module, attributes=tuple(sorted(attrs.items()))), operation.dependencies
    )


def resources_for(*items):
    return {f"{prefix}__{r.name}": r.type for prefix, op in items for r in op.module.resources}


def invoke(builder, operation, prefix="", **arguments):
    mapping = {
        r.name: f"{prefix}__{r.name}" if prefix else r.name for r in operation.module.resources
    }
    builder.call(operation, resources=mapping, **arguments)


@dataclass(frozen=True)
class XorDatabase:
    operation: Operation

    def __post_init__(self):
        regs = {r.name: r.type for r in self.operation.module.registers}
        if set(regs) != {"address", "data"} or any(r.kind != "bits" for r in regs.values()):
            raise ValidationError("XOR database 需要 address/data bits 接口")

    @property
    def address_width(self):
        return next(r.type.width for r in self.operation.module.registers if r.name == "address")

    @property
    def data_width(self):
        return next(r.type.width for r in self.operation.module.registers if r.name == "data")


@dataclass(frozen=True)
class StatePreparation:
    operation: Operation

    def __post_init__(self):
        regs = {r.name: r.type for r in self.operation.module.registers}
        if set(regs) != {"target", "work"} or any(r.kind != "bits" for r in regs.values()):
            raise ValidationError("StatePreparation 需要 target/work bits 接口")

    @property
    def width(self):
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def work_width(self):
        return next(r.type.width for r in self.operation.module.registers if r.name == "work")


@dataclass(frozen=True)
class StateOracle:
    operation: Operation

    @property
    def width(self):
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def signal_qubits(self):
        return next(r.type.width for r in self.operation.module.registers if r.name == "signal")


@dataclass(frozen=True)
class SparseAccess:
    location: Operation
    entry: Operation
    width: int
    value_width: int
    sparsity: int

    def __post_init__(self):
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


def abstract_database(name, address_width, data_width):
    return XorDatabase(
        declare(
            name,
            {"address": Bits(address_width), "data": Bits(data_width)},
            paradigm="database_xor",
        )
    )


def abstract_state_prep(name, width, work_width=0, *, reversible=True):
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


def abstract_block_encoding(name, width, signal_width, alpha):
    return BlockEncoding(
        declare(
            name,
            {"target": Bits(width), "signal": Bits(signal_width)},
            paradigm="block_encoding",
            attributes={"be_alpha": alpha},
        )
    )


def abstract_sparse_access(name, width, value_width, sparsity, work_width=None):
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


def gate_database(address_width, data_width, table, *, name=None):
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


def qram_database(address_width, data_width, *, name=None):
    b = Builder(
        name or f"qram_xor_{address_width}_{data_width}",
        {"address": Bits(address_width), "data": Bits(data_width)},
        {"table": QRAM(address_width, data_width)},
    )
    b.qram("table", b["address"], b["data"])
    return XorDatabase(annotate(b.finish(), "database_xor", implementation="qram"))


def basis_state(width, value=0, *, work_width=0):
    if type(value) is not int or not 0 <= value < 1 << width:
        raise ValidationError("基态值越界")
    b = Builder(
        _name("basis", width, value, work_width), {"target": Bits(width), "work": Bits(work_width)}
    )
    for bit in range(width):
        if (value >> bit) & 1:
            b.x(b["target"][bit])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def uniform_state(width, *, work_width=0):
    b = Builder(
        _name("uniform", width, work_width), {"target": Bits(width), "work": Bits(work_width)}
    )
    b.h(b["target"])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def _state_angles(amplitudes):
    values = tuple(complex(v) for v in amplitudes)
    if not values or len(values) & (len(values) - 1):
        raise ValidationError("态向量长度必须是二的幂")
    if not all(math.isfinite(v.real) and math.isfinite(v.imag) for v in values):
        raise ValidationError("态幅度必须有限")
    if sum(abs(v) ** 2 for v in values) == 0:
        raise ValidationError("零向量不能制备为归一化态")
    n = (len(values) - 1).bit_length()
    nodes = []
    for depth in range(n):
        bit = n - depth - 1
        for prefix in range(1 << depth):
            start, size = prefix << (bit + 1), 1 << bit
            left = sum(abs(v) ** 2 for v in values[start : start + size])
            right = sum(abs(v) ** 2 for v in values[start + size : start + 2 * size])
            nodes.append((depth, prefix, bit, 2 * math.atan2(math.sqrt(right), math.sqrt(left))))
    return values, n, nodes


def gate_state_prep(amplitudes, *, work_width=0, name=None):
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


def qram_state_prep(width, angle_width=8):
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


def qram_state_angles(amplitudes, angle_width=8):
    values, _, nodes = _state_angles(amplitudes)
    if any(v.imag or v.real < 0 for v in values):
        raise ValidationError("当前 QRAM 角表实现接收非负实幅度")
    scale = (1 << angle_width) / (2 * math.pi)
    return {
        ((1 << depth) - 1 + prefix): round(angle * scale) % (1 << angle_width)
        for depth, prefix, _, angle in nodes
    }


def phase_marks(width, marked):
    marked = tuple(sorted(set(marked)))
    b = Builder(_name("phase_marks", width, marked), {"target": Bits(width)})
    for value in marked:
        with b.control(b["target"], value):
            b.global_phase(math.pi)
    return annotate(b.finish(), "phase_oracle")


def _transposition(b, ref, first, second):
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


def sparse_location_gate(width, permutations, *, work_width=None, name=None):
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
        visited = set()
        with b.control(b["column"], column):
            for first in range(dimension):
                if first in visited:
                    continue
                cycle, current = [], first
                while current not in visited:
                    visited.add(current)
                    cycle.append(current)
                    current = permutation[current]
                for second in cycle[1:]:
                    _transposition(b, b["index"], first, second)
    return annotate(b.finish(), "sparse_location_inplace", full_permutation_extension=True)


def sparse_location_qram(width):
    b = Builder(
        f"sparse_position_qram_{width}",
        {"column": Bits(width), "index": Bits(width), "work": Bits(width)},
        {"forward": QRAM(2 * width, width), "inverse": QRAM(2 * width, width)},
    )
    b.qram("forward", fuse(b["column"], b["index"]), b["work"])
    b.swap(b["index"], b["work"])
    b.qram("inverse", fuse(b["column"], b["index"]), b["work"])
    return annotate(b.finish(), "sparse_location_inplace", full_permutation_extension=True)


def sparse_entry(database: XorDatabase, width):
    if database.address_width != 2 * width:
        raise ValidationError("矩阵条目查询需要 row 和 column 两组地址位")
    b = Builder(
        _name("matrix_entry", database.operation),
        {"row": Bits(width), "column": Bits(width), "data": Bits(database.data_width)},
        resources_for(("db", database.operation)),
    )
    invoke(b, database.operation, "db", address=fuse(b["row"], b["column"]), data=b["data"])
    return annotate(b.finish(), "sparse_entry_xor")


def diagonal_block_encoding(database: XorDatabase, *, alpha=1.0, angle_scale=None):
    """按查询字控制信号 Ry，编码由该角表决定的对角矩阵。"""
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


def banked_database(address_width, data_width, *, name="BankedData", abstract=False):
    """将一个宽数据字拆成不超过 64 位的寄存器与 QRAM bank。"""
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
