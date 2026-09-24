"""Select-Swap QROM 数据加载（Low–Kliuchnikov–Schaeffer 2018；Babbush et al. 2018）。

QROM 查找把经典表 ``T`` 实现为 XOR 数据库 ``|address, data> -> |address, data XOR T[address]>``。
基线 ``qrom_lookup`` 是逐地址受控 XOR 的一元迭代；``select_swap_qrom`` 把地址拆为
高位 ``h``（k 位）与低位 ``y``（l 位），``λ = partitions = 2^l`` 个窗口分区共享高位
地址并行加载子表 ``T[h·λ + i]``，再按 ``y == i`` 受控交换归并，以 T 计数约
``4(2^k + λ·b)`` 在查询深度与辅助比特之间权衡。数据表与 gate_database 一样在生成时
展开为门级子数据库，不进入 IR JSON；资源估算由 ``qrom_cost`` 纯经典给出。
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.input_model.contracts import positive_integer
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import XorDatabase, annotate, gate_database, invoke
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, Ref, ValidationError

__all__ = ["QromCost", "qrom_cost", "qrom_lookup", "select_swap_qrom"]


def _normalize_table(
    table: Mapping[int, int] | Sequence[int], data_bits: int | None, path: str
) -> tuple[tuple[tuple[int, int], ...], int, int]:
    """校验查询表并给出 ``(有序条目, 地址位宽, 数据位宽)``；缺失地址按 0 处理。"""
    if isinstance(table, dict):
        pairs = tuple(sorted(table.items()))
    else:
        try:
            pairs = tuple(enumerate(table))
        except TypeError as exc:
            raise ValidationError(f"{path}：查询表必须是字序列或稀疏字典") from exc
    if not pairs:
        raise ValidationError(f"{path}：查询表不能为空")
    for address, value in pairs:
        if type(address) is not int or address < 0:
            raise ValidationError(f"{path}：表地址必须是非负整数")
        if type(value) is not int or value < 0:
            raise ValidationError(f"{path}：表字必须是非负整数")
    address_bits = max(1, max(address for address, _ in pairs).bit_length())
    peak = max(value for _, value in pairs)
    if data_bits is not None:
        positive_integer(data_bits, path + ".data_bits", maximum=64)
        if peak >= 1 << data_bits:
            raise ValidationError(f"{path}：表字超出 data_bits 位宽")
        width = data_bits
    else:
        width = max(1, peak.bit_length())
    return pairs, address_bits, width


def _check_partitions(partitions: int, address_bits: int, path: str) -> None:
    """校验分区数是不超过地址空间且为二的幂的正整数。"""
    positive_integer(partitions, path + ".partitions", maximum=1 << address_bits)
    if partitions & (partitions - 1):
        raise ValidationError(f"{path}：partitions 必须是二的幂")


@dataclass(frozen=True)
class QromCost:
    """Select-Swap QROM 的资源估算快照；Toffoli/T 计数为 compute 单程。"""

    n_addresses: int
    data_bits: int
    partitions: int
    address_bits: int
    select_toffoli: int
    swap_toffoli: int
    work_qubits: int
    fanout_qubits: int

    @property
    def t_count(self) -> int:
        """compute 单程 T 计数，按每个 Toffoli 4 个 T 估计。"""
        return 4 * (self.select_toffoli + self.swap_toffoli)

    @property
    def t_depth(self) -> int:
        """T 深度：分区并行加载时 select 段为串行一元迭代，交换段按分区逐级归并。"""
        return self.select_toffoli + max(0, self.partitions - 1)

    @property
    def round_trip_t_count(self) -> int:
        """相干复净（compute 加伴随 uncompute）的往返 T 计数；测量复净可省掉 uncompute。"""
        return 2 * self.t_count

    @property
    def ancilla_qubits(self) -> int:
        """辅助比特总数：λ 个数据窗口加低位扇出副本（后者可用 dirty qubit）。"""
        return self.work_qubits + self.fanout_qubits

    def to_dict(self) -> dict[str, int]:
        """供目录与后端报告使用的结构化字典。

        Returns:
            dict[str, int]: 各项资源计数与宽度字段的扁平字典。
        """
        return {
            "n_addresses": self.n_addresses,
            "data_bits": self.data_bits,
            "partitions": self.partitions,
            "address_bits": self.address_bits,
            "select_toffoli": self.select_toffoli,
            "swap_toffoli": self.swap_toffoli,
            "t_count": self.t_count,
            "t_depth": self.t_depth,
            "round_trip_t_count": self.round_trip_t_count,
            "work_qubits": self.work_qubits,
            "fanout_qubits": self.fanout_qubits,
            "ancilla_qubits": self.ancilla_qubits,
        }


def qrom_cost(n_addresses: int, data_bits: int, partitions: int = 1) -> QromCost:
    """Select-Swap QROM 的纯经典资源估算：T 计数/深度约 ``4(N/λ + λ·b)``，λ 即 partitions。

    Args:
        n_addresses: 查询表字数 N。
        data_bits: 数据字位宽 b。
        partitions: 分区数 λ，必须是二的幂且不超过 ``2^ceil(log2 N)``。

    Returns:
        QromCost: 该配置下的 Toffoli/T 计数与辅助比特估算快照。
    """
    positive_integer(n_addresses, "qrom_cost.n_addresses")
    positive_integer(data_bits, "qrom_cost.data_bits", maximum=64)
    address_bits = max(1, (n_addresses - 1).bit_length())
    _check_partitions(partitions, address_bits, "qrom_cost")
    low_bits = partitions.bit_length() - 1
    high_bits = address_bits - low_bits
    select_toffoli = (1 << high_bits) - 1
    swap_toffoli = data_bits * (partitions - 1)
    work_qubits = partitions * data_bits
    fanout_qubits = (partitions - 1) * low_bits
    return QromCost(
        n_addresses,
        data_bits,
        partitions,
        address_bits,
        select_toffoli,
        swap_toffoli,
        work_qubits,
        fanout_qubits,
    )


def _cost_attributes(cost: QromCost) -> dict[str, int]:
    """把 ``QromCost`` 快照转成可写入模块属性的扁平字典。"""
    return {
        "qrom_partitions": cost.partitions,
        "qrom_address_bits": cost.address_bits,
        "qrom_data_bits": cost.data_bits,
        "select_toffoli": cost.select_toffoli,
        "swap_toffoli": cost.swap_toffoli,
        "t_count": cost.t_count,
        "round_trip_t_count": cost.round_trip_t_count,
        "work_qubits": cost.work_qubits,
        "dirty_fanout_qubits": cost.fanout_qubits,
    }


def qrom_lookup(
    table: Mapping[int, int] | Sequence[int],
    *,
    data_bits: int | None = None,
    name: str | None = None,
) -> XorDatabase:
    """QROM 基线：逐地址受控 XOR 的一元迭代（即 partitions=1 的 Select-Swap）。

    结构与 gate_database 相同，但标注 qrom_unary_iteration 并附带 qrom_cost 估算，
    作为 Select-Swap 权衡曲线的 λ=1 端点。

    Args:
        table: 字序列或稀疏字典，缺失地址按 0 处理。
        data_bits: 数据位宽，缺省取最大表字的位宽。
        name: 模块名，缺省按表内容生成。

    Returns:
        XorDatabase: 带 qrom_cost 标注的一元迭代数据库视图。
    """
    pairs, address_bits, width = _normalize_table(table, data_bits, "qrom_lookup")
    base = gate_database(address_bits, width, dict(pairs), name=name)
    cost = qrom_cost(1 << address_bits, width, 1)
    return XorDatabase(
        annotate(
            base.operation,
            "database_xor",
            implementation="qrom_unary_iteration",
            **_cost_attributes(cost),
        )
    )


def select_swap_qrom(
    table: Mapping[int, int] | Sequence[int],
    *,
    partitions: int,
    data_bits: int | None = None,
    name: str | None = None,
) -> XorDatabase:
    """Select-Swap QROM：高位地址共享的分区并行加载加低位受控交换归并。

    λ 个窗口寄存器（clean 局部寄存器）各自执行子表查询 ``window_i ^= T[h·λ + i]``，
    全部窗口共享高位地址 h，因此 select 段只支付一次 2^k 一元迭代；随后按
    ``y == i`` 对 ``(window_i, data)`` 施加受控 swap–xor–swap 复合（在 data 任意初值下
    保持 database_xor 的 XOR 语义，窗口内容不被破坏），最后用伴随重放子查询复净全部
    窗口。低位 y 的扇出副本（可用 dirty qubit，计入 qrom_cost 的 fanout_qubits）在
    RIR 中由 Control 原语隐含，留待后端降低时展开。

    Args:
        table: 字序列或稀疏字典，缺失地址按 0 处理。
        partitions: 分区数 λ，必须是二的幂且不超过地址数。
        data_bits: 数据位宽，缺省取最大表字的位宽。
        name: 模块名，缺省按表内容生成。

    Returns:
        XorDatabase: 带 qrom_cost 标注的 Select-Swap 数据库视图。
    """
    pairs, address_bits, width = _normalize_table(table, data_bits, "select_swap_qrom")
    _check_partitions(partitions, address_bits, "select_swap_qrom")
    low_bits = partitions.bit_length() - 1
    high_bits = address_bits - low_bits
    values = dict(pairs)
    sub_tables = tuple(
        {high: values.get((high << low_bits) | index, 0) for high in range(1 << high_bits)}
        for index in range(partitions)
    )
    sub_databases = (
        None
        if high_bits == 0
        else tuple(gate_database(high_bits, width, sub) for sub in sub_tables)
    )
    b = Builder(
        name or _name("select_swap", address_bits, width, partitions, pairs),
        {"address": Bits(address_bits), "data": Bits(width)},
    )
    windows = [b.local(f"window{index}", Bits(width)) for index in range(partitions)]

    def load() -> None:
        """每个窗口并行加载自己分区对应的子表。"""
        for index, window in enumerate(windows):
            if high_bits:
                invoke(
                    b,
                    cast("tuple[XorDatabase, ...]", sub_databases)[index].operation,
                    f"sub{index}",
                    address=b["address"][low_bits:],
                    data=window,
                )
            else:
                for bit in range(width):
                    if (sub_tables[index][0] >> bit) & 1:
                        b.x(window[bit])

    def merge(index: int, window: Ref) -> None:
        """按低位地址比较把窗口内容归并进数据寄存器。"""
        data = b["data"]
        if partitions == 1:
            b.swap(window, data)
            b.xor(data, window)
            b.swap(window, data)
        else:
            with b.control(b["address"][:low_bits], index):
                b.swap(window, data)
                b.xor(data, window)
                b.swap(window, data)

    load()
    for index, window in enumerate(windows):
        merge(index, window)
    with b.adjoint():
        load()
    cost = qrom_cost(1 << address_bits, width, partitions)
    return XorDatabase(
        annotate(
            b.finish(),
            "database_xor",
            implementation="select_swap",
            **_cost_attributes(cost),
        )
    )
