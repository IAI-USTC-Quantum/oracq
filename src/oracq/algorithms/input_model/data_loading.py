"""Select-Swap QROM data loading (Low–Kliuchnikov–Schaeffer 2018; Babbush et al. 2018).

A QROM lookup realizes a classical table ``T`` as the XOR database
``|address, data> -> |address, data XOR T[address]>``. The baseline
``qrom_lookup`` is a unary iteration of per-address controlled XOR;
``select_swap_qrom`` splits the address into high bits ``h`` (k bits) and low
bits ``y`` (l bits); ``λ = partitions = 2^l`` window partitions share the high
address and load the sub-tables ``T[h·λ + i]`` in parallel, then merge by
controlled swaps keyed on ``y == i``, trading query depth against ancilla
qubits at a T count of about ``4(2^k + λ·b)``. Like gate_database, the data
table is expanded into gate-level sub-databases at generation time and does not
enter the IR JSON; resource estimation is given purely classically by
``qrom_cost``.
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
    """Validate the lookup table and return ``(ordered entries, address bit width, data bit width)``; missing addresses are treated as 0."""
    if isinstance(table, dict):
        pairs = tuple(sorted(table.items()))
    else:
        try:
            pairs = tuple(enumerate(table))
        except TypeError as exc:
            raise ValidationError(f"{path}: lookup table must be a sequence of words or a sparse dict") from exc
    if not pairs:
        raise ValidationError(f"{path}: lookup table cannot be empty")
    for address, value in pairs:
        if type(address) is not int or address < 0:
            raise ValidationError(f"{path}: table addresses must be non-negative integers")
        if type(value) is not int or value < 0:
            raise ValidationError(f"{path}: table words must be non-negative integers")
    address_bits = max(1, max(address for address, _ in pairs).bit_length())
    peak = max(value for _, value in pairs)
    if data_bits is not None:
        positive_integer(data_bits, path + ".data_bits", maximum=64)
        if peak >= 1 << data_bits:
            raise ValidationError(f"{path}: table word exceeds the data_bits width")
        width = data_bits
    else:
        width = max(1, peak.bit_length())
    return pairs, address_bits, width


def _check_partitions(partitions: int, address_bits: int, path: str) -> None:
    """Validate that the partition count is a positive integer, a power of two, and no larger than the address space."""
    positive_integer(partitions, path + ".partitions", maximum=1 << address_bits)
    if partitions & (partitions - 1):
        raise ValidationError(f"{path}: partitions must be a power of two")


@dataclass(frozen=True)
class QromCost:
    """Resource-estimation snapshot for Select-Swap QROM; the Toffoli/T counts are for a single compute pass."""

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
        """T count of a single compute pass, estimated at 4 T per Toffoli."""
        return 4 * (self.select_toffoli + self.swap_toffoli)

    @property
    def t_depth(self) -> int:
        """T depth: with partitioned parallel loading the select segment is a serial unary iteration and the swap segment merges level by level per partition."""
        return self.select_toffoli + max(0, self.partitions - 1)

    @property
    def round_trip_t_count(self) -> int:
        """Round-trip T count for coherent uncomputation (compute plus adjoint uncompute); measurement-based uncomputation can drop the uncompute."""
        return 2 * self.t_count

    @property
    def ancilla_qubits(self) -> int:
        """Total ancilla qubits: the λ data windows plus the low-bit fanout copies (the latter may use dirty qubits)."""
        return self.work_qubits + self.fanout_qubits

    def to_dict(self) -> dict[str, int]:
        """Structured dictionary for gallery and backend reports.

        Returns:
            dict[str, int]: Flat dictionary of the resource counts and width fields.
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
    """Purely classical resource estimation for Select-Swap QROM: T count/depth about ``4(N/λ + λ·b)``, with λ being partitions.

    Args:
        n_addresses: Number N of words in the lookup table.
        data_bits: Data word bit width b.
        partitions: Partition count λ; must be a power of two not exceeding ``2^ceil(log2 N)``.

    Returns:
        QromCost: Snapshot of the Toffoli/T counts and ancilla-qubit estimates for this configuration.
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
    """Convert a ``QromCost`` snapshot into a flat dictionary writable as module attributes."""
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
    """QROM baseline: a unary iteration of per-address controlled XOR (i.e. Select-Swap with partitions=1).

    Same structure as gate_database but annotated qrom_unary_iteration and
    carrying a qrom_cost estimate, serving as the λ=1 endpoint of the
    Select-Swap trade-off curve.

    Args:
        table: Sequence of words or a sparse dict; missing addresses are treated as 0.
        data_bits: Data bit width; defaults to the width of the largest table word.
        name: Module name; generated from the table content when omitted.

    Returns:
        XorDatabase: Unary-iteration database view annotated with qrom_cost.
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
    """Select-Swap QROM: partitioned parallel loading sharing the high address bits, merged by low-bit controlled swaps.

    The λ window registers (clean local registers) each perform the sub-table
    query ``window_i ^= T[h·λ + i]``; all windows share the high address h, so
    the select segment pays for only one 2^k unary iteration. Then, keyed on
    ``y == i``, the controlled swap–xor–swap composite is applied to
    ``(window_i, data)`` (preserving the database_xor XOR semantics under
    arbitrary initial data values without destroying the window contents), and
    finally the sub-queries are replayed adjointly to uncompute all windows. The
    fanout copies of the low bits y (which may use dirty qubits and are counted
    in qrom_cost's fanout_qubits) are implied by the Control primitive in the
    RIR and left for expansion during backend lowering.

    Args:
        table: Sequence of words or a sparse dict; missing addresses are treated as 0.
        partitions: Partition count λ; must be a power of two not exceeding the number of addresses.
        data_bits: Data bit width; defaults to the width of the largest table word.
        name: Module name; generated from the table content when omitted.

    Returns:
        XorDatabase: Select-Swap database view annotated with qrom_cost.
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
        """Each window loads the sub-table of its own partition in parallel."""
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
        """Merge the window contents into the data register by comparing the low address bits."""
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
