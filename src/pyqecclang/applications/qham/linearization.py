"按 HAM 阶数和有序张量字惰性生成 QCL；不物化空间矩阵。"

from __future__ import annotations

import json
import math
from dataclasses import dataclass

from pyqecclang.applications.qham.pde import PolynomialPDE
from pyqecclang.infrastructure.ir import ValidationError


def compositions(total, length):
    if length == 0:
        if total == 0:
            yield ()
    elif length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, length - 1):
                yield (first, *rest)


@dataclass(frozen=True)
class HomotopyWeight:
    kind: str = "one"
    power: int = 0

    def evaluate(self, eta):
        if self.kind == "one":
            return 1.0
        if self.kind == "correction":
            return -eta * (1 + eta) ** self.power
        if self.kind == "physical":
            return 1 - (1 + eta) ** self.power
        raise ValidationError("未知同伦权重")

    def formula(self):
        if self.kind == "one":
            return "1"
        if self.kind == "correction":
            return f"-eta*(1+eta)^{self.power}"
        return f"1-(1+eta)^{self.power}"


@dataclass(frozen=True)
class Block:
    kind: str
    orders: tuple[int, ...] = ()

    @property
    def rank(self):
        return 1 if self.kind == "physical" else len(self.orders)

    @property
    def label(self):
        if self.kind == "physical":
            return "u_sum"
        if not self.orders:
            return "one"
        return "Y_" + "_".join(map(str, self.orders))


@dataclass(frozen=True)
class Coupling:
    row: Block
    column: Block
    operator: str
    position: int
    arity: int
    weight: HomotopyWeight = HomotopyWeight()


@dataclass(frozen=True)
class QHAMPlan:
    pde: PolynomialPDE
    order: int
    version: str = "0.1"

    def __post_init__(self):
        self.pde.validate()
        if self.version != "0.1" or type(self.order) is not int or self.order < 0:
            raise ValidationError("QHAM 截断阶必须为非负整数")

    @property
    def grade(self):
        return max(1, self.pde.degree - 1)

    @property
    def max_rank(self):
        return self.grade * self.order + 1

    @property
    def has_forcing(self):
        return any(port.name == "F" for port in self.pde.ports)

    def weight(self, word):
        return self.grade * sum(word) + len(word)

    def contains(self, block):
        if block.kind == "physical":
            return block.orders == ()
        if block.kind != "tensor":
            return False
        if not block.orders:
            return self.has_forcing
        return (
            all(type(a) is int and a >= 0 for a in block.orders)
            and self.weight(block.orders) <= self.max_rank
        )

    def rank_count(self, rank):
        limit = (self.max_rank - rank) // self.grade
        return math.comb(limit + rank, rank) if limit >= 0 and rank >= 1 else 0

    @property
    def block_count(self):
        return (
            1 + int(self.has_forcing) + sum(self.rank_count(k) for k in range(1, self.max_rank + 1))
        )

    def blocks(self, *, max_blocks=None):
        if max_blocks is not None and self.block_count > max_blocks:
            raise ValidationError("显式块生成超过预算；QHAMPlan 仍可惰性查询")
        yield Block("physical")
        for rank in range(1, self.max_rank + 1):
            limit = (self.max_rank - rank) // self.grade
            for total in range(limit + 1):
                for word in compositions(total, rank):
                    yield Block("tensor", word)
        if self.has_forcing:
            yield Block("tensor")

    def row_terms(self, row):
        if not self.contains(row):
            raise ValidationError("QCL 行不属于有限闭包")
        if row.kind == "tensor" and not row.orders:
            return ()
        ports = self.pde.ports
        linear = any(p.name == "L" for p in ports)
        nonlinear = [p for p in ports if p.arity >= 2]
        result = []
        if row.kind == "physical":
            if linear:
                result.append(Coupling(row, row, "L", 0, 1))
            if self.has_forcing:
                result.append(Coupling(row, Block("tensor"), "F", 0, 0))
            for total in range(self.order):
                for port in nonlinear:
                    for word in compositions(total, port.arity):
                        result.append(
                            Coupling(
                                row,
                                Block("tensor", word),
                                port.name,
                                0,
                                port.arity,
                                HomotopyWeight("physical", self.order - total),
                            )
                        )
        else:
            for position, order in enumerate(row.orders):
                if linear:
                    result.append(Coupling(row, row, "L", position, 1))
                if order == 0:
                    if self.has_forcing:
                        word = row.orders[:position] + row.orders[position + 1 :]
                        result.append(Coupling(row, Block("tensor", word), "F", position, 0))
                else:
                    for total in range(order):
                        for port in nonlinear:
                            for word in compositions(total, port.arity):
                                source = row.orders[:position] + word + row.orders[position + 1 :]
                                result.append(
                                    Coupling(
                                        row,
                                        Block("tensor", source),
                                        port.name,
                                        position,
                                        port.arity,
                                        HomotopyWeight("correction", order - 1 - total),
                                    )
                                )
        if any(not self.contains(term.column) for term in result):
            raise AssertionError("QCL 闭包构造错误")
        return tuple(result)

    def raw_dimension(self, dimension):
        return (
            dimension
            + sum(self.rank_count(k) * dimension**k for k in range(1, self.max_rank + 1))
            + int(self.has_forcing)
        )

    def offset(self, block, dimension):
        if not self.contains(block):
            raise ValidationError("未知 QCL 块")
        if block.kind == "physical":
            return 0
        if not block.orders:
            return self.raw_dimension(dimension) - 1
        rank, total = len(block.orders), sum(block.orders)
        offset = dimension + sum(self.rank_count(k) * dimension**k for k in range(1, rank))
        before = math.comb(total + rank - 1, rank) if total else 0
        remainder = total
        for pos, value in enumerate(block.orders[:-1]):
            slots = rank - pos - 1
            for first in range(value):
                before += math.comb(remainder - first + slots - 1, slots - 1)
            remainder -= value
        return offset + before * dimension**rank

    def locate(self, index, dimension):
        if not 0 <= index < self.raw_dimension(dimension):
            raise ValidationError("QCL 原始索引越界")
        if index < dimension:
            return Block("physical"), index
        index -= dimension
        for rank in range(1, self.max_rank + 1):
            size = dimension**rank
            group = self.rank_count(rank) * size
            if index >= group:
                index -= group
                continue
            ordinal, local = divmod(index, size)
            total = 0
            while ordinal >= math.comb(total + rank - 1, rank - 1):
                ordinal -= math.comb(total + rank - 1, rank - 1)
                total += 1
            word = []
            remaining = total
            for pos in range(rank - 1):
                slots = rank - pos - 1
                first = 0
                while ordinal >= math.comb(remaining - first + slots - 1, slots - 1):
                    ordinal -= math.comb(remaining - first + slots - 1, slots - 1)
                    first += 1
                word.append(first)
                remaining -= first
            word.append(remaining)
            return Block("tensor", tuple(word)), local
        return Block("tensor"), 0

    def summary(self):
        return {
            "version": self.version,
            "order": self.order,
            "pde_degree": self.pde.degree,
            "max_tensor_rank": self.max_rank,
            "block_count": self.block_count,
            "homogeneous_constant": self.has_forcing,
            "closure": "(D-1)*sum(orders)+rank <= (D-1)*m+1; D<=1 uses grade 1",
            "operators": [
                {"name": p.name, "input_rank": p.arity, "output_rank": 1} for p in self.pde.ports
            ],
        }

    def dumps(self):
        return (
            json.dumps(
                {"version": self.version, "order": self.order, "pde": json.loads(self.pde.dumps())},
                ensure_ascii=False,
                sort_keys=True,
                indent=2,
            )
            + "\n"
        )

    @classmethod
    def loads(cls, text):
        try:
            raw = json.loads(text)
            if set(raw) != {"version", "order", "pde"}:
                raise ValidationError("未知 QCL plan 字段")
            return cls(PolynomialPDE.loads(json.dumps(raw["pde"])), raw["order"], raw["version"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError("无效 QCL plan：" + str(exc)) from exc
