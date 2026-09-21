"""量子牛顿法（arXiv:2109.08470）的实现：M_F 数据结构与差分 Jacobian oracle。

算法逐轮执行牛顿迭代 A Δx = -F(x)，其中线性求解与数据装载使用量子构件：

- ``NewtonTree`` 是论文的 M_F：叶存 x_i 与 f_i(x)，内部节点存平方部分和，
  角度缓存按 RY 半角约定写入 bank；单点更新只重算依赖该分量的方程。
- ``newton_fd_entry`` 是 O_A2（论文 Eq. 15-17）：给定 (row, slot)，经
  模式库 O_f1 查变量下标、经 x 库 O_M1 取变量值，对被选变量加 delta 后
  两次调用可逆编译的 f_row（O_f2），差商 (f(x+Δe_k)-f(x))/(Δ·scale) 由
  编译算术给出；所有工作位经 XOR 复净。
- ``newton_b_preparation`` 是 O_b：符号残差角树制备 ``|b> = Σ f_i |i> / C_b``。

线性求解（Hermitian 扩张 + 稀疏 QLSS）、l_inf 层析与范数恢复由
``newton_linear_step`` 组装；全部构件的语义验证见 tests/core/test_newton.py。
"""

from __future__ import annotations

import inspect
import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import cast

from pyqecclang.algorithms.common.arithmetic import FixedFormat
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.algorithms.input_model.oracles import (
    StatePreparation,
    annotate,
    invoke,
    qram_database,
    qram_state_prep,
    resources_for,
)
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, Ref, ValidationError, fuse
from pyqecclang.infrastructure.mathfunc import CompiledFunction


def _log2_exact(value: int) -> int:
    """校验 value 是正的 2 的幂并返回其以 2 为底的对数。"""
    if value < 1 or value & (value - 1):
        raise ValidationError("牛顿法的维度必须是 2 的幂")
    return value.bit_length() - 1


@dataclass(frozen=True)
class NewtonFunction:
    """稀疏非线性系统 F(x)=0：每方程一个纯 Python 数值函数与其依赖下标。"""

    functions: tuple[Callable[..., float], ...]
    pattern: tuple[tuple[int, ...], ...]

    @classmethod
    def declare(
        cls, equations: Sequence[tuple[Sequence[int], Callable[..., float]]]
    ) -> NewtonFunction:
        """从 (依赖下标, 数值函数) 序列声明一个稀疏非线性系统。

        Args:
            equations: ``(indices, function)`` 序列；indices 为该方程依赖的变量
                下标（与 function 的位置形参一一对应），function 为纯 Python 数值函数。

        Returns:
            NewtonFunction: 通过维度与稀疏依赖校验后的系统描述。

        Raises:
            ValidationError: 方程数不是 2 的幂，或某方程的下标重复、越界。
            OSError: 某个函数取不到源码（``inspect.getsource``，量子编译需要）。
        """
        pattern = tuple(tuple(indices) for indices, _ in equations)
        functions = tuple(function for _, function in equations)
        size = len(equations)
        _log2_exact(size)
        for row, indices in enumerate(pattern):
            if len(set(indices)) != len(indices) or any(not 0 <= j < size for j in indices):
                raise ValidationError(f"方程 {row} 的稀疏依赖非法")
            inspect.getsource(functions[row])
        return cls(functions, pattern)

    @property
    def size(self) -> int:
        """方程个数；系统为方系统，也等于变量个数。"""
        return len(self.functions)

    @property
    def width(self) -> int:
        """变量（方程）下标寄存器的位宽，即 log2(size)；声明时已保证 size 是 2 的幂。"""
        return _log2_exact(self.size)

    @property
    def sparsity(self) -> int:
        """稀疏度：单个方程依赖变量个数的最大值。"""
        return max(len(row) for row in self.pattern)

    def evaluate(self, x: Sequence[float]) -> list[float]:
        """经典求值 F(x)。

        Args:
            x: 长度为 ``size`` 的实数向量。

        Returns:
            list: 按声明顺序排列的各方程函数值 f_i(x)。
        """
        return [
            function(*(x[j] for j in indices))
            for function, indices in zip(self.functions, self.pattern, strict=True)
        ]

    def jacobian(self, x: Sequence[float], delta: float) -> list[dict[int, float]]:
        """有限差分 Jacobian（经典镜像，作为量子 oracle 的独立参考）。

        Args:
            x: 长度为 ``size`` 的实数向量，差分的基准点。
            delta: 前向差分步长，取非零实数。

        Returns:
            list[dict[int, float]]: 每行一个稀疏字典，键为列号、值为差商偏导。
        """
        base = self.evaluate(x)
        rows: list[dict[int, float]] = []
        for row, indices in enumerate(self.pattern):
            entry: dict[int, float] = {}
            for _k, j in enumerate(indices):
                shifted = list(x)
                shifted[j] = x[j] + delta
                up = self.functions[row](*(shifted[j2] for j2 in indices))
                entry[j] = (up - base[row]) / delta
            rows.append(entry)
        return rows


class NewtonTree:
    """M_F 数据结构（经典侧）。堆为 1-based：根为 1，叶 N..2N-1。"""

    def __init__(
        self,
        function: NewtonFunction,
        x: Sequence[float],
        *,
        fmt: FixedFormat | None = None,
        angle_width: int = 10,
        delta: float = 2.0**-6,
    ) -> None:
        """以系统描述与初值构建 M_F 树，并自底向上填充全部存储表。

        Args:
            function: 被编码的稀疏非线性系统。
            x: 长度与 ``function.size`` 一致的实数初值向量。
            fmt: 数值 bank 使用的定点格式；省略时取 ``FixedFormat(12, 6)``。
            angle_width: RY 角度缓存的位宽。
            delta: 差分步长，必须在 ``fmt`` 下精确可表示且非零。

        Raises:
            ValidationError: delta 不可精确表示，或初值维度与系统不匹配。
        """
        fmt = fmt or FixedFormat(12, 6)
        if fmt.encode(delta) == 0:
            raise ValidationError("delta 必须在定点格式下精确可表示且非零")
        if len(x) != function.size:
            raise ValidationError("初值维度与系统不匹配")
        self.function: NewtonFunction = function
        self.fmt: FixedFormat = fmt
        self.angle_width: int = angle_width
        self.delta: float = delta
        self.x: list[float] = [float(value) for value in x]
        self.f: list[float] = function.evaluate(self.x)
        self.tree: list[float] = [0.0] * (2 * function.size)
        self.banks: dict[str, dict[int, int]] = {
            name: {} for name in ("x", "f", "f_sign", "f_angles")
        }
        self._write_leaves(range(function.size))
        for node in range(function.size - 1, 0, -1):
            self._refresh(node)

    def _refresh(self, node: int) -> None:
        """重算单个内部节点的平方部分和并刷新 RY 半角缓存。"""
        self.tree[node] = self.tree[2 * node] + self.tree[2 * node + 1]
        total = self.tree[node]
        angle = 0.0 if total == 0 else 2 * math.acos(math.sqrt(self.tree[2 * node] / total))
        self.banks["f_angles"][node - 1] = round(
            angle * (1 << self.angle_width) / (2 * math.pi)
        ) % (1 << self.angle_width)

    def _write_leaves(self, rows: Iterable[int]) -> None:
        """把指定方程行的变量字、函数字与符号位写入 bank，并更新叶平方和。"""
        for row in rows:
            self.banks["x"][row] = self.fmt.encode(self.x[row])
            self.banks["f"][row] = self.fmt.encode(self.f[row])
            self.banks["f_sign"][row] = int(self.f[row] < 0)
            self.tree[self.function.size + row] = self.f[row] ** 2

    @property
    def norm_f(self) -> float:
        """当前的 ||F(x)||_2；树根保存各 f_i(x) 的平方和。"""
        return math.sqrt(self.tree[1])

    def update(self, delta_x: Sequence[float]) -> NewtonPatch:
        """x <- x + delta_x；只重算依赖被改分量的方程（论文 III D 的局部更新）。

        Args:
            delta_x: 与 x 等长的增量向量；零分量不触发任何重算。

        Returns:
            NewtonPatch: 记录被改分量、被重算方程与更新后 ||F(x)|| 的更新结果。
        """
        changed = [j for j, value in enumerate(delta_x) if value]
        if not changed:
            return NewtonPatch((), (), self.norm_f)
        rows = [
            row
            for row, indices in enumerate(self.function.pattern)
            if any(j in indices for j in changed)
        ]
        for j in changed:
            self.x[j] += delta_x[j]
        for row in rows:
            self.f[row] = self.function.functions[row](
                *(self.x[j] for j in self.function.pattern[row])
            )
        self._write_leaves(rows)
        affected: set[int] = set()
        for row in rows:
            node = self.function.size + row
            while node > 1:
                node //= 2
                affected.add(node)
        for node in sorted(affected, reverse=True):
            self._refresh(node)
        return NewtonPatch(tuple(changed), tuple(rows), self.norm_f)

    def pattern_memories(self) -> dict[str, dict[int, int]]:
        """O_f1 静态表：每行非零列的全置换扩张及其逆（fuse(row, slot) 位序）。

        Returns:
            dict[str, dict[int, int]]: pattern_forward 与 pattern_inverse 两张
            静态模式表，按 fuse(row, slot) 位序编址。
        """
        size, width = self.function.size, self.function.width
        forward: dict[int, int] = {}
        inverse: dict[int, int] = {}
        for row in range(size):
            columns = list(self.function.pattern[row])
            rest = [c for c in range(size) if c not in columns]
            for index, column in enumerate(columns + rest):
                forward[(index << width) | row] = column
                inverse[(index << width) | column] = row
        return {"pattern_forward": forward, "pattern_inverse": inverse}

    def memories(self) -> dict[str, dict[int, int]]:
        """导出全部 QRAM 存储表。

        Returns:
            dict[str, dict[int, int]]: 四个数据 bank（``x``、``f``、``f_sign``、
            ``f_angles``）的副本，外加 ``pattern_memories`` 生成的两张模式表；
            所有表均为新建副本，修改不影响树的内部状态。
        """
        result = {name: dict(bank) for name, bank in self.banks.items()}
        result.update(self.pattern_memories())
        return result


@dataclass(frozen=True)
class NewtonPatch:
    """一次局部更新的记录：被改分量、被重算方程与更新后的 ||F(x)||。"""

    changed: tuple[int, ...]
    recomputed: tuple[int, ...]
    norm_f: float


def _compile_row(
    function: NewtonFunction, row: int, fmt: FixedFormat
) -> tuple[CompiledFunction, tuple[str, ...]]:
    """把单个方程编译为可逆算术模块，并返回其形参名元组。"""
    from pyqecclang.infrastructure.mathfunc import compile_function

    parameters = tuple(inspect.signature(function.functions[row]).parameters)
    return compile_function(
        function.functions[row],
        inputs={name: "real" for name in parameters},
        fmt=fmt,
    ), parameters


def _compile_difference(fmt: FixedFormat, scale: float) -> CompiledFunction:
    """编译差商函数 (y1 - y2) / scale 的可逆算术模块。"""
    from pyqecclang.infrastructure.mathfunc import compile_function

    source = f"def difference(y1, y2):\n    return (y1 - y2) / {scale!r}\n"
    return compile_function(source, inputs={"y1": "real", "y2": "real"}, fmt=fmt)


def newton_fd_entry(
    function: NewtonFunction, fmt: FixedFormat, *, delta: float, max_scale: float
) -> Operation:
    """O_A2：``|row, slot>|0> → |row, slot>|A_{row,slot}>``，A = F'/(Δ·scale) 语义。

    每个方程按其稀疏依赖静态特化（arity 个变量）；被 slot 选中的变量字加
    delta 后求 f_row，差商经编译算术写入 value；status 传播算术失败。
    所有工作位（变量字、扰动字、上下函数值、模式/x 查询）XOR 复净。

    Args:
        function: 被编码的稀疏非线性系统。
        fmt: 变量字与函数值的定点格式。
        delta: 差分步长，须在 fmt 下精确可表示且非零。
        max_scale: 差商缩放上界；delta 与其乘积也须可精确表示。

    Returns:
        Operation: 按 (row, slot) 查询差商矩阵元素的 QRAM 型 oracle 操作。
    """
    scale = delta * max_scale
    if fmt.encode(delta) == 0 or fmt.encode(scale) == 0:
        raise ValidationError("delta 与 delta*max_scale 必须在定点格式下可精确表示")
    row_modules = {row: _compile_row(function, row, fmt) for row in range(function.size)}
    parameters = {row: names for row, (_, names) in row_modules.items()}
    helper = _compile_difference(fmt, scale)
    width, delta_raw = function.width, fmt.encode(delta)
    slot_width = max(1, (function.sparsity - 1).bit_length())
    x_bank = qram_database(width, fmt.width)
    pattern_bank = qram_database(2 * width, width)
    b = Builder(
        _name("newton_fd_entry", function.pattern, fmt, delta, max_scale),
        {
            "row": Bits(width),
            "slot": Bits(slot_width),
            "value": Bits(fmt.width),
            "status": Bits(2),
        },
        resources_for(
            ("x", x_bank.operation),
            ("pattern", pattern_bank.operation),
            *(("f" + str(i), cast("Operation", op.operation)) for i, (op, _) in sorted(row_modules.items())),
            ("fd", cast("Operation", helper.operation)),
        ),
        attributes={
            "algorithm": "quantum_newton_fd_jacobian",
            "jacobian_source": "finite_difference",
            "delta": delta,
            "scale": max_scale,
            "correctness": "pending",
        },
    )
    row_shadow = b.local("row_shadow", Bits(width))
    b.xor(b["row"], row_shadow)
    for row in range(function.size):
        arity = len(function.pattern[row])
        with b.control(row_shadow, row):
            vars_: list[Ref] = []
            originals: list[Ref] = []
            shifted: list[Ref] = []
            selectors: list[Ref] = []
            for k in range(arity):
                index = b.local(f"r{row}_idx_{k}", Bits(slot_width))
                for bit in range(slot_width):
                    if (k >> bit) & 1:
                        b.x(index[bit])
                var = b.local(f"r{row}_var_{k}", Bits(width))
                invoke(
                    b,
                    pattern_bank.operation,
                    "pattern",
                    address=fuse(b["row"], index),
                    data=var,
                )
                word = b.local(f"r{row}_x_{k}", Bits(fmt.width))
                invoke(b, x_bank.operation, "x", address=var, data=word)
                perturbed = b.local(f"r{row}_s_{k}", Bits(fmt.width))
                b.xor(word, perturbed)
                with b.control(b["slot"], k):
                    b.add_const(perturbed.reinterpret("uint"), delta_raw)
                vars_.append(var)
                originals.append(word)
                shifted.append(perturbed)
                selectors.append(index)
            up = b.local(f"r{row}_up", Bits(fmt.width))
            down = b.local(f"r{row}_down", Bits(fmt.width))
            up_status = b.local(f"r{row}_up_status", Bits(2))
            down_status = b.local(f"r{row}_down_status", Bits(2))
            names = parameters[row]
            shifted_args = {names[k]: shifted[k] for k in range(arity)}
            plain_args = {names[k]: originals[k] for k in range(arity)}
            b.call(cast("Operation", row_modules[row][0].operation), **shifted_args, out=up, status=up_status)  # type: ignore[arg-type]
            b.call(cast("Operation", row_modules[row][0].operation), **plain_args, out=down, status=down_status)  # type: ignore[arg-type]
            b.call(cast("Operation", helper.operation), y1=up, y2=down, out=b["value"], status=b["status"])
            # XOR 语义复净（逆序）：重复调用同参数把 out/status 异或回零。
            b.call(cast("Operation", row_modules[row][0].operation), **plain_args, out=down, status=down_status)  # type: ignore[arg-type]
            b.call(cast("Operation", row_modules[row][0].operation), **shifted_args, out=up, status=up_status)  # type: ignore[arg-type]
            for k in reversed(range(arity)):
                with b.control(b["slot"], k):
                    b.add_const(
                        shifted[k].reinterpret("uint"), (-delta_raw) % (1 << fmt.width)
                    )
                b.xor(originals[k], shifted[k])
                invoke(b, x_bank.operation, "x", address=vars_[k], data=originals[k])
                invoke(
                    b,
                    pattern_bank.operation,
                    "pattern",
                    address=fuse(b["row"], selectors[k]),
                    data=vars_[k],
                )
                for bit in range(slot_width):
                    if (k >> bit) & 1:
                        b.x(selectors[k][bit])
    b.xor(b["row"], row_shadow)
    return b.finish()


def newton_b_preparation(function: NewtonFunction, angle_width: int) -> StatePreparation:
    """O_b（论文 Eq. 13）：``|b> = Σ_i f_i(x)|i>/C_b`` 的符号残差角树制备。

    Args:
        function: 被编码的稀疏非线性系统，决定树宽与符号表编址。
        angle_width: 残差角树的旋转角量化位宽。

    Returns:
        StatePreparation: 符号残差态 ``|b>`` 的制备句柄。
    """
    n = function.width
    prep = qram_state_prep(n, angle_width)
    sign = qram_database(n, 1)
    b = Builder(
        _name("newton_b", n, angle_width),
        {"target": Bits(n + 1), "work": Bits(n + angle_width + 1)},
        resources_for(("prep", prep.operation), ("sign", sign.operation)),
        attributes={"algorithm": "newton_signed_residual_tree", "correctness": "pending"},
    )
    invoke(b, prep.operation, "prep", target=b["target"][:n], work=b["work"][: n + angle_width])
    flag = b["work"][n + angle_width :]
    invoke(b, sign.operation, "sign", address=b["target"][:n], data=flag)
    b.z(flag)
    invoke(b, sign.operation, "sign", address=b["target"][:n], data=flag)
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))
