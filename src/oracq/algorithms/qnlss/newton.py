"""Implementation of the quantum Newton method (arXiv:2109.08470): the M_F data structure and the finite-difference Jacobian oracle.

The algorithm runs Newton iterations A dx = -F(x) round by round, with the linear solve and
data loading using quantum building blocks:

- ``NewtonTree`` is the paper's M_F: leaves store x_i and f_i(x), internal nodes store squared
  partial sums, and the angle cache is written into the bank under the RY half-angle convention;
  a single-point update only recomputes equations depending on that component.
- ``newton_fd_entry`` is O_A2 (paper Eq. 15-17): given (row, slot), it looks up the variable
  index via the pattern database O_f1 and fetches the variable value via the x database O_M1,
  adds delta to the selected variable, and calls the reversibly compiled f_row (O_f2) twice;
  the difference quotient (f(x+delta e_k)-f(x))/(delta*scale) is produced by compiled
  arithmetic; all work bits are cleaned via XOR.
- ``newton_b_preparation`` is O_b: a signed-residual angle tree prepares ``|b> = Σ f_i |i> / C_b``.

The linear solve (Hermitian dilation + sparse QLSS), l_inf tomography, and norm recovery are
assembled by ``newton_linear_step``; semantic validation of all building blocks is in
tests/core/test_newton.py.
"""

from __future__ import annotations

import inspect
import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    annotate,
    invoke,
    qram_database,
    qram_state_prep,
    resources_for,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Ref, ValidationError, fuse
from oracq.infrastructure.mathfunc import CompiledFunction


def _log2_exact(value: int) -> int:
    """Validate that value is a positive power of 2 and return its base-2 logarithm."""
    if value < 1 or value & (value - 1):
        raise ValidationError("Newton method dimensions must be a power of 2")
    return value.bit_length() - 1


@dataclass(frozen=True)
class NewtonFunction:
    """Sparse nonlinear system F(x)=0: one pure-Python numeric function per equation with its dependency indices."""

    functions: tuple[Callable[..., float], ...]
    pattern: tuple[tuple[int, ...], ...]

    @classmethod
    def declare(
        cls, equations: Sequence[tuple[Sequence[int], Callable[..., float]]]
    ) -> NewtonFunction:
        """Declare a sparse nonlinear system from a sequence of (dependency indices, numeric function) pairs.

        Args:
            equations: A sequence of ``(indices, function)`` pairs; indices are the variable
                indices the equation depends on (in one-to-one correspondence with the
                positional parameters of function), and function is a pure-Python numeric
                function.

        Returns:
            NewtonFunction: The system description after dimension and sparse-dependency validation.

        Raises:
            ValidationError: The number of equations is not a power of 2, or some equation has duplicate or out-of-range indices.
            OSError: Source code is unavailable for some function (``inspect.getsource``; required for quantum compilation).
        """
        pattern = tuple(tuple(indices) for indices, _ in equations)
        functions = tuple(function for _, function in equations)
        size = len(equations)
        _log2_exact(size)
        for row, indices in enumerate(pattern):
            if len(set(indices)) != len(indices) or any(not 0 <= j < size for j in indices):
                raise ValidationError(f"Equation {row} has invalid sparse dependencies")
            inspect.getsource(functions[row])
        return cls(functions, pattern)

    @property
    def size(self) -> int:
        """Number of equations; the system is square, so this also equals the number of variables."""
        return len(self.functions)

    @property
    def width(self) -> int:
        """Bit width of the variable (equation) index register, i.e. log2(size); declaration guarantees size is a power of 2."""
        return _log2_exact(self.size)

    @property
    def sparsity(self) -> int:
        """Sparsity: the maximum number of variables a single equation depends on."""
        return max(len(row) for row in self.pattern)

    def evaluate(self, x: Sequence[float]) -> list[float]:
        """Classically evaluate F(x).

        Args:
            x: A real vector of length ``size``.

        Returns:
            list: The equation function values f_i(x) in declaration order.
        """
        return [
            function(*(x[j] for j in indices))
            for function, indices in zip(self.functions, self.pattern, strict=True)
        ]

    def jacobian(self, x: Sequence[float], delta: float) -> list[dict[int, float]]:
        """Finite-difference Jacobian (classical mirror, serving as an independent reference for the quantum oracle).

        Args:
            x: A real vector of length ``size``; the base point of the differences.
            delta: Forward-difference step; a nonzero real number.

        Returns:
            list[dict[int, float]]: One sparse dictionary per row, keyed by column index with difference-quotient partials as values.
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
    """The M_F data structure (classical side). The heap is 1-based: root at 1, leaves N..2N-1."""

    def __init__(
        self,
        function: NewtonFunction,
        x: Sequence[float],
        *,
        fmt: FixedFormat | None = None,
        angle_width: int = 10,
        delta: float = 2.0**-6,
    ) -> None:
        """Build the M_F tree from the system description and initial values, filling all storage tables bottom-up.

        Args:
            function: The sparse nonlinear system being encoded.
            x: Real initial value vector with length matching ``function.size``.
            fmt: Fixed-point format used by the numeric banks; defaults to ``FixedFormat(12, 6)`` when omitted.
            angle_width: Bit width of the RY angle cache.
            delta: Difference step; must be exactly representable and nonzero under ``fmt``.

        Raises:
            ValidationError: delta is not exactly representable, or the initial value dimension does not match the system.
        """
        fmt = fmt or FixedFormat(12, 6)
        if fmt.encode(delta) == 0:
            raise ValidationError("delta must be exactly representable and nonzero in the fixed-point format")
        if len(x) != function.size:
            raise ValidationError("Initial value dimension does not match the system")
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
        """Recompute the squared partial sum of a single internal node and refresh the RY half-angle cache."""
        self.tree[node] = self.tree[2 * node] + self.tree[2 * node + 1]
        total = self.tree[node]
        angle = 0.0 if total == 0 else 2 * math.acos(math.sqrt(self.tree[2 * node] / total))
        self.banks["f_angles"][node - 1] = round(
            angle * (1 << self.angle_width) / (2 * math.pi)
        ) % (1 << self.angle_width)

    def _write_leaves(self, rows: Iterable[int]) -> None:
        """Write the variable word, function word, and sign bit of the given equation rows into the banks, and update the leaf squared sums."""
        for row in rows:
            self.banks["x"][row] = self.fmt.encode(self.x[row])
            self.banks["f"][row] = self.fmt.encode(self.f[row])
            self.banks["f_sign"][row] = int(self.f[row] < 0)
            self.tree[self.function.size + row] = self.f[row] ** 2

    @property
    def norm_f(self) -> float:
        """The current ||F(x)||_2; the tree root stores the sum of squares of the f_i(x)."""
        return math.sqrt(self.tree[1])

    def update(self, delta_x: Sequence[float]) -> NewtonPatch:
        """x <- x + delta_x; only equations depending on changed components are recomputed (the local update of paper III D).

        Args:
            delta_x: Increment vector of the same length as x; zero components trigger no recomputation.

        Returns:
            NewtonPatch: The update result recording changed components, recomputed equations, and the updated ||F(x)||.
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
        """O_f1 static tables: the full-permutation expansion of each row's nonzero columns and its inverse (fuse(row, slot) bit order).

        Returns:
            dict[str, dict[int, int]]: The two static pattern tables pattern_forward and
            pattern_inverse, addressed in fuse(row, slot) bit order.
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
        """Export all QRAM storage tables.

        Returns:
            dict[str, dict[int, int]]: Copies of the four data banks (``x``, ``f``,
            ``f_sign``, ``f_angles``) plus the two pattern tables produced by
            ``pattern_memories``; every table is a fresh copy, and modifying them
            does not affect the tree's internal state.
        """
        result = {name: dict(bank) for name, bank in self.banks.items()}
        result.update(self.pattern_memories())
        return result


@dataclass(frozen=True)
class NewtonPatch:
    """Record of one local update: changed components, recomputed equations, and the updated ||F(x)||."""

    changed: tuple[int, ...]
    recomputed: tuple[int, ...]
    norm_f: float


def _compile_row(
    function: NewtonFunction, row: int, fmt: FixedFormat
) -> tuple[CompiledFunction, tuple[str, ...]]:
    """Compile a single equation into a reversible arithmetic module and return the tuple of its parameter names."""
    from oracq.infrastructure.mathfunc import compile_function

    parameters = tuple(inspect.signature(function.functions[row]).parameters)
    return compile_function(
        function.functions[row],
        inputs={name: "real" for name in parameters},
        fmt=fmt,
    ), parameters


def _compile_difference(fmt: FixedFormat, scale: float) -> CompiledFunction:
    """Compile the difference-quotient function (y1 - y2) / scale into a reversible arithmetic module."""
    from oracq.infrastructure.mathfunc import compile_function

    source = f"def difference(y1, y2):\n    return (y1 - y2) / {scale!r}\n"
    return compile_function(source, inputs={"y1": "real", "y2": "real"}, fmt=fmt)


def newton_fd_entry(
    function: NewtonFunction, fmt: FixedFormat, *, delta: float, max_scale: float
) -> Operation:
    """O_A2: ``|row, slot>|0> -> |row, slot>|A_{row,slot}>``, with A = F'/(delta*scale) semantics.

    Each equation is statically specialized to its sparse dependencies (arity variables); the
    variable word selected by slot gets delta added before evaluating f_row, and the difference
    quotient is written into value via compiled arithmetic; status propagates arithmetic
    failures. All work bits (variable words, perturbed words, up/down function values,
    pattern/x lookups) are cleaned via XOR.

    Args:
        function: The sparse nonlinear system being encoded.
        fmt: Fixed-point format of the variable words and function values.
        delta: Difference step; must be exactly representable and nonzero under fmt.
        max_scale: Upper bound on the difference-quotient scaling; its product with delta must also be exactly representable.

    Returns:
        Operation: A QRAM-style oracle operation returning difference-quotient matrix elements addressed by (row, slot).
    """
    scale = delta * max_scale
    if fmt.encode(delta) == 0 or fmt.encode(scale) == 0:
        raise ValidationError("delta and delta*max_scale must be exactly representable in the fixed-point format")
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
            # XOR-semantic cleanup (reverse order): repeating the call with the same arguments XORs out/status back to zero.
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
    """O_b (paper Eq. 13): signed-residual angle-tree preparation of ``|b> = Σ_i f_i(x)|i>/C_b``.

    Args:
        function: The sparse nonlinear system being encoded; determines the tree width and sign table addressing.
        angle_width: Quantization bit width of the residual angle tree's rotation angles.

    Returns:
        StatePreparation: Preparation handle of the signed residual state ``|b>``.
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
