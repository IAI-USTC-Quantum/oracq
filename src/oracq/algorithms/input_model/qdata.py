"""KP/QFVM style quantum data structures: QVector (qsample) with a squared-norm tree, and the sample-and-query QMatrix.

QVector implements the data structure shared by Kerenidis-Prakash (arXiv:1603.08675
Thm 5.1 + Appendix A.1) and QFVM (arXiv:2102.03557 Eqs. 15–22): a binary tree whose
leaves store squared components rolled up layer by layer, with internal nodes
caching RY rotation angle words; the quantum side prepares the normalized state
layer by layer, with two QRAM queries per layer. Signs are written through an
independent 1-bit bank by phase kickback. On top of this, QMatrix implements
sample-and-query access (arXiv:1704.04992 §I.1): the entries bank supports queries
in arbitrary superposition, the row trees give Ũ|i⟩|0⟩→|i⟩|Ā_i⟩, and the row-norm
root tree gives Ṽ|0⟩|j⟩→|Ã⟩|j⟩.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Sequence
from typing import cast

from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.input_model.oracles import StatePreparation, XorDatabase, annotate
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import QRAM, Bits, Ref, ValidationError
from oracq.infrastructure.qmem import QMem, QPtr


def _angle_word(left_square: float, total_square: float, angle_width: int) -> int:
    """Layer rotation angle word: θ = 2·acos(√(S_left/S_node)), same convention as qram_state_angles."""
    if total_square <= 0:
        return 0
    ratio = math.sqrt(max(0.0, min(1.0, left_square / total_square)))
    return round(2 * math.acos(ratio) * (1 << angle_width) / (2 * math.pi)) % (1 << angle_width)


def _build_tree(
    squares: Sequence[float], angle_width: int
) -> tuple[list[float], dict[int, int]]:
    """Build the 1-based heap tree and the internal-node angle word table from the leaf square array (address = node-1)."""
    n = len(squares)
    tree = [0.0] * (2 * n)
    tree[n:] = [float(square) for square in squares]
    angles = {}
    for node in range(n - 1, 0, -1):
        tree[node] = tree[2 * node] + tree[2 * node + 1]
        angles[node - 1] = _angle_word(tree[2 * node], tree[node], angle_width)
    return tree, angles


def _check_shape(values: Sequence[object], what: str) -> int:
    """Validate that the sequence length is a power of two no smaller than 2 and return that length."""
    n = len(values)
    if n < 2 or n & (n - 1):
        raise ValidationError(f"{what} length must be a power of two no smaller than 2")
    return n


def _tree_layers(
    builder: Builder, cell_of: Callable[[Ref], QPtr], target: Ref, angle: Ref, width: int
) -> None:
    """Walk the tree layer by layer: cell_of(prefix) returns the tree node pointer at that layer, with the prefix taken from the already-fixed high bits of the target register."""
    for depth in range(width):
        bit = width - 1 - depth
        offset = (1 << depth) - 1
        cell = cell_of(target[bit + 1 :]) + offset
        cell.load(angle)
        for k in range(angle.width):
            with builder.control(angle[k]):
                builder.ry(target[bit], 2 * math.pi * (1 << k) / (1 << angle.width))
        cell.load(angle)


class QVector:
    """qsample vector: squared-norm binary tree + angle-word/sign banks, supporting O(log) classical local updates."""

    def __init__(
        self,
        values: Iterable[float],
        *,
        fmt: FixedFormat | None = None,
        angle_width: int = 8,
        name: str = "qvec",
    ) -> None:
        """Initialize the qsample vector and build the squared-norm tree.

        Args:
            values: Sequence of components; the length must be a power of two no smaller than 2.
            fmt: Fixed-point format of the components; defaults to ``FixedFormat(10, 5)``.
            angle_width: Bit width of the rotation angle word, range 1..64.
            name: Naming prefix of the QRAM banks.

        Raises:
            ValidationError: The component length is invalid or angle_width is out of range.
        """
        self.fmt: FixedFormat = fmt or FixedFormat(10, 5)
        if not 1 <= angle_width <= 64:
            raise ValidationError("angle_width must be in the range 1..64")
        values = tuple(float(v) for v in values)
        self.width: int = _check_shape(values, "QVector").bit_length() - 1
        self.angle_width: int = angle_width
        self.name: str = name
        self.words: list[int] = [self.fmt.encode(v) for v in values]
        self._refresh()

    def _decoded(self, index: int) -> float:
        """Decode the fixed-point value of the component at the given index."""
        return self.fmt.decode(self.words[index])

    def _refresh(self) -> None:
        """Rebuild the squared-norm tree and the angle word table from the current word values."""
        self.tree: list[float]
        self.angles: dict[int, int]
        self.tree, self.angles = _build_tree(
            [self._decoded(i) ** 2 for i in range(1 << self.width)], self.angle_width
        )

    def signs(self) -> dict[int, int]:
        """Sign bank content; when all signs are zero the signed preparation can be omitted.

        Returns:
            dict[int, int]: Mapping from indices to sign bits, containing only negative-component entries.
        """
        return {i: int(self._decoded(i) < 0) for i in range(1 << self.width) if self._decoded(i) < 0}

    def update(self, index: int, value: float) -> dict[int, int]:
        """Single-point update: only the sums and angle words along the leaf-to-root path are recomputed; returns the changes to the angle word bank.

        Args:
            index: Index of the component to update, range 0..2^width-1.
            value: New component value, which must be encodable by the fixed-point format.

        Returns:
            dict[int, int]: Changes to the angle word bank, keyed by the address, i.e. the tree node number minus one.
        """
        if not 0 <= index < 1 << self.width:
            raise ValidationError("QVector update address is out of range")
        self.words[index] = self.fmt.encode(float(value))
        node = (1 << self.width) + index
        self.tree[node] = self._decoded(index) ** 2
        changed: dict[int, int] = {}
        while node > 1:
            node //= 2
            self.tree[node] = self.tree[2 * node] + self.tree[2 * node + 1]
            changed[node - 1] = _angle_word(self.tree[2 * node], self.tree[node], self.angle_width)
        self.angles.update(changed)
        return changed

    @property
    def norm(self) -> float:
        """Euclidean norm of the vector: the square root of the squared-norm tree root (computed from the fixed-point decoded values)."""
        return math.sqrt(self.tree[1])

    def amplitudes(self) -> list[float]:
        """Normalized expected amplitudes (including quantization and signs), as the classical reference.

        Returns:
            list[float]: Amplitude table of length 2^width after division by the Euclidean norm.
        """
        norm = self.norm
        if norm == 0:
            raise ValidationError("the zero vector has no normalized state")
        return [self._decoded(i) / norm for i in range(1 << self.width)]

    def snapshot(self, *, signed: bool = False) -> dict[str, dict[int, int]]:
        """Export a snapshot of the QRAM bank contents, for ``simulate`` and the like to bind data by name.

        Args:
            signed: When True, attach the sign bank, containing only negative-component entries.

        Returns:
            dict: Keys are ``{name}_angles`` (address is the tree node number
            minus one), plus the optional ``{name}_sign``.
        """
        banks = {f"{self.name}_angles": dict(self.angles)}
        if signed:
            banks[f"{self.name}_sign"] = self.signs()
        return banks

    def preparation(self, *, signed: bool = False, name: str | None = None) -> StatePreparation:
        """Layer-by-layer tree preparation addressed via QMem; when signed, the phase is written via Load→Z→inverse Load.

        Args:
            signed: When True, attach the phase writing from the sign bank.
            name: Name of the generated operation; derived from the vector name and bit widths by default.

        Returns:
            StatePreparation: View preparing the normalized state via layer-by-layer QRAM queries.
        """
        width, angle_width = self.width, self.angle_width
        angles_name, sign_name = f"{self.name}_angles", f"{self.name}_sign"
        resources = {angles_name: QRAM(width, angle_width)}
        work_width = angle_width + (1 if signed else 0)
        if signed:
            resources[sign_name] = QRAM(width, 1)
        b = Builder(
            name or f"{self.name}_prep_{width}_{angle_width}",
            {"target": Bits(width), "work": Bits(work_width)},
            resources,
        )
        angle = b["work"][:angle_width]
        _tree_layers(b, lambda p: cast("QPtr", QMem(b, angles_name)[p]), b["target"], angle, width)
        if signed:
            flag = b["work"][angle_width:]
            sign = QMem(b, sign_name)
            cast("QPtr", sign[b["target"]]).load(flag)
            b.z(flag)
            cast("QPtr", sign[b["target"]]).load(flag)
        return StatePreparation(
            annotate(
                b.finish(),
                "state_prep_isometry",
                zero_input=True,
                clean_work=True,
                implementation="qram_rotation_tree",
                qram_queries=2 * width + (2 if signed else 0),
            )
        )


class QMatrix:
    """Sample-and-query matrix: entry query + row qsample Ũ + row-norm root tree Ṽ; entries are restricted to nonnegative values."""

    def __init__(
        self,
        matrix: Sequence[Sequence[float]],
        *,
        fmt: FixedFormat | None = None,
        angle_width: int = 8,
        name: str = "qmat",
    ) -> None:
        """Initialize the sample-and-query matrix and build the row trees and the root tree.

        Args:
            matrix: Nonnegative real matrix; both the row count and the column count must be powers of two no smaller than 2.
            fmt: Fixed-point format of the entries; defaults to ``FixedFormat(10, 5)``.
            angle_width: Bit width of the rotation angle word, range 1..64.
            name: Naming prefix of the QRAM banks.

        Raises:
            ValidationError: The matrix shape is invalid, or an entry falls outside the nonnegative range of the fixed-point format.
        """
        self.fmt: FixedFormat = fmt or FixedFormat(10, 5)
        if not 1 <= angle_width <= 64:
            raise ValidationError("angle_width must be in the range 1..64")
        rows = _check_shape(matrix, "QMatrix rows")
        cols = _check_shape(matrix[0], "QMatrix columns")
        if any(len(row) != cols for row in matrix):
            raise ValidationError("QMatrix rows must all have the same length")
        maximum = self.fmt.decode((1 << (self.fmt.width - 1)) - 1)
        for row in matrix:
            for value in row:
                if not 0 <= float(value) <= maximum:
                    raise ValidationError("QMatrix entries must lie within the nonnegative range of the fixed-point format")
        self.rows: int = rows.bit_length() - 1
        self.cols: int = cols.bit_length() - 1
        self.angle_width: int = angle_width
        self.name: str = name
        self.words: list[list[int]] = [[self.fmt.encode(float(v)) for v in row] for row in matrix]
        self._refresh()

    def _row_squares(self, index: int) -> list[float]:
        """Decode the squared values of all entries in the given row."""
        return [self.fmt.decode(w) ** 2 for w in self.words[index]]

    def _refresh(self) -> None:
        """Rebuild the row trees, the row angle word table, and the root tree from the current word values."""
        self.row_trees: list[list[float]]
        self.row_angles: dict[int, int]
        self.row_trees, self.row_angles = [], {}
        for i in range(1 << self.rows):
            tree, angles = _build_tree(self._row_squares(i), self.angle_width)
            self.row_trees.append(tree)
            for address, word in angles.items():
                self.row_angles[i * (1 << self.cols) + address] = word
        squares = [self.row_trees[i][1] for i in range(1 << self.rows)]
        self.root_tree: list[float]
        self.root_angles: dict[int, int]
        self.root_tree, self.root_angles = _build_tree(squares, self.angle_width)

    def update(self, row: int, column: int, value: float) -> dict[str, dict[int, int]]:
        """Single-entry update: recompute the affected row tree and the root tree path; returns the changes to each bank.

        Args:
            row: Row index of the entry to update, range 0..2^rows-1.
            column: Column index of the entry to update, range 0..2^cols-1.
            value: New entry value, which must lie within the nonnegative range of the fixed-point format.

        Returns:
            dict[str, dict[int, int]]: Changes to the three banks ``entries``,
            ``row_angles`` and ``root_angles``.
        """
        if not (0 <= row < 1 << self.rows and 0 <= column < 1 << self.cols):
            raise ValidationError("QMatrix update address is out of range")
        maximum = self.fmt.decode((1 << (self.fmt.width - 1)) - 1)
        if not 0 <= float(value) <= maximum:
            raise ValidationError("QMatrix entries must lie within the nonnegative range of the fixed-point format")
        self.words[row][column] = self.fmt.encode(float(value))
        self.row_trees[row], angles = _build_tree(self._row_squares(row), self.angle_width)
        changed: dict[str, dict[int, int]] = {
            "entries": {row * (1 << self.cols) + column: self.words[row][column]},
            "row_angles": {},
            "root_angles": {},
        }
        base = row * (1 << self.cols)
        for address, word in angles.items():
            if self.row_angles.get(base + address) != word:
                changed["row_angles"][base + address] = word
            self.row_angles[base + address] = word
        squares = [self.row_trees[i][1] for i in range(1 << self.rows)]
        self.root_tree, self.root_angles = _build_tree(squares, self.angle_width)
        node = (1 << self.rows) + row
        while node > 1:
            node //= 2
            changed["root_angles"][node - 1] = self.root_angles[node - 1]
        return changed

    @property
    def frobenius(self) -> float:
        """Frobenius norm of the matrix: the square root of the root tree root, i.e. the square root of the sum of squared entries."""
        return math.sqrt(self.root_tree[1])

    def row_norm(self, index: int) -> float:
        """Return the Euclidean norm of row ``index`` (the square root of the row tree root).

        Args:
            index: Row index, range 0..2^rows-1.

        Returns:
            float: Square root of the sum of squared entries of that row.
        """
        return math.sqrt(self.row_trees[index][1])

    def row_amplitudes(self, index: int) -> list[float]:
        """Normalized expected amplitudes of row index (the reference ground truth for the row qsample).

        Args:
            index: Row index, range 0..2^rows-1.

        Returns:
            list[float]: Amplitude table of that row after division by the row norm.
        """
        return [self.fmt.decode(w) / self.row_norm(index) for w in self.words[index]]

    def user_amplitudes(self) -> list[float]:
        """Normalized expected amplitudes of the root tree (the reference ground truth for the user superposition prepared by Ṽ).

        Returns:
            list[float]: The i-th entry is the ratio of the i-th row norm to the Frobenius norm.
        """
        total = self.root_tree[1]
        return [math.sqrt(self.row_trees[i][1] / total) for i in range(1 << self.rows)]

    def snapshot(self) -> dict[str, dict[int, int]]:
        """Export a snapshot of the entries bank and the row/root angle banks.

        Returns:
            dict: Keys are ``entries`` (row-major flattened addresses to
            fixed-point words), ``row_angles`` and ``root_angles``.
        """
        return {
            "entries": {
                i * (1 << self.cols) + j: self.words[i][j]
                for i in range(1 << self.rows)
                for j in range(1 << self.cols)
            },
            "row_angles": dict(self.row_angles),
            "root_angles": dict(self.root_angles),
        }

    def query(self) -> XorDatabase:
        """Entry XOR query: the word at address (row, column) is XORed into data, addressed by a two-dimensional pointer.

        Returns:
            XorDatabase: Query view backed by the row-major flattened entries bank.
        """
        r, c = self.rows, self.cols
        b = Builder(
            f"{self.name}_query_{r}_{c}",
            {"address": Bits(r + c), "data": Bits(self.fmt.width)},
            {"entries": QRAM(r + c, self.fmt.width)},
        )
        entries = QMem(b, "entries", shape=(1 << r, 1 << c))
        cast("QPtr", entries[b["address"][c:], b["address"][:c]]).load(b["data"])
        return XorDatabase(annotate(b.finish(), "database_xor", implementation="qmem_2d"))

    def row_preparation(self) -> Operation:
        """Ũ: prepare the item register from the zero state into the normalized row vector of row i; the row trees are addressed two-dimensionally by (row, prefix).

        Returns:
            Operation: The unitary operation implementing Ũ.
        """
        r, c, aw = self.rows, self.cols, self.angle_width
        b = Builder(
            f"{self.name}_row_prep_{r}_{c}_{aw}",
            {"row": Bits(r), "item": Bits(c), "work": Bits(aw)},
            {"row_angles": QRAM(r + c, aw)},
        )
        mem = QMem(b, "row_angles", shape=(1 << r, 1 << c))
        _tree_layers(b, lambda p: cast("QPtr", mem[b["row"], p]), b["item"], b["work"], c)
        return annotate(
            b.finish(),
            "unitary",
            implementation="qram_row_trees",
            qram_queries=2 * c,
        )

    def amplitude_preparation(self) -> Operation:
        """Ṽ: prepare the zero-state row register into the user distribution Ã, with the item register passed through.

        Returns:
            Operation: The unitary operation implementing Ṽ.
        """
        r, aw = self.rows, self.angle_width
        b = Builder(
            f"{self.name}_amp_prep_{r}_{aw}",
            {"row": Bits(r), "item": Bits(self.cols), "work": Bits(aw)},
            {"root_angles": QRAM(r, aw)},
        )
        _tree_layers(b, lambda p: cast("QPtr", QMem(b, "root_angles")[p]), b["row"], b["work"], r)
        return annotate(
            b.finish(),
            "unitary",
            implementation="qram_root_tree",
            qram_queries=2 * r,
        )

    def row_state_prep(self, user: int) -> StatePreparation:
        """One-dimensional tree preparation of the classical row user, with the row number folded into a constant address.

        Args:
            user: Classical row number, range 0..2^rows-1.

        Returns:
            StatePreparation: 1D tree preparation view of that row's normalized row vector.
        """
        c, aw = self.cols, self.angle_width
        if not 0 <= user < 1 << self.rows:
            raise ValidationError("row number is out of range")
        b = Builder(
            f"{self.name}_row_{user}_prep_{c}_{aw}",
            {"target": Bits(c), "work": Bits(aw)},
            {"row_angles": QRAM(self.rows + c, aw)},
        )
        mem = QMem(b, "row_angles", shape=(1 << self.rows, 1 << c))
        _tree_layers(b, lambda p: cast("QPtr", mem[user, p]), b["target"], b["work"], c)
        return StatePreparation(
            annotate(
                b.finish(),
                "state_prep_isometry",
                zero_input=True,
                clean_work=True,
                implementation="qram_row_tree",
                qram_queries=2 * c,
            )
        )
