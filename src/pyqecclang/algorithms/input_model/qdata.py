"""KP/QFVM 式量子数据结构：平方范数树的 QVector（qsample）与 sample-and-query 的 QMatrix。

QVector 实现 Kerenidis-Prakash（arXiv:1603.08675 Thm 5.1 + 附录 A.1）与 QFVM
（arXiv:2102.03557 式 15–22）共用的数据结构：叶子存分量平方、逐层上卷的二叉树，
内部节点缓存 RY 旋转角字；量子侧按层制备归一化态，每层两次 QRAM 查询。
符号经独立的 1 位 bank 以相位反冲写入。QMatrix 在此之上实现 sample-and-query
访问（arXiv:1704.04992 §I.1）：条目 bank 支持任意叠加查询，行树给出
Ũ|i⟩|0⟩→|i⟩|Ā_i⟩，行范数根树给出 Ṽ|0⟩|j⟩→|Ã⟩|j⟩。
"""

from __future__ import annotations

import math

from pyqecclang.algorithms.common.arithmetic import FixedFormat
from pyqecclang.algorithms.input_model.oracles import StatePreparation, XorDatabase, annotate
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import QRAM, Bits, ValidationError
from pyqecclang.infrastructure.qmem import QMem


def _angle_word(left_square, total_square, angle_width):
    """层旋转角字：θ = 2·acos(√(S_left/S_node))，与 qram_state_angles 同一约定。"""
    if total_square <= 0:
        return 0
    ratio = math.sqrt(max(0.0, min(1.0, left_square / total_square)))
    return round(2 * math.acos(ratio) * (1 << angle_width) / (2 * math.pi)) % (1 << angle_width)


def _build_tree(squares, angle_width):
    """由叶平方数组构建 1-based 堆树与内部节点角字表（地址 = node-1）。"""
    n = len(squares)
    tree = [0.0] * (2 * n)
    tree[n:] = [float(square) for square in squares]
    angles = {}
    for node in range(n - 1, 0, -1):
        tree[node] = tree[2 * node] + tree[2 * node + 1]
        angles[node - 1] = _angle_word(tree[2 * node], tree[node], angle_width)
    return tree, angles


def _check_shape(values, what):
    n = len(values)
    if n < 2 or n & (n - 1):
        raise ValidationError(f"{what} 长度必须是不少于 2 的二次幂")
    return n


def _tree_layers(builder, cell_of, target, angle, width):
    """按层走树：cell_of(prefix) 返回该层树节点指针，前缀取目标寄存器已定高位。"""
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
    """qsample 向量：平方范数二叉树 + 角字/符号 bank，支持 O(log) 经典局部更新。"""

    def __init__(self, values, *, fmt=None, angle_width=8, name="qvec"):
        self.fmt = fmt or FixedFormat(10, 5)
        if not 1 <= angle_width <= 64:
            raise ValidationError("angle_width 必须是 1..64")
        values = tuple(float(v) for v in values)
        self.width = _check_shape(values, "QVector").bit_length() - 1
        self.angle_width, self.name = angle_width, name
        self.words = [self.fmt.encode(v) for v in values]
        self._refresh()

    def _decoded(self, index):
        return self.fmt.decode(self.words[index])

    def _refresh(self):
        self.tree, self.angles = _build_tree(
            [self._decoded(i) ** 2 for i in range(1 << self.width)], self.angle_width
        )

    def signs(self):
        """符号 bank 内容；全零符号可以省略 signed 制备。"""
        return {i: int(self._decoded(i) < 0) for i in range(1 << self.width) if self._decoded(i) < 0}

    def update(self, index, value):
        """单点更新：只重算叶到根路径上的和与角字；返回角字 bank 的变化。"""
        if not 0 <= index < 1 << self.width:
            raise ValidationError("QVector 更新地址越界")
        self.words[index] = self.fmt.encode(float(value))
        node = (1 << self.width) + index
        self.tree[node] = self._decoded(index) ** 2
        changed = {}
        while node > 1:
            node //= 2
            self.tree[node] = self.tree[2 * node] + self.tree[2 * node + 1]
            changed[node - 1] = _angle_word(self.tree[2 * node], self.tree[node], self.angle_width)
        self.angles.update(changed)
        return changed

    @property
    def norm(self):
        """向量的欧几里得范数：平方范数树根节点开方（按定点解码值计算）。"""
        return math.sqrt(self.tree[1])

    def amplitudes(self):
        """归一化期望振幅（含量化与符号），供经典侧对照。"""
        norm = self.norm
        if norm == 0:
            raise ValidationError("零向量没有归一化态")
        return [self._decoded(i) / norm for i in range(1 << self.width)]

    def snapshot(self, *, signed=False):
        """导出 QRAM bank 内容快照，供 ``simulate`` 等按名绑定数据。

        Args:
            signed: 为 True 时附带符号 bank，只含负分量的条目。

        Returns:
            dict: 键为 ``{name}_angles``（地址为树节点编号减一），
            以及可选的 ``{name}_sign``。
        """
        banks = {f"{self.name}_angles": dict(self.angles)}
        if signed:
            banks[f"{self.name}_sign"] = self.signs()
        return banks

    def preparation(self, *, signed=False, name=None):
        """按层 QMem 寻址的树制备；带符号时经 Load→Z→反 Load 写入相位。"""
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
        _tree_layers(b, lambda p: QMem(b, angles_name)[p], b["target"], angle, width)
        if signed:
            flag = b["work"][angle_width:]
            sign = QMem(b, sign_name)
            sign[b["target"]].load(flag)
            b.z(flag)
            sign[b["target"]].load(flag)
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
    """sample-and-query 矩阵：条目查询 + 行 qsample Ũ + 行范数根树 Ṽ；条目限非负。"""

    def __init__(self, matrix, *, fmt=None, angle_width=8, name="qmat"):
        self.fmt = fmt or FixedFormat(10, 5)
        if not 1 <= angle_width <= 64:
            raise ValidationError("angle_width 必须是 1..64")
        rows = _check_shape(matrix, "QMatrix 行")
        cols = _check_shape(matrix[0], "QMatrix 列")
        if any(len(row) != cols for row in matrix):
            raise ValidationError("QMatrix 每行长度必须一致")
        maximum = self.fmt.decode((1 << (self.fmt.width - 1)) - 1)
        for row in matrix:
            for value in row:
                if not 0 <= float(value) <= maximum:
                    raise ValidationError("QMatrix 条目必须处于定点格式非负值域内")
        self.rows, self.cols = rows.bit_length() - 1, cols.bit_length() - 1
        self.angle_width, self.name = angle_width, name
        self.words = [[self.fmt.encode(float(v)) for v in row] for row in matrix]
        self._refresh()

    def _row_squares(self, index):
        return [self.fmt.decode(w) ** 2 for w in self.words[index]]

    def _refresh(self):
        self.row_trees, self.row_angles = [], {}
        for i in range(1 << self.rows):
            tree, angles = _build_tree(self._row_squares(i), self.angle_width)
            self.row_trees.append(tree)
            for address, word in angles.items():
                self.row_angles[i * (1 << self.cols) + address] = word
        squares = [self.row_trees[i][1] for i in range(1 << self.rows)]
        self.root_tree, self.root_angles = _build_tree(squares, self.angle_width)

    def update(self, row, column, value):
        """单条目更新：重算所在行树与根树路径；返回各 bank 的变化。"""
        if not (0 <= row < 1 << self.rows and 0 <= column < 1 << self.cols):
            raise ValidationError("QMatrix 更新地址越界")
        maximum = self.fmt.decode((1 << (self.fmt.width - 1)) - 1)
        if not 0 <= float(value) <= maximum:
            raise ValidationError("QMatrix 条目必须处于定点格式非负值域内")
        self.words[row][column] = self.fmt.encode(float(value))
        self.row_trees[row], angles = _build_tree(self._row_squares(row), self.angle_width)
        changed = {
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
    def frobenius(self):
        """矩阵的 Frobenius 范数：根树根节点开方，即全部条目平方和的平方根。"""
        return math.sqrt(self.root_tree[1])

    def row_norm(self, index):
        """返回第 ``index`` 行的欧几里得范数（行树根节点开方）。"""
        return math.sqrt(self.row_trees[index][1])

    def row_amplitudes(self, index):
        """第 index 行的归一化期望振幅（行 qsample 的对照真值）。"""
        return [self.fmt.decode(w) / self.row_norm(index) for w in self.words[index]]

    def user_amplitudes(self):
        """根树的归一化期望振幅（Ṽ 制备用户叠加态的对照真值）。"""
        total = self.root_tree[1]
        return [math.sqrt(self.row_trees[i][1] / total) for i in range(1 << self.rows)]

    def snapshot(self):
        """导出条目 bank 与行/根角度 bank 的内容快照。

        Returns:
            dict: 键为 ``entries``（行主序扁平地址到定点字）、
            ``row_angles`` 与 ``root_angles``。
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

    def query(self):
        """条目 XOR 查询：地址 (行, 列) 处的字异或进 data，二维指针寻址。"""
        r, c = self.rows, self.cols
        b = Builder(
            f"{self.name}_query_{r}_{c}",
            {"address": Bits(r + c), "data": Bits(self.fmt.width)},
            {"entries": QRAM(r + c, self.fmt.width)},
        )
        entries = QMem(b, "entries", shape=(1 << r, 1 << c))
        entries[b["address"][c:], b["address"][:c]].load(b["data"])
        return XorDatabase(annotate(b.finish(), "database_xor", implementation="qmem_2d"))

    def row_preparation(self):
        """Ũ：把条目寄存器从零态制备到第 i 行的归一化行向量；行树按 (行, 前缀) 二维寻址。"""
        r, c, aw = self.rows, self.cols, self.angle_width
        b = Builder(
            f"{self.name}_row_prep_{r}_{c}_{aw}",
            {"row": Bits(r), "item": Bits(c), "work": Bits(aw)},
            {"row_angles": QRAM(r + c, aw)},
        )
        mem = QMem(b, "row_angles", shape=(1 << r, 1 << c))
        _tree_layers(b, lambda p: mem[b["row"], p], b["item"], b["work"], c)
        return annotate(
            b.finish(),
            "unitary",
            implementation="qram_row_trees",
            qram_queries=2 * c,
        )

    def amplitude_preparation(self):
        """Ṽ：把零态行寄存器制备到用户分布 Ã，条目寄存器直通。"""
        r, aw = self.rows, self.angle_width
        b = Builder(
            f"{self.name}_amp_prep_{r}_{aw}",
            {"row": Bits(r), "item": Bits(self.cols), "work": Bits(aw)},
            {"root_angles": QRAM(r, aw)},
        )
        _tree_layers(b, lambda p: QMem(b, "root_angles")[p], b["row"], b["work"], r)
        return annotate(
            b.finish(),
            "unitary",
            implementation="qram_root_tree",
            qram_queries=2 * r,
        )

    def row_state_prep(self, user):
        """经典行 user 的一维树制备，行号折叠为常量地址。"""
        c, aw = self.cols, self.angle_width
        if not 0 <= user < 1 << self.rows:
            raise ValidationError("行号越界")
        b = Builder(
            f"{self.name}_row_{user}_prep_{c}_{aw}",
            {"target": Bits(c), "work": Bits(aw)},
            {"row_angles": QRAM(self.rows + c, aw)},
        )
        mem = QMem(b, "row_angles", shape=(1 << self.rows, 1 << c))
        _tree_layers(b, lambda p: mem[user, p], b["target"], b["work"], c)
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
