"经典 Roe Riemann 局部更新与可查询 QRAM 数据结构；不保存矩阵元表。"

from __future__ import annotations

import math
from dataclasses import dataclass

from pyqecclang.algorithms.arithmetic import FixedFormat
from pyqecclang.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class MemoryPatch:
    """一次更新返回的不可变写入记录。

    Attributes:
        version: 应用写入之后的存储版本号；无实际写入时保持不变。
        changes: 实际发生变化的存储项，按 bank 名映射地址到整数值。
        recomputed_faces: 本次更新重算 Riemann 通量的界面编号；存储层直写时为空。
        recomputed_cells: 本次更新重算残差的单元编号；存储层直写时为空。
    """

    version: int
    changes: dict[str, dict[int, int]]
    recomputed_faces: tuple[int, ...] = ()
    recomputed_cells: tuple[int, ...] = ()


class QRAMStore:
    """带版本号的多 bank 经典整数存储，作为可查询 QRAM 的数据源。

    Attributes:
        banks: 按 bank 名组织的存储内容；每个 bank 是地址到整数值的映射。
        version: 版本号；每次产生实际写入的 apply 调用自增一。
    """

    def __init__(self):
        self.banks, self.version = {}, 0

    def apply(self, changes):
        """写入一批更新并返回实际发生的写入记录。

        与当前存储值相同（未写地址按零计）的项不会记入结果；存在实际写入时
        版本号自增。

        Args:
            changes: 写入内容，按 bank 名映射地址到整数值。

        Returns:
            MemoryPatch: 写入后的版本号与实际写入项，不含重算范围。
        """
        actual = {}
        for bank, cells in changes.items():
            target = self.banks.setdefault(bank, {})
            for address, value in cells.items():
                if target.get(address, 0) != value:
                    target[address] = value
                    actual.setdefault(bank, {})[address] = value
        if actual:
            self.version += 1
        return MemoryPatch(self.version, actual)

    def snapshot(self):
        """返回全部 bank 内容的独立副本。

        Returns:
            dict: 每个 bank 的地址到值映射的拷贝，修改快照不影响存储。
        """
        return {k: dict(v) for k, v in self.banks.items()}

    def materialize_changed(self, patch, factory):
        """当前 PySparQ 无写入接口：只重建发生变化的 bank 对象，不宣称物理局部写入。"""
        return {bank: factory(bank, self.banks[bank]) for bank in patch.changes}


def _matmul(a, b):
    return [
        [sum(x * y for x, y in zip(row, col, strict=True)) for col in zip(*b, strict=True)]
        for row in a
    ]


def _inverse3(a):
    aug = [list(row) + [float(i == j) for j in range(3)] for i, row in enumerate(a)]
    for j in range(3):
        pivot = max(range(j, 3), key=lambda i: abs(aug[i][j]))
        aug[j], aug[pivot] = aug[pivot], aug[j]
        if abs(aug[j][j]) < 1e-15:
            raise ValidationError("Roe 特征矩阵奇异")
        scale = aug[j][j]
        aug[j] = [x / scale for x in aug[j]]
        for i in range(3):
            if i != j:
                scale = aug[i][j]
                aug[i] = [x - scale * y for x, y in zip(aug[i], aug[j], strict=True)]
    return [row[3:] for row in aug]


def riemann_flux(left, right, *, gamma=1.4, entropy_delta=0.125):
    """计算单个界面上的经典 Roe 近似 Riemann 通量。

    以左右守恒状态的 Roe 平均构造特征向量矩阵 ``R``，特征值 ``u-c``、``u``、
    ``u+c`` 取熵修正后的绝对值，组装 ``|A| = R*diag(|lambda|)*R**-1`` 参与
    通量公式。

    Args:
        left: 左侧单元的守恒变量三元组，依次为密度、动量、能量。
        right: 右侧单元的守恒变量三元组，分量顺序同 ``left``。
        gamma: 比热比。
        entropy_delta: Harten 熵修正阈值；绝对值小于它的特征值改用平滑值。

    Returns:
        tuple: 三个守恒分量的数值通量 ``0.5*(f_l+f_r) - 0.5*|A|*(u_r-u_l)``。

    Raises:
        ValidationError: 任一侧密度非正，或 Roe 平均的声速平方非正。
    """
    def primitive(state):
        rho, m, e = state
        if rho <= 0:
            raise ValidationError("经典流场密度必须为正")
        u = m / rho
        p = (gamma - 1) * (e - 0.5 * m * u)
        return u, p, (e + p) / rho, [m, m * u + p, u * (e + p)]

    ul, _, hl, fl = primitive(left)
    ur, _, hr, fr = primitive(right)
    wl, wr = math.sqrt(left[0]), math.sqrt(right[0])
    u, h = (wl * ul + wr * ur) / (wl + wr), (wl * hl + wr * hr) / (wl + wr)
    c2 = (gamma - 1) * (h - 0.5 * u * u)
    if c2 <= 0:
        raise ValidationError("经典 Roe 声速平方必须为正")
    c = math.sqrt(c2)
    r = [[1, 1, 1], [u - c, u, u + c], [h - u * c, 0.5 * u * u, h + u * c]]
    vals = [
        abs(x) if abs(x) >= entropy_delta else (x * x + entropy_delta**2) / (2 * entropy_delta)
        for x in (u - c, u, u + c)
    ]
    abs_a = _matmul([[x * vals[j] for j, x in enumerate(row)] for row in r], _inverse3(r))
    jump = [y - x for x, y in zip(left, right, strict=True)]
    return tuple(
        0.5 * (fl[i] + fr[i]) - 0.5 * sum(abs_a[i][j] * jump[j] for j in range(3)) for i in range(3)
    )


class RoeFlowData:
    """周期一维网格，三分量补齐到四分量；符号叶子+平方范数/角度二叉树。"""

    def __init__(
        self,
        states,
        *,
        fmt=None,
        gamma=1.4,
        entropy_delta=0.125,
        angle_width=10,
        dx=1.0,
    ):
        fmt = fmt or FixedFormat(10, 5)
        self.states = [tuple(map(float, x)) for x in states]
        n = len(states)
        if n < 4 or n & (n - 1) or any(len(x) != 3 for x in states) or dx <= 0:
            raise ValidationError("流场需要至少四个、数量为二次幂的三分量单元和正 dx")
        self.fmt, self.gamma, self.delta = fmt, gamma, entropy_delta
        self.angle_width, self.dx, self.n = angle_width, dx, n
        self.store = QRAMStore()
        self.fluxes, self.residuals = {}, {}
        self.tree = [0.0] * (8 * n)
        self.last_patch = self.update(dict(enumerate(self.states)), initialize=True)

    def update(self, changed, *, initialize=False):
        """用给定的单元新值执行一次局部更新并写入 QRAM 存储。

        只重算受影响界面（含周期邻居）的 Riemann 通量及其相邻单元的残差；
        守恒量、符号、右端平方范数二叉树与角度各 bank 同步维护。

        Args:
            changed: 单元编号到新守恒变量三元组的映射。
            initialize: 为 True 时残差对全部单元重算，仅用于构造时的首次填充。

        Returns:
            MemoryPatch: 本次实际写入与重算范围，同时保存为 ``last_patch``。

        Raises:
            ValidationError: 单元编号越界或分量数不是三。
        """
        for cell, values in changed.items():
            if not 0 <= cell < self.n or len(values) != 3:
                raise ValidationError("流场局部更新的地址/分量无效")
        changed = {cell: tuple(map(float, values)) for cell, values in changed.items()}

        def state_at(cell):
            return changed.get(cell, self.states[cell])

        faces = sorted({face for cell in changed for face in ((cell - 1) % self.n, cell)})
        updates = {}
        for face in faces:
            updates[face] = riemann_flux(
                state_at(face),
                state_at((face + 1) % self.n),
                gamma=self.gamma,
                entropy_delta=self.delta,
            )

        def flux_at(face):
            return updates[face] if face in updates else self.fluxes[face]

        cells = sorted({cell for face in faces for cell in (face, (face + 1) % self.n)})
        if initialize:
            cells = list(range(self.n))
        changes = {
            bank: {}
            for bank in ("rho", "momentum", "energy", "rhs_values", "rhs_sign", "rhs_angles")
        }
        for cell in changed:
            for j, bank in enumerate(("rho", "momentum", "energy")):
                changes[bank][cell] = self.fmt.encode(changed[cell][j])
        ancestors = set()
        for cell in cells:
            residual = tuple(
                self.fmt.decode(
                    self.fmt.encode((flux_at((cell - 1) % self.n)[j] - flux_at(cell)[j]) / self.dx)
                )
                for j in range(3)
            )
            self.residuals[cell] = residual
            for j, value in enumerate((*residual, 0.0)):
                address = 4 * cell + j
                changes["rhs_values"][address] = self.fmt.encode(value)
                changes["rhs_sign"][address] = int(value < 0)
                node = 4 * self.n + address
                self.tree[node] = value * value
                while node > 1:
                    node //= 2
                    ancestors.add(node)
        for node in sorted(ancestors, reverse=True):
            self.tree[node] = self.tree[2 * node] + self.tree[2 * node + 1]
            angle = (
                0
                if self.tree[node] == 0
                else 2 * math.acos(math.sqrt(self.tree[2 * node] / self.tree[node]))
            )
            changes["rhs_angles"][node - 1] = round(
                angle * (1 << self.angle_width) / (2 * math.pi)
            ) % (1 << self.angle_width)
        for cell, values in changed.items():
            self.states[cell] = values
        self.fluxes.update(updates)
        patch = self.store.apply(changes)
        self.last_patch = MemoryPatch(patch.version, patch.changes, tuple(faces), tuple(cells))
        return self.last_patch

    @property
    def rhs_norm(self):
        """残差场的 L2 范数，即平方和二叉树根节点值开平方。"""
        return math.sqrt(self.tree[1])
