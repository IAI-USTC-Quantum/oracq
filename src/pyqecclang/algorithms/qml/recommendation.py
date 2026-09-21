"""Kerenidis-Prakash 量子推荐系统（arXiv:1603.08675）：QSVE + 阈值投影 + 采样。

数据面是 QMatrix 的 sample-and-query 结构（algorithms/qdata.py）。按论文
Lemma 5.3 构造 W = U·V，其中 U = Ũ R₁ Ũ⁻¹、V = Ṽ R₀ Ṽ⁻¹：Ũ 由行树把
(i, 0) 映到 (i, 行 i 的归一化向量)，Ṽ 由行范数根树把 (0, j) 映到
(用户分布 Ã, j)，R₀/R₁ 是对行/条目寄存器零基矢的反射。对 W 做相位估计
得到 θ，满足 cos(θ_i/2) = σ_i/‖A‖_F；逐相位字估计 σ̂ 并按阈值翻转 flag
（§5.3 Alg 2 的确定性投影版），逆相位估计后测 (flag, item)：flag=1 分支的
条目寄存器分布即推荐采样分布（§6）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from pyqecclang.algorithms.common.estimation import phase_estimation
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.algorithms.input_model.qdata import QMatrix
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import QRAM, Bits, ValidationError


def sigma_from_phase(value, precision, frobenius):
    """相位读数 t → 奇异值估计 σ̂ = ‖A‖_F 乘 cos(π t/2^precision) 的绝对值。

    W 的特征值成对出现 e^{±iθ}（同一 σ 的两个旋转方向），镜像相位
    2^precision−t 必须映射到同一奇异值，因此取绝对值。
    """
    if not 0 <= value < 1 << precision:
        raise ValidationError("相位读数越界")
    return frobenius * abs(math.cos(math.pi * value / (1 << precision)))


def _reflect_zero(builder, register):
    """对零基矢的反射（2 倍投影减恒等）：X 全翻 + 多控 Z + X 全翻；单比特即 Z。"""
    if register.width == 1:
        builder.z(register)
        return
    for bit in range(register.width):
        builder.x(register[bit])
    with builder.control(register[1:], (1 << (register.width - 1)) - 1):
        builder.z(register[0])
    for bit in range(register.width):
        builder.x(register[bit])


@dataclass(frozen=True)
class KPRecommendationConfig:
    """precision 为相位寄存器位数；sigma 为奇异值阈值，缺省取 0.5·‖A‖_F。"""

    precision: int = 4
    sigma: float | None = None


@dataclass(frozen=True)
class RecommendationResult:
    """KP 推荐采样电路及其读出契约。

    Attributes:
        operation: 推荐采样电路；寄存器 row、item、phase、flag，QRAM 资源为
            row_angles 与 root_angles。
        matrix: 构造电路所用的 QMatrix 输入。
        user: 目标用户编号。
        precision: 相位寄存器位数。
        sigma: 实际采用的奇异值阈值；config 未指定时取 0.5·frobenius。
        frobenius: 矩阵的 Frobenius 范数，sigma_from_phase 解码时使用。
    """

    operation: Operation
    matrix: QMatrix
    user: int
    precision: int
    sigma: float
    frobenius: float

    def memories(self):
        """提取电路所需两座 QRAM 角度库的初值。

        Returns:
            dict: 键为 row_angles（行树）与 root_angles（行范数根树），值为
            地址到角度字的映射，供执行入口绑定电路声明的 QRAM 资源。"""
        snapshot = self.matrix.snapshot()
        return {"row_angles": snapshot["row_angles"], "root_angles": snapshot["root_angles"]}

    def readout(self, state):
        """把模拟/执行结果折算为 (成功概率, 条件推荐分布)。"""
        registers = self.operation.module.registers
        index = {r.name: i for i, r in enumerate(registers)}
        success, items = 0.0, {}
        for key, amplitude in state.amplitudes.items():
            probability = abs(amplitude) ** 2
            if key[index["flag"]]:
                success += probability
                item = key[index["item"]]
                items[item] = items.get(item, 0.0) + probability
        distribution = {item: p / success for item, p in items.items()} if success else {}
        return success, distribution


def _walk_unitary(matrix: QMatrix):
    """W = Ũ R₁ Ũ⁻¹ · Ṽ R₀ Ṽ⁻¹；寄存器 row/item，资源与 QMatrix bank 同名。"""
    r, c, aw = matrix.rows, matrix.cols, matrix.angle_width
    b = Builder(
        _name("kp_walk", matrix.rows, matrix.cols, matrix.angle_width),
        {"row": Bits(r), "item": Bits(c)},
        {"row_angles": QRAM(r + c, aw), "root_angles": QRAM(r, aw)},
        attributes={"algorithm": "kp_recommendation_walk"},
    )
    amp = matrix.amplitude_preparation()
    row = matrix.row_preparation()
    amp_work = b.local("amp_work", Bits(aw))
    row_work = b.local("row_work", Bits(aw))

    def sandwich(prep, register, work, resource):
        """prep 共轭的零基矢反射：伴随先行，反射居中，正向收尾。"""
        with b.adjoint():
            b.call(prep, row=b["row"], item=b["item"], work=work, resources=resource)
        _reflect_zero(b, register)
        b.call(prep, row=b["row"], item=b["item"], work=work, resources=resource)

    # 先 V = Ṽ R₀ Ṽ⁻¹，后 U = Ũ R₁ Ũ⁻¹；算符乘积 U·V。
    sandwich(amp, b["row"], amp_work, {"root_angles": "root_angles"})
    sandwich(row, b["item"], row_work, {"row_angles": "row_angles"})
    return b.finish()


def kp_recommendation(matrix: QMatrix, user: int, config: KPRecommendationConfig | None = None):
    """对用户 user 生成推荐采样电路；返回含读出契约的 RecommendationResult。"""
    config = config or KPRecommendationConfig()
    if not isinstance(matrix, QMatrix):
        raise ValidationError("kp_recommendation 需要 QMatrix 输入")
    if not 1 <= config.precision <= 12:
        raise ValidationError("推荐相位寄存器位数必须为 1..12")
    r, c, aw = matrix.rows, matrix.cols, matrix.angle_width
    sigma = config.sigma
    if sigma is None:
        sigma = 0.5 * matrix.frobenius
    if not 0 < sigma <= matrix.frobenius:
        raise ValidationError("奇异值阈值必须处于 (0, ‖A‖_F]")
    walk = _walk_unitary(matrix)
    qpe = phase_estimation(walk, precision=config.precision)
    prep = matrix.row_state_prep(user)
    b = Builder(
        _name("kp_recommendation", matrix.rows, matrix.cols, aw, user, config.precision, sigma),
        {
            "row": Bits(r),
            "item": Bits(c),
            "phase": Bits(config.precision),
            "flag": Bits(1),
        },
        {"row_angles": QRAM(r + c, aw), "root_angles": QRAM(r, aw)},
        attributes={
            "algorithm": "kp_recommendation",
            "correctness": "pending",
            "user": user,
            "sigma": sigma,
            "precision": config.precision,
            "frobenius": matrix.frobenius,
            "reference": "arXiv:1603.08675",
            "assumptions": "; ".join(
                (
                    "nonnegative fixed-point entries",
                    "rotation-angle quantization at angle_width bits",
                    "deterministic flag projection; no amplitude amplification",
                )
            ),
        },
    )
    qpe_resources = {f"u__{name}": name for name in ("row_angles", "root_angles")}
    amp = matrix.amplitude_preparation()
    prep_work = b.local("prep_work", Bits(aw))
    amp_work = b.local("amp_work", Bits(aw))
    b.call(prep.operation, target=b["item"], work=prep_work, resources={"row_angles": "row_angles"})
    b.call(amp, row=b["row"], item=b["item"], work=amp_work, resources={"root_angles": "root_angles"})
    b.call(qpe, row=b["row"], item=b["item"], phase=b["phase"], resources=qpe_resources)
    for value in range(1 << config.precision):
        if sigma_from_phase(value, config.precision, matrix.frobenius) >= sigma:
            with b.control(b["phase"], value):
                b.x(b["flag"])
    with b.adjoint():
        b.call(qpe, row=b["row"], item=b["item"], phase=b["phase"], resources=qpe_resources)
    with b.adjoint():
        b.call(
            amp,
            row=b["row"],
            item=b["item"],
            work=amp_work,
            resources={"root_angles": "root_angles"},
        )
    return RecommendationResult(
        b.finish(),
        matrix,
        user,
        config.precision,
        sigma,
        matrix.frobenius,
    )
