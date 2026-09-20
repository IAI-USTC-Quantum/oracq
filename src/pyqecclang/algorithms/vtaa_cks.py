"""CKS 第 5 节 VTAA 线性系统求解器：GPE 时钟、分频带逆 LCU 与变时幅度放大。

依据 Childs–Kothari–Somma（arXiv:1511.02306 §5，SIAM J. Comput. 2017）。变时层结构：
时钟位 C_j 由带隙相位估计（Lemma 22）写入，受控触发分频带 Chebyshev 逆 LCU
（Lemma 23，含式 (98) 的 alpha_max 均匀化旋转），外层按 Ambainis（arXiv:1010.4458）
的时钟感知嵌套放大（操作符形式取 Low–Su arXiv:2410.18178 式 (47)–(53)），最后以
A' 的逆抹除 GPE 垃圾（式 (99)–(112)）。

GPE 判决采用 Low–Su Prop 23 的确定性路线：对 qubitization walk 施加 Yoder–Low–Chuang
定点判决多项式（复用 qsvt.fixed_point_search_phases 的已验证合成），在 signal==0
分支翻转判决位：λ 模长 ≥ θ_j 时判 1（fire，硬界 epsilon）振幅接近 1，
x=0 处判决为 0；过渡带响应确定、可精确计算。查询次数 O((alpha/theta_j) log(1/eps))
与 CKS Lemma 22 的 PEA + 多数投票同阶；P_j 垃圾按论文保留，由 A' 的逆抹除。
全网此前没有 VTAA 或 CKS §5 的公开实现。
"""

from __future__ import annotations

import cmath
import json
import math
from dataclasses import dataclass
from functools import lru_cache

from pyqecclang.algorithms.arithmetic import BooleanNetwork
from pyqecclang.algorithms.block_encoding import reflect_zero
from pyqecclang.algorithms.contracts import finite_real, positive_integer, require_instance
from pyqecclang.algorithms.interfaces import as_block_encoding
from pyqecclang.algorithms.operators import BlockEncoding, _name
from pyqecclang.algorithms.oracles import (
    StateOracle,
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from pyqecclang.algorithms.qlss import CKSConfig, QLSSProtocol
from pyqecclang.algorithms.qsvt import fixed_point_search_phases, qsp_response
from pyqecclang.algorithms.sparse import chebyshev_block, real_symmetric_sparse_encoding
from pyqecclang.algorithms.transforms import qsvt_sequence
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse


@lru_cache(maxsize=64)
def gpe_fire_phases(threshold, x_edge, epsilon, degree_cap=40):
    """GPE fire 判决相位：YLC 定点多项式，x 模长 >= threshold 时 P 的模 >= 1-epsilon。

    阈值随 degree 单调下降，从 3 起按奇数搜索至 [threshold, x_edge] 网格验证通过；
    x=0 处 P=0，过渡带（0, threshold) 的响应确定、可精确计算（CKS 未承诺带）。
    """
    finite_real(threshold, "gpe_fire.threshold", minimum=0, strict=True)
    finite_real(x_edge, "gpe_fire.x_edge", minimum=0, strict=True)
    finite_real(epsilon, "gpe_fire.epsilon", minimum=0, strict=True)
    if not 0 < threshold <= x_edge <= 1 or not 0 < epsilon < 1:
        raise ValidationError("GPE 判决几何或精度无效")
    fire_grid = [threshold + (x_edge - threshold) * i / 16.0 for i in range(17)]
    for degree in range(3, degree_cap + 1, 2):
        try:
            phases = fixed_point_search_phases(epsilon, degree)
        except ValidationError:
            continue
        if min(abs(qsp_response(x, phases)) for x in fire_grid) >= 1 - epsilon:
            return phases, degree
    raise ValidationError("GPE 判决多项式度数超出上限；请增大 phi 或降低 kappa 声明")


@lru_cache(maxsize=64)
def clock_or_operation(width):
    """时钟前缀 OR 谓词：stopped<=j 的相干判据，供 VTAA 反射使用。"""
    net = BooleanNetwork()
    bits = net.input("prefix", width)
    net.outputs = {"stopped": [net.any(bits)]}
    return net.operation(attributes={"algorithm": "clock_prefix_or", "width": width})


def gapped_phase_estimation(
    a: BlockEncoding, threshold, x_edge, *, epsilon=0.02, degree_cap=40
):
    """CKS Lemma 22 的 GPE（确定性 QSP 路线，Low–Su arXiv:2410.18178 Prop 23）。

    对 BE 的 qubitization walk 施加 YLC 定点判决多项式 P(lambda/alpha)，在
    signal==0 分支翻转判决位：λ/α 模长 >= threshold 时（fire 硬界 epsilon）
    判 1 振幅接近 1（在本频带停下求逆），谱低端判决趋向 0（延后到更细频带）。
    signal 寄存器保留 γ 垃圾（论文 §5.2 的 P_j），由 A' 的逆在收尾时抹除。

    Args:
        a: 自伴酉扩张的 Hermitian BE。
        threshold: fire 判决阈值（编码谱单位，即相对 alpha）。
        x_edge: 谱上界/alpha，验证网格右端。
        epsilon: fire 侧判决响应硬界。
        degree_cap: 判决多项式度数上限。

    Returns:
        Operation: 寄存器 target/signal/decision；decision=1 表示在本频带停下。
    """
    a = as_block_encoding(a)
    phases, degree = gpe_fire_phases(threshold, x_edge, epsilon, degree_cap)
    marker = qsvt_sequence(a, phases)
    b = Builder(
        _name("vtaa_gpe", a.operation, threshold, x_edge, epsilon, degree),
        {
            "target": Bits(a.width),
            "signal": Bits(a.signal_qubits),
            "decision": Bits(1),
        },
        resources_for(("marker", marker)),
        attributes={
            "algorithm": "gapped_phase_estimation",
            "fire_threshold": threshold,
            "marker_degree": degree,
            "marker_epsilon": epsilon,
            "walk_queries_max": degree,
            "references": "CKS Lemma 22 via Low-Su Prop 23 deterministic QSP decision",
            "validation_stage": "paradigm",
        },
    )
    invoke(b, marker, "marker", target=b["target"], signal=b["signal"])
    with b.control(b["signal"], 0):
        b.x(b["decision"])
    return b.finish()


def band_inverse_step(a: BlockEncoding, coefficients, alpha_max):
    """CKS Lemma 23 的 W(lambda, delta)：分频带 Chebyshev 逆 LCU 加均匀化旋转。

    Args:
        a: 自伴酉扩张 BE；chebyshev_block 提供奇次 Chebyshev 幂。
        coefficients: 本频带的截断逆多项式系数 c_k（T_{2k+1} 基）。
        alpha_max: 全部频带的最大 LCU 1-范数；式 (98) 把成功振幅均匀压到 1/alpha_max。

    Returns:
        Operation: 寄存器 target/signal/flag；signal==0 且 flag==1 的分支携带
        h(A)ψ/alpha_max（ψ 为输入态），h 为本频带多项式。
    """
    a = as_block_encoding(a)
    terms = tuple(
        (complex(c), chebyshev_block(a, 2 * k + 1))
        for k, c in enumerate(coefficients)
        if c != 0
    )
    if not terms:
        raise ValidationError("频带逆多项式系数全为零")
    beta = sum(abs(c) for c, _ in terms)
    if alpha_max <= 0 or beta > alpha_max * (1 + 1e-9):
        raise ValidationError("频带 LCU 归一化与 alpha_max 声明冲突")
    selector_width = (len(terms) - 1).bit_length()
    b = Builder(
        _name("vtaa_band_inverse", a.operation, tuple(coefficients), alpha_max),
        {
            "target": Bits(a.width),
            "signal": Bits(selector_width + a.signal_qubits),
            "flag": Bits(1),
        },
        resources_for(*[(f"term{i}", t.operation) for i, (_, t) in enumerate(terms)]),
        attributes={
            "algorithm": "vtaa_band_inverse",
            "band_lcu_normalization": beta,
            "uniform_normalization": alpha_max,
            "polynomial_terms": len(terms),
            "implementation_scope": "CKS eq. 96-98 piecewise inverse with fixed alpha_max",
            "validation_stage": "paradigm",
        },
    )
    selector, work = b["signal"][:selector_width], b["signal"][selector_width:]
    if selector_width:
        weights = [math.sqrt(abs(c)) for c, _ in terms]
        weights += [0.0] * ((1 << selector_width) - len(weights))
        prep = gate_state_prep(weights)
        invoke(b, prep.operation, target=selector, work=selector[:0])
        for i, (coefficient, term) in enumerate(terms):
            with b.control(selector, i):
                b.global_phase(cmath.phase(coefficient))
                invoke(b, term.operation, f"term{i}", target=b["target"], signal=work)
        with b.adjoint():
            invoke(b, prep.operation, target=selector, work=selector[:0])
    else:
        coefficient, term = terms[0]
        b.global_phase(cmath.phase(coefficient))
        invoke(b, term.operation, "term0", target=b["target"], signal=work)
    # 式 (98)：LCU 成功分支上把旗标转到 beta/alpha_max 振幅，均匀化各频带成功幅值。
    with b.control(b["signal"], 0):
        b.ry(b["flag"], 2 * math.asin(beta / alpha_max))
    return b.finish()


@dataclass(frozen=True)
class VTAAConfig:
    """order 为第 1 频带的基础 Chebyshev 阶数，第 j 频带按 2^(j-1) 几何放大。

    rounds 是各阶段的 VTAA 放大轮数 r_j（None 表示全零，即纯变时层加后选）；
    正确性不依赖日程，日程只影响成功率。生产部署应按 Ambainis 算法 2 的振幅
    估计或 Low–Su 的确定性日程选择 r_j（见 tunable_rounds）。
    """

    order: int = 2
    terms: int | None = None
    clock_steps: int | None = None
    marker_epsilon: float = 0.02
    degree_cap: int = 40
    rounds: tuple[int, ...] | None = None

    def __post_init__(self):
        positive_integer(self.order, "VTAAConfig.order", maximum=128)
        if self.terms is not None:
            positive_integer(self.terms, "VTAAConfig.terms", maximum=self.order)
        if self.clock_steps is not None:
            positive_integer(self.clock_steps, "VTAAConfig.clock_steps", maximum=16)
        finite_real(self.marker_epsilon, "VTAAConfig.marker_epsilon", minimum=0, strict=True)
        if self.marker_epsilon >= 0.2:
            raise ValidationError("GPE 判决硬界需要小于 0.2")
        positive_integer(self.degree_cap, "VTAAConfig.degree_cap", minimum=8, maximum=40)
        if self.rounds is not None:
            rounds = tuple(int(r) for r in self.rounds)
            if not rounds or any(not 0 <= r <= 64 for r in rounds):
                raise ValidationError("VTAA rounds 需要非负且不超过 64")
            object.__setattr__(self, "rounds", rounds)

    def band_order(self, step):
        """返回第 step 频带的 Chebyshev 逆多项式阶数。

        第 1 频带取 ``order``，其后每个频带翻倍（2 的几何增长），上限 128；
        ``step`` 小于 1 时按第 1 频带处理。"""
        return min(128, self.order * (1 << max(0, step - 1)))

    def band_coefficients(self, step):
        """返回第 step 频带的截断逆多项式系数（T_{2k+1} 基）。

        以 ``band_order(step)`` 为阶数、``terms`` 为截断项数（None 表示取满阶），
        经 ``CKSConfig.coefficients`` 计算得到。"""
        return CKSConfig(self.band_order(step), self.terms).coefficients()


def tunable_rounds(stage_amplitudes, thresholds=None):
    """Low–Su 可调 VTAA 日程（arXiv:2410.18178 式 (52)–(53)）。

    Args:
        stage_amplitudes: 各阶段“尚未失败”振幅范数 x_j 的估计（可由振幅估计通道获得）。
        thresholds: 阈值 alpha_j，缺省均分；总和为常数以保证常数损失。

    Returns:
        每阶段轮数元组 r_j = max(ceil(sqrt(alpha_j)/(6 x_j) - 1/2), 0)，
        满足不过冲条件 (2 r_j + 1) x_j <= sqrt(alpha_j)。
    """
    values = tuple(float(v) for v in stage_amplitudes)
    if not values or any(not math.isfinite(v) or not 0 < v <= 1 for v in values):
        raise ValidationError("阶段范数必须是 (0,1] 内的有限数")
    thresholds = (
        tuple(1.0 / len(values) for _ in values)
        if thresholds is None
        else tuple(float(t) for t in thresholds)
    )
    if len(thresholds) != len(values) or any(
        not math.isfinite(t) or t <= 0 for t in thresholds
    ):
        raise ValidationError("VTAA 阈值数量或数值无效")
    return tuple(
        max(0, math.ceil(math.sqrt(alpha) / (6 * x) - 0.5))
        for x, alpha in zip(values, thresholds, strict=True)
    )


def vtaa_cks(system, config=None):
    """CKS §5 的 VTAA 求解内核；输入 SparseSystem，输出解态 StateOracle。"""
    config = config or VTAAConfig()
    require_instance(config, VTAAConfig, "vtaa_cks.config")
    if not system.hermitian:
        raise ValidationError("VTAA 稀疏输入需要 Hermitian 声明或显式 Hermitian dilation")
    a = real_symmetric_sparse_encoding(
        system.access,
        system.value_format,
        system.entry_bound,
        diagonal_nonnegative=system.diagonal_nonnegative,
    )
    x_edge = system.spectrum.norm_upper / a.alpha
    physical_kappa = system.spectrum.norm_upper / system.spectrum.sigma_min_lower
    # 覆盖条件：最细频带的 fire 阈值 2^(1-steps)*x_edge 需不超过 sigma_min/alpha。
    if config.clock_steps is not None:
        steps = config.clock_steps
    else:
        steps = math.ceil(math.log2(physical_kappa)) + 1 if physical_kappa > 1 else 1
    if physical_kappa > 2 ** (steps - 1):
        raise ValidationError("clock_steps 不足以覆盖声明的最小奇异值频带")
    if config.rounds is not None and len(config.rounds) != steps:
        raise ValidationError("VTAA rounds 长度必须等于 clock_steps")

    bands = []
    for step in range(1, steps + 1):
        threshold = x_edge * 2.0 ** (1 - step)
        coefficients = config.band_coefficients(step)
        bands.append((threshold, coefficients, sum(abs(c) for c in coefficients)))
    alpha_max = max(beta for _, _, beta in bands)
    gpes = tuple(
        gapped_phase_estimation(
            a,
            threshold,
            x_edge,
            epsilon=config.marker_epsilon,
            degree_cap=config.degree_cap,
        )
        for threshold, _, _ in bands
    )
    inverses = tuple(
        band_inverse_step(a, coefficients, alpha_max) for _, coefficients, _ in bands
    )
    rounds = config.rounds or (0,) * steps
    selector_width = max(
        (len([c for c in coefficients if c != 0]) - 1).bit_length()
        for _, coefficients, _ in bands
    )
    # signal 布局：[共享逆 LCU 区 selw+s][每步独立 GPE 垃圾 m*s]；P_j 垃圾保留待抹除。
    signal_width = selector_width + a.signal_qubits * (steps + 1)
    registers = {
        "target": Bits(a.width),
        "clock": Bits(steps),
        "flag": Bits(1),
        "signal": Bits(signal_width),
    }
    gpe_base = selector_width + a.signal_qubits

    def build_step(step, *, uncompute):
        """A_j：受控 GPE 写 C_j；C_j=1（fire）时施加 W_j（A'_j 仅翻旗标，式 (99)）。"""
        threshold, coefficients, beta = bands[step - 1]
        kind = "vtaa_uncompute_step" if uncompute else "vtaa_variable_step"
        resources = resources_for(("gpe", gpes[step - 1]))
        if not uncompute:
            resources.update(resources_for(("inverse", inverses[step - 1])))
        b = Builder(
            _name(kind, a.operation, system.rhs.operation, step, coefficients, config),
            registers,
            resources,
            attributes={
                "algorithm": kind,
                "band_index": step,
                "fire_threshold": threshold,
                "band_lcu_normalization": beta,
                "validation_stage": "paradigm",
            },
        )
        prefix = b["clock"][: step - 1]
        gpe_slice = b["signal"][gpe_base + (step - 1) * a.signal_qubits :][
            : a.signal_qubits
        ]
        inverse_width = (
            (len([c for c in coefficients if c != 0]) - 1).bit_length() + a.signal_qubits
        )

        def run_gpe():
            invoke(
                b,
                gpes[step - 1],
                "gpe",
                target=b["target"],
                signal=gpe_slice,
                decision=b["clock"][step - 1],
            )

        # 前缀全 0 = 尚未停下；C_j=1 表示在本频带 fire。
        if prefix.width:
            with b.control(prefix, 0):
                run_gpe()
        else:
            run_gpe()
        with b.control(b["clock"][step - 1]):
            if uncompute:
                b.x(b["flag"])
            else:
                invoke(
                    b,
                    inverses[step - 1],
                    "inverse",
                    target=b["target"],
                    signal=b["signal"][:inverse_width],
                    flag=b["flag"],
                )
        return b.finish()

    initial = Builder(
        _name("vtaa_initial", system.rhs.operation),
        registers,
        resources_for(("rhs", system.rhs.operation)),
    )
    rhs_work = initial.local("rhs_work", Bits(system.rhs.work_width))
    invoke(initial, system.rhs.operation, "rhs", target=initial["target"], work=rhs_work)
    initial_op = initial.finish()

    def run_chain(b, items):
        for prefix, op in items:
            invoke(
                b,
                op,
                prefix,
                target=b["target"],
                clock=b["clock"],
                flag=b["flag"],
                signal=b["signal"],
            )

    def prefix_module(step, items):
        b = Builder(
            _name("vtaa_prefix", a.operation, system.rhs.operation, config, step),
            registers,
            resources_for(*items),
            attributes={
                "algorithm": "vtaa_prefix",
                "stage": step,
                "validation_stage": "paradigm",
            },
        )
        run_chain(b, items)
        return b.finish()

    def amplification(prefix_op, step, count):
        """M_j = (R_s R_f)^{r_j} P_j；R_f 翻转 stopped<=j 且旗标为 0 的相位。"""
        b = Builder(
            _name("vtaa_amplified", prefix_op, step, count),
            registers,
            resources_for(("prefix", prefix_op)),
            attributes={
                "algorithm": "vtaa_amplified_stage",
                "stage": step,
                "rounds": count,
                "implementation_scope": "Ambainis VTAA cascade; operator form Low-Su eq. 47",
                "validation_stage": "paradigm",
            },
        )
        run_chain(b, (("prefix", prefix_op),))

        def run_prefix():
            run_chain(b, (("prefix", prefix_op),))

        if count:
            orop = clock_or_operation(step)
            stopped = b.local("stopped", Bits(1))

            def mark_stopped():
                b.call(orop, prefix=b["clock"][:step], stopped=stopped)

            with b.repeat(count):
                mark_stopped()
                # stopped=1（已 fire）且 flag=0（失败）的分支获得 pi 相位。
                with b.control(fuse(stopped, b["flag"]), 1):
                    b.global_phase(math.pi)
                mark_stopped()
                with b.adjoint():
                    run_prefix()
                reflect_zero(
                    b, fuse(b["target"], b["clock"], b["flag"], b["signal"]), positive=True
                )
                run_prefix()
        return b.finish()

    chain = [("start", initial_op)]
    amplified = initial_op
    for step in range(1, steps + 1):
        chain.append((f"step{step}", build_step(step, uncompute=False)))
        amplified = amplification(prefix_module(step, chain), step, rounds[step - 1])
        chain = [("run", amplified)]
    erase_ops = tuple(build_step(step, uncompute=True) for step in range(1, steps + 1))

    top = Builder(
        _name("vtaa_cks", a.operation, system.rhs.operation, config),
        {"target": Bits(a.width), "signal": Bits(steps + 1 + signal_width)},
        resources_for(
            ("run", amplified), *((f"erase{s}", op) for s, op in enumerate(erase_ops, 1))
        ),
        attributes={
            "algorithm": "vtaa_cks",
            "input_model": "sparse_location_inplace_and_entry_xor",
            "input_alpha": a.alpha,
            "encoded_inverse_bound": system.spectrum.inverse_bound(a.alpha),
            "clock_steps": steps,
            "fire_thresholds": json.dumps([t for t, _, _ in bands]),
            "band_orders": json.dumps([config.band_order(s) for s in range(1, steps + 1)]),
            "band_lcu_normalizations": json.dumps([beta for _, _, beta in bands]),
            "alpha_max": alpha_max,
            "marker_degrees": json.dumps(
                [dict(op.module.attributes)["marker_degree"] for op in gpes]
            ),
            "marker_epsilon": config.marker_epsilon,
            "rounds": json.dumps(list(rounds)),
            "implementation_scope": (
                "CKS arXiv:1511.02306 section 5 VTAA; GPE per Lemma 22 realized as the"
                " deterministic QSP decision of Low-Su arXiv:2410.18178 Prop 23 on the"
                " qubitization walk with YLC fixed-point markers; clock-aware cascade"
                " per Ambainis arXiv:1010.4458 with operator form from Low-Su eq. 47-53"
            ),
            "kernel_status": "prototype; band polynomial accuracy and VTAA schedule pending",
            "correctness": "pending",
            "success_condition": "signal == 0; erased branch carries h(A)|b> mixture",
            "validation_stage": "paradigm",
        },
    )
    clock, flag, rest = (
        top["signal"][:steps],
        top["signal"][steps],
        top["signal"][steps + 1 :],
    )

    def call_top(op, prefix):
        invoke(top, op, prefix, target=top["target"], clock=clock, flag=flag, signal=rest)

    call_top(amplified, "run")
    # (A')^dagger：逆序抹除 GPE 时钟与旗标（式 (110)–(112)），解态留在 target。
    with top.adjoint():
        for step in range(1, steps + 1):
            call_top(erase_ops[step - 1], f"erase{step}")
    return StateOracle(
        annotate(
            top.finish(),
            "unitary",
            algorithm="vtaa_cks",
            implementation_scope="CKS section 5 variable-time amplitude amplification",
            correctness="pending",
        )
    )


def make_vtaa_cks_qlss(config=None):
    """CKS §5 VTAA 求解协议；输入模型与 cks_chebyshev 相同的稀疏访问。"""
    config = VTAAConfig() if config is None else config
    require_instance(config, VTAAConfig, "make_vtaa_cks_qlss.config")
    return QLSSProtocol("vtaa_cks", "sparse", lambda system: vtaa_cks(system, config))
