"""量子奇异值变换（QSVT）标准变换库（Gilyén et al. 2019, arXiv:1806.01838）。

约定（与 transforms.qsvt_sequence 逐步对齐，并经参考模拟逐点验证）：
``qsvt_sequence(a, Φ)`` 在奇异值 x 的二维不变子空间上实现

    p(x) = [S(φ_0) W(x) S(φ_1) W(x) … W(x) S(φ_d)]_00,
    W(x) = [[x, s], [s, -x]],  s = √(1 − x²),  S(φ) = diag(e^{iφ}, e^{−iφ}),

相位按时间顺序排列（φ_0 最先作用），共 d 次 BE 调用、 d+1 个相位。
可实现对 (P, Q) 的充要条件：deg P ≤ d、deg Q ≤ d−1、P 的奇偶性为 d mod 2、
Q 的奇偶性为 (d−1) mod 2，且多项式恒等式 P P̄ + (1−x²) Q Q̄ ≡ 1 成立；
特别地必有 ``|P(±1)| = 1`` （端点饱和）。

相位合成采用 Gilyén 论文的补多项式求根 + 逐层剥离（layer stripping）：
给定实目标 f 与可选虚部补全 h（均为实系数多项式），令 P = f + i·h，
从 R = (1 − f² − h²)/(1−x²) 的根构造 Q，再逆向递推恢复相位。
端点不饱和的实目标（如 1/x 的截断近似）必须借助非零虚部补全；
此时用 (U_Φ + U_{−Φ})/2 的 LCU 组合提取实部（−Φ 恰好实现 P̄），
块编码即为 f(A/α)。

数值边界：求根用纯 Python Durand–Kerner 迭代，剥离用双精度复数；
合成度数限制为 40，且每次合成后用 qsp_response 做往返自检，
超限抛出 ValidationError。度数几十以内、补多项式根分离良好时
往返误差通常在 1e-9 量级。
"""

from __future__ import annotations

import cmath
import math

from pyqecclang.algorithms.contracts import require_instance
from pyqecclang.algorithms.operators import BlockEncoding, identity, linear_combination
from pyqecclang.algorithms.oracles import annotate
from pyqecclang.algorithms.transforms import qsvt_sequence
from pyqecclang.infrastructure.ir import ValidationError

__all__ = [
    "eigenstate_filter",
    "fixed_point_search",
    "fixed_point_search_phases",
    "qsp_phases",
    "qsp_response",
    "qsvt_hamiltonian_simulation",
    "qsvt_matrix_inversion",
]

_MAX_DEGREE = 40
_BOUND_TOL = 1e-9
_GRID = 4096
_STRIP_TOL = 1e-5


# ---------------------------------------------------------------------------
# 多项式工具：升幂系数（常数项在前），实系数用 float、复系数用 complex。
# ---------------------------------------------------------------------------


def _trim(p, tol=1e-12):
    p = list(p)
    while len(p) > 1 and abs(p[-1]) <= tol * max(1.0, max(abs(c) for c in p)):
        p.pop()
    return tuple(p)


def _add(a, b):
    n = max(len(a), len(b))
    return tuple((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n))


def _sub(a, b):
    return _add(a, tuple(-v for v in b))


def _mul(a, b):
    out = [0j] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def _scale(c, p):
    return tuple(c * v for v in p)


def _eval(p, x):
    v = 0j
    for c in reversed(p):
        v = v * x + c
    return v


def _conj(p):
    return tuple(complex(v).conjugate() for v in p)


def _realify(p, *, tol=1e-10):
    scale = max(1.0, max(abs(c) for c in p))
    if any(abs(complex(c).imag) > tol * scale for c in p):
        raise ValidationError("内部多项式应为实系数")
    return tuple(float(complex(c).real) for c in p)


def _deg(p):
    return len(p) - 1


def _grid(n=_GRID):
    return tuple(math.cos(math.pi * j / n) for j in range(n + 1))


def _sup_norm(p):
    return max(abs(_eval(p, x)) for x in _grid())


def _roots(coeffs, *, iters=4000, tol=1e-30):
    """Durand–Kerner 同时求根；输入升幂首一化前的任意实/复系数。"""
    coeffs = _trim(tuple(complex(c) for c in coeffs))
    n = _deg(coeffs)
    if n <= 0:
        return []
    lead = coeffs[n]
    a = tuple(c / lead for c in coeffs)
    roots = [cmath.exp(2j * math.pi * (k + 0.318) / n) * (1.0 + 0.4j) for k in range(n)]
    for _ in range(iters):
        worst = 0.0
        for i in range(n):
            denom = 1 + 0j
            for j in range(n):
                if i != j:
                    denom *= roots[i] - roots[j]
            step = _eval(a, roots[i]) / denom
            roots[i] -= step
            worst = max(worst, abs(step))
        if worst < tol:
            break
    return roots


def _cluster(roots, *, rel=2e-5):
    """按相对距离把数值根聚成重根簇，返回 [(质心, 重数), …]。"""
    clusters = []
    for r in sorted(roots, key=lambda z: (abs(z), z.real, z.imag)):
        for cluster in clusters:
            if abs(r - cluster[0]) <= rel * max(1.0, abs(r)):
                cluster[1] += 1
                break
        else:
            clusters.append([r, 1])
    return [(c, m) for c, m in clusters]


def _chebyshev_t(n):
    """T_n 的升幂单项式系数。"""
    if n == 0:
        return (1.0,)
    if n == 1:
        return (0.0, 1.0)
    a, b = (1.0,), (0.0, 1.0)
    for _ in range(2, n + 1):
        a, b = b, _sub(_scale(2.0, (0.0,) + b), a)
    return b


# ---------------------------------------------------------------------------
# QSP 响应与相位合成（reflection 约定，与 qsvt_sequence 一致）。
# ---------------------------------------------------------------------------


def qsp_response(x, phases):
    """相位序列 Φ 在反射约定下实现的顶层左块 p(x)（x ∈ [−1, 1]）。"""
    phases = tuple(float(p) for p in phases)
    s = math.sqrt(max(0.0, 1.0 - x * x))
    m00, m01, m10, m11 = 1 + 0j, 0j, 0j, 1 + 0j
    for i, phi in enumerate(phases):
        if i:  # W(x) 左乘（时间上先于本步相位）
            m00, m01, m10, m11 = (
                x * m00 + s * m10,
                x * m01 + s * m11,
                s * m00 - x * m10,
                s * m01 - x * m11,
            )
        e, f = cmath.exp(1j * phi), cmath.exp(-1j * phi)
        m00, m01, m10, m11 = e * m00, e * m01, f * m10, f * m11
    return m00


def _check_real_poly(coeffs, label):
    values = []
    for c in coeffs:
        z = complex(c)
        if not math.isfinite(z.real) or abs(z.imag) > 1e-12:
            raise ValidationError(f"{label} 必须是有限实系数多项式")
        values.append(z.real)
    if not values:
        raise ValidationError(f"{label} 不能为空")
    return _trim(values)


def _check_parity(coeffs, d, label):
    scale = max(1.0, max(abs(c) for c in coeffs))
    for k, c in enumerate(coeffs):
        if (d - k) % 2 and abs(c) > 1e-9 * scale:
            raise ValidationError(f"{label} 的奇偶性与度数 d mod 2 不符")


def _div_1mx2(dpoly, *, tol=1e-8):
    """计算 R = D/(1−x²)；要求 D(±1) = 0，否则抛出 ValidationError。"""
    scale = max(1.0, max(abs(c) for c in dpoly))
    n = len(dpoly)
    r = [0j] * (n + 2)
    for j in range(n + 2):
        r[j] = (dpoly[j] if j < n else 0) + (r[j - 2] if j >= 2 else 0)
    if max(abs(r[n - 2]), abs(r[n - 1])) > tol * scale:
        raise ValidationError("补多项式条件不满足：1 − f² − h² 不能被 1−x² 整除（端点未饱和）")
    return _realify(_trim(r[: max(1, n - 2)]))


def _q_from_roots(rpoly, d):
    """由 R = (1−f²−h²)/(1−x²) 的根构造奇偶性为 (d−1) mod 2 的复系数 Q，使 Q Q̄ = R。"""
    rpoly = _trim(rpoly)
    if _deg(rpoly) <= 0:
        if rpoly[0] <= 0:
            raise ValidationError("补多项式 R 恒为非正，无法谱分解")
        return (math.sqrt(rpoly[0]),) if d % 2 == 1 else None
    if _deg(rpoly) != 2 * d - 2:
        raise ValidationError("补多项式度数与目标不匹配")
    if rpoly[-1] <= 0:
        raise ValidationError("补多项式首项系数必须为正")
    clusters = _cluster(_roots(rpoly))
    factors = []  # 每个因子均为偶多项式；零根单独给出奇偶性
    zero_mult = 0
    used = set()
    ctol = 1e-5
    for i, (c, m) in enumerate(clusters):
        if i in used:
            continue
        if abs(c) <= ctol:
            zero_mult += m
            used.add(i)
            continue
        if abs(c.imag) <= ctol * abs(c):  # 实根 r 与 −r
            partner = next(
                (j for j, (c2, m2) in enumerate(clusters) if j > i and j not in used
                 and abs(c2 + c) <= ctol * max(1.0, abs(c)) and m2 == m),
                None,
            )
            if partner is None or m % 2:
                raise ValidationError("补多项式实根必须成对且为偶数重")
            rr = (abs(c.real) + abs(clusters[partner][0].real)) / 2
            for _ in range(m // 2):
                factors.append((-(rr * rr), 0.0, 1.0))
            used.add(partner)
        elif abs(c.real) <= ctol * abs(c):  # 纯虚根 i b 与 −i b（互为共轭）
            partner = next(
                (j for j, (c2, m2) in enumerate(clusters) if j > i and j not in used
                 and abs(c2 + c) <= ctol * max(1.0, abs(c)) and m2 == m),
                None,
            )
            if partner is None or m % 2:
                raise ValidationError("补多项式纯虚根重数必须为偶数")
            bb = (abs(c.imag) + abs(clusters[partner][0].imag)) / 2
            for _ in range(m // 2):
                factors.append((bb * bb, 0.0, 1.0))  # x² − (ib)² = x² + b²
            used.add(partner)
        else:  # 一般复根：四元组 {±ρ, ±ρ̄} 各 m 重；用对称化质心降低数值偏差
            sib = [
                j
                for j, (c2, m2) in enumerate(clusters)
                if j > i
                and j not in used
                and m2 == m
                and (abs(c2 - c.conjugate()) <= ctol * max(1.0, abs(c))
                     or abs(c2 + c) <= ctol * max(1.0, abs(c))
                     or abs(c2 + c.conjugate()) <= ctol * max(1.0, abs(c)))
            ]
            if len(sib) != 3:
                raise ValidationError("补多项式复根必须成共轭-反号四元组")
            members = [c] + [clusters[j][0] for j in sib]
            ra = sum(abs(v.real) for v in members) / 4
            rb = sum(abs(v.imag) for v in members) / 4
            rho = complex(ra, rb)
            for _ in range(m):
                factors.append((-rho * rho, 0.0, 1.0))  # (x−ρ)(x+ρ) = x² − ρ²
            used.update(sib)
        used.add(i)
    if len(used) != len(clusters):
        raise ValidationError("补多项式根结构不完整，无法谱分解")
    if zero_mult % 2:
        raise ValidationError("补多项式零根重数必须为偶数")
    m0 = zero_mult // 2
    if m0 % 2 != (d - 1) % 2:
        raise ValidationError("补多项式谱因子的奇偶性与度数要求不符")
    q = (math.sqrt(float(rpoly[-1])),)
    for f in factors:
        q = _mul(q, f)
    q = _mul(q, (0.0, 1.0)) if m0 else q
    q = _trim(tuple(0.0 if abs(v) < 1e-12 else v for v in _scale(1.0, q)))
    # 奇偶性：只保留与 d−1 同奇偶的幂次（数值噪声清零）
    q = tuple(v if (d - 1 - k) % 2 == 0 else 0.0 for k, v in enumerate(q))
    resid = _sub(_trim(_mul(q, _conj(q))), rpoly)
    if max((abs(v) for v in resid), default=0.0) > 1e-5 * max(1.0, max(abs(c) for c in rpoly)):
        raise ValidationError("补多项式谱分解自检失败（数值精度不足）")
    return q


def _strip(ppoly, qpoly, d):
    """layer stripping：由 (P, Q) 逐层恢复相位，时间正序返回。"""
    p, q = list(ppoly), list(qpoly)
    phases = []
    for k in range(d, 0, -1):
        if abs(q[k - 1]) < 1e-13:
            raise ValidationError("逐层剥离退化：补多项式首项过小，相位数值不稳定")
        phi = cmath.phase(p[k] / q[k - 1]) / 2
        phases.append(phi)
        ei, ej = cmath.exp(-1j * phi), cmath.exp(1j * phi)
        pn = _add(_scale(ei, (0j,) + tuple(p)), _scale(ej, _sub(q, (0j, 0j) + tuple(q))))
        qn = _sub(_scale(ei, p), _scale(ej, (0j,) + tuple(q)))
        # 理论上 x^k、x^{k+1} 与 x^{k−1} 的首部系数应精确相消；记录残差由自检兜底
        p, q = list(_trim(pn[:k])), list(_trim(qn[: max(1, k - 1)]))
        if len(p) < k:
            p += [0j] * (k - len(p))
        if len(q) < max(1, k - 1):
            q += [0j] * (max(1, k - 1) - len(q))
    if abs(abs(p[0]) - 1.0) > 1e-6:
        raise ValidationError("逐层剥离自检失败：零层相位模长偏离 1")
    phases.append(cmath.phase(p[0]))
    return tuple(reversed(phases))


def qsp_phases(coeffs, imag=None):
    """由实系数目标多项式合成 QSP 相位序列（时间正序，长度 d+1）。

    coeffs 为升幂实系数（常数项在前），目标为 P = f（imag 为 None）或
    P = f + i·h（imag 为 h 的升幂实系数）。可实现条件：f 与 h 的奇偶性均为
    d mod 2、在 [−1,1] 上 f² + h² ≤ 1、端点饱和 f(±1)² + h(±1)² = 1，
    且 R = (1 − f² − h²)/(1−x²) 非负并满足谱分解的根重数条件。
    不提供 imag 时即为纯实目标，此时必须有 ``|f(±1)| = 1``。
    """
    f = _check_real_poly(coeffs, "目标多项式")
    d = _deg(f)
    if d > _MAX_DEGREE:
        raise ValidationError(f"目标多项式度数超过合成上限 {_MAX_DEGREE}")
    if d == 0:
        if abs(abs(f[0]) - 1.0) > 1e-9:
            raise ValidationError("零次目标必须是单位模常数")
        return (cmath.phase(f[0]),)
    _check_parity(f, d, "目标多项式")
    h = _check_real_poly(imag, "虚部补全") if imag is not None else (0.0,)
    if _trim(h) != (0.0,):
        if _deg(h) > d:
            raise ValidationError("虚部补全度数不能超过目标度数")
        _check_parity(h + (0.0,) * (d + 1 - len(h)), d, "虚部补全")
    if _sup_norm(f) > 1.0 + 1e-6 and imag is None:
        raise ValidationError("目标多项式在 [−1,1] 上超过上界 1")
    bound = max(abs(complex(_eval(f, x), _eval(h, x))) for x in _grid())
    if bound > 1.0 + 1e-6:
        raise ValidationError("P = f + i·h 在 [−1,1] 上超过单位圆盘")
    sat = max(abs(_eval(f, 1.0).real ** 2 + _eval(h, 1.0).real ** 2 - 1.0),
              abs(_eval(f, -1.0).real ** 2 + _eval(h, -1.0).real ** 2 - 1.0))
    if sat > 1e-6:
        raise ValidationError("端点未饱和：需要 f(±1)² + h(±1)² = 1（可传入虚部补全）")
    dpoly = _trim(_sub((1.0,), _add(_mul(f, f), _mul(h, h))))
    rpoly = _div_1mx2(dpoly)
    if min((_eval(rpoly, x).real for x in _grid()), default=0.0) < -1e-9:
        raise ValidationError("补多项式 R 在 [−1,1] 上取负值，谱分解不存在")
    q = _q_from_roots(rpoly, d)
    if q is None:
        raise ValidationError("补多项式谱因子的奇偶性与度数要求不符")
    ppoly = _add(f, tuple(1j * v for v in h))
    phases = _strip(ppoly, q, d)
    err = max(abs(qsp_response(x, phases) - _eval(ppoly, x)) for x in _grid(512))
    if err > _STRIP_TOL:
        raise ValidationError(f"相位合成往返自检失败（误差 {err:.2e}）：度数过高或补多项式病态")
    return phases


# ---------------------------------------------------------------------------
# 组装辅助：相位序列 → BE，以及实部提取 LCU。
# ---------------------------------------------------------------------------


def _wrap_qsvt_be(a, phases):
    return BlockEncoding(
        annotate(qsvt_sequence(a, phases), "block_encoding", be_alpha=1.0)
    )


def _real_qsvt_be(a, phases):
    """块编码 (P + P̄)(A/α)/2 = f(A/α)：−Φ 恰好实现 P̄，经 LCU 各半提取实部。"""
    plus = _wrap_qsvt_be(a, phases)
    minus = _wrap_qsvt_be(a, tuple(-p for p in phases))
    return linear_combination(0.5, plus, 0.5, minus)


def _finish(be, algorithm, **attributes):
    return BlockEncoding(
        annotate(
            be.operation,
            "block_encoding",
            be_alpha=be.alpha,
            algorithm=algorithm,
            **attributes,
        )
    )


def _synthesize_with_imag(f, imag_candidates, label):
    for h in imag_candidates:
        try:
            return qsp_phases(f, imag=h)
        except ValidationError:
            continue
    raise ValidationError(f"{label} 的虚部补全失败：请降低目标度数或放宽参数")


# ---------------------------------------------------------------------------
# 标准变换族。
# ---------------------------------------------------------------------------


def qsvt_matrix_inversion(a, kappa, *, error=0.05):
    """近似 A⁻¹ 的 QSVT 块编码（奇扩展多项式 J_b(x) = (1−(1−x²)^b)/x 的缩放）。

    目标多项式 f(x) = c·J_b(x)（升幂系数 ( −1)^m C(b, m+1) 解析给出），
    在 ``|x| ≥ 1/κ`` 上相对误差不超过 error 地逼近 c/x，``||f||∞ ≤ 1/3``。
    返回的 BE 的零信号块约为 inverse_scale · A⁻¹，inverse_scale = c·α。
    """
    require_instance(a, BlockEncoding, "qsvt_matrix_inversion.a")
    if not (math.isfinite(kappa) and kappa >= 1):
        raise ValidationError("条件数 κ 必须是不小于 1 的有限数")
    if not (0 < error < 1):
        raise ValidationError("近似误差 error 必须在 (0,1) 内")
    if kappa == 1:
        b = 1
    else:
        b = max(1, math.ceil(math.log(1 / error) / -math.log(1 - 1 / kappa**2)))
    d = 2 * b - 1
    if d > _MAX_DEGREE:
        raise ValidationError(f"κ={kappa} 需要度数 {d}，超过合成上限 {_MAX_DEGREE}；请放宽 error")
    f = tuple(
        ((-1.0) ** m) * math.comb(b, m + 1) if i % 2 == 1 else 0.0
        for i in range(d + 1)
        for m in [i // 2]
    )
    f = _trim(f)
    norm = _sup_norm(f)
    c_scale = 1.0 / (3.0 * norm)
    f = _scale(c_scale, f)
    sat = math.sqrt(max(0.0, 1.0 - _eval(f, 1.0).real ** 2))
    phases = _synthesize_with_imag(f, [(0.0, sat)], "矩阵求逆多项式")
    be = _real_qsvt_be(a, phases)
    return _finish(
        be,
        "qsvt_matrix_inversion",
        kappa=float(kappa),
        error=float(error),
        qsp_degree=d,
        inverse_scale=c_scale * a.alpha,
    )


def eigenstate_filter(a, gap, degree, *, center=0.0):
    """特征态过滤：块编码在 ``|x−center| ≤ gap`` 外被压到 1/T_d(r) 以下的尖峰多项式。

    目标为 Lin–Tong 型过滤多项式 f(x) = T_d(g(x²))/T_d(r)，
    g(y) = 2(y−Δ²)/(1−Δ²) − 1，r = (1+Δ²)/(1−Δ²)；f 在 x=0 处饱和（``|f(0)|=1``），
    在 ``|x| ≥ Δ`` 上 ``|f| ≤ 1/T_d(r)``。center 非零时先经 BE 线性组合平移谱。
    """
    require_instance(a, BlockEncoding, "eigenstate_filter.a")
    if not (0 < gap < 1):
        raise ValidationError("过滤宽度 gap 必须在 (0,1) 内")
    if type(degree) is not int or degree < 1:
        raise ValidationError("Chebyshev 度数必须为正整数")
    d2 = 2 * degree
    if d2 > _MAX_DEGREE:
        raise ValidationError(f"合成度数 {d2} 超过上限 {_MAX_DEGREE}")
    shifted = a
    if center != 0.0:
        if not math.isfinite(center) or abs(center) >= 1:
            raise ValidationError("过滤中心 center 必须在 (−1,1) 内")
        shifted = linear_combination(1.0, a, -complex(center), identity(a.width))
    r = (1 + gap**2) / (1 - gap**2)
    norm_d = math.cosh(degree * math.acosh(r))  # T_d(r)
    u = (-1.0 - 2 * gap**2 / (1 - gap**2), 0.0, 2.0 / (1 - gap**2))  # g(x²)
    t0, t1 = (1.0,), u
    for _ in range(2, degree + 1):
        t0, t1 = t1, _sub(_scale(2.0, _mul(u, t1)), t0)
    f = _scale(1.0 / norm_d, t1 if degree >= 1 else t0)
    sat = math.sqrt(max(0.0, 1.0 - _eval(f, 1.0).real ** 2))
    phases = _synthesize_with_imag(f, [(0.0, 0.0, sat)], "特征态过滤多项式")
    be = _real_qsvt_be(shifted, phases)
    return _finish(
        be,
        "eigenstate_filter",
        gap=float(gap),
        filter_degree=degree,
        center=float(center),
        qsp_degree=d2,
        suppression=1.0 / norm_d,
    )


def _bessel_j(n, x):
    """第一类 Bessel 函数 J_n(x)，幂级数纯 Python 实现。"""
    term = (x / 2) ** n / math.factorial(n)
    total = term
    m = 0
    while abs(term) > 1e-18 * max(1.0, abs(total)) and m < 100000:
        term *= -((x / 2) ** 2) / ((m + 1) * (m + n + 1))
        total += term
        m += 1
    return total


def _jacobi_anger(t, error):
    """e^{itx} 的 Jacobi–Anger 截断：返回 (偶支 cos 系数, 奇支 sin 系数, 截断度数 K)。"""
    kmax = min(_MAX_DEGREE, int(math.ceil(abs(t))) + 8 * int(math.ceil(math.log10(4 / error))) + 8)
    js = [_bessel_j(k, abs(t)) for k in range(kmax + 2)]
    suffix = [0.0] * (kmax + 3)
    for j in range(kmax + 1, -1, -1):
        suffix[j] = suffix[j + 1] + 2.0 * abs(js[j])
    k = next((j for j in range(kmax + 1) if 2 * suffix[j + 1] <= error / 4), kmax)
    sign = 1.0 if t >= 0 else -1.0
    fc = [0.0] * (k + 1)
    fs = [0.0] * (k + 1)
    for j in range(0, k + 1):
        tk = _chebyshev_t(j)
        if j % 2 == 0:
            coef = (1.0 if j == 0 else 2.0) * (-1.0) ** (j // 2) * js[j]
            for i, v in enumerate(tk):
                fc[i] += coef * v
        else:
            coef = 2.0 * (-1.0) ** ((j - 1) // 2) * js[j] * sign
            for i, v in enumerate(tk):
                fs[i] += coef * v
    return _trim(fc), _trim(fs), k


def qsvt_hamiltonian_simulation(a, t, *, error=0.01):
    """e^{itA/α} 的 QSVT 块编码：Jacobi–Anger 偶/奇两支分别合成，再经 LCU 组合。

    偶支近似 cos(tx)、奇支近似 sin(tx)，统一缩放 s 使两支均留出虚部补全余量；
    每支用 (U_Φ + U_{−Φ})/2 提取实部，最后按 1 与 i 做 LCU。
    返回 BE 的零信号块约为 e^{itA/α}/sim_scale，sim_scale = 2s。
    """
    require_instance(a, BlockEncoding, "qsvt_hamiltonian_simulation.a")
    if not (math.isfinite(t) and t != 0):
        raise ValidationError("演化时间 t 必须是非零有限实数")
    if not (0 < error < 1):
        raise ValidationError("近似误差 error 必须在 (0,1) 内")
    fc, fs, k = _jacobi_anger(t, error)
    s = 1.5 * max(_sup_norm(fc), _sup_norm(fs), 1e-3)
    fc, fs = _scale(1.0 / s, fc), _scale(1.0 / s, fs)
    dc, ds = _deg(fc), _deg(fs)
    if dc % 2 or ds % 2 == 0:
        raise ValidationError("Jacobi–Anger 分支奇偶性异常")

    def cos_imags():
        a0 = math.sqrt(max(0.0, 1.0 - _eval(fc, 0.0).real ** 2))
        a1 = math.sqrt(max(0.0, 1.0 - _eval(fc, 1.0).real ** 2))
        for m in range(1, dc // 2 + 1):
            yield tuple(
                (a0 if i == 0 else 0.0) + ((a1 - a0) if i == 2 * m else 0.0)
                for i in range(dc + 1)
            )

    def sin_imags():
        a1 = math.sqrt(max(0.0, 1.0 - _eval(fs, 1.0).real ** 2))
        for m in range(0, (ds - 1) // 2 + 1):
            yield tuple(a1 * (1.0 if i == 2 * m + 1 else 0.0) for i in range(ds + 1))

    phases_c = _synthesize_with_imag(fc, cos_imags(), "哈密顿模拟 cos 支")
    phases_s = _synthesize_with_imag(fs, sin_imags(), "哈密顿模拟 sin 支")
    uc = _real_qsvt_be(a, phases_c)
    us = _real_qsvt_be(a, phases_s)
    be = linear_combination(1.0, uc, 1j, us)
    return _finish(
        be,
        "qsvt_hamiltonian_simulation",
        time=float(t),
        error=float(error),
        qsp_degree=k,
        sim_scale=2.0 * s,
    )


def fixed_point_search_phases(delta, degree):
    """Yoder–Low–Chuang 定点振幅放大的相位序列（时间正序，长度 degree+1）。

    构造依据 YLC 闭式补多项式：记 L = degree（BE 调用次数，必须为奇数）、
    c = T_{1/L}(1/δ)，则 Q(x) = δc·R(c²(1−x²))，R(u) = T_L(√u)/√u 为解析多项式；
    P 由 1 − (1−x²)Q² 的求根谱分解得到。实现的成功概率恰为
    P_S(x) = 1 − δ² T_L²(c√(1−x²))：``|x| ≥ √(1−1/c²)`` 时 P_S ≥ 1 − δ²，
    且阈值随 L 单调下降趋于 0（不动点性质）。
    """
    if not (0 < delta < 1):
        raise ValidationError("定点搜索误差 δ 必须在 (0,1) 内")
    if type(degree) is not int or degree < 1 or degree % 2 == 0:
        raise ValidationError("定点搜索度数（BE 调用次数）必须为正奇数")
    if degree > _MAX_DEGREE // 2:
        raise ValidationError(f"定点搜索度数超过上限 {_MAX_DEGREE // 2}")
    L = degree
    c = math.cosh(math.acosh(1.0 / delta) / L)
    tcoeff = _chebyshev_t(L)
    qpoly = (0.0,)
    for m in range((L + 1) // 2):
        term = (tcoeff[2 * m + 1] * c ** (2 * m),)
        for _ in range(m):
            term = _mul(term, (1.0, 0.0, -1.0))
        qpoly = _add(qpoly, term)
    qpoly = _trim(_scale(delta * c, qpoly))
    fpoly = _trim(_sub((1.0,), _mul((1.0, 0.0, -1.0), _mul(qpoly, qpoly))))
    if abs(fpoly[0]) > 1e-8 or abs(fpoly[1]) > 1e-8:
        raise ValidationError("YLC 补多项式构造异常：零根缺失")
    ft = _trim(fpoly[2:])  # 除以 x²（零点二重根）
    clusters = _cluster(_roots(ft))
    factors = []
    used = set()
    ctol = 1e-5
    for i, (rt, m) in enumerate(clusters):
        if i in used:
            continue
        sib = [
            j
            for j, (c2, m2) in enumerate(clusters)
            if j > i
            and j not in used
            and m2 == m
            and (abs(c2 - rt.conjugate()) <= ctol * max(1.0, abs(rt))
                 or abs(c2 + rt) <= ctol * max(1.0, abs(rt))
                 or abs(c2 + rt.conjugate()) <= ctol * max(1.0, abs(rt)))
        ]
        if len(sib) != 3:
            raise ValidationError("YLC 谱分解根结构异常")
        for _ in range(m):
            factors.append((-rt * rt, 0.0, 1.0))
        used.update(sib)
        used.add(i)
    if 1 + 2 * len(factors) != L:
        raise ValidationError("YLC 谱分解因子计数异常")
    ppoly = (0.0, 1.0)
    for fct in factors:
        ppoly = _mul(ppoly, fct)
    ppoly = _trim(_scale(abs(qpoly[-1]), ppoly))
    phases = _strip(ppoly, qpoly, L)
    err = max(
        abs(abs(qsp_response(x, phases)) ** 2 - (1.0 - (1 - x * x) * _eval(qpoly, x) ** 2))
        for x in _grid(512)
    )
    if err > _STRIP_TOL:
        raise ValidationError(f"定点搜索相位自检失败（误差 {err:.2e}）")
    return phases


def fixed_point_search(a, delta, degree):
    """定点振幅放大的 QSVT 组装：对 BE 应用 fixed_point_search_phases 的序列。

    零信号块为复多项式 P(A/α)，成功概率 ``|P(x)|²`` 满足 YLC 不动点保证；
    threshold 属性给出 √(1−1/c²) 的放大阈值。
    """
    require_instance(a, BlockEncoding, "fixed_point_search.a")
    phases = fixed_point_search_phases(delta, degree)
    L = degree
    gamma = 1.0 / math.cosh(math.acosh(1.0 / delta) / L)
    be = _wrap_qsvt_be(a, phases)
    return _finish(
        be,
        "fixed_point_search",
        delta=float(delta),
        qsp_degree=L,
        threshold=math.sqrt(max(0.0, 1.0 - gamma * gamma)),
    )
