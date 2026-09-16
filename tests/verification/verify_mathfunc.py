"""数学函数廊与 frozen Roe 面的论文级数值验证。

覆盖 examples/math_functions.py 的 5 个函数（pressure / roe_speed / phase_response /
guarded_reciprocal / polynomial）与 applications/roe_formulas.py 的 frozen_roe_face，
全部经 pyqecclang.compile_function 编译为定点可逆模块后在真实后端上执行：

- rir-pysparq（PySparQ 原生 RIR 解释器）为主力路径，叠加态一次穷举输入域；
- reference / adapter-pysparq 在代表性程序上做振幅级三方对拍；
- OriginIR-ext 态向量路径因编译函数工作区远超 24 量子比特预算而不适用
  （workspace_table 实测 10^2–10^3 比特，各 case 的 parameters 记录实测值）。

oracle 独立性：期望值由本脚本内用 float64（math/cmath）按数学定义独立转写的公式
给出，不导入被测实现；phase_response 另把模块属性 math_approximation 中的
Chebyshev 系数作为"实现参考"、真函数作为"方法参考"，实现误差与方法误差分开报告。
status 旗标语义（位 0 = 定义域失效、位 1 = 值域/字长越界）逐分支核对。

运行：PYTHONPATH=src <含 pysparq+uniqc 的 python> tests/verification/verify_mathfunc.py
环境变量 VERIFY_WORKERS 控制并行进程数（默认 10；本组门级仿真单分支最重约 2.5 s，
靠多进程把总墙钟压进 5 分钟）。所有执行均为叠加态穷举，不做逐基态循环。
"""

from __future__ import annotations

import cmath
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

from harness import (
    ROOT,
    Report,
    adapter_pysparq,
    amplitude_error,
    reference,
    rir_pysparq,
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))  # 使 worker 进程可导入 examples/ 下的被测函数

from pyqecclang import Builder, FixedFormat
from pyqecclang.infrastructure.mathfunc import Index, MathConfig, compile_function

FMTS = {"6.2": (6, 2), "8.3": (8, 3)}

# phase_response 的 Chebyshev 阶数：默认 degree=6 的核门数对门级仿真过重，取
# degree=3 控制预算；方法误差（多项式 vs 真函数）单独报告，不混入实现误差。
PHASE_DEGREE = 3

# 各函数量化容差（单位：量子），按截断算子（mul/div/sqrt/核求值）数量与传播放宽。
TOLERANCE = {
    "polynomial": 0.0,
    "guarded_reciprocal": 1.0,
    "pressure": 0.0,
    "roe_speed": 0.0,
    "phase_response": 0.0,
    "frozen_roe_face": 0.0,
}

# ---------------------------------------------------------------------------
# 被测函数的编译规格（worker 进程内按需编译；编译确定且廉价）
# ---------------------------------------------------------------------------


def _compile(fn_key, fmt):
    if fn_key == "pressure":
        from examples.math_functions import pressure

        return compile_function(pressure, fmt=fmt, constants={"gamma": 1.4})
    if fn_key == "roe_speed":
        from examples.math_functions import roe_speed

        return compile_function(roe_speed, fmt=fmt)
    if fn_key == "guarded_reciprocal":
        from examples.math_functions import guarded_reciprocal

        return compile_function(guarded_reciprocal, fmt=fmt)
    if fn_key == "polynomial":
        from examples.math_functions import polynomial

        return compile_function(polynomial, fmt=fmt, constants={"order": 3})
    if fn_key == "phase_response":
        from examples.math_functions import phase_response

        return compile_function(
            phase_response, fmt=fmt, config=MathConfig(degree=PHASE_DEGREE)
        )
    if fn_key == "frozen_roe_face":
        from pyqecclang.applications.roe_formulas import frozen_roe_face

        return compile_function(
            frozen_roe_face,
            fmt=fmt,
            inputs={
                **{k: "real" for k in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r")},
                "row": Index(2),
                "col": Index(2),
            },
            constants={"gamma": 1.4, "entropy_delta": 0.125},
            output_names=("left", "right"),
        )
    raise KeyError(fn_key)


def _workspace(fn_key, fmt):
    """编译程序的工作区峰值（locals 叠调用链最大值），用于 OriginIR 预算说明。"""
    from pyqecclang.infrastructure.layout import workspace_table

    program = _compile(fn_key, fmt).program()
    return workspace_table(program)[program.entry]


# ---------------------------------------------------------------------------
# 独立 float64 oracle（按数学定义独立转写）
# 返回 (期望输出|None, 参考输出|None, 期望 status|None)：
#   期望输出 None —— 该分支不比对数值（已旗标或灰区）；
#   期望 status：位 0 = 定义域失效，位 1 = 值域/字长越界；None = 阈值灰区不比对；
#   参考输出仅 phase_response 使用（真函数值，用于方法误差）。
# ---------------------------------------------------------------------------


def _limit(fmt):
    return 2.0 ** (fmt.width - fmt.fraction - 1)


def _quantize(fmt, value):
    """生成期常量的定点编码-解码往返（encode 向零截断）。"""
    return fmt.decode(fmt.encode(value))


class _Fx:
    """fixed_arithmetic 的逐比特仿真（raw 整数域）：
    mul/div/sqrt 幅度向零截断，add/sub 模 wrap，status 位 0 = 定义域、位 1 = 越界，
    select 只传播被选中分支的旗标。值用 (raw, flag) 二元组沿路径携带。"""

    def __init__(self, fmt):
        self.w, self.f = fmt.width, fmt.fraction
        self.sign = 1 << (self.w - 1)
        self.mask = (1 << self.w) - 1

    def wrap(self, raw):
        raw &= self.mask
        return raw - (1 << self.w) if raw & self.sign else raw

    def const(self, x):
        # 与 fmt.encode 一致：int() 向零截断后取模；此处常量均在字长范围内
        return self.wrap(int(x * (1 << self.f))), 0

    def input(self, x):
        # 网格值（decode 过的输入）回 raw；浮点往返精确
        return int(round(x * (1 << self.f))), 0

    def _ovf(self, mag):
        return 2 if mag >= self.sign else 0

    def mul(self, a, b):
        ra, fa = a
        rb, fb = b
        mag = (abs(ra) * abs(rb)) >> self.f
        out = mag if (ra < 0) == (rb < 0) else -mag
        return self.wrap(out), fa | fb | self._ovf(mag)

    def div(self, a, b):
        ra, fa = a
        rb, fb = b
        if rb == 0:
            return 0, fa | fb | 1
        mag = (abs(ra) << self.f) // abs(rb)
        out = mag if (ra < 0) == (rb < 0) else -mag
        return self.wrap(out), fa | fb | self._ovf(mag)

    def add(self, a, b):
        ra, fa = a
        rb, fb = b
        out = self.wrap(ra + rb)
        overflow = ((ra >= 0) == (rb >= 0)) and ((out >= 0) != (ra >= 0))
        return out, fa | fb | (2 if overflow else 0)

    def sub(self, a, b):
        ra, fa = a
        rb, fb = b
        out = self.wrap(ra - rb)
        overflow = ((ra >= 0) != (rb >= 0)) and ((out >= 0) != (ra >= 0))
        return out, fa | fb | (2 if overflow else 0)

    def neg(self, a):
        ra, fa = a
        return self.wrap(-ra), fa | (2 if ra == -self.sign else 0)

    def abs_(self, a):
        ra, fa = a
        return self.wrap(abs(ra)), fa | (2 if ra == -self.sign else 0)

    def sqrt(self, a):
        ra, fa = a
        if ra < 0:
            return 0, fa | 1
        mag = math.isqrt(ra << self.f)
        return self.wrap(mag), fa | self._ovf(mag)

    def ge(self, a, b):
        # 布尔比较（电路为 lt 取反；比较本身无旗标）
        return (1 if a[0] >= b[0] else 0), a[1] | b[1]

    def select(self, test, yes, no):
        # 两分支都被计算，但被选中分支的（值, 旗标）才传播
        chosen = yes if test[0] else no
        return chosen[0], test[1] | chosen[1]

    def decode(self, a):
        return a[0] / (1 << self.f)


def oracle_polynomial(v, fmt):
    # 截断链精确仿真：x*x 幅度向零截断至量子网格，加法无截断；溢出按字长判定。
    x = v["x"]
    q = 1.0 / (1 << fmt.fraction)
    xx = math.trunc(x * x / q) * q
    val = xx + x + 1.0
    lim = _limit(fmt)
    flag = 2 if (x * x >= lim or not -lim <= val < lim) else 0
    return ({"out": val} if flag == 0 else None), None, flag


def oracle_guarded_reciprocal(v, fmt):
    x = v["x"]
    if x == 0:
        return {"out": 0.0}, None, 0  # 守护分支吸收奇点，未选分支除零被掩码
    val = 1.0 / x
    flag = 2 if abs(val) >= _limit(fmt) else 0
    return ({"out": val} if flag == 0 else None), None, flag


def oracle_pressure(v, fmt):
    # 截断链精确仿真（与 polynomial 同法）：div/mul 幅度向零截断一个量子，add/sub
    # 无截断；gamma-1 经前端常量折叠为 0.4 再编码（6.2 → 0.25，8.3 → 0.375）。
    rho, m, e = v["rho"], v["momentum"], v["energy"]
    if rho == 0:
        return None, None, 1  # 除零 → 位 0
    q = 1.0 / (1 << fmt.fraction)

    def trunc(x):
        return math.trunc(x / q) * q

    g1 = _quantize(fmt, 0.4)
    vel = trunc(m / rho)
    kin = trunc(trunc(0.5 * m) * vel)
    inner = e - kin
    val = trunc(g1 * inner)
    lim = _limit(fmt)
    flag = (
        2
        if (
            abs(m / rho) >= lim
            or abs(trunc(0.5 * m) * vel) >= lim
            or not -lim <= inner < lim
            or abs(g1 * inner) >= lim
        )
        else 0
    )
    exact = 0.4 * (e - 0.5 * m * (m / rho))
    return ({"out": val} if flag == 0 else None), {"out": exact}, flag


def oracle_roe_speed(v, fmt):
    # _Fx 逐比特仿真：sqrt 为 isqrt(raw<<f)，mul/div 向零截断，旗标语义同电路。
    F = _Fx(fmt)
    rl, rr = v["rho_l"], v["rho_r"]
    if rl <= 0 or rr <= 0:
        return None, None, 1  # rho<0 平方根定义域 / rho=0 除零（电路同置位 0）
    left = F.sqrt(F.input(rl))
    right = F.sqrt(F.input(rr))
    b = F.div(F.mul(left, F.input(v["momentum_l"])), F.input(rl))
    d = F.div(F.mul(right, F.input(v["momentum_r"])), F.input(rr))
    val = F.div(F.add(b, d), F.add(left, right))
    status = val[1]
    # 方法参考：float64 原式（含 sqrt 精确值）
    lw, rw = math.sqrt(rl), math.sqrt(rr)
    exact = (lw * v["momentum_l"] / rl + rw * v["momentum_r"] / rr) / (lw + rw)
    return ({"out": F.decode(val)} if status == 0 else None), {"out": exact}, status


def _pick3(index, first, second, third):
    return (first, second, third, 0.0)[index] if 0 <= index <= 3 else 0.0


def _entropy_abs(lam, delta):
    value = abs(lam)
    if delta <= 0 or value >= delta:
        return value
    return (lam * lam + delta * delta) / (2 * delta)


def _roe_face_formula(v, consts):
    """Roe 面公式主体；consts = (g1, gm, gm3, tmg, delta) 分别为
    (gamma-1, gamma, gamma-3, 3-gamma, entropy_delta) 的取值。"""
    g1, gm, gm3, tmg, delta = consts
    rl, ml, el = v["rho_l"], v["m_l"], v["e_l"]
    rr, mr, er = v["rho_r"], v["m_r"], v["e_r"]
    row, col = v["row"], v["col"]
    ul = ml / rl
    pl = g1 * (el - 0.5 * ml * ul)
    hl = (el + pl) / rl
    ur = mr / rr
    pr = g1 * (er - 0.5 * mr * ur)
    hr = (er + pr) / rr
    wl, wr = math.sqrt(rl), math.sqrt(rr)
    total = wl + wr
    u = (wl * ul + wr * ur) / total
    h = (wl * hl + wr * hr) / total
    u2 = u * u
    c2 = g1 * (h - 0.5 * u2)
    return rl, rr, u, h, u2, c2, ul, hl, ur, hr, pl, pr, row, col


def _roe_face_outputs(parts, consts):
    g1, gm, gm3, tmg, delta = consts
    rl, rr, u, h, u2, c2, ul, hl, ur, hr, pl, pr, row, col = parts
    c = math.sqrt(c2)
    beta = g1 / c2
    inverse_c = 0.5 / c
    bu, bu2 = beta * u, beta * u2
    r0 = _pick3(row, 1.0, u - c, h - u * c)
    r1 = _pick3(row, 1.0, u, 0.5 * u2)
    r2 = _pick3(row, 1.0, u + c, h + u * c)
    i0 = _pick3(col, 0.25 * bu2 + u * inverse_c, -0.5 * bu - inverse_c, 0.5 * beta)
    i1 = _pick3(col, 1 - 0.5 * bu2, bu, -beta)
    i2 = _pick3(col, 0.25 * bu2 - u * inverse_c, -0.5 * bu + inverse_c, 0.5 * beta)
    a0 = _entropy_abs(u - c, delta)
    a1 = _entropy_abs(u, delta)
    a2 = _entropy_abs(u + c, delta)
    absolute = r0 * a0 * i0 + r1 * a1 * i1 + r2 * a2 * i2

    def euler(velocity, enthalpy):
        square = velocity * velocity
        first = _pick3(col, 0.0, 1.0, 0.0)
        second = _pick3(col, 0.5 * gm3 * square, tmg * velocity, g1)
        third = _pick3(
            col,
            velocity * (0.5 * g1 * square - enthalpy),
            enthalpy - g1 * square,
            gm * velocity,
        )
        return _pick3(row, first, second, third)

    left = 0.5 * (euler(ul, hl) + absolute)
    right = 0.5 * (euler(ur, hr) - absolute)
    mags = [
        abs(x)
        for x in (
            ul, pl, hl, ur, pr, hr, u, h, u2, c2, beta, inverse_c, bu, bu2,
            r0, r1, r2, i0, i1, i2, a0, a1, a2, absolute, left, right,
        )
    ]
    return {"left": left, "right": right}, mags


def _roe_face_vm(v, fmt):
    """frozen_roe_face 的 _Fx 逐比特仿真；按 roe_formulas.py 源码顺序转写。

    旗标语义的两大要点（对照 lowering.py / numeric.py）：
    - helper 调用把**全部参数**的 status 与模块内 status 合并（combined），
      因此 pick3/euler_entry 等 helper 的未选中参数分支旗标不被掩码；
    - helper 形参（gamma、delta）在 helper 模块内是动态输入，不参与前端常量折叠
      （如 euler_entry 内 3-gamma = 3-1.25 = 1.75，而非折叠值 1.5）。
    返回 (left, right, flag)。"""
    F = _Fx(fmt)
    rl, ml, el = (F.input(v[k]) for k in ("rho_l", "m_l", "e_l"))
    rr, mr, er = (F.input(v[k]) for k in ("rho_r", "m_r", "e_r"))
    row, col = v["row"], v["col"]
    half = F.const(0.5)
    gm = F.const(1.4)  # gamma 作为 helper 实参被物化为常量寄存器
    g1 = F.const(0.4)  # 顶层函数体内 gamma-1 经前端折叠后编码

    def primitive(rho, m, e):
        gm1 = F.sub(gm, F.const(1.0))  # helper 内动态计算 gamma-1
        u = F.div(m, rho)
        p = F.mul(gm1, F.sub(e, F.mul(F.mul(half, m), u)))
        h = F.div(F.add(e, p), rho)
        merged = rho[1] | m[1] | e[1] | gm[1] | u[1] | p[1] | h[1]
        return (u[0], merged), (p[0], merged), (h[0], merged)

    def pick3(index, first, second, third):
        chosen = (first, second, third, F.const(0.0))[index] if 0 <= index <= 3 else F.const(0.0)
        return chosen[0], first[1] | second[1] | third[1]

    def entropy_abs(lam, delta):
        value = F.abs_(lam)
        smooth = F.div(
            F.add(F.mul(lam, lam), F.mul(delta, delta)), F.mul(F.const(2.0), delta)
        )
        rest = F.select(F.ge(value, delta), value, smooth)
        # `delta <= 0` ⟺ `0 >= delta`；6.2 下 delta raw = 0 → 恒取 value 分支
        out = F.select(F.ge(F.const(0.0), delta), value, rest)
        return out[0], lam[1] | delta[1] | out[1]

    def euler(vel, ent):
        gm1 = F.sub(gm, F.const(1.0))
        square = F.mul(vel, vel)
        first = pick3(col, F.const(0.0), F.const(1.0), F.const(0.0))
        second = pick3(
            col,
            F.mul(F.mul(half, F.sub(gm, F.const(3.0))), square),
            F.mul(F.sub(F.const(3.0), gm), vel),
            gm1,
        )
        third = pick3(
            col,
            F.mul(vel, F.sub(F.mul(F.mul(half, gm1), square), ent)),
            F.sub(ent, F.mul(gm1, square)),
            F.mul(gm, vel),
        )
        out = pick3(row, first, second, third)
        return out[0], vel[1] | ent[1] | gm[1] | out[1]

    ul, pl, hl = primitive(rl, ml, el)
    ur, pr, hr = primitive(rr, mr, er)
    wl = F.sqrt(rl)
    wr = F.sqrt(rr)
    total = F.add(wl, wr)
    u = F.div(F.add(F.mul(wl, ul), F.mul(wr, ur)), total)
    h = F.div(F.add(F.mul(wl, hl), F.mul(wr, hr)), total)
    u2 = F.mul(u, u)
    c2 = F.mul(g1, F.sub(h, F.mul(half, u2)))
    c = F.sqrt(c2)
    beta = F.div(g1, c2)
    inverse_c = F.div(half, c)
    bu = F.mul(beta, u)
    bu2 = F.mul(beta, u2)
    r0 = pick3(row, F.const(1.0), F.sub(u, c), F.sub(h, F.mul(u, c)))
    r1 = pick3(row, F.const(1.0), u, F.mul(half, u2))
    r2 = pick3(row, F.const(1.0), F.add(u, c), F.add(h, F.mul(u, c)))
    i0 = pick3(
        col,
        F.add(F.mul(F.const(0.25), bu2), F.mul(u, inverse_c)),
        F.sub(F.mul(F.const(-0.5), bu), inverse_c),
        F.mul(half, beta),
    )
    i1 = pick3(col, F.sub(F.const(1.0), F.mul(half, bu2)), bu, F.neg(beta))
    i2 = pick3(
        col,
        F.sub(F.mul(F.const(0.25), bu2), F.mul(u, inverse_c)),
        F.add(F.mul(F.const(-0.5), bu), inverse_c),
        F.mul(half, beta),
    )
    delta = F.const(0.125)
    a0 = entropy_abs(F.sub(u, c), delta)
    a1 = entropy_abs(u, delta)
    a2 = entropy_abs(F.add(u, c), delta)
    absolute = F.add(
        F.add(F.mul(F.mul(r0, a0), i0), F.mul(F.mul(r1, a1), i1)),
        F.mul(F.mul(r2, a2), i2),
    )
    left = F.mul(half, F.add(euler(ul, hl), absolute))
    right = F.mul(half, F.sub(euler(ur, hr), absolute))
    flag = left[1] | right[1]
    return left, right, flag


def oracle_frozen_roe_face(v, fmt):
    if v["rho_l"] <= 0 or v["rho_r"] <= 0:
        return None, None, 1  # rho<0 平方根定义域 / rho=0 除零
    left, right, flag = _roe_face_vm(v, fmt)
    if flag != 0:
        return None, None, flag
    outputs = {
        "left": left[0] / (1 << fmt.fraction),
        "right": right[0] / (1 << fmt.fraction),
    }
    try:
        exact = _roe_face_outputs(
            _roe_face_formula(v, (0.4, 1.4, -1.6, 1.6, 0.125)),
            (0.4, 1.4, -1.6, 1.6, 0.125),
        )[0]
    except (ValueError, ZeroDivisionError):
        exact = None  # 精确常量公式在该点奇异（如 c2≤0），方法误差点不计
    return outputs, exact, 0


def _kernel_fx(F, coeffs, kind, a):
    """初等核模块的 _Fx 仿真：区间常量经编码截断，Clenshaw 递推按 numeric.py 的
    乘加顺序逐比特复现；越出区间置位 1（对照模块属性 math_approximation 的系数）。"""
    (lo, hi), cs = coeffs[kind]
    t = F.mul(F.sub(a, F.const((lo + hi) / 2)), F.const(1 / ((hi - lo) / 2)))
    nxt = F.const(0.0)
    nx2 = F.const(0.0)
    for coefficient in reversed(cs[1:]):
        nxt, nx2 = (
            F.add(F.sub(F.mul(F.const(2.0), F.mul(t, nxt)), nx2), F.const(coefficient)),
            nxt,
        )
    value = F.add(F.sub(F.mul(t, nxt), nx2), F.const(cs[0]))
    outside = 2 if (a[0] < F.const(lo)[0] or F.const(hi)[0] < a[0]) else 0
    return value[0], value[1] | outside


def _phase_vm(v, fmt, coeffs):
    """phase_response 的 _Fx 逐比特仿真（含三个初等核）。"""
    F = _Fx(fmt)
    x = F.input(v["x"])
    y = F.input(v["y"])
    # w = e^{i z}：1j·z = (−y, x)（乘 0/1 精确）
    mag = _kernel_fx(F, coeffs, "exp", F.sub(F.const(0.0), y))
    wre = F.mul(mag, _kernel_fx(F, coeffs, "cos", x))
    wim = F.mul(mag, _kernel_fx(F, coeffs, "sin", x))
    zzre = F.sub(F.mul(x, x), F.mul(y, y))
    zzim = F.add(F.mul(x, y), F.mul(y, x))
    bre = F.add(zzre, F.const(1.0))
    den = F.add(F.mul(bre, bre), F.mul(zzim, zzim))
    ore = F.div(F.add(F.mul(wre, bre), F.mul(wim, zzim)), den)
    oim = F.div(F.sub(F.mul(wim, bre), F.mul(wre, zzim)), den)
    return ore, oim, ore[1] | oim[1]


def make_phase_oracle(coeffs):
    """phase_response oracle：_Fx 仿真（实现参考，对照系数配方）+ 真函数（方法参考）。"""

    def oracle(v, fmt):
        ore, oim, flag = _phase_vm(v, fmt, coeffs)
        if flag != 0:
            return None, None, flag
        scale = 1 << fmt.fraction
        outputs = {"out_real": ore[0] / scale, "out_imag": oim[0] / scale}
        z = complex(v["x"], v["y"])
        true = cmath.exp(1j * z) / (1 + z * z)
        return outputs, {"out_real": true.real, "out_imag": true.imag}, 0

    return oracle


def _phase_coefficients(fmt):
    coeffs = {}
    for module in _compile("phase_response", fmt).program().modules:
        for key, value in module.attributes:
            if key == "math_approximation":
                data = json.loads(value)
                coeffs[data["function"]] = (tuple(data["interval"]), tuple(data["coefficients"]))
    if set(coeffs) != {"exp", "sin", "cos"}:
        raise AssertionError(f"phase_response 核系数缺失：{sorted(coeffs)}")
    return coeffs


# ---------------------------------------------------------------------------
# 执行 worker（顶层函数，可 pickled；每个任务 = 一次后端执行）
# ---------------------------------------------------------------------------


def execute(task):
    fmt = FixedFormat(*task["fmt"])
    compiled = _compile(task["fn"], fmt)
    operation = compiled.operation
    safe = "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in task["task"])
    builder = Builder("verify_" + safe, {r.name: r.type for r in operation.module.registers})
    nbits = 0
    for name, bits in task["sweeps"].items():
        for i in bits:
            builder.h(builder[name][i])
            nbits += 1
    for name, raw in task["fixed"].items():
        for i in range(builder[name].width):
            if (raw >> i) & 1:
                builder.x(builder[name][i])
    builder.call(operation, **{r.name: builder[r.name] for r in operation.module.registers})
    program = builder.finish().program()
    runner = {"rir": rir_pysparq, "reference": reference, "adapter": adapter_pysparq}[
        task["path"]
    ]
    kwargs = {"max_states": task["max_states"]}
    if task.get("max_steps"):
        # reference/adapter 的静态展开估计（受控区按位宽加权）高于实际事件数，
        # phase/roe_face 的静态估计约 1.2e6/2.5e6，默认 1e6 预算不足，故显式放宽。
        kwargs["max_steps"] = task["max_steps"]
    state = runner(program, **kwargs)
    return task["task"], nbits, {k: complex(v) for k, v in state.items()}


# ---------------------------------------------------------------------------
# 任务清单
# ---------------------------------------------------------------------------


def _task(name, fn, fmt_key, sweeps, fixed, path="rir", max_steps=None):
    nbits = sum(len(bits) for bits in sweeps.values())
    return {
        "task": name,
        "fn": fn,
        "fmt": FMTS[fmt_key],
        "sweeps": sweeps,
        "fixed": fixed,
        "path": path,
        "max_states": max(64, 1 << (nbits + 2)),
        "max_steps": max_steps,
    }


def build_tasks():
    fmt62, fmt83 = FixedFormat(6, 2), FixedFormat(8, 3)
    e62, e83 = fmt62.encode, fmt83.encode
    tasks = []

    def add(name, fn, fmt_key, sweeps, fixed, path="rir", max_steps=None):
        tasks.append(_task(name, fn, fmt_key, sweeps, fixed, path, max_steps))

    # polynomial / guarded_reciprocal：单输入全域穷举（64 / 256 分支）+ 三方对拍
    for fn in ("polynomial", "guarded_reciprocal"):
        add(f"{fn}-full-6.2", fn, "6.2", {"x": list(range(6))}, {})
        add(f"{fn}-full-8.3", fn, "8.3", {"x": list(range(8))}, {})
        add(f"{fn}-cross-ref", fn, "6.2", {"x": list(range(6))}, {}, "reference")
        add(f"{fn}-cross-adp", fn, "6.2", {"x": list(range(6))}, {}, "adapter")

    # pressure：三轴纤维穷举 + 低位联合立方体（4×4×4）
    pressure_fix = {
        "rho": ({"momentum": e62(0.5), "energy": e62(2.5)}, {"momentum": e83(0.5), "energy": e83(2.5)}),
        "momentum": ({"rho": e62(1.0), "energy": e62(2.5)}, {"rho": e83(1.0), "energy": e83(2.5)}),
        "energy": ({"rho": e62(1.0), "momentum": e62(0.5)}, {"rho": e83(1.0), "momentum": e83(0.5)}),
    }
    for axis, (f62, f83) in pressure_fix.items():
        add(f"pressure-fiber-{axis}-6.2", "pressure", "6.2", {axis: list(range(6))}, f62)
        add(f"pressure-fiber-{axis}-8.3", "pressure", "8.3", {axis: list(range(8))}, f83)
    cube = {"rho": [0, 1], "momentum": [0, 1], "energy": [0, 1]}
    base62 = {"rho": e62(1.0), "momentum": 0, "energy": e62(2.0)}
    base83 = {"rho": e83(1.0), "momentum": 0, "energy": e83(2.0)}
    add("pressure-cube-6.2", "pressure", "6.2", cube, base62)
    add("pressure-cube-8.3", "pressure", "8.3", cube, base83)
    add("pressure-cross-ref", "pressure", "6.2", cube, base62, "reference")
    add("pressure-cross-adp", "pressure", "6.2", cube, base62, "adapter")

    # roe_speed：四轴纤维穷举 + 联合立方体（2^4）
    roe_fix62 = {"rho_l": e62(1.0), "momentum_l": e62(0.5), "rho_r": e62(1.5), "momentum_r": e62(-0.25)}
    roe_fix83 = {"rho_l": e83(1.0), "momentum_l": e83(0.5), "rho_r": e83(1.5), "momentum_r": e83(-0.25)}
    for axis in ("rho_l", "momentum_l", "rho_r", "momentum_r"):
        f62 = {k: v for k, v in roe_fix62.items() if k != axis}
        f83 = {k: v for k, v in roe_fix83.items() if k != axis}
        add(f"roe_speed-fiber-{axis}-6.2", "roe_speed", "6.2", {axis: list(range(6))}, f62)
        add(f"roe_speed-fiber-{axis}-8.3", "roe_speed", "8.3", {axis: list(range(8))}, f83)
    rs_cube = {"rho_l": [0], "momentum_l": [0], "rho_r": [0], "momentum_r": [0]}
    rs_base62 = {"rho_l": e62(1.0), "momentum_l": 0, "rho_r": e62(1.0), "momentum_r": 0}
    rs_base83 = {"rho_l": e83(1.0), "momentum_l": 0, "rho_r": e83(1.0), "momentum_r": 0}
    add("roe_speed-cube-6.2", "roe_speed", "6.2", rs_cube, rs_base62)
    add("roe_speed-cube-8.3", "roe_speed", "8.3", rs_cube, rs_base83)
    add("roe_speed-cross-ref", "roe_speed", "6.2", rs_cube, rs_base62, "reference")
    add("roe_speed-cross-adp", "roe_speed", "6.2", rs_cube, rs_base62, "adapter")

    # phase_response：实/虚轴纤维穷举 + 奇异线 y=1（6.2）+ 联合立方体（4×4）
    add("phase-fiber-re-6.2", "phase_response", "6.2", {"z_real": list(range(6))}, {"z_imag": e62(0.5)})
    add("phase-fiber-im-6.2", "phase_response", "6.2", {"z_imag": list(range(6))}, {"z_real": e62(0.5)})
    add("phase-fiber-sg-6.2", "phase_response", "6.2", {"z_real": list(range(6))}, {"z_imag": e62(1.0)})
    add("phase-fiber-re-8.3", "phase_response", "8.3", {"z_real": list(range(8))}, {"z_imag": e83(0.5)})
    add("phase-fiber-im-8.3", "phase_response", "8.3", {"z_imag": list(range(8))}, {"z_real": e83(0.5)})
    ph_cube = {"z_real": [0, 1], "z_imag": [0, 1]}
    ph_base62 = {"z_real": e62(0.5), "z_imag": e62(0.5)}
    ph_base83 = {"z_real": e83(0.5), "z_imag": e83(0.5)}
    add("phase-cube-6.2", "phase_response", "6.2", ph_cube, ph_base62)
    add("phase-cube-8.3", "phase_response", "8.3", ph_cube, ph_base83)
    add("phase-cross-ref", "phase_response", "6.2", ph_cube, ph_base62, "reference", 8_000_000)
    add("phase-cross-adp", "phase_response", "6.2", ph_cube, ph_base62, "adapter", 8_000_000)

    # frozen_roe_face：Index 全域联合（row×col 16 组）+ 六实轴全幅值网格纤维
    roe_face_fix62 = {
        "rho_l": e62(1.0), "m_l": e62(0.5), "e_l": e62(2.5),
        "rho_r": e62(1.5), "m_r": e62(-0.25), "e_r": e62(3.0), "row": 1, "col": 2,
    }
    roe_face_fix83 = {
        "rho_l": e83(1.0), "m_l": e83(0.5), "e_l": e83(2.5),
        "rho_r": e83(1.5), "m_r": e83(-0.25), "e_r": e83(3.0), "row": 1, "col": 2,
    }
    idx_fix62 = {k: v for k, v in roe_face_fix62.items() if k not in ("row", "col")}
    idx_fix83 = {k: v for k, v in roe_face_fix83.items() if k not in ("row", "col")}
    add("roe-face-index-6.2", "frozen_roe_face", "6.2", {"row": [0, 1], "col": [0, 1]}, idx_fix62)
    add("roe-face-index-8.3", "frozen_roe_face", "8.3", {"row": [0, 1], "col": [0, 1]}, idx_fix83)
    # 跨后端对拍用 8 分支缩小版（reference/adapter 较慢）
    add("roe-face-cross-rir", "frozen_roe_face", "6.2", {"row": [0, 1], "col": [0]}, idx_fix62)
    add("roe-face-cross-ref", "frozen_roe_face", "6.2", {"row": [0, 1], "col": [0]}, idx_fix62, "reference", 8_000_000)
    add("roe-face-cross-adp", "frozen_roe_face", "6.2", {"row": [0, 1], "col": [0]}, idx_fix62, "adapter", 8_000_000)
    for axis in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r"):
        f62 = {k: v for k, v in roe_face_fix62.items() if k != axis}
        f83 = {k: v for k, v in roe_face_fix83.items() if k != axis}
        # 全幅值均匀网格：6.2 扫 bits[2..6)（步长 1.0），8.3 扫 bits[4..8)（步长 2.0）
        add(f"roe-face-fiber-{axis}-6.2", "frozen_roe_face", "6.2", {axis: [2, 3, 4, 5]}, f62)
        add(f"roe-face-fiber-{axis}-8.3", "frozen_roe_face", "8.3", {axis: [4, 5, 6, 7]}, f83)

    # 基态确定性抽点（单基态 → 单基态，输入保持、输出符合参考）
    basis_points = {
        "polynomial": {"x": e62(1.5)},
        "guarded_reciprocal": {"x": e62(0.5)},
        "pressure": {"rho": e62(1.0), "momentum": e62(0.5), "energy": e62(2.5)},
        "roe_speed": {"rho_l": e62(1.0), "momentum_l": e62(0.5), "rho_r": e62(1.5), "momentum_r": e62(-0.25)},
        "phase_response": {"z_real": e62(0.5), "z_imag": e62(0.25)},
        "frozen_roe_face": {
            "rho_l": e62(1.0), "m_l": e62(0.5), "e_l": e62(2.5),
            "rho_r": e62(1.5), "m_r": e62(-0.25), "e_r": e62(3.0), "row": 1, "col": 2,
        },
    }
    for fn, fixed in basis_points.items():
        add(f"basis-{fn}", fn, "6.2", {}, fixed)
    return tasks


# ---------------------------------------------------------------------------
# 结果评估
# ---------------------------------------------------------------------------


def _decoded_inputs(fn, fmt, values):
    if fn in ("polynomial", "guarded_reciprocal"):
        return {"x": fmt.decode(values["x"])}
    if fn == "pressure":
        return {k: fmt.decode(values[k]) for k in ("rho", "momentum", "energy")}
    if fn == "roe_speed":
        return {k: fmt.decode(values[k]) for k in ("rho_l", "momentum_l", "rho_r", "momentum_r")}
    if fn == "phase_response":
        return {"x": fmt.decode(values["z_real"]), "y": fmt.decode(values["z_imag"])}
    if fn == "frozen_roe_face":
        out = {k: fmt.decode(values[k]) for k in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r")}
        out["row"], out["col"] = values["row"], values["col"]
        return out
    raise KeyError(fn)


def _decoded_outputs(fn, fmt, values):
    if fn == "phase_response":
        return {
            "out_real": fmt.decode(values["out_real"]),
            "out_imag": fmt.decode(values["out_imag"]),
        }
    if fn == "frozen_roe_face":
        return {"left": fmt.decode(values["left"]), "right": fmt.decode(values["right"])}
    return {"out": fmt.decode(values["out"])}


def _register_names(fn, fmt):
    return [r.name for r in _compile(fn, fmt).operation.module.registers]


class Sweep:
    """单 case 的逐分支累积评估：误差、旗标一致性、叠加均匀性、方法误差。"""

    def __init__(self, fn, fmt, tolerance_quanta):
        self.fn, self.fmt = fn, fmt
        self.quantum = 1.0 / (1 << fmt.fraction)
        self.tolerance = tolerance_quanta * self.quantum
        self.max_error = 0.0
        self.per_output = {}
        self.method_error = 0.0
        self.exact = 0
        self.compared = 0
        self.flagged = 0
        self.misflagged = 0
        self.gray = 0
        self.branches = 0
        self.uniform_dev = 0.0

    def absorb(self, state, names, nbits, oracle):
        if len(state) != 1 << nbits:
            raise AssertionError(f"{self.fn}: 分支数 {len(state)} ≠ 2^{nbits}")
        self.branches += len(state)
        uniform = 2.0 ** (-nbits / 2)
        for basis, amplitude in state.items():
            self.uniform_dev = max(self.uniform_dev, abs(abs(amplitude) - uniform))
            values = dict(zip(names, basis, strict=True))
            inputs = _decoded_inputs(self.fn, self.fmt, values)
            expected, reference_out, exp_status = oracle(inputs, self.fmt)
            status = values["status"]
            if exp_status is None:
                self.gray += 1
                continue
            if status != exp_status:
                self.misflagged += 1
                continue
            if status:
                self.flagged += 1
                continue
            if expected is None:
                continue
            actual = _decoded_outputs(self.fn, self.fmt, values)
            for key, want in expected.items():
                err = abs(actual[key] - want)
                self.per_output[key] = max(self.per_output.get(key, 0.0), err)
                self.max_error = max(self.max_error, err)
                if reference_out is not None:
                    self.method_error = max(
                        self.method_error, abs(want - reference_out[key])
                    )
            self.compared += 1
            self.exact += int(
                all(actual[key] == want for key, want in expected.items())
            )

    @property
    def passed(self):
        return (
            self.max_error <= self.tolerance + 1e-12
            and self.misflagged == 0
            and self.uniform_dev < 1e-9
        )

    def metrics(self):
        out = {
            "max_error": round(self.max_error, 10),
            "max_error_quanta": round(self.max_error / self.quantum, 4),
            "exact_fraction": round(self.exact / max(1, self.compared), 6),
            "compared": self.compared,
            "flagged": self.flagged,
            "misflagged": self.misflagged,
            "gray": self.gray,
            "branches": self.branches,
            "uniform_dev": self.uniform_dev,
        }
        for key, value in sorted(self.per_output.items()):
            out[f"err_{key}"] = round(value, 10)
        if self.method_error:
            out["method_error"] = round(self.method_error, 10)
        return out


SIMPLE_ORACLES = {
    "polynomial": oracle_polynomial,
    "guarded_reciprocal": oracle_guarded_reciprocal,
    "pressure": oracle_pressure,
    "roe_speed": oracle_roe_speed,
    "frozen_roe_face": oracle_frozen_roe_face,
}


def run():
    report = Report(
        "mathfunc",
        "数学函数廊（5 函数）与 frozen Roe 面在 FixedFormat(6,2)/(8,3) 下的叠加穷举"
        "数值验证；rir-pysparq 主路径 + reference/adapter 对拍；实现/方法误差分离，"
        "status 定义域语义逐分支核对。",
    )
    tasks = build_tasks()
    workers = int(os.environ.get("VERIFY_WORKERS", "10"))
    results = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for name, nbits, state in pool.map(execute, tasks):
            results[name] = (nbits, state)

    # 各编译函数的工作区峰值（说明 OriginIR-ext 24 比特预算为何不适用的实测值）
    workspaces = {
        f"{fn}/{fmt_key}": _workspace(fn, FixedFormat(*FMTS[fmt_key]))
        for fn in TOLERANCE
        for fmt_key in FMTS
    }

    def sweep_case(case, fn, fmt_key, task_names, criterion, oracle=None, extra=None):
        fmt = FixedFormat(*FMTS[fmt_key])
        sweep = Sweep(fn, fmt, TOLERANCE[fn])
        names = _register_names(fn, fmt)
        for task_name in task_names:
            nbits, state = results[task_name]
            sweep.absorb(state, names, nbits, oracle or SIMPLE_ORACLES[fn])
        report.case(
            case,
            paths=["rir-pysparq"],
            parameters={
                "function": fn,
                "format": fmt_key,
                "sweeps": len(task_names),
                "workspace_qubits": workspaces[f"{fn}/{fmt_key}"],
                "originir_ext": "不适用：工作区远超 24 比特预算",
                **(extra or {}),
            },
            metrics=sweep.metrics(),
            criterion=criterion,
            passed=sweep.passed,
        )
        return sweep

    def cross_case(case, fn, fmt_key, rir_name, ref_name, adp_name, extra=None):
        _, rir = results[rir_name]
        _, ref = results[ref_name]
        _, adp = results[adp_name]
        deviation = max(
            amplitude_error(rir, ref),
            amplitude_error(rir, adp),
            amplitude_error(ref, adp),
        )
        report.case(
            case,
            paths=["rir-pysparq", "reference", "adapter-pysparq"],
            parameters={
                "function": fn,
                "format": fmt_key,
                "branches": len(rir),
                "workspace_qubits": workspaces[f"{fn}/{fmt_key}"],
                **(extra or {}),
            },
            metrics={"max_pairwise_deviation": deviation},
            criterion="同一叠加程序三路径振幅两两一致（< 1e-9）",
            passed=deviation < 1e-9,
        )

    # --- 单输入函数：全域穷举 ---
    for fmt_key, bits in (("6.2", 6), ("8.3", 8)):
        sweep_case(
            f"polynomial-exhaustive-{fmt_key}",
            "polynomial",
            fmt_key,
            [f"polynomial-full-{fmt_key}"],
            f"无旗标分支与截断链参考逐比特一致（max_error == 0）；"
            f"旗标集合与值域越界参考完全一致（misflagged == 0）；2^{bits} 分支均匀",
        )
        sweep_case(
            f"guarded-reciprocal-exhaustive-{fmt_key}",
            "guarded_reciprocal",
            fmt_key,
            [f"guarded_reciprocal-full-{fmt_key}"],
            f"全域无旗标（守护分支吸收 x=0 奇点）；max_error ≤ 1 量子；2^{bits} 分支均匀",
        )

    # --- pressure ---
    for fmt_key in FMTS:
        sweep_case(
            f"pressure-fibers-{fmt_key}",
            "pressure",
            fmt_key,
            [f"pressure-fiber-{axis}-{fmt_key}" for axis in ("rho", "momentum", "energy")],
            "三轴逐轴穷举：与截断链参考逐比特一致（max_error == 0）；rho=0 定义域旗标（位 0）"
            "与越界旗标（位 1）逐分支符合参考（misflagged == 0）",
        )
        sweep_case(
            f"pressure-joint-cube-{fmt_key}",
            "pressure",
            fmt_key,
            [f"pressure-cube-{fmt_key}"],
            "4×4×4 联合立方体：max_error == 0；misflagged == 0",
        )

    # --- roe_speed ---
    for fmt_key in FMTS:
        sweep_case(
            f"roe-speed-fibers-{fmt_key}",
            "roe_speed",
            fmt_key,
            [f"roe_speed-fiber-{axis}-{fmt_key}" for axis in ("rho_l", "momentum_l", "rho_r", "momentum_r")],
            "四轴逐轴穷举：与定点语义仿真逐比特一致（max_error == 0）；rho≤0 定义域旗标逐分支符合参考",
        )
        sweep_case(
            f"roe-speed-joint-cube-{fmt_key}",
            "roe_speed",
            fmt_key,
            [f"roe_speed-cube-{fmt_key}"],
            "2^4 联合立方体：max_error == 0；misflagged == 0",
        )

    # --- phase_response（复函数，按实/虚部分开评估）---
    phase_oracle = make_phase_oracle(_phase_coefficients(FixedFormat(6, 2)))
    for fmt_key in FMTS:
        fiber_names = (
            ["phase-fiber-re-6.2", "phase-fiber-im-6.2", "phase-fiber-sg-6.2"]
            if fmt_key == "6.2"
            else ["phase-fiber-re-8.3", "phase-fiber-im-8.3"]
        )
        sweep_case(
            f"phase-response-fibers-{fmt_key}",
            "phase_response",
            fmt_key,
            fiber_names,
            "实/虚轴穷举（6.2 含奇异线 y=1）：与核配方 _Fx 仿真逐比特一致（max_error == 0，"
            "对照 math_approximation Chebyshev 系数）；核区间越界与 z=±i 定义域旗标逐分支符合",
            oracle=phase_oracle,
            extra={"degree": PHASE_DEGREE, "singular_line": fmt_key == "6.2"},
        )
        sweep_case(
            f"phase-response-cube-{fmt_key}",
            "phase_response",
            fmt_key,
            [f"phase-cube-{fmt_key}"],
            "4×4 联合立方体：max_error == 0；misflagged == 0",
            oracle=phase_oracle,
            extra={"degree": PHASE_DEGREE},
        )

    # 核方法误差（纯经典稠密网格；系数取自模块属性，不跑后端）
    coeffs = _phase_coefficients(FixedFormat(6, 2))
    method = {}
    for kind, fn_ref in (("exp", math.exp), ("sin", math.sin), ("cos", math.cos)):
        (lo, hi), cs = coeffs[kind]

        def cheb(x, cs=cs, lo=lo, hi=hi):
            t = (x - (lo + hi) / 2) / ((hi - lo) / 2)
            b1 = b2 = 0.0
            for a in reversed(cs[1:]):
                b1, b2 = 2 * t * b1 - b2 + a, b1
            return t * b1 - b2 + cs[0]

        method[f"method_error_{kind}"] = max(
            abs(cheb(lo + (hi - lo) * i / 4000) - fn_ref(lo + (hi - lo) * i / 4000))
            for i in range(4001)
        )
    report.case(
        "phase-response-kernel-method",
        paths=["classic-dense-grid"],
        parameters={
            "degree": PHASE_DEGREE,
            "intervals": {k: list(v[0]) for k, v in coeffs.items()},
            "note": "方法误差为 Chebyshev 阶数固有近似误差，信息性指标，不计入实现判据",
        },
        metrics={k: round(v, 8) for k, v in method.items()},
        criterion="信息性：报告 degree=3 系数在 4001 点稠密网格上与真函数的最大偏差",
        passed=True,
    )

    # --- frozen_roe_face ---
    for fmt_key in FMTS:
        sweep_case(
            f"roe-face-index-joint-{fmt_key}",
            "frozen_roe_face",
            fmt_key,
            [f"roe-face-index-{fmt_key}"],
            "row×col 全 16 组（含越界索引 3 → 0.0）：max_error == 0；misflagged == 0",
            extra={"physical_state": "(1.0,0.5,2.5)/(1.5,-0.25,3.0)", "outputs": ["left", "right"]},
        )
        sweep_case(
            f"roe-face-fibers-{fmt_key}",
            "frozen_roe_face",
            fmt_key,
            [f"roe-face-fiber-{axis}-{fmt_key}" for axis in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r")],
            "六实轴全幅值 16 点网格：max_error == 0；rho≤0 / c2≤0 定义域旗标逐分支符合参考",
            extra={"grid": "全幅值均匀网格（6.2 步长 1.0，8.3 步长 2.0）"},
        )

    # --- 跨后端三方对拍 ---
    cross_case(
        "cross-backend-polynomial-6.2", "polynomial", "6.2",
        "polynomial-full-6.2", "polynomial-cross-ref", "polynomial-cross-adp",
    )
    cross_case(
        "cross-backend-guarded-reciprocal-6.2", "guarded_reciprocal", "6.2",
        "guarded_reciprocal-full-6.2", "guarded_reciprocal-cross-ref", "guarded_reciprocal-cross-adp",
    )
    cross_case(
        "cross-backend-pressure-cube-6.2", "pressure", "6.2",
        "pressure-cube-6.2", "pressure-cross-ref", "pressure-cross-adp",
    )
    cross_case(
        "cross-backend-roe-speed-cube-6.2", "roe_speed", "6.2",
        "roe_speed-cube-6.2", "roe_speed-cross-ref", "roe_speed-cross-adp",
    )
    cross_case(
        "cross-backend-phase-cube-6.2", "phase_response", "6.2",
        "phase-cube-6.2", "phase-cross-ref", "phase-cross-adp",
    )
    cross_case(
        "cross-backend-roe-face-index-6.2", "frozen_roe_face", "6.2",
        "roe-face-cross-rir", "roe-face-cross-ref", "roe-face-cross-adp",
        extra={"sweeps": "row 2 bit × col 1 bit（8 分支缩小版）"},
    )

    # --- 基态确定性抽点 ---
    basis_metrics = {}
    basis_failures = 0
    phase_oracle_b = phase_oracle
    for fn in TOLERANCE:
        nbits, state = results[f"basis-{fn}"]
        if len(state) != 1:
            basis_failures += 1
            basis_metrics[f"{fn}_error"] = math.inf
            continue
        basis, amplitude = next(iter(state.items()))
        fmt = FixedFormat(6, 2)
        names = _register_names(fn, fmt)
        values = dict(zip(names, basis, strict=True))
        inputs = _decoded_inputs(fn, fmt, values)
        oracle = phase_oracle_b if fn == "phase_response" else SIMPLE_ORACLES[fn]
        expected, _, exp_status = oracle(inputs, fmt)
        err = math.inf
        if (
            abs(amplitude - 1) < 1e-12
            and expected is not None
            and values["status"] == exp_status
        ):
            actual = _decoded_outputs(fn, fmt, values)
            err = max(abs(actual[k] - expected[k]) for k in expected)
        if not (
            abs(amplitude - 1) < 1e-12
            and err <= TOLERANCE[fn] * (1.0 / (1 << fmt.fraction)) + 1e-12
        ):
            basis_failures += 1
        basis_metrics[f"{fn}_error"] = round(err, 10) if math.isfinite(err) else "n/a"
    basis_metrics["failures"] = basis_failures
    report.case(
        "basis-determinism-6.2",
        paths=["rir-pysparq"],
        parameters={"points": {fn: "物理代表点" for fn in TOLERANCE}, "format": "6.2"},
        metrics=basis_metrics,
        criterion="单基态输入 → 单基态输出（|振幅|=1，输入保持），输出与参考误差在容差内，旗标一致",
        passed=basis_failures == 0,
    )

    report.write()
    return report


if __name__ == "__main__":
    run()
