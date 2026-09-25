"""Publication-grade numerical validation of the math-function gallery and the frozen Roe face.

Covers the 5 functions of examples/math_functions.py (pressure / roe_speed /
phase_response / guarded_reciprocal / polynomial) and the frozen_roe_face of
applications/roe_formulas.py; all are compiled via oracq.compile_function into
fixed-point reversible modules and executed on real backends:

- rir-pysparq (the PySparQ native RIR interpreter) is the primary path; a
  superposition run exhausts the input domain at once;
- reference / adapter-pysparq provide amplitude-level three-way cross-checks
  on representative programs;
- the OriginIR-ext state-vector path is not applicable because the workspace
  of compiled functions far exceeds the 24-qubit budget (workspace_table
  measures 10^2-10^3 bits; each case's parameters record the measured value).

Oracle independence: expected values come from formulas independently
transcribed in this script from their mathematical definitions using float64
(math/cmath), without importing the implementation under test;
phase_response additionally treats the Chebyshev coefficients in the module
attribute math_approximation as the "implementation reference" and the true
function as the "method reference", reporting implementation error and method
error separately. The status flag semantics (bit 0 = domain failure, bit 1 =
range/word-length overflow) are checked branch by branch.

Run: PYTHONPATH=src <python with pysparq+uniqc> tests/verification/verify_mathfunc.py
The environment variable VERIFY_WORKERS controls the number of parallel
processes (default 10; the heaviest gate-level simulation of a single branch
in this group takes about 2.5 s, and multiprocessing keeps the total wall
clock under 5 minutes). All executions are superposition-exhaustive; no
per-basis-state loops are used.
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
    sys.path.insert(0, str(ROOT))  # lets worker processes import the functions under test from examples/

from oracq import Builder, FixedFormat
from oracq.infrastructure.mathfunc import Index, MathConfig, compile_function

FMTS = {"6.2": (6, 2), "8.3": (8, 3)}

# Chebyshev degree for phase_response: the kernel gate count of the default
# degree=6 is too heavy for gate-level simulation, so degree=3 is used to keep
# the budget; the method error (polynomial vs true function) is reported
# separately and never mixed into the implementation error.
PHASE_DEGREE = 3

# Quantization tolerance per function (unit: quanta), relaxed according to the
# number of truncating operators (mul/div/sqrt/kernel evaluation) and propagation.
TOLERANCE = {
    "polynomial": 0.0,
    "guarded_reciprocal": 1.0,
    "pressure": 0.0,
    "roe_speed": 0.0,
    "phase_response": 0.0,
    "frozen_roe_face": 0.0,
}

# ---------------------------------------------------------------------------
# Compilation specs of the functions under test (compiled on demand inside worker
# processes; compilation is deterministic and cheap)
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
        from oracq.applications.roe_formulas import frozen_roe_face

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
    """Workspace peak of the compiled program (max over locals stacked along the call chain), used for the OriginIR budget note."""
    from oracq.infrastructure.layout import workspace_table

    program = _compile(fn_key, fmt).program()
    return workspace_table(program)[program.entry]


# ---------------------------------------------------------------------------
# Independent float64 oracle (independently transcribed from the mathematical definitions)
# Returns (expected output|None, reference output|None, expected status|None):
#   expected output None -- that branch skips the numeric comparison (flagged or gray zone);
#   expected status: bit 0 = domain failure, bit 1 = range/word-length overflow; None = threshold gray zone, not compared;
#   reference output used only by phase_response (true function value, for the method error).
# ---------------------------------------------------------------------------


def _limit(fmt):
    return 2.0 ** (fmt.width - fmt.fraction - 1)


def _quantize(fmt, value):
    """Fixed-point encode-decode round trip of a generation-time constant (encode truncates toward zero)."""
    return fmt.decode(fmt.encode(value))


class _Fx:
    """Bit-by-bit simulation of fixed_arithmetic on the raw integer domain:
    mul/div/sqrt truncate magnitudes toward zero, add/sub wrap modulo the word,
    status bit 0 = domain, bit 1 = overflow; select propagates only the chosen
    branch's flag. Values travel along paths as (raw, flag) tuples."""

    def __init__(self, fmt):
        self.w, self.f = fmt.width, fmt.fraction
        self.sign = 1 << (self.w - 1)
        self.mask = (1 << self.w) - 1

    def wrap(self, raw):
        raw &= self.mask
        return raw - (1 << self.w) if raw & self.sign else raw

    def const(self, x):
        # Consistent with fmt.encode: int() truncates toward zero then wraps; constants here all fit the word length
        return self.wrap(int(x * (1 << self.f))), 0

    def input(self, x):
        # Grid values (decoded inputs) back to raw; the float round trip is exact
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
        # Boolean comparison (the circuit computes lt and inverts; the comparison itself carries no flag)
        return (1 if a[0] >= b[0] else 0), a[1] | b[1]

    def select(self, test, yes, no):
        # Both branches are computed, but only the chosen branch's (value, flag) propagates
        chosen = yes if test[0] else no
        return chosen[0], test[1] | chosen[1]

    def decode(self, a):
        return a[0] / (1 << self.f)


def oracle_polynomial(v, fmt):
    # Exact simulation of the truncation chain: the x*x magnitude truncates toward zero onto the quantum grid, addition has no truncation; overflow judged by word length.
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
        return {"out": 0.0}, None, 0  # the guard branch absorbs the singularity; the unselected branch's division by zero is masked
    val = 1.0 / x
    flag = 2 if abs(val) >= _limit(fmt) else 0
    return ({"out": val} if flag == 0 else None), None, flag


def oracle_pressure(v, fmt):
    # Exact simulation of the truncation chain (same method as polynomial): div/mul truncate the magnitude by one quantum,
    # add/sub do not truncate; gamma-1 is front-end constant-folded to 0.4 before encoding (6.2 -> 0.25, 8.3 -> 0.375).
    rho, m, e = v["rho"], v["momentum"], v["energy"]
    if rho == 0:
        return None, None, 1  # division by zero -> bit 0
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
    # _Fx bit-by-bit simulation: sqrt is isqrt(raw<<f), mul/div truncate toward zero, flag semantics as in the circuit.
    F = _Fx(fmt)
    rl, rr = v["rho_l"], v["rho_r"]
    if rl <= 0 or rr <= 0:
        return None, None, 1  # rho<0 sqrt domain / rho=0 division by zero (the circuit also sets bit 0)
    left = F.sqrt(F.input(rl))
    right = F.sqrt(F.input(rr))
    b = F.div(F.mul(left, F.input(v["momentum_l"])), F.input(rl))
    d = F.div(F.mul(right, F.input(v["momentum_r"])), F.input(rr))
    val = F.div(F.add(b, d), F.add(left, right))
    status = val[1]
    # Method reference: the float64 original formula (with exact sqrt values)
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
    """Main body of the Roe face formula; consts = (g1, gm, gm3, tmg, delta) are the
    values of (gamma-1, gamma, gamma-3, 3-gamma, entropy_delta)."""
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
    """_Fx bit-by-bit simulation of frozen_roe_face, transcribed in the source order of roe_formulas.py.

    Two key points of the flag semantics (against lowering.py / numeric.py):
    - a helper call merges the status of **all arguments** into the module's
      internal status (combined), so flags of unselected argument branches of
      helpers like pick3/euler_entry are not masked;
    - helper formal parameters (gamma, delta) are dynamic inputs inside the
      helper module and do not participate in front-end constant folding
      (e.g. 3-gamma = 3-1.25 = 1.75 inside euler_entry, not the folded 1.5).
    Returns (left, right, flag)."""
    F = _Fx(fmt)
    rl, ml, el = (F.input(v[k]) for k in ("rho_l", "m_l", "e_l"))
    rr, mr, er = (F.input(v[k]) for k in ("rho_r", "m_r", "e_r"))
    row, col = v["row"], v["col"]
    half = F.const(0.5)
    gm = F.const(1.4)  # gamma materialized as a constant register for the helper argument
    g1 = F.const(0.4)  # gamma-1 folded by the front end inside the top-level function body, then encoded

    def primitive(rho, m, e):
        gm1 = F.sub(gm, F.const(1.0))  # gamma-1 computed dynamically inside the helper
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
        # `delta <= 0` is equivalent to `0 >= delta`; under 6.2, delta raw = 0 -> the value branch is always taken
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
        return None, None, 1  # rho<0 sqrt domain / rho=0 division by zero
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
        exact = None  # the exact-constant formula is singular at this point (e.g. c2<=0); the method-error point is skipped
    return outputs, exact, 0


def _kernel_fx(F, coeffs, kind, a):
    """_Fx simulation of an elementary kernel module: interval constants are encode-truncated, and the Clenshaw recursion
    reproduces the multiply-add order of numeric.py bit by bit; leaving the interval sets bit 1 (against the coefficients
    of the module attribute math_approximation)."""
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
    """_Fx bit-by-bit simulation of phase_response (three elementary kernels included)."""
    F = _Fx(fmt)
    x = F.input(v["x"])
    y = F.input(v["y"])
    # w = e^{i z}: 1j*z = (-y, x) (multiplication by 0/1 is exact)
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
    """phase_response oracle: _Fx simulation (implementation reference, against the coefficient recipe) + the true function (method reference)."""

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
        raise AssertionError(f"phase_response kernel coefficients missing: {sorted(coeffs)}")
    return coeffs


# ---------------------------------------------------------------------------
# Execution workers (top-level functions, picklable; each task = one backend execution)
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
        # The reference/adapter static expansion estimate (controlled regions weighted by bit width) exceeds the actual
        # event count; the static estimates for phase/roe_face are about 1.2e6/2.5e6, above the default 1e6 budget, so
        # it is explicitly relaxed here.
        kwargs["max_steps"] = task["max_steps"]
    state = runner(program, **kwargs)
    return task["task"], nbits, {k: complex(v) for k, v in state.items()}


# ---------------------------------------------------------------------------
# Task list
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

    # polynomial / guarded_reciprocal: full-domain sweep on the single input (64 / 256 branches) + three-way cross-check
    for fn in ("polynomial", "guarded_reciprocal"):
        add(f"{fn}-full-6.2", fn, "6.2", {"x": list(range(6))}, {})
        add(f"{fn}-full-8.3", fn, "8.3", {"x": list(range(8))}, {})
        add(f"{fn}-cross-ref", fn, "6.2", {"x": list(range(6))}, {}, "reference")
        add(f"{fn}-cross-adp", fn, "6.2", {"x": list(range(6))}, {}, "adapter")

    # pressure: three-axis fiber sweeps + low-bit joint cube (4x4x4)
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

    # roe_speed: four-axis fiber sweeps + joint cube (2^4)
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

    # phase_response: real/imaginary-axis fiber sweeps + singular line y=1 (6.2) + joint cube (4x4)
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

    # frozen_roe_face: Index full joint (row x col, 16 groups) + six real-axis full-magnitude grid fibers
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
    # Reduced 8-branch version for the cross-backend check (reference/adapter are slower)
    add("roe-face-cross-rir", "frozen_roe_face", "6.2", {"row": [0, 1], "col": [0]}, idx_fix62)
    add("roe-face-cross-ref", "frozen_roe_face", "6.2", {"row": [0, 1], "col": [0]}, idx_fix62, "reference", 8_000_000)
    add("roe-face-cross-adp", "frozen_roe_face", "6.2", {"row": [0, 1], "col": [0]}, idx_fix62, "adapter", 8_000_000)
    for axis in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r"):
        f62 = {k: v for k, v in roe_face_fix62.items() if k != axis}
        f83 = {k: v for k, v in roe_face_fix83.items() if k != axis}
        # Full-magnitude uniform grid: 6.2 sweeps bits[2..6) (step 1.0), 8.3 sweeps bits[4..8) (step 2.0)
        add(f"roe-face-fiber-{axis}-6.2", "frozen_roe_face", "6.2", {axis: [2, 3, 4, 5]}, f62)
        add(f"roe-face-fiber-{axis}-8.3", "frozen_roe_face", "8.3", {axis: [4, 5, 6, 7]}, f83)

    # Basis-state determinism spot checks (single basis state -> single basis state, inputs preserved, outputs match the reference)
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
# Result evaluation
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
    """Per-branch accumulating evaluation for a single case: errors, flag consistency, superposition uniformity, method error."""

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
            raise AssertionError(f"{self.fn}: branch count {len(state)} != 2^{nbits}")
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
        "Superposition-exhaustive numerical validation of the math-function gallery (5 functions) and the frozen Roe face "
        "under FixedFormat(6,2)/(8,3); rir-pysparq primary path + reference/adapter cross-checks; implementation/method "
        "errors separated, status domain semantics checked branch by branch.",
    )
    tasks = build_tasks()
    workers = int(os.environ.get("VERIFY_WORKERS", "10"))
    results = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for name, nbits, state in pool.map(execute, tasks):
            results[name] = (nbits, state)

    # Workspace peak of each compiled function (measured values explaining why the 24-bit OriginIR-ext budget does not apply)
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
                "originir_ext": "not applicable: workspace far exceeds the 24-bit budget",
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
            criterion="the same superposition program agrees pairwise across the three paths (< 1e-9)",
            passed=deviation < 1e-9,
        )

    # --- Single-input functions: full-domain sweep ---
    for fmt_key, bits in (("6.2", 6), ("8.3", 8)):
        sweep_case(
            f"polynomial-exhaustive-{fmt_key}",
            "polynomial",
            fmt_key,
            [f"polynomial-full-{fmt_key}"],
            f"unflagged branches agree bit by bit with the truncation-chain reference (max_error == 0); "
            f"the flag set fully matches the range-overflow reference (misflagged == 0); 2^{bits} uniform branches",
        )
        sweep_case(
            f"guarded-reciprocal-exhaustive-{fmt_key}",
            "guarded_reciprocal",
            fmt_key,
            [f"guarded_reciprocal-full-{fmt_key}"],
            f"no flags over the whole domain (the guard branch absorbs the x=0 singularity); max_error <= 1 quantum; 2^{bits} uniform branches",
        )

    # --- pressure ---
    for fmt_key in FMTS:
        sweep_case(
            f"pressure-fibers-{fmt_key}",
            "pressure",
            fmt_key,
            [f"pressure-fiber-{axis}-{fmt_key}" for axis in ("rho", "momentum", "energy")],
            "three per-axis sweeps: agrees bit by bit with the truncation-chain reference (max_error == 0); the rho=0 domain flag (bit 0)"
            " and the overflow flag (bit 1) match the reference branch by branch (misflagged == 0)",
        )
        sweep_case(
            f"pressure-joint-cube-{fmt_key}",
            "pressure",
            fmt_key,
            [f"pressure-cube-{fmt_key}"],
            "4x4x4 joint cube: max_error == 0; misflagged == 0",
        )

    # --- roe_speed ---
    for fmt_key in FMTS:
        sweep_case(
            f"roe-speed-fibers-{fmt_key}",
            "roe_speed",
            fmt_key,
            [f"roe_speed-fiber-{axis}-{fmt_key}" for axis in ("rho_l", "momentum_l", "rho_r", "momentum_r")],
            "four per-axis sweeps: agrees bit by bit with the fixed-point semantics simulation (max_error == 0); rho<=0 domain flags match the reference branch by branch",
        )
        sweep_case(
            f"roe-speed-joint-cube-{fmt_key}",
            "roe_speed",
            fmt_key,
            [f"roe_speed-cube-{fmt_key}"],
            "2^4 joint cube: max_error == 0; misflagged == 0",
        )

    # --- phase_response (complex function, real/imaginary parts evaluated separately) ---
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
            "real/imaginary-axis sweeps (6.2 includes the singular line y=1): agrees bit by bit with the kernel-recipe _Fx simulation "
            "(max_error == 0, against the math_approximation Chebyshev coefficients); kernel-interval overflow and z=+/-i domain flags match branch by branch",
            oracle=phase_oracle,
            extra={"degree": PHASE_DEGREE, "singular_line": fmt_key == "6.2"},
        )
        sweep_case(
            f"phase-response-cube-{fmt_key}",
            "phase_response",
            fmt_key,
            [f"phase-cube-{fmt_key}"],
            "4x4 joint cube: max_error == 0; misflagged == 0",
            oracle=phase_oracle,
            extra={"degree": PHASE_DEGREE},
        )

    # Kernel method error (purely classical dense grid; coefficients taken from the module attribute, no backend run)
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
            "note": "the method error is the inherent approximation error of the Chebyshev degree; an informative metric excluded from the implementation criterion",
        },
        metrics={k: round(v, 8) for k, v in method.items()},
        criterion="informative: reports the maximum deviation of the degree=3 coefficients from the true function on a 4001-point dense grid",
        passed=True,
    )

    # --- frozen_roe_face ---
    for fmt_key in FMTS:
        sweep_case(
            f"roe-face-index-joint-{fmt_key}",
            "frozen_roe_face",
            fmt_key,
            [f"roe-face-index-{fmt_key}"],
            "all 16 row x col groups (out-of-range index 3 -> 0.0 included): max_error == 0; misflagged == 0",
            extra={"physical_state": "(1.0,0.5,2.5)/(1.5,-0.25,3.0)", "outputs": ["left", "right"]},
        )
        sweep_case(
            f"roe-face-fibers-{fmt_key}",
            "frozen_roe_face",
            fmt_key,
            [f"roe-face-fiber-{axis}-{fmt_key}" for axis in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r")],
            "six real-axis full-magnitude 16-point grids: max_error == 0; rho<=0 / c2<=0 domain flags match the reference branch by branch",
            extra={"grid": "full-magnitude uniform grid (6.2 step 1.0, 8.3 step 2.0)"},
        )

    # --- Three-way cross-backend check ---
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
        extra={"sweeps": "row 2 bit x col 1 bit (reduced 8-branch version)"},
    )

    # --- Basis-state determinism spot checks ---
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
        parameters={"points": {fn: "physical representative point" for fn in TOLERANCE}, "format": "6.2"},
        metrics=basis_metrics,
        criterion="single basis-state input -> single basis-state output (|amplitude| = 1, inputs preserved), output within tolerance of the reference, flags consistent",
        passed=basis_failures == 0,
    )

    report.write()
    return report


if __name__ == "__main__":
    run()
