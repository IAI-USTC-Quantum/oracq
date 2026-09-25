"""Publication-grade numerical validation of QHAM and QFVM.

QHAM side (applications/qham + algorithms/qham.py + algorithms/pde.py):
- The lifted-generator BE retrieves the full matrix in one run through the
  aux register, compared element-wise against the independent classical
  matrix of reference.py (two input models: the structured stencil port and
  the spectral-embedding Pauli port);
- Generator consistency of the same problem under the three input models
  stencil / spectral / QRAM angle table, plus the pointwise truth of the QRAM
  coefficient angle encoding (diagonal_block_encoding + angle table);
- Lifted initial-state preparation vs a directly constructed tensor-word
  weight vector; the finite-Taylor QODE end-to-end chain vs the classical
  (I+tG)Y_in; matrix-level validation of the explicit dissipative shift
  G - mu*I;
- Homotopy weights, generation-row rules vs the independent chain rule
  (reference.py), and contraction convergence of the HAM partial sums
  against the exact Riccati solution (measured contraction factor vs
  |1+eta|).

QFVM side (applications/qfvm.py + roe.py + roe_formulas.py + flow_data.py):
- The compiled fixed-point roe_face circuit compared bit-level against an
  independent fixed-point simulation oracle (re-implemented in this script
  per the documented fixed_arithmetic semantics toward_zero/modular_wrap);
  the method error (vs the float64 roe_formulas) is reported separately;
- The sparse-entry oracle compared element-wise under an 8-branch
  superposition against the same fixed-point simulation matrix (west/center/
  east bands and zero structure included), plus an independent case for the
  padding diagonal padding_value;
- The location oracle validated with all-column superposition: full
  permutation per column + slot mapping vs the independent geometric
  semantics;
- RHS residual-state preparation (QRAM angle tree + sign) vs the independent
  residual computation;
- flow_data single-step update identities: F*(L,R) = left.L + right.R vs the
  matrix-inversion implementation, M.u - mass.u = residual, local update
  patches identical to a full recomputation, and the ptheta pointwise truth.

Format note: roe_face's compile_function requires the row/column Index(2)
inputs to be representable, so the minimal usable fixed point is
FixedFormat(5,2); entropy_delta is taken as 0.5 so that 2*delta, delta^2 and
delta are all exactly representable in that format. If the fraction bits
cannot represent 2*delta (e.g. FixedFormat(4,1) with the default
delta=0.125, where 2*delta=0.25 truncates to raw 0), the entropy-correction
branch divides by zero and the entry silently zeroes per the documented
totalize behavior; the existing integration tests cover only the padding
diagonal, not the non-zero Roe entries.

Run: PYTHONPATH=src <python with pysparq+uniqc> tests/verification/verify_qham_qfvm.py
"""

from __future__ import annotations

import math

from harness import (
    Report,
    adapter_pysparq,
    amplitude_error,
    amplitudes_to_statevector,
    originir_ext,
    reference,
    rir_pysparq,
    statevector_error,
    superposition_program,
)

from oracq import Binding, Bits, Builder, FixedFormat, bind, unresolved
from oracq.algorithms.input_model.block_encoding import matrix_pauli_encoding
from oracq.algorithms.input_model.operators import scale
from oracq.algorithms.input_model.oracles import gate_state_prep, invoke, qram_database
from oracq.algorithms.input_model.qham import gate_bindings, qham_input_model, taylor_qode
from oracq.algorithms.qpde.pde import DiscretePDE, PDEInput, make_qpde, qpde_solver
from oracq.applications.flow_data import RoeFlowData
from oracq.applications.qfvm import (
    bind_qfvm,
    geometry_cells,
    ptheta_cells,
    qfvm_memories,
    qfvm_sparse_access,
    roe_qfvm_inputs,
)
from oracq.applications.qham import (
    Discretization,
    Field,
    Grid,
    HomotopyWeight,
    Known,
    PolynomialPDE,
    QHAMPlan,
    qram_coefficient_encoding,
    qram_coefficient_memory,
    structured_fd_bindings,
)
from oracq.applications.qham.linearization import compositions
from oracq.applications.roe import roe_face
from oracq.applications.roe_formulas import frozen_roe_face
from oracq.infrastructure.layout import workspace_table

# ---------------------------------------------------------------------------
# Independent fixed-point simulation oracle: re-implemented per the documented
# fixed_arithmetic semantics (toward_zero magnitude truncation, modular wrap,
# status[0] = domain, status[1] = word-length overflow)
# ---------------------------------------------------------------------------


class Fx:
    """Fixed-point word simulation: raw is an n-bit unsigned word, decoded as two's complement."""

    def __init__(self, fmt, raw, flags):
        self.fmt = fmt
        self.raw = raw & ((1 << fmt.width) - 1)
        self.flags = flags

    @classmethod
    def const(cls, fmt, value, flags):
        return cls(fmt, int(value * (1 << fmt.fraction)) & ((1 << fmt.width) - 1), flags)

    def decoded(self):
        value = self.raw
        if self.fmt.signed and value >> (self.fmt.width - 1):
            value -= 1 << self.fmt.width
        return value / (1 << self.fmt.fraction)

    def _mag_sign(self):
        n = self.fmt.width
        sign = (self.raw >> (n - 1)) & 1 if self.fmt.signed else 0
        mag = (-self.raw) & ((1 << n) - 1) if sign else self.raw
        return mag, sign

    def _from_mag_sign(self, mag, sign):
        raw = (-mag) & ((1 << self.fmt.width) - 1) if sign else mag & ((1 << self.fmt.width) - 1)
        return Fx(self.fmt, raw, self.flags)

    def add(self, other):
        n = self.fmt.width
        out = (self.raw + other.raw) & ((1 << n) - 1)
        if self.fmt.signed:
            sa, sb, so = self.raw >> (n - 1), other.raw >> (n - 1), out >> (n - 1)
            if (not (sa ^ sb)) and (sa ^ so):
                self.flags[1] = 1
        return Fx(self.fmt, out, self.flags)

    def sub(self, other):
        return self.add(other.neg())

    def neg(self):
        n = self.fmt.width
        if self.fmt.signed and self.raw == (1 << (n - 1)):
            self.flags[1] = 1
        return Fx(self.fmt, (-self.raw) & ((1 << n) - 1), self.flags)

    def abs(self):
        mag, sign = self._mag_sign()
        if self.fmt.signed and sign and mag == (1 << (self.fmt.width - 1)):
            self.flags[1] = 1
        return Fx(self.fmt, mag, self.flags)

    def mul(self, other):
        n, f = self.fmt.width, self.fmt.fraction
        ma, sa = self._mag_sign()
        mb, sb = other._mag_sign()
        mag = (ma * mb) >> f
        if mag >> (n - int(self.fmt.signed)):
            self.flags[1] = 1
        return self._from_mag_sign(mag, sa ^ sb)

    def div(self, other):
        n, f = self.fmt.width, self.fmt.fraction
        ma, sa = self._mag_sign()
        mb, sb = other._mag_sign()
        if mb == 0:
            self.flags[0] = 1
            return Fx(self.fmt, 0, self.flags)
        mag = (ma << f) // mb
        if mag >> (n - int(self.fmt.signed)):
            self.flags[1] = 1
        return self._from_mag_sign(mag, sa ^ sb)

    def sqrt(self):
        n, f = self.fmt.width, self.fmt.fraction
        _, sign = self._mag_sign()
        if sign:
            self.flags[0] = 1
            return Fx(self.fmt, 0, self.flags)
        mag = math.isqrt(self.raw << f)
        if mag >> (n - int(self.fmt.signed)):
            self.flags[1] = 1
        return Fx(self.fmt, mag, self.flags)


def _pick3(index, first, second, third):
    return (first, second, third, Fx.const(first.fmt, 0.0, first.flags))[index]


def roe_face_fx(fmt, rho_l, m_l, e_l, rho_r, m_r, e_r, row, col, *, gamma, entropy_delta):
    """Simulate the fixed-point pipeline in the evaluation order of the frozen_roe_face source; inputs are actual numeric values.

    Constant handling matches the compiler: purely constant subexpressions
    are evaluated in float first, then truncated into the fixed-point domain
    via fmt.encode. Returns (left, right, (invalid, overflow)).
    """
    flags = [0, 0]

    def C(value):
        return Fx.const(fmt, value, flags)

    def S(value):
        return Fx.const(fmt, value, flags)

    def conserved(rho, momentum, energy):
        velocity = S(momentum).div(S(rho))
        pressure = C(gamma - 1).mul(S(energy).sub(C(0.5).mul(S(momentum)).mul(velocity)))
        enthalpy = S(energy).add(pressure).div(S(rho))
        return velocity, pressure, enthalpy

    def euler_entry(velocity, enthalpy, r, c):
        square = velocity.mul(velocity)
        first = _pick3(c, C(0.0), C(1.0), C(0.0))
        second = _pick3(
            c, C(0.5 * (gamma - 3)).mul(square), C(3 - gamma).mul(velocity), C(gamma - 1)
        )
        third = _pick3(
            c,
            velocity.mul(C(0.5 * (gamma - 1)).mul(square).sub(enthalpy)),
            enthalpy.sub(C(gamma - 1).mul(square)),
            C(gamma).mul(velocity),
        )
        return _pick3(r, first, second, third)

    def entropy_abs(eigenvalue):
        value = eigenvalue.abs()
        branch = eigenvalue.mul(eigenvalue).add(C(entropy_delta * entropy_delta)).div(
            C(2 * entropy_delta)
        )
        return value if value.decoded() >= entropy_delta else branch

    ul, _, hl = conserved(rho_l, m_l, e_l)
    ur, _, hr = conserved(rho_r, m_r, e_r)
    wl, wr = S(rho_l).sqrt(), S(rho_r).sqrt()
    total = wl.add(wr)
    u = wl.mul(ul).add(wr.mul(ur)).div(total)
    h = wl.mul(hl).add(wr.mul(hr)).div(total)
    u2 = u.mul(u)
    c2 = C(gamma - 1).mul(h.sub(C(0.5).mul(u2)))
    c = c2.sqrt()
    beta = C(gamma - 1).div(c2)
    inverse_c = C(0.5).div(c)
    bu, bu2 = beta.mul(u), beta.mul(u2)
    r0 = _pick3(row, C(1.0), u.sub(c), h.sub(u.mul(c)))
    r1 = _pick3(row, C(1.0), u, C(0.5).mul(u2))
    r2 = _pick3(row, C(1.0), u.add(c), h.add(u.mul(c)))
    inv0 = _pick3(
        col,
        C(0.25).mul(bu2).add(u.mul(inverse_c)),
        C(-0.5).mul(bu).sub(inverse_c),
        C(0.5).mul(beta),
    )
    inv1 = _pick3(col, C(1.0).sub(C(0.5).mul(bu2)), bu, beta.neg())
    inv2 = _pick3(
        col,
        C(0.25).mul(bu2).sub(u.mul(inverse_c)),
        C(-0.5).mul(bu).add(inverse_c),
        C(0.5).mul(beta),
    )
    a0, a1, a2 = entropy_abs(u.sub(c)), entropy_abs(u), entropy_abs(u.add(c))
    absolute = r0.mul(a0).mul(inv0).add(r1.mul(a1).mul(inv1)).add(r2.mul(a2).mul(inv2))
    left = C(0.5).mul(euler_entry(ul, hl, row, col).add(absolute))
    right = C(0.5).mul(euler_entry(ur, hr, row, col).sub(absolute))
    return left, right, tuple(flags)


def qfvm_entry_fx(fmt, cells, source, rowvar, colvar, band, *, mass, dx, gamma, entropy_delta):
    """Simulate the band selection and mass term of roe_entry; cells[i] is an (rho, m, e) tuple of actual values."""
    left0, right0, flags0 = roe_face_fx(
        fmt, *cells[(source - 1) % len(cells)], *cells[source % len(cells)],
        rowvar, colvar, gamma=gamma, entropy_delta=entropy_delta,
    )
    left1, right1, flags1 = roe_face_fx(
        fmt, *cells[source % len(cells)], *cells[(source + 1) % len(cells)],
        rowvar, colvar, gamma=gamma, entropy_delta=entropy_delta,
    )
    flags = [flags0[0] | flags1[0], flags0[1] | flags1[1]]
    dx_word = Fx.const(fmt, dx, flags)
    west = left0.neg().div(dx_word)
    center = left1.sub(right0).div(dx_word)
    if rowvar == colvar:
        center = center.add(Fx.const(fmt, mass, flags))
    east = right1.div(dx_word)
    value = _pick3(band, west, center, east)
    if flags != [0, 0]:
        return 0, tuple(flags)
    return value.raw, (0, 0)


# Unified QFVM flow field (fmt=(5,2) exactly representable; all intermediates within [-4, 3.75])
QFVM_FMT = FixedFormat(5, 2)
QFVM_STATES = [(1.0, 0.0, 1.0), (1.0, 0.25, 1.0), (1.25, 0.0, 1.0), (1.25, -0.25, 1.0)]
QFVM_DELTA = 0.5
QFVM_MEMORY = {
    "rho": {0: 4, 1: 4, 2: 5, 3: 5},
    "momentum": {1: 1, 3: 31},
    "energy": {i: 4 for i in range(4)},
}


def _qfvm_matrices(fmt, cells, *, mass, dx, gamma, entropy_delta):
    """Independently rebuild D_float (float64 Roe formulas) and D_emu (fixed-point simulation, raw words) per the geometric ABI.

    D = [[0, M], [M^T, 0]] acts on (half, cell*4+var) coordinates; the rows
    and columns of the padding variable var=3 receive padding_value only on
    the exact diagonal. Returns (D_float, D_emu_raw).
    """
    n = len(cells)
    size = 2 * 4 * n

    def slot_params(column, band, other_var):
        cell, var, half = (column >> 2) % n, column % 4, column >> 4
        neighbor_cell = (cell + band - 1) % n
        row = other_var + 4 * neighbor_cell + ((1 - half) << 4)
        if half == 0:
            return row, cell, var, other_var, band
        return row, neighbor_cell, other_var, var, 2 - band

    d_float = [[0.0] * size for _ in range(size)]
    d_emu = [[0] * size for _ in range(size)]
    for column in range(size):
        var = column % 4
        if var == 3:
            continue
        for slot in range(9):
            band, other_var = divmod(slot, 3)
            row, source, rowvar, colvar, sband = slot_params(column, band, other_var)
            left0 = [frozen_roe_face(*cells[(source - 1) % n], *cells[source % n],
                                     rowvar, colvar, gamma=gamma, entropy_delta=entropy_delta)]
            left1 = [frozen_roe_face(*cells[source % n], *cells[(source + 1) % n],
                                     rowvar, colvar, gamma=gamma, entropy_delta=entropy_delta)]
            west = -left0[0][0] / dx
            center = (left1[0][0] - left0[0][1]) / dx + mass * (rowvar == colvar)
            east = left1[0][1] / dx
            d_float[row][column] = (west, center, east, 0.0)[sband]
            raw, _ = qfvm_entry_fx(
                fmt, cells, source, rowvar, colvar, sband,
                mass=mass, dx=dx, gamma=gamma, entropy_delta=entropy_delta,
            )
            d_emu[row][column] = raw
    for i in range(size):
        if i % 4 == 3:
            d_float[i][i] = 1.0
            d_emu[i][i] = fmt.encode(1.0)
    return d_float, d_emu


def _face_block_float(cell_l, cell_r, *, gamma, entropy_delta):
    left = [[0.0] * 3 for _ in range(3)]
    right = [[0.0] * 3 for _ in range(3)]
    for r in range(3):
        for c in range(3):
            left[r][c], right[r][c] = frozen_roe_face(
                *cell_l, *cell_r, r, c, gamma=gamma, entropy_delta=entropy_delta
            )
    return left, right


def _script_flux(cell_l, cell_r, *, gamma, entropy_delta):
    """Script-side Roe numerical flux: F* = left.U_L + right.U_R (Euler is degree-1 homogeneous, F = A.U)."""
    left, right = _face_block_float(cell_l, cell_r, gamma=gamma, entropy_delta=entropy_delta)
    return tuple(
        sum(left[i][j] * cell_l[j] for j in range(3))
        + sum(right[i][j] * cell_r[j] for j in range(3))
        for i in range(3)
    )


# ---------------------------------------------------------------------------
# QHAM shared constructions
# ---------------------------------------------------------------------------


def _quadratic_problem():
    u = Field("u")
    pde = PolynomialPDE.from_equations({"u": -0.2 * u + 0.1 * u * u})
    plan = QHAMPlan(pde, 2)
    disc = Discretization(pde, Grid(("x",), (2,), (1.0,)))
    return plan, disc


def _aux_matrix(operation, width, signal_qubits, alpha, raw_dimension, matrix, *, name):
    """aux-register superposition trick: retrieve the full BE matrix (signal==0 block) in one run.

    Under the |c>_aux |c>_target superposition, the amplitude of (c, r, 0) is
    exactly M[r][c]/alpha/sqrt(2^w). Returns (max element deviation, missing
    support count, out-of-range support count).
    """
    b = Builder(
        name,
        {"aux": Bits(width), "target": Bits(width), "signal": Bits(signal_qubits)},
    )
    b.h(b["aux"])
    b.xor(b["aux"], b["target"])
    b.call(operation, target=b["target"], signal=b["signal"])
    state = rir_pysparq(b.finish().program())
    uniform = 1 / math.sqrt(1 << width)
    expected_support = {
        (r, c) for r in range(raw_dimension) for c in range(raw_dimension) if matrix[r][c]
    }
    seen = set()
    worst = 0.0
    for (aux, target, signal), amplitude in state.items():
        if signal != 0:
            continue
        seen.add((target, aux))
        expected = (
            matrix[target][aux] / alpha * uniform
            if target < raw_dimension and aux < raw_dimension
            else 0j
        )
        worst = max(worst, abs(amplitude - expected))
    missing = len(expected_support - seen)
    pad_leak = len({(r, c) for (r, c) in seen if r >= raw_dimension or c >= raw_dimension})
    return worst, missing, pad_leak, len(expected_support)


def verify_qham_generator_stencil(report, context):
    plan, disc = _quadratic_problem()
    matrix = disc.matrix(plan, -0.4)
    raw = plan.raw_dimension(disc.dimension)
    model = qham_input_model(plan, structured_fd_bindings(disc, [0.2, 0]), eta=-0.4)
    be = model.generator
    worst, missing, pad, support = _aux_matrix(
        be.operation, be.width, be.signal_qubits, be.alpha, raw, matrix, name="aux_stencil"
    )
    context["stencil_matrix"] = matrix
    context["stencil_model"] = model
    report.case(
        "qham-generator-matrix-stencil",
        paths=["rir-pysparq"],
        parameters={"raw_dimension": raw, "width": be.width, "alpha": be.alpha},
        metrics={"max_error": worst, "missing_support": missing, "pad_leakage": pad},
        criterion="full matrix element-wise < 1e-9 (against the independent reference matrix), support exact and padding block zero",
        passed=worst < 1e-9 and missing == 0 and pad == 0,
    )


def verify_qham_generator_spectral(report, context):
    plan, disc = _quadratic_problem()
    matrix = context["stencil_matrix"]
    raw = plan.raw_dimension(disc.dimension)
    spectral = qham_input_model(plan, gate_bindings(disc, [0.2, 0]), eta=-0.4)
    be = spectral.generator
    worst, missing, pad, _ = _aux_matrix(
        be.operation, be.width, be.signal_qubits, be.alpha, raw, matrix, name="aux_spectral"
    )
    context["spectral_model"] = spectral
    report.case(
        "qham-generator-matrix-spectral",
        paths=["rir-pysparq"],
        parameters={"raw_dimension": raw, "width": be.width, "alpha": be.alpha},
        metrics={"max_error": worst, "missing_support": missing, "pad_leakage": pad},
        criterion="the spectral-embedding port's full matrix matches the reference matrix element-wise < 1e-9 (both input models give the same matrix)",
        passed=worst < 1e-9 and missing == 0 and pad == 0,
    )


def verify_qham_generator_cross_paths(report, context):
    model = context["stencil_model"]
    be = model.generator
    program = superposition_program(be.operation, ["target"], name="cross_generator")
    total_qubits = sum(r.type.width for r in program.main.registers) + workspace_table(program)[
        program.entry
    ]
    ref = reference(program)
    rir = rir_pysparq(program)
    adapter = adapter_pysparq(program)
    origin = originir_ext(program)
    widths = [r.type.width for r in program.main.registers]
    deviation = max(
        amplitude_error(ref, rir),
        amplitude_error(ref, adapter),
        statevector_error(origin, amplitudes_to_statevector(ref, widths)),
    )
    report.case(
        "qham-generator-cross-paths",
        paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        parameters={"total_qubits": total_qubits, "budget": 24},
        metrics={"max_pairwise_deviation": deviation},
        criterion="pairwise amplitude deviation across the four real backend paths < 1e-9",
        passed=deviation < 1e-9,
    )


def verify_qham_input_models(report):
    """Consistency of the three port implementations for a linear PDE (spatially varying coefficient nu)."""
    u = Field("u")
    pde = PolynomialPDE.from_equations({"u": 0.5 * Known("nu") * u.d("x", 1) - 0.2 * u})
    plan = QHAMPlan(pde, 1)
    grid = Grid(("x",), (4,), (1.0,))
    disc = Discretization(pde, grid, {"nu": [1.0, 0.5, 1.5, 0.25]})
    initial = disc.encode_fields({"u": [0.3, 0.1, 0.2, 0.15]})
    matrix = disc.matrix(plan, -0.4)
    raw = plan.raw_dimension(disc.dimension)
    measured = {}
    alphas = {}
    for tag, bindings in (
        ("stencil", structured_fd_bindings(disc, initial)),
        ("spectral", gate_bindings(disc, initial)),
    ):
        model = qham_input_model(plan, bindings, eta=-0.4)
        be = model.generator
        worst, missing, pad, _ = _aux_matrix(
            be.operation, be.width, be.signal_qubits, be.alpha, raw, matrix,
            name="aux_model_" + tag,
        )
        measured[tag] = (worst, missing, pad)
        alphas[tag] = be.alpha
    # QRAM angle-table coefficients: the same program swaps in qram_database + angle-table memory at binding time
    from functools import partial

    bindings = structured_fd_bindings(
        disc, initial, coefficient_encoder=partial(qram_coefficient_encoding, angle_width=8)
    )
    model = qham_input_model(plan, bindings, eta=-0.4)
    be = model.generator
    b = Builder(
        "aux_model_qram",
        {"aux": Bits(be.width), "target": Bits(be.width), "signal": Bits(be.signal_qubits)},
    )
    b.h(b["aux"])
    b.xor(b["aux"], b["target"])
    b.call(be.operation, target=b["target"], signal=b["signal"])
    opened = b.finish().program()
    (slot,) = [r.name for r in unresolved(opened)]
    (monomial,) = (
        term.monomial for port in disc.pde.ports for term in port.terms if term.monomial.known
    )
    table = qram_coefficient_memory(disc, monomial, angle_width=8)
    program = bind(
        opened, {slot: Binding(qram_database(2, 8).operation, {"table": "nu_angles"})}
    )
    state = rir_pysparq(program, {"nu_angles": table})
    uniform = 1 / math.sqrt(1 << be.width)
    qram_worst = 0.0
    for (aux, target, signal), amplitude in state.items():
        if signal != 0:
            continue
        expected = (
            matrix[target][aux] / be.alpha * uniform
            if target < raw and aux < raw
            else 0j
        )
        qram_worst = max(qram_worst, abs(amplitude - expected))
    angle_bound = 0.75 * math.pi / (1 << 8)
    report.case(
        "qham-input-models-consistency",
        paths=["rir-pysparq", "rir-pysparq+qram"],
        parameters={
            "raw_dimension": raw,
            "alphas": alphas | {"qram": be.alpha},
            "angle_width": 8,
        },
        metrics={
            "stencil_max_error": measured["stencil"][0],
            "spectral_max_error": measured["spectral"][0],
            "qram_max_error": qram_worst,
            "qram_angle_bound": angle_bound,
        },
        criterion=(
            "stencil/spectral full matrices < 1e-9; the QRAM angle-table deviation does not exceed alpha*pi/2^angle_width"
        ),
        passed=measured["stencil"][0] < 1e-9
        and measured["spectral"][0] < 1e-9
        and measured["stencil"][1] == measured["spectral"][1] == 0
        and measured["stencil"][2] == measured["spectral"][2] == 0
        and qram_worst < angle_bound,
    )


def verify_qram_coefficient_pointwise(report):
    """QRAM coefficient angle-encoding pointwise truth: the diagonal BE's amplitude at every address vs v/alpha."""
    u = Field("u")
    pde = PolynomialPDE.from_equations({"u": 0.5 * Known("nu") * u.d("x", 1) - 0.2 * u})
    grid = Grid(("x",), (4,), (1.0,))
    disc = Discretization(pde, grid, {"nu": [1.0, 0.5, 1.5, 0.25]})
    (monomial,) = (
        term.monomial for port in disc.pde.ports for term in port.terms if term.monomial.known
    )
    be = qram_coefficient_encoding(disc, monomial, angle_width=8)
    table = qram_coefficient_memory(disc, monomial, angle_width=8)
    b = Builder(
        "coeff_drv",
        {r.name: r.type for r in be.operation.module.registers},
    )
    b.h(b["target"])
    b.call(be.operation, target=b["target"], signal=b["signal"])
    opened = b.finish().program()
    (slot,) = [r.name for r in unresolved(opened)]
    program = bind(opened, {slot: Binding(qram_database(2, 8).operation, {"table": "nu"})})
    state = rir_pysparq(program, {"nu": table})
    origin = originir_ext(program, {"nu": table})
    values = [disc.known_product(monomial, row).real for row in range(4)]
    uniform = 0.5
    worst = 0.0
    for (address, signal), amplitude in state.items():
        if signal != 0:
            continue
        worst = max(worst, abs(amplitude - values[address] / be.alpha * uniform))
    widths = [r.type.width for r in program.main.registers]
    cross = statevector_error(origin, amplitudes_to_statevector(state, widths))
    bound = math.pi / (1 << 8) / 2
    report.case(
        "qram-coefficient-pointwise",
        paths=["rir-pysparq", "originir-ext"],
        parameters={"alpha": be.alpha, "angle_width": 8, "values": values},
        metrics={
            "max_amplitude_error": worst,
            "quantization_bound": bound,
            "cross_path_deviation": cross,
        },
        criterion="per-address amplitude deviation from v/alpha <= pi/2^(aw+1) (the angle quantization bound), and both paths agree",
        passed=worst <= bound and cross < 1e-9,
    )


def verify_qham_lifted_initial(report, context):
    plan, disc = _quadratic_problem()
    model = context["stencil_model"]
    u_in = [0.2, 0]
    lifted = disc.lift(plan, [u_in, [0, 0], [0, 0]])
    norm = math.sqrt(sum(abs(v) ** 2 for v in lifted))
    b = Builder(
        "initial_drv", {r.name: r.type for r in model.initial.operation.module.registers}
    )
    b.call(model.initial.operation, target=b["target"], work=b["work"])
    program = b.finish().program()
    total_qubits = sum(r.type.width for r in program.main.registers) + workspace_table(program)[
        program.entry
    ]
    state = rir_pysparq(program)
    worst = 0.0
    dirty_work = 0
    for (target, work), amplitude in state.items():
        if work != 0:
            dirty_work += 1
        worst = max(worst, abs(amplitude - lifted[target] / norm))
    report.case(
        "qham-lifted-initial",
        paths=["rir-pysparq"],
        parameters={
            "raw_dimension": len(lifted),
            "total_qubits": total_qubits,
            "note": "71 qubits exceeds the UniQC budget of 24; only the pysparq path is used",
        },
        metrics={
            "max_error": worst,
            "log_initial_norm": model.log_initial_norm,
            "classical_log_norm": math.log(norm),
            "dirty_work_states": dirty_work,
        },
        criterion="lifted initial state amplitude by amplitude < 1e-9 (vs the direct tensor-word weight construction), log norms agree, work returned clean",
        passed=worst < 1e-9
        and abs(model.log_initial_norm - math.log(norm)) < 1e-12
        and dirty_work == 0,
    )


def verify_qham_taylor_and_shift(report, context):
    plan, disc = _quadratic_problem()
    matrix = context["stencil_matrix"]
    raw = plan.raw_dimension(disc.dimension)
    time = 0.01
    initial = disc.lift(plan, [[0.2, 0], [0, 0], [0, 0]])
    expected = [
        initial[row] + time * sum(a * b for a, b in zip(matrix[row], initial, strict=True))
        for row in range(2)
    ]
    results = {}
    for tag in ("stencil_model", "spectral_model"):
        model = context[tag]
        state = model.solve(lambda g, i, t: taylor_qode(g, i, t, degree=1), time)
        amplitudes = rir_pysparq(state.operation.program())
        normalization = math.exp(model.log_initial_norm) * (1 + time * model.generator.alpha)
        worst = max(
            abs(amplitudes.get((row, 0), 0) - expected[row] / normalization) for row in range(2)
        )
        results[tag] = worst
    pairwise = abs(results["stencil_model"] - results["spectral_model"])
    # Explicit dissipative shift: G - mu*I acts on the full 2^w space (the padding subspace gets the -mu diagonal)
    model = context["stencil_model"]
    shifted = model.dissipative_shift()
    be = shifted.generator
    mu = model.generator.alpha
    size = 1 << be.width
    shifted_matrix = [
        [
            (matrix[r][c] if r < raw and c < raw else 0.0) - (mu if r == c else 0.0)
            for c in range(size)
        ]
        for r in range(size)
    ]
    worst_shift, missing, pad, _ = _aux_matrix(
        be.operation, be.width, be.signal_qubits, be.alpha, size, shifted_matrix, name="aux_shift"
    )
    report.case(
        "qham-taylor-solve-and-shift",
        paths=["rir-pysparq"],
        parameters={"time": time, "degree": 1, "shift": mu, "shifted_alpha": be.alpha},
        metrics={
            "taylor_stencil_error": results["stencil_model"],
            "taylor_spectral_error": results["spectral_model"],
            "solve_pairwise": pairwise,
            "shifted_matrix_error": worst_shift,
            "missing_support": missing,
            "growth_shift": shifted.growth_shift,
        },
        criterion=(
            "both input models' Taylor solutions match the classical (I+tG)Y_in < 1e-9; "
            "the shifted generator's full matrix on the complete 2^w space (padding subspace -mu diagonal included) < 1e-9"
        ),
        passed=max(results.values()) < 1e-9 and pairwise < 1e-9 and worst_shift < 1e-9
        and missing == 0 and pad == 0 and shifted.growth_shift == mu,
    )


def _exact_riccati(t, u0=0.2, a=-0.2, b=0.1):
    e = math.exp(a * t)
    return u0 * e / (1 + (b * u0 / a) * (1 - e))


def _ham_partial_sum(disc, order, eta, u0, final_time, dt):
    """Synchronous RK4 integration of the HAM recursion: U'_k = L U_k + prev_k (prev recursed by the homotopy-step algebra)."""
    dim = disc.dimension
    port_l = [[disc.entry("L", r, c) for c in range(dim)] for r in range(dim)]
    U = [list(u0)] + [[0.0] * dim for _ in range(order)]

    def derivs(state):
        out = [
            [
                sum(port_l[r][c] * state[0][c] for c in range(dim)).real
                for r in range(dim)
            ]
        ]
        previous = [0.0] * dim
        for k in range(1, len(state)):
            nonlinear = [0.0] * dim
            for indices in compositions(k - 1, 2):
                image = disc.apply_port("B_0", [list(state[i]) for i in indices])
                nonlinear = [x + y.real for x, y in zip(nonlinear, image, strict=True)]
            previous = [(1 + eta) * p - eta * n for p, n in zip(previous, nonlinear, strict=True)]
            out.append(
                [
                    sum(port_l[r][c] * state[k][c] for c in range(dim)).real + previous[r]
                    for r in range(dim)
                ]
            )
        return out

    for _ in range(round(final_time / dt)):
        k1 = derivs(U)
        s2 = [[u + 0.5 * dt * d for u, d in zip(uk, dk, strict=True)] for uk, dk in zip(U, k1, strict=True)]
        k2 = derivs(s2)
        s3 = [[u + 0.5 * dt * d for u, d in zip(uk, dk, strict=True)] for uk, dk in zip(U, k2, strict=True)]
        k3 = derivs(s3)
        s4 = [[u + dt * d for u, d in zip(uk, dk, strict=True)] for uk, dk in zip(U, k3, strict=True)]
        k4 = derivs(s4)
        U = [
            [u + dt / 6 * (a + 2 * b + 2 * c + d) for u, a, b, c, d in zip(uk, d1, d2, d3, d4, strict=True)]
            for uk, d1, d2, d3, d4 in zip(U, k1, k2, k3, k4, strict=True)
        ]
    return [sum(x) for x in zip(*U, strict=True)]


def verify_qham_homotopy_contraction(report):
    # (a) Homotopy weights against the explicit formulas
    weight_worst = 0.0
    for eta in (-1.0, -0.4, 0.2):
        for power in range(5):
            weight_worst = max(
                weight_worst,
                abs(HomotopyWeight("correction", power).evaluate(eta) - (-eta * (1 + eta) ** power)),
                abs(HomotopyWeight("physical", power).evaluate(eta) - (1 - (1 + eta) ** power)),
            )
    weight_worst = max(weight_worst, abs(HomotopyWeight().evaluate(-0.4) - 1.0))
    # (b) Generation-row rules (linear_action) vs the independent chain rule (chain_rule)
    rule_residual = 0.0
    plan, disc = _quadratic_problem()
    for eta in (-1.0, -0.4, 0.2):
        values = [[0.2, -0.1], [0.05, 0.07], [-0.03, 0.02]]
        lifted = disc.lift(plan, values)
        via_rules = disc.linear_action(plan, eta, lifted)
        via_chain = disc.chain_rule(plan, values, eta)
        rule_residual = max(
            rule_residual,
            max(abs(a - b) for a, b in zip(via_rules, via_chain, strict=True)),
        )
    u = Field("u")
    burgers = PolynomialPDE.from_equations(
        {"u": 0.1 * u.d("x", 2) - u * u.d("x") + Known("f")}, label="forced_burgers"
    )
    plan_b = QHAMPlan(burgers, 2)
    disc_b = Discretization(burgers, Grid(("x",), (4,), (1.0,)), {"f": [0.05, 0, -0.05, 0]})
    values_b = [[0.1, 0.2, 0, -0.1], [0.02, -0.01, 0.03, 0], [0.01, 0, -0.02, 0.01]]
    lifted_b = disc_b.lift(plan_b, values_b)
    residual_b = max(
        abs(a - b)
        for a, b in zip(
            disc_b.linear_action(plan_b, -0.4, lifted_b),
            disc_b.chain_rule(plan_b, values_b, -0.4),
            strict=True,
        )
    )
    rule_residual = max(rule_residual, residual_b)
    # (c) Contraction convergence of the HAM partial sums vs the exact Riccati solution; measured contraction factor vs |1+eta|
    final_time = 0.5
    convergence = {}
    ratio_ok = True
    for eta in (-0.8, -0.4):
        errors = []
        for order in range(1, 7):
            summed = _ham_partial_sum(disc, order, eta, [0.2, 0.0], final_time, 1e-3)
            errors.append(abs(summed[0] - _exact_riccati(final_time)))
        ratios = [errors[i + 1] / errors[i] for i in range(len(errors) - 1)]
        mean_ratio = sum(ratios) / len(ratios)
        convergence[str(eta)] = {"errors": errors, "mean_ratio": mean_ratio}
        ratio_ok = ratio_ok and abs(mean_ratio - abs(1 + eta)) < 0.05
    monotone = all(
        convergence[e]["errors"][-1] < convergence[e]["errors"][0]
        and convergence[e]["errors"][-1] < 1e-3
        for e in convergence
    )
    report.case(
        "qham-homotopy-contraction",
        paths=["classical-independent"],
        parameters={"final_time": final_time, "orders": [1, 6], "etas": [-0.8, -0.4]},
        metrics={
            "weight_max_error": weight_worst,
            "rule_vs_chain_residual": rule_residual,
            "errors_eta_-0.4": convergence["-0.4"]["errors"],
            "mean_ratio_-0.4": convergence["-0.4"]["mean_ratio"],
            "errors_eta_-0.8": convergence["-0.8"]["errors"],
            "mean_ratio_-0.8": convergence["-0.8"]["mean_ratio"],
        },
        criterion=(
            "weights exact; row rules vs chain rule < 1e-9; HAM errors contract with order and the mean contraction factor is approximately |1+eta|"
        ),
        passed=weight_worst < 1e-12 and rule_residual < 1e-9 and ratio_ok and monotone,
    )


def verify_pde_wrappers(report):
    """Numerical pass-through of DiscretePDE / make_qpde / qpde_solver from algorithms/pde.py."""
    laplacian = [[2.0 if r == c else (-1.0 if (r - c) % 4 in (1, 3) else 0.0) for c in range(4)] for r in range(4)]
    be = scale(-0.1, matrix_pauli_encoding(laplacian))
    prep = gate_state_prep([1.0, 0.0, 0.0, 0.0])
    discrete = DiscretePDE(be, prep, label="ring_heat")
    time = 0.05
    g_matrix = [[-0.1 * v for v in row] for row in laplacian]
    psi = [1.0, 0.0, 0.0, 0.0]
    alpha_e = 1 + time * be.alpha
    expected = [
        (psi[r] + time * sum(g_matrix[r][c] * psi[c] for c in range(4))) / alpha_e
        for r in range(4)
    ]
    via_make = make_qpde(lambda g, i, t: taylor_qode(g, i, t, degree=1))(discrete, time)
    via_qpde = qpde_solver(lambda d, t: taylor_qode(d.generator, d.initial, t, degree=1))(
        PDEInput(discrete), time
    )
    worst = 0.0
    for state in (via_make, via_qpde):
        amplitudes = rir_pysparq(state.operation.program())
        worst = max(
            worst,
            max(abs(amplitudes.get((row, 0), 0) - expected[row]) for row in range(4)),
        )
    report.case(
        "qham-pde-wrappers",
        paths=["rir-pysparq"],
        parameters={"time": time, "alpha_evolution": alpha_e},
        metrics={"max_error": worst},
        criterion="both wrapper chains' states match the classical (I+tG)psi/alpha < 1e-9",
        passed=worst < 1e-9,
    )


# ---------------------------------------------------------------------------
# QFVM cases
# ---------------------------------------------------------------------------


def verify_roe_face_pointwise(report):
    fmt = QFVM_FMT
    enc = fmt.encode
    slices = [
        (
            "rowcol",
            {"rho_l": enc(1.0), "m_l": enc(0.25), "e_l": enc(1.0),
             "rho_r": enc(1.25), "m_r": enc(0.0), "e_r": enc(1.0)},
            ["row", "col"],
        ),
        (
            "rho",
            {"m_l": 0, "m_r": 0, "e_l": enc(1.0), "e_r": enc(1.0), "row": 1, "col": 1,
             "rho_l": enc(1.0), "rho_r": enc(1.0)},
            ["rho_l_bit01", "rho_r_bit0"],
        ),
        (
            "momentum",
            {"rho_l": enc(1.25), "rho_r": enc(1.25), "e_l": enc(1.0), "e_r": enc(1.0),
             "row": 1, "col": 1, "m_l": 0, "m_r": 0},
            ["m_l_bit01", "m_r_bit0"],
        ),
    ]
    method_worst = 0.0
    total = 0
    for tag, base, hspec in slices:
        h_registers = []
        for spec in hspec:
            if spec.endswith("_bit01"):
                h_registers.append((spec[: -len("_bit01")], (0, 1)))
            elif spec.endswith("_bit0"):
                h_registers.append((spec[: -len("_bit0")], (0,)))
            else:
                h_registers.append((spec, (0, 1)))
        face = roe_face(fmt=fmt, entropy_delta=QFVM_DELTA)
        names = [r.name for r in face.module.registers]
        b = Builder("face_" + tag, {r.name: r.type for r in face.module.registers})
        for key, value in base.items():
            register = b[key]
            for bit in range(register.width):
                if (value >> bit) & 1:
                    b.x(register[bit])
        for key, bits in h_registers:
            register = b[key]
            for bit in bits:
                b.h(register[bit])
        b.call(face, **{r.name: b[r.name] for r in face.module.registers})
        state = rir_pysparq(b.finish().program(), max_steps=8_000_000)
        for basis, amplitude in state.items():
            if abs(abs(amplitude) - 1 / math.sqrt(len(state))) > 1e-9:
                raise AssertionError(f"roe_face/{tag}: non-uniform branch amplitude {amplitude}")
            values = dict(zip(names, basis, strict=True))
            args = [fmt.decode(values[k]) for k in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r")]
            left, right, flags = roe_face_fx(
                fmt, *args, values["row"], values["col"], gamma=1.4, entropy_delta=QFVM_DELTA
            )
            if values["left"] != left.raw or values["right"] != right.raw:
                raise AssertionError(
                    f"roe_face/{tag}: bit-level mismatch {values} vs emu ({left.raw}, {right.raw})"
                )
            if values["status"] != flags[0] + 2 * flags[1]:
                raise AssertionError(f"roe_face/{tag}: status mismatch {values['status']} vs {flags}")
            fl, fr = frozen_roe_face(*args, values["row"], values["col"], entropy_delta=QFVM_DELTA)
            method_worst = max(
                method_worst, abs(fmt.decode(values["left"]) - fl), abs(fmt.decode(values["right"]) - fr)
            )
            total += 1
    report.case(
        "roe-face-pointwise",
        paths=["rir-pysparq"],
        parameters={"format": "5.2", "entropy_delta": QFVM_DELTA, "branches": total},
        metrics={"exact_raw_match_branches": total, "method_max_error": method_worst},
        criterion=(
            "all superposition branches (16+8+8=32) have outputs raw and status bit-identical to the independent fixed-point simulation (any mismatch raises); "
            "the method error (vs the float64 Roe formulas) is reported separately"
        ),
        passed=True,
    )


def _qfvm_entry_driver(inputs, entry, column_x_bits, row_x_bits, row_h_bits, name):
    b = Builder(name, {r.name: r.type for r in entry.module.registers})
    for bit in column_x_bits:
        b.x(b["column"][bit])
    for bit in row_x_bits:
        b.x(b["row"][bit])
    for bit in row_h_bits:
        b.h(b["row"][bit])
    b.call(entry, row=b["row"], column=b["column"], data=b["data"])
    return bind_qfvm(b.finish().program(), inputs)


def verify_qfvm_entry_matrix(report):
    fmt = QFVM_FMT
    inputs = roe_qfvm_inputs(fmt=fmt)
    entry = qfvm_sparse_access(inputs, entropy_delta=QFVM_DELTA).entry
    d_float, d_emu = _qfvm_matrices(
        fmt, QFVM_STATES, mass=1.0, dx=1.0, gamma=1.4, entropy_delta=QFVM_DELTA
    )
    memory = dict(QFVM_MEMORY)
    memory["geometry"] = geometry_cells(inputs)
    # Column 5 (cell1/var1/half0), rows take half1's cells{0..3} x vars{0,1}: 8 branches
    program = _qfvm_entry_driver(inputs, entry, [0, 2], [4], [0, 2, 3], "entry_sample")
    state = rir_pysparq(program, memory, max_steps=60_000_000)
    impl_worst = 0.0
    method_worst = 0.0
    zero_ok = True
    for (row, column, data), amplitude in state.items():
        if abs(amplitude - 1 / math.sqrt(8)) > 1e-9:
            raise AssertionError(f"entry_sample: anomalous branch amplitude {amplitude}")
        emulated = d_emu[row][column]
        if data != emulated:
            raise AssertionError(
                f"entry_sample: ({row},{column}) raw {data} != emulation {emulated}"
            )
        impl_worst = max(impl_worst, abs(fmt.decode(data) - fmt.decode(emulated)))
        method_worst = max(method_worst, abs(fmt.decode(emulated) - d_float[row][column]))
        if emulated == 0 and data != 0:
            zero_ok = False
    report.case(
        "qfvm-entry-matrix-sample",
        paths=["rir-pysparq"],
        parameters={
            "column": 5,
            "rows": "half=1, cells 0-3, vars 0-1 (8-branch superposition)",
            "format": "5.2",
            "note": "the entry circuit has about 1.6e7 expansion gates; the superposition sweep covers this column's three non-zero bands and the zero structure",
        },
        metrics={
            "impl_max_error": impl_worst,
            "method_max_error": method_worst,
            "zero_structure_ok": zero_ok,
        },
        criterion="sampled entry raws are bit-identical to the independent fixed-point simulation; the method error vs the float64 matrix is reported separately",
        passed=zero_ok,
    )


def verify_qfvm_entry_padding(report):
    fmt = QFVM_FMT
    inputs = roe_qfvm_inputs(fmt=fmt)
    entry = qfvm_sparse_access(inputs, entropy_delta=QFVM_DELTA).entry
    memory = dict(QFVM_MEMORY)
    memory["geometry"] = geometry_cells(inputs)
    # Column 7 (cell1/var3/half0), rows {7, 23}: the exact diagonal receives padding, outside the half-block is zero
    program = _qfvm_entry_driver(inputs, entry, [0, 1, 2], [0, 1, 2], [4], "entry_padding")
    state = rir_pysparq(program, memory, max_steps=60_000_000)
    padding_raw = fmt.encode(1.0)
    matched = {}
    for (row, column, data), _amplitude in state.items():
        matched[(row, column)] = data
    diag_ok = matched.get((7, 7)) == padding_raw
    off_ok = matched.get((23, 7)) == 0
    report.case(
        "qfvm-entry-padding",
        paths=["rir-pysparq"],
        parameters={"column": 7, "rows": [7, 23], "padding_value": 1.0},
        metrics={"diagonal_raw": matched.get((7, 7)), "off_block_raw": matched.get((23, 7))},
        criterion="the padding variable's exact diagonal holds the padding_value raw, and outside the expanded block is zero",
        passed=diag_ok and off_ok and len(matched) == 2,
    )


def _geometry_expected(inputs):
    """Independently recompute the geometry table per the documented ABI (without calling geometry_cells)."""
    cw, w = inputs.cell_width, inputs.width
    n, table = 1 << cw, {}
    for row in range(1 << w):
        cell, var, half = (row >> 2) % n, row % 4, row >> (cw + 2)
        for slot in range(16):
            band, other_var = divmod(slot, 3)
            valid = int(slot < 9 and var < 3)
            neighbor_cell = (cell + band - 1) % n
            neighbor = other_var + 4 * neighbor_cell + ((1 - half) << (cw + 2))
            reverse = (2 - band) * 3 + var if valid else 0
            source = cell if half == 0 else neighbor_cell
            rowvar, colvar, source_band = (
                (var, other_var, band) if half == 0 else (other_var, var, 2 - band)
            )
            packed = 0
            offset = 0
            for value, width in (
                (neighbor, w), (reverse, 4), (source, cw), (rowvar, 2),
                (colvar, 2), (source_band, 2), (valid, 1),
            ):
                packed |= (value & ((1 << width) - 1)) << offset
                offset += width
            table[row + (slot << w)] = packed
    return table


def verify_qfvm_location(report):
    inputs = roe_qfvm_inputs(fmt=QFVM_FMT)
    geometry = geometry_cells(inputs)
    pointwise = sum(
        1 for key, value in _geometry_expected(inputs).items() if geometry.get(key) != value
    )
    location = qfvm_sparse_access(inputs, entropy_delta=QFVM_DELTA).location
    b = Builder(
        "loc_drv",
        {"column": Bits(5), "index": Bits(5), "record": Bits(5), "work": Bits(0)},
    )
    b.h(b["column"])
    b.h(b["index"])
    b.xor(b["index"], b["record"])
    b.call(location, column=b["column"], index=b["index"], work=b["work"])
    program = bind_qfvm(b.finish().program(), inputs)
    state = rir_pysparq(program, {"geometry": geometry}, max_states=2048)
    rows = {}
    dirty_work = 0
    for (column, index, record, work), _amplitude in state.items():
        if work != 0:
            dirty_work += 1
        rows.setdefault(column, {})[record] = index
    permutation_ok = all(set(mapping.values()) == set(range(32)) for mapping in rows.values())
    n = 4
    mismatches = 0
    for column, mapping in rows.items():
        cell, var, half = (column >> 2) % n, column % 4, column >> 4
        for slot in range(9):
            band, other_var = divmod(slot, 3)
            if var < 3:
                expected = other_var + 4 * ((cell + band - 1) % n) + ((1 - half) << 4)
            else:
                expected = (column + slot) % 32
            if mapping.get(slot) != expected:
                mismatches += 1
    report.case(
        "qfvm-location-permutation",
        paths=["rir-pysparq"],
        parameters={"columns": len(rows), "branches": len(state)},
        metrics={
            "geometry_pointwise_mismatch": pointwise,
            "permutation_ok": permutation_ok,
            "slot_mismatches": mismatches,
            "dirty_work_states": dirty_work,
        },
        criterion=(
            "geometry table pointwise identical; under the 32-column full superposition every column is a complete permutation, and the 9 slot mappings match the independent geometric semantics"
        ),
        passed=pointwise == 0
        and permutation_ok
        and mismatches == 0
        and dirty_work == 0
        and len(rows) == 32,
    )


def verify_qfvm_rhs_preparation(report):
    fmt, aw = QFVM_FMT, 8
    inputs = roe_qfvm_inputs(fmt=fmt, angle_width=aw)
    flow = RoeFlowData(QFVM_STATES, fmt=fmt, entropy_delta=QFVM_DELTA, angle_width=aw)
    memories = qfvm_memories(inputs, flow)
    # Independent residual: script-side fluxes + fmt quantization round trip
    fluxes = [
        _script_flux(QFVM_STATES[i], QFVM_STATES[(i + 1) % 4], gamma=1.4, entropy_delta=QFVM_DELTA)
        for i in range(4)
    ]
    residual = []
    for cell in range(4):
        for j in range(3):
            raw = fmt.encode((fluxes[(cell - 1) % 4][j] - fluxes[cell][j]) / 1.0)
            residual.append(fmt.decode(raw))
        residual.append(0.0)
    norm = math.sqrt(sum(v * v for v in residual))
    b = Builder("rhs_drv", {"target": Bits(5), "work": Bits(13)})
    invoke(b, inputs.rhs.operation, "rhs", target=b["target"], work=b["work"])
    program = bind_qfvm(b.finish().program(), inputs)
    state = rir_pysparq(
        program, {"rhs_angles": memories["rhs_angles"], "rhs_sign": memories["rhs_sign"]}
    )
    amplitude_worst = 0.0
    sign_bad = 0
    works = set()
    for (target, work), amplitude in state.items():
        works.add(work)
        expected = residual[target] / norm if target < 16 else 0.0
        amplitude_worst = max(amplitude_worst, abs(abs(amplitude) - abs(expected)))
        if expected != 0 and abs(amplitude) > 1e-6 and (amplitude.real > 0) != (expected > 0):
            sign_bad += 1
    total_prob = sum(abs(v) ** 2 for v in state.values())
    # Angle/sign banks pointwise truth (against an independent sum-of-squares tree)
    tree = [0.0] * (8 * 4)
    for address, value in enumerate(residual):
        tree[16 + address] = value * value
    for node in range(15, 0, -1):
        tree[node] = tree[2 * node] + tree[2 * node + 1]
    step = 2 * math.pi / (1 << aw)
    angle_bad = 0
    for node in range(1, 16):
        angle = 0.0 if tree[node] == 0 else 2 * math.acos(math.sqrt(tree[2 * node] / tree[node]))
        expected_word = round(angle / step) % (1 << aw)
        if memories["rhs_angles"].get(node - 1, 0) != expected_word:
            angle_bad += 1
    sign_bad_bank = sum(
        1
        for address, value in enumerate(residual)
        if memories["rhs_sign"].get(address, 0) != int(value < 0)
    )
    report.case(
        "qfvm-rhs-preparation",
        paths=["rir-pysparq"],
        parameters={"angle_width": aw, "rhs_norm": norm},
        metrics={
            "amplitude_max_error": amplitude_worst,
            "sign_mismatches": sign_bad,
            "total_probability": total_prob,
            "work_patterns": len(works),
            "angle_word_mismatches": angle_bad,
            "sign_bank_mismatches": sign_bad_bank,
            "flow_norm_vs_tree": abs(flow.rhs_norm - norm),
        },
        criterion=(
            "prepared amplitudes match the independent residual/norm (angle quantization tolerance 0.01), signs exact, "
            "angle tree/sign bank pointwise truth"
        ),
        passed=amplitude_worst < 0.01
        and sign_bad == 0
        and abs(total_prob - 1) < 1e-9
        and angle_bad == 0
        and sign_bad_bank == 0
        and abs(flow.rhs_norm - norm) < 1e-12,
    )


def verify_qfvm_single_step(report):
    """Classical identities of the flow_data single-step update and the Roe flux (vs an independent matrix-inversion implementation)."""
    from oracq.applications.flow_data import riemann_flux

    gamma, delta = 1.4, 0.125
    flux_worst = 0.0
    consistency_worst = 0.0
    samples = [
        ((1.0, 0.2, 2.0), (1.1, -0.1, 2.2)),
        ((0.8, 0.0, 1.5), (1.2, 0.3, 2.5)),
        ((1.3, -0.4, 3.0), (0.9, 0.1, 1.8)),
        ((1.0, 0.5, 2.0), (1.0, 0.5, 2.0)),
    ]
    for left_state, right_state in samples:
        via_inverse = riemann_flux(left_state, right_state, gamma=gamma, entropy_delta=delta)
        via_blocks = _script_flux(left_state, right_state, gamma=gamma, entropy_delta=delta)
        flux_worst = max(
            flux_worst,
            max(abs(a - b) for a, b in zip(via_inverse, via_blocks, strict=True)),
        )
        if left_state == right_state:
            rho, m, e = left_state
            u, p = m / rho, (gamma - 1) * (e - 0.5 * m * m / rho)
            physical = (m, m * u + p, u * (e + p))
            consistency_worst = max(
                consistency_worst,
                max(abs(a - b) for a, b in zip(via_blocks, physical, strict=True)),
            )
    # M.u - mass.u == -residual: the matrix encodes mass*I minus the residual flux difference (the standard sign
    # convention of the implicit FVM linearization; residual is defined by flow_data as F*_{i-1} - F*_i, and M is
    # built from the frozen_roe_face three bands)
    d_float, _ = _qfvm_matrices(
        QFVM_FMT, QFVM_STATES, mass=1.0, dx=1.0, gamma=1.4, entropy_delta=QFVM_DELTA
    )
    m_matrix = [row[16:] for row in d_float[:16]]
    u_vec = [component for cell in QFVM_STATES for component in (*cell, 0.0)]
    fluxes = [
        _script_flux(QFVM_STATES[i], QFVM_STATES[(i + 1) % 4], gamma=1.4, entropy_delta=QFVM_DELTA)
        for i in range(4)
    ]
    residual_worst = 0.0
    for i in range(16):
        cell, var = i // 4, i % 4
        mu = sum(m_matrix[i][j] * u_vec[j] for j in range(16)) - 1.0 * u_vec[i]
        expected = (
            (fluxes[cell][var] - fluxes[(cell - 1) % 4][var]) / 1.0 if var < 3 else 0.0 - u_vec[i]
        )
        residual_worst = max(residual_worst, abs(mu - expected))
    # The RoeFlowData residual bank matches the independent quantized residual; local update patches are exact
    flow = RoeFlowData(QFVM_STATES, fmt=QFVM_FMT, entropy_delta=QFVM_DELTA, angle_width=8)
    quantized_worst = 0.0
    for cell in range(4):
        for j in range(3):
            raw = QFVM_FMT.encode((fluxes[(cell - 1) % 4][j] - fluxes[cell][j]) / 1.0)
            quantized_worst = max(
                quantized_worst, abs(flow.residuals[cell][j] - QFVM_FMT.decode(raw))
            )
    changed = {2: (1.0, 0.0, 1.25)}
    patch = flow.update(changed)
    fresh = RoeFlowData(
        [changed.get(i, s) for i, s in enumerate(QFVM_STATES)],
        fmt=QFVM_FMT, entropy_delta=QFVM_DELTA, angle_width=8,
    )
    bank_diff = sum(
        1
        for bank, cells in flow.store.snapshot().items()
        for address, value in cells.items()
        if fresh.store.snapshot().get(bank, {}).get(address, 0) != value
    )
    locality_ok = (
        set(patch.recomputed_faces) == {1, 2} and set(patch.recomputed_cells) == {1, 2, 3}
    )
    # ptheta pointwise truth
    table = ptheta_cells(QFVM_FMT, 8, 8.0)
    ptheta_bad = sum(
        1
        for raw, word in table.items()
        if word
        != round(
            2 * math.acos(min(1, abs(QFVM_FMT.decode(raw)) / 8.0)) * (1 << 8) / (2 * math.pi)
        )
        % (1 << 8)
    )
    report.case(
        "qfvm-single-step-identities",
        paths=["classical-independent"],
        parameters={"cells": 4, "entropy_delta": delta},
        metrics={
            "flux_block_vs_inverse": flux_worst,
            "flux_consistency": consistency_worst,
            "matrix_vector_vs_residual": residual_worst,
            "flow_bank_vs_independent": quantized_worst,
            "local_patch_bank_diff": bank_diff,
            "locality_ok": locality_ok,
            "ptheta_mismatches": ptheta_bad,
        },
        criterion=(
            "F* = left.L + right.R vs the inversion implementation < 1e-9; M.u - mass.u = -residual (implicit sign convention) < 1e-9; "
            "local updates match a full recomputation; ptheta exact pointwise"
        ),
        passed=flux_worst < 1e-9
        and consistency_worst < 1e-9
        and residual_worst < 1e-9
        and quantized_worst == 0.0
        and bank_diff == 0
        and locality_ok
        and ptheta_bad == 0,
    )


def run():
    report = Report(
        "qham_qfvm",
        "Publication-grade numerical validation of the QHAM generator/homotopy steps/initial state/solve chain and the QFVM Roe fluxes/entries/locations/RHS.",
    )
    context = {}
    verify_qham_generator_stencil(report, context)
    verify_qham_generator_spectral(report, context)
    verify_qham_generator_cross_paths(report, context)
    verify_qham_input_models(report)
    verify_qram_coefficient_pointwise(report)
    verify_qham_lifted_initial(report, context)
    verify_qham_taylor_and_shift(report, context)
    verify_qham_homotopy_contraction(report)
    verify_pde_wrappers(report)
    verify_roe_face_pointwise(report)
    verify_qfvm_entry_matrix(report)
    verify_qfvm_entry_padding(report)
    verify_qfvm_location(report)
    verify_qfvm_rhs_preparation(report)
    verify_qfvm_single_step(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
