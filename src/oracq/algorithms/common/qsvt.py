"""Standard transformation library for quantum singular value transformation, QSVT
(Gilyén et al. 2019, arXiv:1806.01838).

Conventions, aligned step by step with transforms.qsvt_sequence and verified
pointwise against the reference simulator: on the two-dimensional invariant
subspace of a singular value x, ``qsvt_sequence(a, Φ)`` realizes

    p(x) = [S(φ_0) W(x) S(φ_1) W(x) … W(x) S(φ_d)]_00,
    W(x) = [[x, s], [s, -x]],  s = √(1 − x²),  S(φ) = diag(e^{iφ}, e^{−iφ}),

with the phases ordered in time, φ_0 acting first, for d block-encoding calls
and d+1 phases in total. The pair, P and Q, is realizable if and only if
deg P ≤ d, deg Q ≤ d−1, the parity of P equals d mod 2, the parity of Q equals
(d−1) mod 2, and the polynomial identity P P̄ + (1−x²) Q Q̄ ≡ 1 holds; in
particular ``|P(±1)| = 1`` must hold, the endpoint saturation.

Phase synthesis follows the Gilyén paper's complementary-polynomial root
finding plus layer stripping: given a real target f and an optional imaginary
completion h, both real-coefficient polynomials, set P = f + i·h, construct
Q from the roots of R = (1 − f² − h²)/(1−x²), then recover the phases by
backward recursion. Real targets that do not saturate at the endpoints, such
as truncated approximations of 1/x, must use a nonzero imaginary completion;
the real part is then extracted by the LCU combination of (U_Φ + U_{−Φ})/2,
where −Φ realizes exactly P̄, and the block encoding equals f(A/α).

Numerical boundaries: root finding uses pure-Python Durand–Kerner iteration
and stripping uses double-precision complex arithmetic; the synthesis degree
is capped at 40, and every synthesis is followed by a round-trip self-check
with qsp_response that raises ValidationError beyond the limits. For degrees
within a few tens and well-separated complementary-polynomial roots, the
round-trip error is typically on the order of 1e-9.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Iterable, Iterator, Sequence
from typing import cast

from oracq.algorithms.common.transforms import qsvt_sequence
from oracq.algorithms.input_model.contracts import require_instance
from oracq.algorithms.input_model.operators import BlockEncoding, identity, linear_combination
from oracq.algorithms.input_model.oracles import annotate
from oracq.infrastructure.ir import ValidationError

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
# Polynomial utilities: ascending coefficients with the constant term first;
# real coefficients use float and complex coefficients use complex.
# ---------------------------------------------------------------------------


def _trim(p: Sequence[float | complex], tol: float = 1e-12) -> tuple[float | complex, ...]:
    """Drop trailing coefficients below the tolerance and tighten the polynomial into a tuple."""
    p = list(p)
    while len(p) > 1 and abs(p[-1]) <= tol * max(1.0, max(abs(c) for c in p)):
        p.pop()
    return tuple(p)


def _add(
    a: Sequence[float | complex], b: Sequence[float | complex]
) -> tuple[float | complex, ...]:
    """Add two polynomials coefficient by coefficient, zero-padding the shorter one."""
    n = max(len(a), len(b))
    return tuple((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n))


def _sub(
    a: Sequence[float | complex], b: Sequence[float | complex]
) -> tuple[float | complex, ...]:
    """Subtract two polynomials coefficient by coefficient."""
    return _add(a, tuple(-v for v in b))


def _mul(
    a: Sequence[float | complex], b: Sequence[float | complex]
) -> tuple[float | complex, ...]:
    """Multiply two polynomials by convolution."""
    out = [0j] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def _scale(c: float | complex, p: Sequence[float | complex]) -> tuple[float | complex, ...]:
    """Multiply every coefficient of the polynomial by the same scalar."""
    return tuple(c * v for v in p)


def _eval(p: Sequence[float | complex], x: float | complex) -> complex:
    """Evaluate the polynomial at ``x`` with Horner's method."""
    v = 0j
    for c in reversed(p):
        v = v * x + c
    return v


def _conj(p: Sequence[float | complex]) -> tuple[complex, ...]:
    """Take the complex conjugate coefficient by coefficient."""
    return tuple(complex(v).conjugate() for v in p)


def _realify(p: Sequence[float | complex], *, tol: float = 1e-10) -> tuple[float, ...]:
    """Verify the imaginary parts are near zero and tighten the coefficients into a ``float`` tuple."""
    scale = max(1.0, max(abs(c) for c in p))
    if any(abs(complex(c).imag) > tol * scale for c in p):
        raise ValidationError("An internal polynomial must have real coefficients")
    return tuple(float(complex(c).real) for c in p)


def _deg(p: Sequence[float | complex]) -> int:
    """Return the degree of the polynomial given by ascending coefficients."""
    return len(p) - 1


def _grid(n: int = _GRID) -> tuple[float, ...]:
    """Return the Chebyshev cosine-node grid on [−1, 1]."""
    return tuple(math.cos(math.pi * j / n) for j in range(n + 1))


def _sup_norm(p: Sequence[float | complex]) -> float:
    """Return the maximum magnitude of the polynomial on the Chebyshev grid."""
    return max(abs(_eval(p, x)) for x in _grid())


def _roots(
    coeffs: Iterable[float | complex], *, iters: int = 4000, tol: float = 1e-30
) -> list[complex]:
    """Simultaneous Durand–Kerner root finding; accepts arbitrary real or complex ascending
    coefficients before monic normalization."""
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


def _cluster(roots: Iterable[complex], *, rel: float = 2e-5) -> list[tuple[complex, int]]:
    """Cluster numerical roots into multiple-root groups by relative distance, returning
    [(centroid, multiplicity), …]."""
    clusters: list[list[int | complex]] = []
    for r in sorted(roots, key=lambda z: (abs(z), z.real, z.imag)):
        for cluster in clusters:
            if abs(r - cluster[0]) <= rel * max(1.0, abs(r)):
                cluster[1] += 1
                break
        else:
            clusters.append([r, 1])
    return [(cast("complex", c), cast("int", m)) for c, m in clusters]


def _chebyshev_t(n: int) -> tuple[float, ...]:
    """Ascending monomial coefficients of T_n."""
    if n == 0:
        return (1.0,)
    if n == 1:
        return (0.0, 1.0)
    a: tuple[float, ...]
    b: tuple[float, ...]
    a, b = (1.0,), (0.0, 1.0)
    for _ in range(2, n + 1):
        a, b = b, cast("tuple[float, ...]", _sub(_scale(2.0, (0.0,) + b), a))
    return b


# ---------------------------------------------------------------------------
# QSP response and phase synthesis, in the reflection convention shared with
# qsvt_sequence.
# ---------------------------------------------------------------------------


def qsp_response(x: float, phases: Iterable[float]) -> complex:
    """Top-left block p(x) realized by the phase sequence Φ in the reflection
    convention, for x ∈ [−1, 1].

    Args:
        x: Evaluation point, a real number in [−1, 1].
        phases: QSP phase sequence Φ, in radians.

    Returns:
        complex: Value of the response polynomial p(x) realized by the phase sequence.
    """
    phases = tuple(float(p) for p in phases)
    s = math.sqrt(max(0.0, 1.0 - x * x))
    m00, m01, m10, m11 = 1 + 0j, 0j, 0j, 1 + 0j
    for i, phi in enumerate(phases):
        if i:  # W(x) applied on the left, earlier in time than this step's phase
            m00, m01, m10, m11 = (
                x * m00 + s * m10,
                x * m01 + s * m11,
                s * m00 - x * m10,
                s * m01 - x * m11,
            )
        e, f = cmath.exp(1j * phi), cmath.exp(-1j * phi)
        m00, m01, m10, m11 = e * m00, e * m01, f * m10, f * m11
    return m00


def _check_real_poly(coeffs: Iterable[float | complex], label: str) -> tuple[float, ...]:
    """Validate finite real coefficients and tighten them into a ``float`` tuple."""
    values: list[float] = []
    for c in coeffs:
        z = complex(c)
        if not math.isfinite(z.real) or abs(z.imag) > 1e-12:
            raise ValidationError(f"{label} must be a polynomial with finite real coefficients")
        values.append(z.real)
    if not values:
        raise ValidationError(f"{label} must not be empty")
    return cast("tuple[float, ...]", _trim(values))


def _check_parity(coeffs: Sequence[float | complex], d: int, label: str) -> None:
    """Verify that every nonzero power of the polynomial has the same parity as ``d``."""
    scale = max(1.0, max(abs(c) for c in coeffs))
    for k, c in enumerate(coeffs):
        if (d - k) % 2 and abs(c) > 1e-9 * scale:
            raise ValidationError(f"{label} has parity inconsistent with the degree d mod 2")


def _div_1mx2(dpoly: Sequence[float | complex], *, tol: float = 1e-8) -> tuple[float, ...]:
    """Compute R = D/(1−x²); requires D(±1) = 0, otherwise a ValidationError is raised."""
    scale = max(1.0, max(abs(c) for c in dpoly))
    n = len(dpoly)
    r = [0j] * (n + 2)
    for j in range(n + 2):
        r[j] = (dpoly[j] if j < n else 0) + (r[j - 2] if j >= 2 else 0)
    if max(abs(r[n - 2]), abs(r[n - 1])) > tol * scale:
        raise ValidationError(
            "Complementary polynomial condition failed: 1 − f² − h² is not divisible by 1−x²"
            " because the endpoints are not saturated"
        )
    return _realify(_trim(r[: max(1, n - 2)]))


def _q_from_roots(rpoly: Sequence[float], d: int) -> tuple[float | complex, ...] | None:
    """Construct the complex-coefficient Q with parity (d−1) mod 2 from the roots of
    R = (1−f²−h²)/(1−x²) such that Q Q̄ = R."""
    rpoly = cast("tuple[float, ...]", _trim(rpoly))
    if _deg(rpoly) <= 0:
        if rpoly[0] <= 0:
            raise ValidationError("The complementary polynomial R is nonpositive everywhere,"
                                  " so no spectral decomposition exists")
        return (math.sqrt(rpoly[0]),) if d % 2 == 1 else None
    if _deg(rpoly) != 2 * d - 2:
        raise ValidationError("The complementary polynomial degree does not match the target")
    if rpoly[-1] <= 0:
        raise ValidationError("The leading coefficient of the complementary polynomial must be"
                              " positive")
    clusters = _cluster(_roots(rpoly))
    factors: list[tuple[float | complex, ...]] = []  # even polynomial factors; zero roots alone carry the parity
    zero_mult = 0
    used: set[int] = set()
    ctol = 1e-5
    for i, (c, m) in enumerate(clusters):
        if i in used:
            continue
        if abs(c) <= ctol:
            zero_mult += m
            used.add(i)
            continue
        if abs(c.imag) <= ctol * abs(c):  # real root r together with −r
            partner = next(
                (j for j, (c2, m2) in enumerate(clusters) if j > i and j not in used
                 and abs(c2 + c) <= ctol * max(1.0, abs(c)) and m2 == m),
                None,
            )
            if partner is None or m % 2:
                raise ValidationError("Real roots of the complementary polynomial must come in"
                                      " pairs and have even multiplicity")
            rr = (abs(c.real) + abs(clusters[partner][0].real)) / 2
            for _ in range(m // 2):
                factors.append((-(rr * rr), 0.0, 1.0))
            used.add(partner)
        elif abs(c.real) <= ctol * abs(c):  # purely imaginary root i b with −i b, mutual conjugates
            partner = next(
                (j for j, (c2, m2) in enumerate(clusters) if j > i and j not in used
                 and abs(c2 + c) <= ctol * max(1.0, abs(c)) and m2 == m),
                None,
            )
            if partner is None or m % 2:
                raise ValidationError("Purely imaginary roots of the complementary polynomial"
                                      " must have even multiplicity")
            bb = (abs(c.imag) + abs(clusters[partner][0].imag)) / 2
            for _ in range(m // 2):
                factors.append((bb * bb, 0.0, 1.0))  # x² − (ib)² = x² + b²
            used.add(partner)
        else:  # generic complex root: the quadruple {±ρ, ±ρ̄} each with multiplicity m; the symmetrized centroid reduces numerical bias
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
                raise ValidationError("Complex roots of the complementary polynomial must form"
                                      " conjugate-negated quadruples")
            members = [c] + [clusters[j][0] for j in sib]
            ra = sum(abs(v.real) for v in members) / 4
            rb = sum(abs(v.imag) for v in members) / 4
            rho = complex(ra, rb)
            for _ in range(m):
                factors.append((-rho * rho, 0.0, 1.0))  # (x−ρ)(x+ρ) = x² − ρ²
            used.update(sib)
        used.add(i)
    if len(used) != len(clusters):
        raise ValidationError("The complementary polynomial root structure is incomplete,"
                              " so no spectral decomposition exists")
    if zero_mult % 2:
        raise ValidationError("Zero roots of the complementary polynomial must have even"
                              " multiplicity")
    m0 = zero_mult // 2
    if m0 % 2 != (d - 1) % 2:
        raise ValidationError("The parity of the complementary polynomial spectral factor"
                              " contradicts the degree requirement")
    q: tuple[float | complex, ...] = (math.sqrt(float(rpoly[-1])),)
    for f in factors:
        q = _mul(q, f)
    q = _mul(q, (0.0, 1.0)) if m0 else q
    q = _trim(tuple(0.0 if abs(v) < 1e-12 else v for v in _scale(1.0, q)))
    # Parity: keep only the powers with the same parity as d−1, zeroing numerical noise
    q = tuple(v if (d - 1 - k) % 2 == 0 else 0.0 for k, v in enumerate(q))
    resid = _sub(_trim(_mul(q, _conj(q))), rpoly)
    if max((abs(v) for v in resid), default=0.0) > 1e-5 * max(1.0, max(abs(c) for c in rpoly)):
        raise ValidationError("The complementary polynomial spectral decomposition failed its"
                              " self-check because of insufficient numerical precision")
    return q


def _strip(
    ppoly: Sequence[float | complex], qpoly: Sequence[float | complex], d: int
) -> tuple[float, ...]:
    """Layer stripping: recover the phases layer by layer from (P, Q), returned in time order."""
    p, q = list(ppoly), list(qpoly)
    phases: list[float] = []
    for k in range(d, 0, -1):
        if abs(q[k - 1]) < 1e-13:
            raise ValidationError("Layer stripping degenerated: the leading coefficient of the"
                                  " complementary polynomial is too small and the phases are"
                                  " numerically unstable")
        phi = cmath.phase(p[k] / q[k - 1]) / 2
        phases.append(phi)
        ei, ej = cmath.exp(-1j * phi), cmath.exp(1j * phi)
        pn = _add(_scale(ei, (0j,) + tuple(p)), _scale(ej, _sub(q, (0j, 0j) + tuple(q))))
        qn = _sub(_scale(ei, p), _scale(ej, (0j,) + tuple(q)))
        # In theory the leading coefficients of x^k, x^{k+1} and x^{k−1} cancel exactly; residual errors are caught by the self-check
        p, q = list(_trim(pn[:k])), list(_trim(qn[: max(1, k - 1)]))
        if len(p) < k:
            p += [0j] * (k - len(p))
        if len(q) < max(1, k - 1):
            q += [0j] * (max(1, k - 1) - len(q))
    if abs(abs(p[0]) - 1.0) > 1e-6:
        raise ValidationError("Layer stripping self-check failed: the magnitude of the zero-level"
                              " phase deviates from 1")
    phases.append(cmath.phase(p[0]))
    return tuple(reversed(phases))


def qsp_phases(coeffs: Iterable[float], imag: Iterable[float] | None = None) -> tuple[float, ...]:
    """Synthesize a QSP phase sequence from a real-coefficient target polynomial, in
    time order with length d+1.

    coeffs holds ascending real coefficients with the constant term first; the
    target is P = f when imag is None, or P = f + i·h when imag holds the
    ascending real coefficients of h. The realizability conditions are: f and h
    both have parity d mod 2; f² + h² ≤ 1 on [−1,1]; endpoint saturation
    f(±1)² + h(±1)² = 1; and R = (1 − f² − h²)/(1−x²) is nonnegative and
    satisfies the root multiplicity conditions of a spectral decomposition.
    Omitting imag means a purely real target, which then requires
    ``|f(±1)| = 1``.

    Args:
        coeffs: Ascending real coefficients of the target real part f, constant term first.
        imag: Ascending real coefficients of the imaginary completion h; omit for a purely real target.

    Returns:
        tuple[float, ...]: QSP phase sequence in time order, with length equal to the target degree plus one.
    """
    f = _check_real_poly(coeffs, "target polynomial")
    d = _deg(f)
    if d > _MAX_DEGREE:
        raise ValidationError(f"Target polynomial degree exceeds the synthesis limit {_MAX_DEGREE}")
    if d == 0:
        if abs(abs(f[0]) - 1.0) > 1e-9:
            raise ValidationError("A degree-zero target must be a unit-modulus constant")
        return (cmath.phase(f[0]),)
    _check_parity(f, d, "target polynomial")
    h = _check_real_poly(imag, "imaginary completion") if imag is not None else (0.0,)
    if _trim(h) != (0.0,):
        if _deg(h) > d:
            raise ValidationError("The imaginary completion degree must not exceed the target"
                                  " degree")
        _check_parity(h + (0.0,) * (d + 1 - len(h)), d, "imaginary completion")
    if _sup_norm(f) > 1.0 + 1e-6 and imag is None:
        raise ValidationError("The target polynomial exceeds the upper bound 1 on the interval"
                              " from −1 to 1")
    bound = max(abs(complex(_eval(f, x), _eval(h, x))) for x in _grid())
    if bound > 1.0 + 1e-6:
        raise ValidationError("P = f + i·h leaves the unit disk on the interval from −1 to 1")
    sat = max(abs(_eval(f, 1.0).real ** 2 + _eval(h, 1.0).real ** 2 - 1.0),
              abs(_eval(f, -1.0).real ** 2 + _eval(h, -1.0).real ** 2 - 1.0))
    if sat > 1e-6:
        raise ValidationError("Endpoints are not saturated: f and h must satisfy the saturation"
                              " identity at x = 1 and x = -1; provide an imaginary completion")
    dpoly = _trim(_sub((1.0,), _add(_mul(f, f), _mul(h, h))))
    rpoly = _div_1mx2(dpoly)
    if min((_eval(rpoly, x).real for x in _grid()), default=0.0) < -1e-9:
        raise ValidationError("The complementary polynomial R takes negative values on the"
                              " interval from −1 to 1, so no spectral decomposition exists")
    q = _q_from_roots(rpoly, d)
    if q is None:
        raise ValidationError("The parity of the complementary polynomial spectral factor"
                              " contradicts the degree requirement")
    ppoly = _add(f, tuple(1j * v for v in h))
    phases = _strip(ppoly, q, d)
    err = max(abs(qsp_response(x, phases) - _eval(ppoly, x)) for x in _grid(512))
    if err > _STRIP_TOL:
        raise ValidationError(f"Phase synthesis round-trip self-check failed with error {err:.2e}:"
                              " degree too high or ill-conditioned complementary polynomial")
    return phases


# ---------------------------------------------------------------------------
# Assembly helpers: phase sequence to block encoding, and the real-part
# extraction LCU.
# ---------------------------------------------------------------------------


def _wrap_qsvt_be(a: BlockEncoding, phases: Iterable[float]) -> BlockEncoding:
    """Assemble the QSVT operation from the phase sequence and wrap it as a block encoding."""
    return BlockEncoding(
        annotate(qsvt_sequence(a, phases), "block_encoding", be_alpha=1.0)
    )


def _real_qsvt_be(a: BlockEncoding, phases: Iterable[float]) -> BlockEncoding:
    """Block encoding of (P + P̄)(A/α)/2 = f(A/α): −Φ realizes exactly P̄, and the
    equal-weight LCU extracts the real part."""
    plus = _wrap_qsvt_be(a, phases)
    minus = _wrap_qsvt_be(a, tuple(-p for p in phases))
    return linear_combination(0.5, plus, 0.5, minus)


def _finish(be: BlockEncoding, algorithm: str, **attributes: float) -> BlockEncoding:
    """Record the algorithm name and quantitative attributes on the final block encoding."""
    return BlockEncoding(
        annotate(
            be.operation,
            "block_encoding",
            be_alpha=be.alpha,
            algorithm=algorithm,
            **attributes,
        )
    )


def _synthesize_with_imag(
    f: Sequence[float], imag_candidates: Iterable[Iterable[float]], label: str
) -> tuple[float, ...]:
    """Try the imaginary completion candidates in turn; raise ValidationError when all fail."""
    for h in imag_candidates:
        try:
            return qsp_phases(f, imag=h)
        except ValidationError:
            continue
    raise ValidationError(f"Imaginary completion failed for {label}: reduce the target degree"
                          " or relax the parameters")


# ---------------------------------------------------------------------------
# Standard transformation family.
# ---------------------------------------------------------------------------


def qsvt_matrix_inversion(a: BlockEncoding, kappa: float, *, error: float = 0.05) -> BlockEncoding:
    """QSVT block encoding approximating A⁻¹, a scaling of the odd extension polynomial
    J_b(x) = (1−(1−x²)^b)/x.

    The target polynomial is f(x) = c·J_b(x), with ascending coefficients
    (−1)^m C(b, m+1) given analytically; it approximates c/x on ``|x| ≥ 1/κ``
    with relative error at most error, and ``||f||∞ ≤ 1/3``. The zero-signal
    block of the returned block encoding approximates inverse_scale · A⁻¹, with
    inverse_scale = c·α.

    Args:
        a: Block encoding of the matrix A to invert.
        kappa: Condition number κ, a finite number not less than 1.
        error: Relative approximation error, in (0,1).

    Returns:
        BlockEncoding: Matrix inversion block encoding whose zero-signal block approximates
        inverse_scale·A⁻¹.
    """
    require_instance(a, BlockEncoding, "qsvt_matrix_inversion.a")
    if not (math.isfinite(kappa) and kappa >= 1):
        raise ValidationError("The condition number κ must be a finite number not less than 1")
    if not (0 < error < 1):
        raise ValidationError("The approximation error must lie strictly between 0 and 1")
    if kappa == 1:
        b = 1
    else:
        b = max(1, math.ceil(math.log(1 / error) / -math.log(1 - 1 / kappa**2)))
    d = 2 * b - 1
    if d > _MAX_DEGREE:
        raise ValidationError(f"κ={kappa} requires degree {d} which exceeds the synthesis limit"
                              f" {_MAX_DEGREE}; please relax error")
    f = tuple(
        ((-1.0) ** m) * math.comb(b, m + 1) if i % 2 == 1 else 0.0
        for i in range(d + 1)
        for m in [i // 2]
    )
    f = cast("tuple[float, ...]", _trim(f))
    norm = _sup_norm(f)
    c_scale = 1.0 / (3.0 * norm)
    f = cast("tuple[float, ...]", _scale(c_scale, f))
    sat = math.sqrt(max(0.0, 1.0 - _eval(f, 1.0).real ** 2))
    phases = _synthesize_with_imag(f, [(0.0, sat)], "matrix inversion polynomial")
    be = _real_qsvt_be(a, phases)
    return _finish(
        be,
        "qsvt_matrix_inversion",
        kappa=float(kappa),
        error=float(error),
        qsp_degree=d,
        inverse_scale=c_scale * a.alpha,
    )


def eigenstate_filter(
    a: BlockEncoding, gap: float, degree: int, *, center: float = 0.0
) -> BlockEncoding:
    """Eigenstate filtering: a peaked polynomial block encoding that suppresses the
    block below 1/T_d(r) outside ``|x−center| ≤ gap``.

    The target is the Lin–Tong filter polynomial f(x) = T_d(g(x²))/T_d(r),
    with g(y) = 2(y−Δ²)/(1−Δ²) − 1 and r = (1+Δ²)/(1−Δ²); f saturates at
    x=0, ``|f(0)|=1``, and on ``|x| ≥ Δ`` it obeys ``|f| ≤ 1/T_d(r)``. A
    nonzero center first shifts the spectrum via a block-encoding linear
    combination.

    Args:
        a: Block encoding of the input operator.
        gap: Filter half-width Δ, in (0,1); spectral components within Δ of the center are kept.
        degree: Chebyshev filter degree, a positive integer; the actual synthesis degree is twice it.
        center: Spectral position of the filter center, in (−1,1); 0 shifts no spectrum.

    Returns:
        BlockEncoding: Peaked filter block encoding, with the suppression attribute 1/T_d(r).
    """
    require_instance(a, BlockEncoding, "eigenstate_filter.a")
    if not (0 < gap < 1):
        raise ValidationError("The filter width gap must lie strictly between 0 and 1")
    if type(degree) is not int or degree < 1:
        raise ValidationError("The Chebyshev degree must be a positive integer")
    d2 = 2 * degree
    if d2 > _MAX_DEGREE:
        raise ValidationError(f"The synthesis degree {d2} exceeds the limit {_MAX_DEGREE}")
    shifted = a
    if center != 0.0:
        if not math.isfinite(center) or abs(center) >= 1:
            raise ValidationError("The filter center must lie strictly between −1 and 1")
        shifted = linear_combination(1.0, a, -complex(center), identity(a.width))
    r = (1 + gap**2) / (1 - gap**2)
    norm_d = math.cosh(degree * math.acosh(r))  # T_d(r)
    u = (-1.0 - 2 * gap**2 / (1 - gap**2), 0.0, 2.0 / (1 - gap**2))  # g(x²)
    t0: tuple[float | complex, ...]
    t1: tuple[float | complex, ...]
    t0, t1 = (1.0,), u
    for _ in range(2, degree + 1):
        t0, t1 = t1, _sub(_scale(2.0, _mul(u, t1)), t0)
    f = cast("tuple[float, ...]", _scale(1.0 / norm_d, t1 if degree >= 1 else t0))
    sat = math.sqrt(max(0.0, 1.0 - _eval(f, 1.0).real ** 2))
    phases = _synthesize_with_imag(f, [(0.0, 0.0, sat)], "eigenstate filter polynomial")
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


def _bessel_j(n: int, x: float) -> float:
    """Bessel function of the first kind J_n(x), a pure-Python power series."""
    term = (x / 2) ** n / math.factorial(n)
    total = term
    m = 0
    while abs(term) > 1e-18 * max(1.0, abs(total)) and m < 100000:
        term *= -((x / 2) ** 2) / ((m + 1) * (m + n + 1))
        total += term
        m += 1
    return total


def _jacobi_anger(t: float, error: float) -> tuple[tuple[float, ...], tuple[float, ...], int]:
    """Jacobi–Anger truncation of e^{itx}: returns the even-branch cos coefficients, the
    odd-branch sin coefficients, and the truncation degree K."""
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
    return cast("tuple[float, ...]", _trim(fc)), cast("tuple[float, ...]", _trim(fs)), k


def qsvt_hamiltonian_simulation(
    a: BlockEncoding, t: float, *, error: float = 0.01
) -> BlockEncoding:
    """QSVT block encoding of e^{itA/α}: the Jacobi–Anger even and odd branches are
    synthesized separately and then combined by an LCU.

    The even branch approximates cos(tx) and the odd branch approximates
    sin(tx); a common scale s leaves headroom for the imaginary completion of
    both branches. Each branch extracts the real part via (U_Φ + U_{−Φ})/2,
    and the two are finally combined by an LCU with weights 1 and i. The
    zero-signal block of the returned block encoding approximates
    e^{itA/α}/sim_scale, with sim_scale = 2s.

    Args:
        a: Block encoding of the evolution generator A.
        t: Evolution time, a nonzero finite real number.
        error: Jacobi–Anger truncation and synthesis error, in (0,1).

    Returns:
        BlockEncoding: Block encoding whose zero-signal block approximates e^{itA/α}/sim_scale.
    """
    require_instance(a, BlockEncoding, "qsvt_hamiltonian_simulation.a")
    if not (math.isfinite(t) and t != 0):
        raise ValidationError("The evolution time t must be a nonzero finite real number")
    if not (0 < error < 1):
        raise ValidationError("The approximation error must lie strictly between 0 and 1")
    fc, fs, k = _jacobi_anger(t, error)
    s = 1.5 * max(_sup_norm(fc), _sup_norm(fs), 1e-3)
    fc, fs = (
        cast("tuple[float, ...]", _scale(1.0 / s, fc)),
        cast("tuple[float, ...]", _scale(1.0 / s, fs)),
    )
    dc, ds = _deg(fc), _deg(fs)
    if dc % 2 or ds % 2 == 0:
        raise ValidationError("Unexpected parity in the Jacobi–Anger branches")

    def cos_imags() -> Iterator[tuple[float, ...]]:
        """Yield the imaginary completion polynomials to try for the cos branch, one at a time."""
        a0 = math.sqrt(max(0.0, 1.0 - _eval(fc, 0.0).real ** 2))
        a1 = math.sqrt(max(0.0, 1.0 - _eval(fc, 1.0).real ** 2))
        for m in range(1, dc // 2 + 1):
            yield tuple(
                (a0 if i == 0 else 0.0) + ((a1 - a0) if i == 2 * m else 0.0)
                for i in range(dc + 1)
            )

    def sin_imags() -> Iterator[tuple[float, ...]]:
        """Yield the imaginary completion polynomials to try for the sin branch, one at a time."""
        a1 = math.sqrt(max(0.0, 1.0 - _eval(fs, 1.0).real ** 2))
        for m in range(0, (ds - 1) // 2 + 1):
            yield tuple(a1 * (1.0 if i == 2 * m + 1 else 0.0) for i in range(ds + 1))

    phases_c = _synthesize_with_imag(fc, cos_imags(), "Hamiltonian simulation cos branch")
    phases_s = _synthesize_with_imag(fs, sin_imags(), "Hamiltonian simulation sin branch")
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


def fixed_point_search_phases(delta: float, degree: int) -> tuple[float, ...]:
    """Phase sequence for Yoder–Low–Chuang fixed-point amplitude amplification, in time
    order with length degree+1.

    The construction follows the closed-form YLC complementary polynomial: let
    L = degree, the number of block-encoding calls, which must be odd, and
    c = T_{1/L}(1/δ); then Q(x) = δc·R(c²(1−x²)) where R(u) = T_L(√u)/√u is
    an analytic polynomial, and P is obtained from the root-based spectral
    decomposition of 1 − (1−x²)Q². The realized success probability is
    exactly P_S(x) = 1 − δ² T_L²(c√(1−x²)): whenever ``|x| ≥ √(1−1/c²)``
    we have P_S ≥ 1 − δ², and the threshold decreases monotonically in L
    toward 0, the fixed-point property.

    Args:
        delta: Failure probability bound δ, in (0,1).
        degree: Amplification degree L, the number of block-encoding calls, a positive odd
            integer not exceeding half the synthesis limit.

    Returns:
        tuple[float, ...]: Phase sequence in time order, with length degree+1.
    """
    if not (0 < delta < 1):
        raise ValidationError("The fixed-point search error δ must lie strictly between 0 and 1")
    if type(degree) is not int or degree < 1 or degree % 2 == 0:
        raise ValidationError("The fixed-point search degree, the number of block-encoding calls,"
                              " must be a positive odd integer")
    if degree > _MAX_DEGREE // 2:
        raise ValidationError(f"The fixed-point search degree exceeds the limit"
                              f" {_MAX_DEGREE // 2}")
    L = degree
    c = math.cosh(math.acosh(1.0 / delta) / L)
    tcoeff = _chebyshev_t(L)
    qpoly: tuple[float | complex, ...] = (0.0,)
    for m in range((L + 1) // 2):
        term: tuple[float | complex, ...] = (tcoeff[2 * m + 1] * c ** (2 * m),)
        for _ in range(m):
            term = _mul(term, (1.0, 0.0, -1.0))
        qpoly = _add(qpoly, term)
    qpoly = _trim(_scale(delta * c, qpoly))
    fpoly = _trim(_sub((1.0,), _mul((1.0, 0.0, -1.0), _mul(qpoly, qpoly))))
    if abs(fpoly[0]) > 1e-8 or abs(fpoly[1]) > 1e-8:
        raise ValidationError("Abnormal YLC complementary polynomial construction: the zero root"
                              " is missing")
    ft = _trim(fpoly[2:])  # divide by x², the double root at zero
    clusters = _cluster(_roots(ft))
    factors: list[tuple[float | complex, ...]] = []
    used: set[int] = set()
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
            raise ValidationError("Abnormal root structure in the YLC spectral decomposition")
        for _ in range(m):
            factors.append((-rt * rt, 0.0, 1.0))
        used.update(sib)
        used.add(i)
    if 1 + 2 * len(factors) != L:
        raise ValidationError("Abnormal factor count in the YLC spectral decomposition")
    ppoly: tuple[float | complex, ...] = (0.0, 1.0)
    for fct in factors:
        ppoly = _mul(ppoly, fct)
    ppoly = _trim(_scale(abs(qpoly[-1]), ppoly))
    phases = _strip(ppoly, qpoly, L)
    err = max(
        abs(abs(qsp_response(x, phases)) ** 2 - (1.0 - (1 - x * x) * _eval(qpoly, x) ** 2))
        for x in _grid(512)
    )
    if err > _STRIP_TOL:
        raise ValidationError(f"Fixed-point search phase self-check failed with error {err:.2e}")
    return phases


def fixed_point_search(a: BlockEncoding, delta: float, degree: int) -> BlockEncoding:
    """QSVT assembly for fixed-point amplitude amplification: applies the phase sequence
    from fixed_point_search_phases to a block encoding.

    The zero-signal block is the complex polynomial P(A/α), whose success
    probability ``|P(x)|²`` satisfies the YLC fixed-point guarantee; the
    threshold attribute gives the amplification threshold √(1−1/c²).

    Args:
        a: The block encoding to amplify.
        delta: Failure probability bound δ, in (0,1).
        degree: Amplification degree, the number of block-encoding calls, a positive odd integer.

    Returns:
        BlockEncoding: Fixed-point amplification block encoding whose zero-signal block is P(A/α).
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
