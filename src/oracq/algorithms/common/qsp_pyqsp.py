"""pyqsp-backed QSP phase synthesis: an optional PhaseSynthesizer component.

This module adapts pyqsp (https://github.com/ichuang/pyqsp, MIT-licensed, a
pure numerical package) to the PhaseSynthesizer contract of
oracq.algorithms.common.qsvt. It is imported only at the usage entry: nothing
in the core imports it, and pyqsp is an optional dependency (extra ``pyqsp``).

Convention mapping, verified pointwise against qsp_response on Chebyshev
grids. pyqsp's symmetric solver returns phases Θ = [θ_0, …, θ_d] for the Wx
signal convention, building U = e^{iθ_0 Z} ∏ W_x(x) e^{iθ_k Z} with θ_0 the
last-applied rotation and the target polynomial appearing as Im U_00. The two
signal operators obey W_x(x) = S(−π/4) R(x) (i·S(−π/4)), so with

    φ_j = θ_{d−j} − π/4  (j = 0 or j = d),   φ_j = θ_{d−j} − π/2  (middle),

the oracq reflection-convention product realizes p(x) = i^{−d}·U_00(x). To pin
Re p = f the adapter additionally negates the target when d ≡ 3 (mod 4) and
adds γ = −π/2 (d ≡ 0), +π/2 (d ≡ 2), or 0 (odd d) to the first-applied phase
φ_0, which multiplies p by e^{iγ}.

Scope: the adapter supports the free-completion mode only (Re p = f, the
complementary part is pyqsp's choice); a pinned imaginary completion is
rejected with ValidationError. Input coefficients follow the library contract
(ascending monomial coefficients) and are converted to the Chebyshev basis
internally; like all coefficient-level pipelines, the conversion and the
library-side grid validation are reliable only in the same moderate-degree
domain as the bundled route. For high-degree work, construct the target
natively in the Chebyshev basis and call chebyshev_phases_via_pyqsp directly.
"""

from __future__ import annotations

import contextlib
import io
import math
from collections.abc import Sequence

from oracq.infrastructure.ir import ValidationError

__all__ = ["PyqspSynthesizer", "chebyshev_phases_via_pyqsp"]


def _trim_coeffs(coeffs: Sequence[float]) -> tuple[float, ...]:
    """Drop numerically negligible trailing coefficients."""
    c = [float(v) for v in coeffs]
    while len(c) > 1 and c[-1] == 0.0:
        c.pop()
    return tuple(c)


def _wx_phases(cheb_coeffs: Sequence[float]) -> tuple[float, ...]:
    """Run pyqsp's symmetric-QSP Newton solver on ascending Chebyshev coefficients."""
    try:
        from pyqsp.angle_sequence import QuantumSignalProcessingPhases
    except ImportError as exc:
        raise ImportError(
            "qsp_pyqsp requires the optional pyqsp dependency: pip install oracq[pyqsp]"
        ) from exc
    with contextlib.redirect_stdout(io.StringIO()):  # silence the solver's progress log
        full_phases, _, _ = QuantumSignalProcessingPhases(
            list(cheb_coeffs), method="sym_qsp", chebyshev_basis=True
        )
    return tuple(float(v) for v in full_phases)


def chebyshev_phases_via_pyqsp(cheb_coeffs: Sequence[float]) -> tuple[float, ...]:
    """Synthesize oracq-convention phases for a real target in ascending Chebyshev coefficients.

    The target must have definite parity (its degree's parity) and satisfy
    ``|f| <= 1`` on ``[-1, 1]`` with enough headroom for a complementary
    polynomial.
    The returned sequence Φ (time order, length d+1) realizes Re p = f under
    qsp_response; the imaginary completion is pyqsp's choice.

    Args:
        cheb_coeffs: Ascending Chebyshev coefficients of the real target, constant term first.

    Returns:
        tuple[float, ...]: Phase sequence in time order in the reflection convention.
    """
    f = _trim_coeffs(cheb_coeffs)
    d = len(f) - 1
    if d < 1:
        raise ValidationError("The Chebyshev target must have degree at least 1")
    target = tuple(-v for v in f) if d % 4 == 3 else f
    thetas = _wx_phases(target)
    if len(thetas) != d + 1:
        raise ValidationError(f"pyqsp returned {len(thetas)} phases for a degree-{d} target")
    gamma = (-math.pi / 2) if d % 4 == 0 else (math.pi / 2 if d % 4 == 2 else 0.0)
    phases = [
        thetas[d - j] - (math.pi / 4 if j in (0, d) else math.pi / 2) for j in range(d + 1)
    ]
    phases[0] += gamma
    return tuple(phases)


class PyqspSynthesizer:
    """PhaseSynthesizer backed by pyqsp's symmetric-QSP Newton solver (free-completion mode).

    Follows the PhaseSynthesizer contract of oracq.algorithms.common.qsvt:
    accepts ascending monomial coefficients (converted to the Chebyshev basis
    internally) and returns time-ordered phases in the reflection convention.
    A pinned imaginary completion (``imag`` given) is not supported, since
    pyqsp chooses the complementary part itself.
    """

    def __call__(
        self, coeffs: Sequence[float], imag: Sequence[float] | None = None
    ) -> tuple[float, ...]:
        """Synthesize phases realizing Re p = f, with f given in monomial coefficients."""
        if imag is not None:
            raise ValidationError(
                "The pyqsp synthesizer supports the free-completion mode only:"
                " the imag argument must be omitted"
            )
        from numpy.polynomial import Chebyshev, Polynomial

        f = _trim_coeffs(coeffs)
        cheb = Polynomial(f).convert(kind=Chebyshev)
        return chebyshev_phases_via_pyqsp(tuple(cheb.coef))
