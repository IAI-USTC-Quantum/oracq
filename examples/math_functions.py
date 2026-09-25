"""Directly callable classical pure functions; these functions do not depend on the quantum Builder."""

from __future__ import annotations

import cmath
import math


def pressure(rho: float, momentum: float, energy: float, gamma: float = 1.4) -> float:
    """Compute the ideal-gas pressure from density, momentum, and total energy."""
    velocity = momentum / rho
    return (gamma - 1) * (energy - 0.5 * momentum * velocity)


def roe_speed(rho_l: float, momentum_l: float, rho_r: float, momentum_r: float) -> float:
    """Compute the Roe-averaged advective velocity of the left and right states."""
    left = math.sqrt(rho_l)
    right = math.sqrt(rho_r)
    return (left * momentum_l / rho_l + right * momentum_r / rho_r) / (left + right)


def phase_response(z: complex) -> complex:
    """Return the phase response exp(i z)/(1+z^2)."""
    return cmath.exp(1j * z) / (1 + z * z)


def guarded_reciprocal(x: float) -> float:
    """Guarded reciprocal; returns 0 for zero input."""
    if x == 0:
        return 0.0
    return 1.0 / x


def polynomial(x: float, order: int = 3) -> float:
    """Compute the sum of x raised to the powers 0 through order-1."""
    result = 0.0
    for k in range(order):
        result = result + x**k
    return result
