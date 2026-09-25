"Plain side-effect-free Roe math functions; callable classically or compilable by compile_function."

from __future__ import annotations

import math


def pick3(index: int, first: float, second: float, third: float) -> float:
    """Select one of three values by integer index; returns ``0.0`` when the index is outside 0..2.

    Args:
        index: Selection index in 0..2; out of range selects no branch.
        first: Value selected for index 0.
        second: Value selected for index 1.
        third: Value selected for index 2.

    Returns:
        float: The selected value; ``0.0`` when the index is out of range.
    """
    if index == 0:
        return first
    if index == 1:
        return second
    if index == 2:
        return third
    return 0.0


def conserved_to_primitive(
    rho: float, momentum: float, energy: float, gamma: float
) -> tuple[float, float, float]:
    """Convert one-dimensional Euler conserved variables to primitive variables.

    Args:
        rho: Density.
        momentum: Momentum density.
        energy: Total energy density.
        gamma: Ratio of specific heats.

    Returns:
        tuple: ``(velocity, pressure, enthalpy)``, with enthalpy defined as total energy plus pressure divided by density.
    """
    velocity = momentum / rho
    pressure = (gamma - 1) * (energy - 0.5 * momentum * velocity)
    enthalpy = (energy + pressure) / rho
    return velocity, pressure, enthalpy


def euler_entry(velocity: float, enthalpy: float, row: int, col: int, gamma: float) -> float:
    """Return the requested element of the one-dimensional Euler flux Jacobian matrix.

    Elements are expressed using only flow velocity and enthalpy; row and column
    indices follow the ordering of the conserved variables ``(rho, momentum, energy)``,
    numbered from zero.

    Args:
        velocity: Flow velocity.
        enthalpy: Specific enthalpy, defined as in the return value of ``conserved_to_primitive``.
        row: Row index, in 0..2.
        col: Column index, in 0..2.
        gamma: Ratio of specific heats.

    Returns:
        float: The element at row ``row`` and column ``col`` of the Jacobian matrix.
    """
    square = velocity * velocity
    first = pick3(col, 0.0, 1.0, 0.0)
    second = pick3(col, 0.5 * (gamma - 3) * square, (3 - gamma) * velocity, gamma - 1)
    third = pick3(
        col,
        velocity * (0.5 * (gamma - 1) * square - enthalpy),
        enthalpy - (gamma - 1) * square,
        gamma * velocity,
    )
    return pick3(row, first, second, third)


def entropy_absolute(eigenvalue: float, delta: float) -> float:
    """Absolute value of an eigenvalue with the Harten entropy fix.

    When ``delta`` is positive and the absolute value of the eigenvalue is smaller,
    the smoothed value ``(eigenvalue**2+delta**2)/(2*delta)`` is returned;
    otherwise the plain absolute value is returned.

    Args:
        eigenvalue: The eigenvalue.
        delta: Entropy fix threshold; the fix is disabled when non-positive.

    Returns:
        float: The entropy-fixed absolute value of the eigenvalue.
    """
    value = abs(eigenvalue)
    if delta <= 0:
        return value
    if value >= delta:
        return value
    return (eigenvalue * eigenvalue + delta * delta) / (2 * delta)


def frozen_roe_face(
    rho_l: float,
    m_l: float,
    e_l: float,
    rho_r: float,
    m_r: float,
    e_r: float,
    row: int,
    col: int,
    gamma: float = 1.4,
    entropy_delta: float = 0.125,
) -> tuple[float, float]:
    """Compute the element at the given row and column of the frozen-Roe left and right coefficient matrices.

    Roe averaging of the left and right conserved states yields the mean flow
    velocity and enthalpy; the eigenvalues ``u-c``, ``u`` and ``u+c`` are made
    absolute with the entropy fix, and the elements of ``|A|`` are assembled as
    ``R*diag(a)*R**-1``. Combined with the flux Jacobian elements on both sides
    given by ``euler_entry``, this returns the values of ``0.5*(A_l+|A|)`` and
    ``0.5*(A_r-|A|)`` at ``row`` and ``col``, as the Roe flux share coefficients
    applied to the left and right states.

    Args:
        rho_l: Left-side density.
        m_l: Left-side momentum density.
        e_l: Left-side total energy density.
        rho_r: Right-side density.
        m_r: Right-side momentum density.
        e_r: Right-side total energy density.
        row: Row index of the matrix element, in 0..2.
        col: Column index of the matrix element, in 0..2.
        gamma: Ratio of specific heats.
        entropy_delta: Harten entropy fix threshold.

    Returns:
        tuple: ``(left, right)``, the elements of the left and right coefficient matrices at the given row and column.
    """
    ul, pl, hl = conserved_to_primitive(rho_l, m_l, e_l, gamma)
    ur, pr, hr = conserved_to_primitive(rho_r, m_r, e_r, gamma)
    wl = math.sqrt(rho_l)
    wr = math.sqrt(rho_r)
    total = wl + wr
    u = (wl * ul + wr * ur) / total
    h = (wl * hl + wr * hr) / total
    u2 = u * u
    c2 = (gamma - 1) * (h - 0.5 * u2)
    c = math.sqrt(c2)
    beta = (gamma - 1) / c2
    inverse_c = 0.5 / c
    bu = beta * u
    bu2 = beta * u2
    r0 = pick3(row, 1.0, u - c, h - u * c)
    r1 = pick3(row, 1.0, u, 0.5 * u2)
    r2 = pick3(row, 1.0, u + c, h + u * c)
    inv0 = pick3(col, 0.25 * bu2 + u * inverse_c, -0.5 * bu - inverse_c, 0.5 * beta)
    inv1 = pick3(col, 1 - 0.5 * bu2, bu, -beta)
    inv2 = pick3(col, 0.25 * bu2 - u * inverse_c, -0.5 * bu + inverse_c, 0.5 * beta)
    a0 = entropy_absolute(u - c, entropy_delta)
    a1 = entropy_absolute(u, entropy_delta)
    a2 = entropy_absolute(u + c, entropy_delta)
    absolute = r0 * a0 * inv0 + r1 * a1 * inv1 + r2 * a2 * inv2
    left = 0.5 * (euler_entry(ul, hl, row, col, gamma) + absolute)
    right = 0.5 * (euler_entry(ur, hr, row, col, gamma) - absolute)
    return left, right
