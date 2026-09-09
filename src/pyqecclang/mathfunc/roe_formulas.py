"""普通、无副作用的 Roe 数学函数；可经典调用，也可由 compile_function 编译。"""

import math


def pick3(index, first, second, third):
    if index == 0:
        return first
    if index == 1:
        return second
    if index == 2:
        return third
    return 0.0


def conserved_to_primitive(rho, momentum, energy, gamma):
    velocity = momentum / rho
    pressure = (gamma - 1) * (energy - 0.5 * momentum * velocity)
    enthalpy = (energy + pressure) / rho
    return velocity, pressure, enthalpy


def euler_entry(velocity, enthalpy, row, col, gamma):
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


def entropy_absolute(eigenvalue, delta):
    value = abs(eigenvalue)
    if delta <= 0:
        return value
    if value >= delta:
        return value
    return (eigenvalue * eigenvalue + delta * delta) / (2 * delta)


def frozen_roe_face(rho_l, m_l, e_l, rho_r, m_r, e_r, row, col, gamma=1.4, entropy_delta=0.125):
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
