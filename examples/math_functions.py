"""可直接调用的经典纯函数；这些函数不依赖量子 Builder。"""

import cmath
import math


def pressure(rho, momentum, energy, gamma=1.4):
    velocity = momentum / rho
    return (gamma - 1) * (energy - 0.5 * momentum * velocity)


def roe_speed(rho_l, momentum_l, rho_r, momentum_r):
    left = math.sqrt(rho_l)
    right = math.sqrt(rho_r)
    return (left * momentum_l / rho_l + right * momentum_r / rho_r) / (left + right)


def phase_response(z: complex):
    return cmath.exp(1j * z) / (1 + z * z)


def guarded_reciprocal(x):
    if x == 0:
        return 0.0
    return 1.0 / x


def polynomial(x, order=3):
    result = 0.0
    for k in range(order):
        result = result + x**k
    return result
