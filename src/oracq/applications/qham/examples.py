"用于自动推导的 PDE 案例；参数与空间离散保持分离。"

from oracq.applications.qham.pde import Field, Known, PolynomialPDE


def example_pde(name: str) -> PolynomialPDE:
    """按名字返回内置的示例 PDE。

    可用案例为 ``burgers``、``kdv``、``reaction``、``coupled`` 和
    ``vector_burgers_2d``；除 ``vector_burgers_2d`` 使用 ``("x", "y")``
    空间轴外，其余均为一维。

    Args:
        name: 案例名。

    Returns:
        PolynomialPDE: 对应案例的 PDE，``label`` 即案例名。

    Raises:
        ValueError: 案例名未知。
    """
    u, v = Field("u"), Field("v")
    if name == "burgers":
        return PolynomialPDE.from_equations(
            {"u": 0.1 * u.d("x", 2) - u * u.d("x") + Known("f")}, label=name
        )
    if name == "kdv":
        return PolynomialPDE.from_equations({"u": -u.d("x", 3) - 6 * u * u.d("x")}, label=name)
    if name == "reaction":
        return PolynomialPDE.from_equations(
            {"u": 0.1 * u.d("x", 2) + u - u**3 + Known("f")}, label=name
        )
    if name == "coupled":
        return PolynomialPDE.from_equations(
            {
                "u": -0.2 * u + 0.1 * v + 0.2 * u * v + Known("f"),
                "v": 0.1 * u - 0.3 * v - 0.1 * v * v,
            },
            label=name,
        )
    if name == "vector_burgers_2d":
        return PolynomialPDE.from_equations(
            {
                "u": 0.1 * (u.d("x", 2) + u.d("y", 2)) - u * u.d("x") - v * u.d("y"),
                "v": 0.1 * (v.d("x", 2) + v.d("y", 2)) - u * v.d("x") - v * v.d("y"),
            },
            axes=("x", "y"),
            label=name,
        )
    raise ValueError("未知 PDE 案例：" + name)
